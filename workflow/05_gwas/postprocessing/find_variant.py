#!/usr/bin/env python3
"""
Map kmers to M. tuberculosis reference and call SNPs/indels.
Annotates each SNP as synonymous / nonsynonymous / stop_gained /
stop_lost / start_lost / intergenic using a GFF3 annotation file.

Libraries used:
  - gffutils   : parses and queries GFF/GFF3 files via SQLite
  - BioPython  : bacterial translation table 11
  - pyfaidx    : fast reference sequence lookup
  - pysam      : SAM/BAM parsing

Usage:
  python snp_tb.py kmers.fa ref.fasta annotation.gff out.tsv [--threads 4]

Install deps:
  pip install gffutils biopython pyfaidx pysam
"""

import argparse
import subprocess
import tempfile
import os
import sys

import pysam
from pyfaidx import Fasta
import gffutils
from Bio.Data import CodonTable

# ── Bacterial (translation table 11) codon table via BioPython ─────────────
_BACT_TABLE = CodonTable.unambiguous_dna_by_id[11]

# Alternative start codons for table 11 (GTG, TTG, CTG → fMet in vivo)
_ALT_START_CODONS = set(_BACT_TABLE.start_codons)  # {'TTG', 'CTG', 'GTG', 'ATG'}


def translate_codon(codon: str, is_start: bool = False) -> str:
    """
    Translate a 3-nt codon using bacterial table 11.
    If is_start=True, any valid start codon → 'M' (fMet).
    Returns 'Stop' for stop codons, single-letter AA, or '?' if unknown.
    """
    codon = codon.upper()
    if codon in _BACT_TABLE.stop_codons:
        return "Stop"
    if is_start and codon in _ALT_START_CODONS:
        return "M"
    return _BACT_TABLE.forward_table.get(codon, "?")


def reverse_complement(seq: str) -> str:
    comp = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(comp)[::-1]


# ── GFF loading ───────────────────────────────────────────────────────────
def load_gff(gff_path: str) -> gffutils.FeatureDB:
    """
    Build (or reuse) a gffutils SQLite database from the GFF file.
    """
    db_path = gff_path + ".db"
    if not os.path.exists(db_path):
        print(f"Indexing GFF → {db_path}  (one-time, may take a minute) …")
        gffutils.create_db(
            gff_path,
            dbfn=db_path,
            force=False,
            keep_order=True,
            merge_strategy="merge",
            sort_attribute_values=True,
            disable_infer_genes=False,
            disable_infer_transcripts=False,
        )
    db = gffutils.FeatureDB(db_path, keep_order=True)
    return db


# ── Core: classify a SNP ─────────────────────────────────────────────────
def classify_snp(
    chrom: str,
    pos1: int,
    ref_base: str,
    alt_base: str,
    db: gffutils.FeatureDB,
    ref_fasta: Fasta,
) -> dict:
    """
    Given a SNP (1-based position), look up overlapping CDS in the GFF,
    reconstruct the affected codon, and return an annotation dict.

    Handles:
      - Strand (reverse-complement for minus-strand genes)
      - GFF3 phase/frame field for reading-frame offset
      - Start-codon variants (start_lost)
    """
    empty = dict(
        region="intergenic",
        effect="intergenic",
        ref_codon="",
        alt_codon="",
        ref_aa="",
        alt_aa="",
        gene="",
        locus_tag="",
        codon_pos="",
        aa_pos="",
    )

    # Query gffutils for CDS features overlapping this 1-based position
    try:
        hits = list(
            db.region(
                seqid=chrom,
                start=pos1,
                end=pos1,
                featuretype="CDS",
                completely_within=False,
            )
        )
    except Exception:
        return empty

    if not hits:
        return empty

    # Use the first CDS hit (TB genes are mostly non-overlapping)
    cds = hits[0]
    strand = cds.strand

    # Pull gene / locus_tag from attributes (NCBI GFF3 style)
    gene_name = cds.attributes.get("gene", [""])[0]
    locus_tag = cds.attributes.get("locus_tag", [""])[0]
    if not gene_name:
        gene_name = locus_tag

    # ── GFF3 phase field ─────────────────────────────────────────────────
    # phase (0, 1, or 2) = number of bases to skip at the START of this
    # CDS feature to reach the first complete codon.
    # For single-exon prokaryotic genes this is almost always 0.
    try:
        phase = int(cds.frame) if cds.frame not in (".", None) else 0
    except (ValueError, TypeError):
        phase = 0

    # ── Get the CDS sequence from the reference ─────────────────────────
    # pyfaidx is 0-based half-open; gffutils start/end are 1-based inclusive
    cds_seq = ref_fasta[chrom][cds.start - 1 : cds.end].seq.upper()

    # Position of our SNP within the CDS (0-based offset from CDS start)
    snp_cds_pos = pos1 - cds.start  # 0-based within cds_seq

    if strand == "-":
        cds_seq = reverse_complement(cds_seq)
        snp_cds_pos = len(cds_seq) - 1 - snp_cds_pos
        ref_base_in_cds = reverse_complement(ref_base)
        alt_base_in_cds = reverse_complement(alt_base)
    else:
        ref_base_in_cds = ref_base.upper()
        alt_base_in_cds = alt_base.upper()

    # Bounds check
    if not (0 <= snp_cds_pos < len(cds_seq)):
        return dict(
            region="CDS", effect="out_of_bounds",
            ref_codon="", alt_codon="", ref_aa="", alt_aa="",
            gene=gene_name, locus_tag=locus_tag, codon_pos="", aa_pos="",
        )

    # Sanity check: does the reference base match?
    if cds_seq[snp_cds_pos] != ref_base_in_cds:
        print(
            f"  WARNING: ref mismatch at {chrom}:{pos1} — "
            f"expected '{ref_base_in_cds}' in CDS but found '{cds_seq[snp_cds_pos]}' "
            f"(gene={gene_name}, strand={strand})",
            file=sys.stderr,
        )

    # ── Codon arithmetic (accounting for phase) ─────────────────────────
    # Skip 'phase' bases at the start to reach the first in-frame codon
    coding_offset = snp_cds_pos - phase
    if coding_offset < 0:
        # SNP falls in the partial-codon leader — can't classify
        return dict(
            region="CDS", effect="partial_codon",
            ref_codon="", alt_codon="", ref_aa="", alt_aa="",
            gene=gene_name, locus_tag=locus_tag, codon_pos="", aa_pos="",
        )

    codon_index = coding_offset // 3
    pos_in_codon = coding_offset % 3
    codon_start = phase + codon_index * 3

    if codon_start + 3 > len(cds_seq):
        return dict(
            region="CDS", effect="incomplete_codon",
            ref_codon="", alt_codon="", ref_aa="", alt_aa="",
            gene=gene_name, locus_tag=locus_tag, codon_pos="", aa_pos="",
        )

    ref_codon = cds_seq[codon_start : codon_start + 3]
    alt_codon = (
        ref_codon[:pos_in_codon] + alt_base_in_cds + ref_codon[pos_in_codon + 1 :]
    )

    # Is this the first codon? (codon_index == 0 means start codon)
    is_start = codon_index == 0

    ref_aa = translate_codon(ref_codon, is_start=is_start)
    alt_aa = translate_codon(alt_codon, is_start=is_start)

    # ── Classify the effect ──────────────────────────────────────────────
    if is_start and alt_codon.upper() not in _ALT_START_CODONS:
        effect = "start_lost"
    elif ref_aa == alt_aa:
        effect = "synonymous"
    elif alt_aa == "Stop" and ref_aa != "Stop":
        effect = "stop_gained"
    elif ref_aa == "Stop" and alt_aa != "Stop":
        effect = "stop_lost"
    else:
        effect = "nonsynonymous"

    return dict(
        region="CDS",
        effect=effect,
        ref_codon=ref_codon,
        alt_codon=alt_codon,
        ref_aa=ref_aa,
        alt_aa=alt_aa,
        gene=gene_name,
        locus_tag=locus_tag,
        codon_pos=str(pos_in_codon + 1),      # 1-based within codon
        aa_pos=str(codon_index + 1),           # 1-based amino acid position
    )


# ── Main ──────────────────────────────────────────────────────────────────
def main():
    ap = argparse.ArgumentParser(
        description="Kmer → BWA alignment → SNP calling + annotation "
        "for M. tuberculosis (or any prokaryote with a GFF3)."
    )
    ap.add_argument("kmers_fa", help="Kmer sequences (FASTA)")
    ap.add_argument("ref", help="Reference genome FASTA")
    ap.add_argument("gff", help="GFF3 annotation file (e.g. from NCBI)")
    ap.add_argument("out", help="Output TSV file")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument(
        "--min_mapq",
        type=int,
        default=20,
        help="Minimum mapping quality (default: 20; use 0 to include multi-mappers)",
    )
    args = ap.parse_args()

    # 1. Index reference if needed
    if not os.path.exists(f"{args.ref}.bwt"):
        print("Indexing reference with BWA …")
        subprocess.run(["bwa", "index", args.ref], check=True)

    # 2. Load reference FASTA
    ref = Fasta(args.ref, rebuild=False)

    # 3. Load GFF
    print("Loading GFF annotation …")
    db = load_gff(args.gff)
    print("GFF loaded.")

    # 4. Read kmer sequences
    kmer_fa = Fasta(args.kmers_fa)
    kmers = {name: str(kmer_fa[name][:].seq) for name in kmer_fa.keys()}

    # 5. Align + call variants
    with tempfile.TemporaryDirectory() as tmp:
        sam = os.path.join(tmp, "k.sam")

        # Fix: properly close the file handle for stdout redirection
        with open(sam, "w") as sam_fh:
            subprocess.run(
                ["bwa", "mem", "-t", str(args.threads), "-a", args.ref, args.kmers_fa],
                stdout=sam_fh,
                stderr=subprocess.DEVNULL,
                check=True,
            )

        header_cols = [
            "kmer_id", "kmer", "chrom", "pos", "type", "ref", "alt",
            "strand", "mapq", "cigar", "region", "effect",
            "ref_codon", "alt_codon", "ref_aa", "alt_aa",
            "gene", "locus_tag", "codon_pos", "aa_pos",
        ]

        seen_variants = set()  # deduplicate across secondary alignments

        with pysam.AlignmentFile(sam) as s, open(args.out, "w") as out:
            out.write("\t".join(header_cols) + "\n")

            for aln in s:
                # ── Filter: unmapped, secondary, supplementary ───────────
                if aln.is_unmapped:
                    continue
                if aln.is_secondary or aln.is_supplementary:
                    continue
                if aln.mapping_quality < args.min_mapq:
                    continue

                kid = aln.query_name
                kseq = kmers.get(kid, aln.query_sequence or "")
                chrom = aln.reference_name
                strand = "-" if aln.is_reverse else "+"
                mapq = aln.mapping_quality
                cigar = aln.cigarstring

                rpos = aln.reference_start  # 0-based
                qpos = 0
                qseq = aln.query_sequence.upper()
                has_variant = False

                def emit(pos1, vtype, rb, ab):
                    """Write one variant row, with annotation for SNPs."""
                    # Deduplicate: same kmer + same position + same variant
                    dedup_key = (kid, chrom, pos1, vtype, rb, ab)
                    if dedup_key in seen_variants:
                        return
                    seen_variants.add(dedup_key)

                    ann = dict(
                        region="", effect="", ref_codon="", alt_codon="",
                        ref_aa="", alt_aa="", gene="", locus_tag="",
                        codon_pos="", aa_pos="",
                    )
                    if vtype == "SNP":
                        ann = classify_snp(chrom, pos1, rb, ab, db, ref)

                    fields = [
                        kid, kseq, chrom, str(pos1), vtype, rb, ab,
                        strand, str(mapq), cigar,
                        ann["region"], ann["effect"],
                        ann["ref_codon"], ann["alt_codon"],
                        ann["ref_aa"], ann["alt_aa"],
                        ann["gene"], ann["locus_tag"],
                        ann.get("codon_pos", ""), ann.get("aa_pos", ""),
                    ]
                    out.write("\t".join(fields) + "\n")

                for op, ln in aln.cigartuples or []:
                    if op == 4:  # soft clip
                        qpos += ln
                    elif op == 5:  # hard clip — consumes neither
                        pass
                    elif op in (0, 7, 8):  # M / = / X
                        for i in range(ln):
                            rb = ref[chrom][rpos : rpos + 1].seq.upper()
                            qb = qseq[qpos]
                            if rb != qb and rb != "N" and qb != "N":
                                has_variant = True
                                # 1-based position for output
                                emit(rpos + 1, "SNP", rb, qb)
                            rpos += 1
                            qpos += 1
                    elif op == 1:  # insertion
                        has_variant = True
                        ins = qseq[qpos : qpos + ln]
                        anc = ref[chrom][rpos - 1 : rpos].seq.upper()
                        # 1-based: anchor base is at rpos (0-based) → rpos+0 in
                        # 0-based = rpos in 1-based? No.
                        # The anchor is at 0-based (rpos-1), so 1-based = rpos.
                        emit(rpos, "INS", anc, anc + ins)
                        qpos += ln
                    elif op == 2:  # deletion
                        has_variant = True
                        dels = ref[chrom][rpos : rpos + ln].seq.upper()
                        anc = ref[chrom][rpos - 1 : rpos].seq.upper()
                        # Anchor at 0-based (rpos-1) → 1-based = rpos
                        emit(rpos, "DEL", anc + dels, anc)
                        rpos += ln

                if not has_variant:
                    fields = [
                        kid, kseq, chrom, str(aln.reference_start + 1),
                        "MATCH", ".", ".", strand, str(mapq), cigar,
                        "", "", "", "", "", "", "", "", "", "",
                    ]
                    out.write("\t".join(fields) + "\n")

    print(f"Done. Results written to {args.out}")
    print(f"  Total kmers processed: {len(kmers)}", file=sys.stderr)


if __name__ == "__main__":
    main()
