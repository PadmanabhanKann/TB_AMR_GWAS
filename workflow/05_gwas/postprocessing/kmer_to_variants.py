#!/usr/bin/env python3
"""
Map kmers (from fasta) to reference and call SNPs/indels.
Usage: python kmer_to_variants.py kmers.fa ref.fasta out.tsv [--threads 8]
"""
import argparse, subprocess, tempfile, os
import pysam
from pyfaidx import Fasta

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("kmers_fa", help="Kmer sequences in fasta format")
    ap.add_argument("ref", help="Reference fasta")
    ap.add_argument("out", help="Output TSV")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--min_mapq", type=int, default=0)
    args = ap.parse_args()

    # Index if needed
    if not os.path.exists(f"{args.ref}.bwt"):
        subprocess.run(["bwa", "index", args.ref], check=True)
    
    ref = Fasta(args.ref, rebuild=False)
    
    # Read kmer fasta
    kmer_fa = Fasta(args.kmers_fa)
    kmers = {name: str(kmer_fa[name][:].seq) for name in kmer_fa.keys()}
    
    with tempfile.TemporaryDirectory() as tmp:
        # Map
        sam = f"{tmp}/k.sam"
        subprocess.run(["bwa", "mem", "-t", str(args.threads), "-a", args.ref, args.kmers_fa],
                       stdout=open(sam, "w"), check=True)
        
        # Parse variants
        with pysam.AlignmentFile(sam) as s, open(args.out, "w") as out:
            out.write("kmer_id\tkmer\tchrom\tpos\ttype\tref\talt\tstrand\tmapq\tcigar\n")
            
            for aln in s:
                if aln.is_unmapped or aln.mapping_quality < args.min_mapq:
                    continue
                
                kid = aln.query_name
                kseq = kmers[kid]
                chrom = aln.reference_name
                strand = "-" if aln.is_reverse else "+"
                mapq = aln.mapping_quality
                cigar = aln.cigarstring
                
                rpos = aln.reference_start
                qpos = 0
                qseq = aln.query_sequence.upper()
                
                has_variant = False
                
                for op, ln in aln.cigartuples or []:
                    if op == 4:  # soft clip
                        qpos += ln
                    elif op in (0, 7, 8):  # M/=/X
                        for i in range(ln):
                            rb = ref[chrom][rpos:rpos+1].seq.upper()
                            qb = qseq[qpos]
                            if rb != qb and rb != "N" and qb != "N":
                                has_variant = True
                                out.write(f"{kid}\t{kseq}\t{chrom}\t{rpos+1}\tSNP\t{rb}\t{qb}\t{strand}\t{mapq}\t{cigar}\n")
                            rpos += 1
                            qpos += 1
                    elif op == 1:  # insertion
                        has_variant = True
                        ins = qseq[qpos:qpos+ln]
                        anc = ref[chrom][rpos-1:rpos].seq.upper()
                        out.write(f"{kid}\t{kseq}\t{chrom}\t{rpos}\tINS\t{anc}\t{anc}{ins}\t{strand}\t{mapq}\t{cigar}\n")
                        qpos += ln
                    elif op == 2:  # deletion
                        has_variant = True
                        dels = ref[chrom][rpos:rpos+ln].seq.upper()
                        anc = ref[chrom][rpos-1:rpos].seq.upper()
                        out.write(f"{kid}\t{kseq}\t{chrom}\t{rpos}\tDEL\t{anc}{dels}\t{anc}\t{strand}\t{mapq}\t{cigar}\n")
                        rpos += ln
                
                # If no variants found, output as MATCH
                if not has_variant:
                    out.write(f"{kid}\t{kseq}\t{chrom}\t{aln.reference_start+1}\tMATCH\t.\t.\t{strand}\t{mapq}\t{cigar}\n")

if __name__ == "__main__":
    main()
