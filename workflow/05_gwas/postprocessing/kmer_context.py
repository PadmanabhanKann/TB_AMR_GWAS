#!/usr/bin/env python3
"""
Analyze genomic context of mutually exclusive kmers to find structural variants.
Usage: python check_kmer_context.py kmers.fa ref.fasta
"""
import argparse
from pyfaidx import Fasta
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

def reverse_complement(seq):
    comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join(comp.get(b, b) for b in reversed(seq.upper()))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("kmers_fa", help="Kmer fasta file")
    ap.add_argument("ref", help="Reference fasta")
    ap.add_argument("--context", type=int, default=100, 
                    help="Bases of context to extract around each kmer (default: 100)")
    args = ap.parse_args()
    
    # Load reference and kmers
    ref = Fasta(args.ref)
    kmers = Fasta(args.kmers_fa)
    
    print("="*80)
    print("KMER CONTEXT ANALYSIS")
    print("="*80)
    
    # Extract each kmer and its genomic context
    contexts = {}
    
    for kmer_id in kmers.keys():
        kmer_seq = str(kmers[kmer_id][:].seq).upper()
        print(f"\n{'='*80}")
        print(f"Kmer: {kmer_id}")
        print(f"Sequence: {kmer_seq}")
        print(f"Length: {len(kmer_seq)}")
        
        # Search in reference
        found = False
        for chrom in ref.keys():
            ref_seq = str(ref[chrom][:].seq).upper()
            
            # Search forward strand
            pos = ref_seq.find(kmer_seq)
            if pos != -1:
                found = True
                strand = "+"
                start = max(0, pos - args.context)
                end = min(len(ref_seq), pos + len(kmer_seq) + args.context)
                context = ref_seq[start:end]
                
                print(f"\nFound on {chrom} (forward strand)")
                print(f"Position: {pos+1} (1-based)")
                print(f"Context region: {start+1}-{end+1}")
                
                # Store context
                contexts[kmer_id] = {
                    'chrom': chrom,
                    'pos': pos,
                    'strand': strand,
                    'context': context,
                    'kmer_in_context': pos - start
                }
                
                # Print context with kmer highlighted
                kmer_start_in_context = pos - start
                kmer_end_in_context = kmer_start_in_context + len(kmer_seq)
                
                print(f"\nGenomic Context ({args.context}bp flanking):")
                print(f"5' flank: {context[:kmer_start_in_context]}")
                print(f"KMER:     {context[kmer_start_in_context:kmer_end_in_context]}")
                print(f"3' flank: {context[kmer_end_in_context:]}")
                
                break
            
            # Search reverse strand
            kmer_rc = reverse_complement(kmer_seq)
            pos = ref_seq.find(kmer_rc)
            if pos != -1:
                found = True
                strand = "-"
                start = max(0, pos - args.context)
                end = min(len(ref_seq), pos + len(kmer_seq) + args.context)
                context = ref_seq[start:end]
                
                print(f"\nFound on {chrom} (reverse strand)")
                print(f"Position: {pos+1} (1-based, forward coordinates)")
                print(f"Context region: {start+1}-{end+1}")
                
                contexts[kmer_id] = {
                    'chrom': chrom,
                    'pos': pos,
                    'strand': strand,
                    'context': context,
                    'kmer_in_context': pos - start
                }
                
                # Print context
                kmer_start_in_context = pos - start
                kmer_end_in_context = kmer_start_in_context + len(kmer_seq)
                
                print(f"\nGenomic Context ({args.context}bp flanking):")
                print(f"5' flank: {context[:kmer_start_in_context]}")
                print(f"KMER(rc): {context[kmer_start_in_context:kmer_end_in_context]}")
                print(f"3' flank: {context[kmer_end_in_context:]}")
                
                break
        
        if not found:
            print(f"NOT FOUND in reference (possible structural variant)")
    
    # Compare contexts between kmers
    if len(contexts) >= 2:
        print(f"\n\n{'='*80}")
        print("CONTEXT COMPARISON")
        print("="*80)
        
        kmer_ids = list(contexts.keys())
        
        for i in range(len(kmer_ids)):
            for j in range(i+1, len(kmer_ids)):
                k1 = kmer_ids[i]
                k2 = kmer_ids[j]
                
                print(f"\n{'-'*80}")
                print(f"Comparing {k1} vs {k2}")
                print(f"{'-'*80}")
                
                c1 = contexts[k1]
                c2 = contexts[k2]
                
                # Check if same chromosome
                if c1['chrom'] != c2['chrom']:
                    print(f"Different chromosomes: {c1['chrom']} vs {c2['chrom']}")
                    continue
                
                # Calculate distance
                dist = abs(c1['pos'] - c2['pos'])
                print(f"Distance: {dist} bp")
                print(f"Strand: {c1['strand']} vs {c2['strand']}")
                
                # Check for overlap
                k1_end = c1['pos'] + len(str(kmers[k1][:].seq))
                k2_end = c2['pos'] + len(str(kmers[k2][:].seq))
                
                if c1['pos'] <= c2['pos'] < k1_end or c2['pos'] <= c1['pos'] < k2_end:
                    print("⚠️  OVERLAPPING KMERS (should be mutually exclusive!)")
                    
                    # Align the contexts to see differences
                    print("\nAligning contexts to find variants:")
                    alignments = pairwise2.align.globalxx(c1['context'], c2['context'])
                    if alignments:
                        print("\nTop alignment:")
                        print(format_alignment(*alignments[0], full_sequences=True))
                        
                        # Find mismatches
                        aln = alignments[0]
                        seq1_aln = aln.seqA
                        seq2_aln = aln.seqB
                        
                        mismatches = []
                        indels = []
                        
                        for pos, (b1, b2) in enumerate(zip(seq1_aln, seq2_aln)):
                            if b1 != b2:
                                if b1 == '-' or b2 == '-':
                                    indels.append(pos)
                                else:
                                    mismatches.append((pos, b1, b2))
                        
                        print(f"\nSummary:")
                        print(f"  Mismatches: {len(mismatches)}")
                        print(f"  Indels: {len(indels)}")
                        
                        if mismatches:
                            print(f"\n  First 10 mismatches:")
                            for pos, b1, b2 in mismatches[:10]:
                                print(f"    Position {pos}: {b1} vs {b2}")
                        
                        if indels:
                            print(f"\n  Indel regions detected at positions: {indels[:10]}")
                
                else:
                    print("Non-overlapping kmers")
    
    print(f"\n{'='*80}")
    print("ANALYSIS COMPLETE")
    print("="*80)

if __name__ == "__main__":
    main()
