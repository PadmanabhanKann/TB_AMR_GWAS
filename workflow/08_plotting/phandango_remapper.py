#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Reposition k-mer coordinates across multiple contigs using GFF data.

Usage:
    python reposition_kmers.py --txt saureus_kmers_GCF_000144955.plot \
                               --gff /path/to/GCF_000144955.gff \
                               --output repositioned_kmers.plot
"""

import pandas as pd
import argparse
import os

def get_contig_lengths(file_path):
    contigs = {}
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('##sequence-region'):
                parts = line.strip().split()
                contig_name = parts[1]
                length = int(parts[-1])
                contigs[contig_name] = length
    return contigs

def add_cumulative_length(txt_file_path, gff_file_path, output_file_path=None):
    contigs = get_contig_lengths(gff_file_path)

    # Calculate cumulative contig lengths
    cumulative_lengths = {}
    total_length = 0
    for contig in contigs:
        cumulative_lengths[contig] = total_length
        total_length += contigs[contig]

    # Read the txt file as DataFrame
    df = pd.read_csv(txt_file_path, sep='\t', header=None, skiprows=1, index_col=0)

    # Adjust coordinates
    def adjust_bp(row):
        contig = row.name
        if contig in cumulative_lengths:
            start, end = map(int, row[2].split('..'))
            start += cumulative_lengths[contig]
            end += cumulative_lengths[contig]
            return f"{start}..{end}"
        return row[2]

    df[2] = df.apply(adjust_bp, axis=1)

    if not output_file_path:
        base = os.path.basename(txt_file_path)
        output_file_path = f"repositioned_{base}"

    df.to_csv(output_file_path, sep='\t')
    print(f"Output written to: {output_file_path}")
    return df

def main():
    parser = argparse.ArgumentParser(description="Reposition k-mer plot coordinates based on contig lengths.")
    parser.add_argument('--txt', required=True, help='Path to the .plot file (k-mer locations)')
    parser.add_argument('--gff', required=True, help='Path to the GFF file with contig length info')
    parser.add_argument('--output', help='Optional output file path')

    args = parser.parse_args()

    add_cumulative_length(args.txt, args.gff, args.output)

if __name__ == '__main__':
    main()
