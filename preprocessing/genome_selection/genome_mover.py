import os
import shutil
import argparse

#to move the top genomes after quast analysis

# --- Argument parser ---
parser = argparse.ArgumentParser(description="Copy genome files listed in a text file to a target folder.")
parser.add_argument("-i", "--input", required=True, help="Path to genome list file (e.g., top_250_genomes.txt)")
parser.add_argument("-s", "--source", required=True, help="Path to source folder containing .fna files")
parser.add_argument("-o", "--output", required=True, help="Path to output folder where files will be copied")

args = parser.parse_args()
genome_list_file = args.input
source_dir = args.source
output_dir = args.output

# --- Ensure output directory exists ---
os.makedirs(output_dir, exist_ok=True)

# --- Read genome IDs ---
with open(genome_list_file, 'r') as f:
    genome_ids = [line.strip() for line in f if line.strip()]

# --- Copy files ---
missing = []
copied = 0

for genome_id in genome_ids:
    src = os.path.join(source_dir, f"{genome_id}.fna")
    dest = os.path.join(output_dir, f"{genome_id}.fna")
    if os.path.isfile(src):
        shutil.copy(src, dest)
        copied += 1
        print(f"Copied: {genome_id}.fna")
    else:
        print(f"Missing: {genome_id}.fna")
        missing.append(genome_id)

# --- Summary ---
print(f"\n✅ Done. {copied} genomes copied to {output_dir}.")
if missing:
    print(f"⚠️ {len(missing)} genomes were missing:")
    for m in missing:
        print(f" - {m}")
