#!/bin/bash
# genomesaver.sh
# Usage: ./genomesaver.sh <input_file_with_genome_ids> <output_directory>

set -euo pipefail

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_file_with_genome_ids> <output_directory>"
    exit 1
fi

GENOME_LIST="$1"
OUTDIR="$2"

mkdir -p "$OUTDIR"
FAILED="$OUTDIR/failed_ids.log"
: > "$FAILED"   # reset log

# Track duplicates
declare -A SEEN

while IFS= read -r GENOME_ID; do
    # Skip empty or comment lines
    [[ -z "${GENOME_ID// }" || "$GENOME_ID" =~ ^[[:space:]]*# ]] && continue

    # Skip duplicates in the list
    if [[ -n "${SEEN[$GENOME_ID]:-}" ]]; then
        echo "[SKIP] $GENOME_ID (duplicate in list)"
        continue
    fi
    SEEN[$GENOME_ID]=1

    # OPTIONAL: only process N. gonorrhoeae (taxon 485). Uncomment to enforce.
    # [[ "$GENOME_ID" =~ ^485[_\.] ]] || continue

    # Resume-safe: skip if final file already exists and is non-empty
    if [ -s "$OUTDIR/$GENOME_ID.fna" ]; then
        echo "[SKIP] $GENOME_ID (exists)"
        continue
    fi

    echo "[GET ] $GENOME_ID"

    # BV-BRC expects dots; keep underscores for local filenames
    CANON_ID="${GENOME_ID//_/.}"

    URL_HTTPS="https://www.bvbrc.org/api/genome_sequence/?eq(genome_id,${CANON_ID})&http_accept=text/x-fasta"
    URL_FTP="ftp://ftp.bvbrc.org/genomes/${CANON_ID}/${CANON_ID}.fna"

    TMP="$OUTDIR/$GENOME_ID.fna.part"
    rm -f "$TMP"

    # 1) Try HTTPS first (best through firewalls/HPC)
    curl -fsSL --retry 3 --retry-delay 2 "$URL_HTTPS" -o "$TMP" >/dev/null 2>&1 || true

    # 2) If empty, fall back to FTP (passive mode is default)
    if [ ! -s "$TMP" ]; then
        wget -q --tries=3 --read-timeout=30 --timeout=30 -O "$TMP" "$URL_FTP" || true
    fi

    # Finalize or record failure
    if [ -s "$TMP" ]; then
        mv -f "$TMP" "$OUTDIR/$GENOME_ID.fna"
        echo "[OK  ] $GENOME_ID"
    else
        rm -f "$TMP" "$OUTDIR/$GENOME_ID.fna"
        echo "[FAIL] $GENOME_ID (empty after download)" | tee -a "$FAILED"
    fi
done < "$GENOME_LIST"

echo "Download complete. Files saved in $OUTDIR/"
if [ -s "$FAILED" ]; then
    echo "Some downloads failed. See: $FAILED"
else
    echo "No failures recorded."
fi
