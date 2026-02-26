#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  barcode_histogram.sh [barcodes.txt]

Description:
  Reads barcodes (one per line) and outputs a TSV histogram:
    barcode<TAB>count

  Uses an in-memory hash table (awk associative array).
  Much faster than sort|uniq for large inputs.

  If no file is provided (or "-" is given), reads from stdin.
EOF
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  usage
  exit 0
fi

INPUT="${1:-"-"}"

awk '
{
  counts[$0]++
}
END {
  for (barcode in counts)
    printf "%s\t%d\n", barcode, counts[barcode]
}
' "$INPUT"
