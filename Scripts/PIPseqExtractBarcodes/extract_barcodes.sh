#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  extract_barcodes_from_fastq.sh [-p "0_7+11_16+20_25+31_38"] [fastq]

Description:
  Extracts barcode substrings from the sequence line (line 2 of each FASTQ record),
  using 0-based inclusive ranges like "start_end" joined by '+'.

Arguments:
  fastq      Input FASTQ file. If omitted or "-", reads from stdin.

Options:
  -p POS     Barcode positions string (default: 0_7+11_16+20_25+31_38)
  -h         Show this help.

Examples:
  ./extract_barcodes_from_fastq.sh reads.fastq > barcodes.txt
  zcat reads.fastq.gz | ./extract_barcodes_from_fastq.sh - > barcodes.txt
  ./extract_barcodes_from_fastq.sh -p "0_7+11_16" reads.fastq
EOF
}

POS="0_7+11_16+20_25+31_38"

while getopts ":p:h" opt; do
  case "$opt" in
    p) POS="$OPTARG" ;;
    h) usage; exit 0 ;;
    \?) echo "Error: Invalid option -$OPTARG" >&2; usage; exit 2 ;;
    :)  echo "Error: Option -$OPTARG requires an argument" >&2; usage; exit 2 ;;
  esac
done
shift $((OPTIND - 1))

INPUT="${1:-"-"}"

awk -v POS="$POS" '
BEGIN{
  n = split(POS, parts, /\+/)
  for (i=1; i<=n; i++) {
    split(parts[i], a, /_/)
    s[i] = a[1] + 1              # awk substr() is 1-based
    len[i] = a[2] - a[1] + 1
  }
}
NR % 4 == 2 {
  sub(/\r$/, "", $0)             # tolerate CRLF
  out = ""
  for (i=1; i<=n; i++) out = out substr($0, s[i], len[i])
  print out
}
' "$INPUT"
