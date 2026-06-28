#!/usr/bin/env python3
# Extract FASTA records whose first header token is listed in a CSV column.
# パスウェイに関連する遺伝子のFASTA配列をCSVに基づいて抽出するスクリプト
import argparse
import csv
import sys


def fasta_records(path):
    header = None
    seq_lines = []
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, seq_lines
                header = line
                seq_lines = []
            else:
                seq_lines.append(line)
    if header is not None:
        yield header, seq_lines


def main():
    parser = argparse.ArgumentParser(
        description="Extract FASTA records whose first header token is listed in a CSV column."
    )
    parser.add_argument("--csv", required=True, help="CSV file with target IDs")
    parser.add_argument("--id-column", default="transcript_id", help="CSV column containing FASTA IDs")
    parser.add_argument("--fasta", required=True, help="Input FASTA")
    parser.add_argument("--out", required=True, help="Output FASTA")
    args = parser.parse_args()

    with open(args.csv, "r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if args.id_column not in reader.fieldnames:
            raise SystemExit(f"ID column not found: {args.id_column}")
        target_ids = [row[args.id_column].strip() for row in reader if row[args.id_column].strip()]

    wanted = set(target_ids)
    found = set()

    with open(args.out, "w", encoding="utf-8") as out:
        for header, seq_lines in fasta_records(args.fasta):
            record_id = header[1:].split()[0]
            if record_id in wanted:
                found.add(record_id)
                out.write(header + "\n")
                out.write("\n".join(seq_lines) + "\n")

    missing = [record_id for record_id in target_ids if record_id not in found]
    if missing:
        print("Missing FASTA records:", ", ".join(missing), file=sys.stderr)
        return 1

    print(f"Extracted {len(found)} records to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
