#!/usr/bin/env python3
# Coffea arabica parental origin classification from BLAST hits　
import argparse
import csv
import math

import pandas as pd


BLAST_COLS = [
    "qseqid", "sseqid", "pident", "length", "evalue", "bitscore",
    "qlen", "slen", "qstart", "qend", "sstart", "send",
]


def load_targets(path):
    with open(path, "r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    return rows


def best_row(df, parent):
    parent_df = df[df["parent"] == parent]
    if parent_df.empty:
        return None
    return parent_df.sort_values(
        ["bitscore", "pident", "cov_q"], ascending=[False, False, False]
    ).iloc[0]


def value(row, column):
    if row is None:
        return ""
    val = row[column]
    if isinstance(val, float) and math.isnan(val):
        return ""
    return val


def classify(canephora, eugenioides, min_delta_bitscore, min_ratio):
    if canephora is None and eugenioides is None:
        return "no_hit"
    if canephora is None:
        return "eugenioides-like"
    if eugenioides is None:
        return "canephora-like"

    c_score = float(canephora["bitscore"])
    e_score = float(eugenioides["bitscore"])
    delta = c_score - e_score

    if delta >= min_delta_bitscore and c_score / e_score >= min_ratio:
        return "canephora-like"
    if -delta >= min_delta_bitscore and e_score / c_score >= min_ratio:
        return "eugenioides-like"
    return "ambiguous"


def main():
    parser = argparse.ArgumentParser(
        description="Assign C. arabica transcripts to the closest parental species from BLAST hits."
    )
    parser.add_argument("--targets", required=True, help="Target CSV with transcript_id and annotations")
    parser.add_argument("--hits", required=True, help="BLAST outfmt 6 with 12 columns")
    parser.add_argument("--out", required=True, help="Output TSV")
    parser.add_argument("--min-cov", type=float, default=0.7, help="Minimum query coverage")
    parser.add_argument("--min-pident", type=float, default=85.0, help="Minimum percent identity")
    parser.add_argument("--min-delta-bitscore", type=float, default=20.0)
    parser.add_argument("--min-ratio", type=float, default=1.02, help="Minimum best-score ratio")
    args = parser.parse_args()

    targets = load_targets(args.targets)
    hits = pd.read_csv(args.hits, sep="\t", header=None, names=BLAST_COLS)

    numeric_cols = ["pident", "length", "evalue", "bitscore", "qlen", "slen"]
    for col in numeric_cols:
        hits[col] = pd.to_numeric(hits[col], errors="coerce")

    hits["cov_q"] = hits["length"] / hits["qlen"]
    hits["parent"] = hits["sseqid"].str.split("|", regex=False).str[0]
    hits["parent_hit_id"] = hits["sseqid"].str.split("|", regex=False).str[1]
    hits = hits[(hits["cov_q"] >= args.min_cov) & (hits["pident"] >= args.min_pident)].copy()

    rows = []
    for target in targets:
        query_id = target["transcript_id"].strip()
        query_hits = hits[hits["qseqid"] == query_id]
        c_best = best_row(query_hits, "canephora")
        e_best = best_row(query_hits, "eugenioides")
        call = classify(c_best, e_best, args.min_delta_bitscore, args.min_ratio)

        c_score = value(c_best, "bitscore")
        e_score = value(e_best, "bitscore")
        delta = ""
        if c_score != "" and e_score != "":
            delta = float(c_score) - float(e_score)

        rows.append({
            "transcript_id": query_id,
            "gene_synonym": target.get("gene_synonym", ""),
            "description": target.get("description", ""),
            "xref_id": target.get("xref_id", ""),
            "call": call,
            "delta_bitscore_canephora_minus_eugenioides": delta,
            "canephora_hit": value(c_best, "parent_hit_id"),
            "canephora_bitscore": c_score,
            "canephora_pident": value(c_best, "pident"),
            "canephora_cov_q": value(c_best, "cov_q"),
            "canephora_evalue": value(c_best, "evalue"),
            "eugenioides_hit": value(e_best, "parent_hit_id"),
            "eugenioides_bitscore": e_score,
            "eugenioides_pident": value(e_best, "pident"),
            "eugenioides_cov_q": value(e_best, "cov_q"),
            "eugenioides_evalue": value(e_best, "evalue"),
        })

    pd.DataFrame(rows).to_csv(args.out, sep="\t", index=False)


if __name__ == "__main__":
    main()
