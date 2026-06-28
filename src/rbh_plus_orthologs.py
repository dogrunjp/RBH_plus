#!/usr/bin/env python3
# rbh_inparalog.py
# Author: ChatGPT
# Description:
#   Compute Reciprocal Best Hits (RBH) between species pairs, then expand
#   co-orthologs by adding in-paralogs using a score ratio threshold (tau)
#   and coverage/length sanity checks. Emits orthogroups per pair and
#   mapping tables to a target species (e.g., Arabidopsis/Saccharomyces -> Nannochloropsis).

import argparse
import os
import sys
from collections import defaultdict
import pandas as pd
import numpy as np

COLS = [
    "qseqid","sseqid","pident","length","evalue","bitscore",
    "qlen","slen","qstart","qend","sstart","send"
]

def read_id_list(path):
    ids = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            t = line.strip()
            if t:
                ids.append(t)
    return set(ids)

def load_hits_table(path, min_cols=len(COLS)):
    # auto-detect columns count (but expect outfmt 6 with these 12 fields)
    df = pd.read_csv(path, sep="\t", header=None, comment="#", low_memory=False)
    if df.shape[1] < min_cols:
        raise ValueError(f"{path}: expected >= {min_cols} columns (outfmt 6 with qlen/slen/coords). Got {df.shape[1]}")
    df = df.iloc[:, :len(COLS)]
    df.columns = COLS
    # types
    df["length"]   = pd.to_numeric(df["length"], errors="coerce")
    df["evalue"]   = pd.to_numeric(df["evalue"], errors="coerce")
    df["bitscore"] = pd.to_numeric(df["bitscore"], errors="coerce")
    df["qlen"]     = pd.to_numeric(df["qlen"], errors="coerce")
    df["slen"]     = pd.to_numeric(df["slen"], errors="coerce")
    # coverage & sanity
    df["cov_q"] = df["length"] / df["qlen"].replace(0, np.nan)
    df["cov_s"] = df["length"] / df["slen"].replace(0, np.nan)
    df["cov_both"] = df[["cov_q","cov_s"]].min(axis=1)
    # length ratio sanity (min(qlen,slen)/max(qlen,slen))
    lr = (df[["qlen","slen"]].min(axis=1) / df[["qlen","slen"]].max(axis=1))
    df["len_ratio"] = lr
    return df

def subset_cross_species(df, ids_q, ids_s):
    df2 = df[df["qseqid"].isin(ids_q) & df["sseqid"].isin(ids_s)].copy()
    return df2

def best_hits(df_cross, metric="bitscore"):
    # For each query, pick subject with maximum metric
    # returns dict: q -> (s, score)
    idx = df_cross.groupby("qseqid")[metric].idxmax()
    best_df = df_cross.loc[idx, ["qseqid","sseqid",metric]].copy()
    best_df = best_df.rename(columns={metric:"score"})
    best = dict(zip(best_df["qseqid"], zip(best_df["sseqid"], best_df["score"])))
    return best

def rbh_pairs(bestAB, bestBA):
    anchors = []
    for a,(b,score_ab) in bestAB.items():
        if b in bestBA:
            a2,score_ba = bestBA[b]
            if a2 == a:
                anchors.append((a,b,float(score_ab)))
    return anchors

def expand_inparalogs(df_all, anchor, tau, min_cov, len_ratio, min_align, evalue):
    """
    anchor: (a, b, score_ab)
    For A-side expansion: pick all a' such that hit(a', b) exists and passes thresholds
                          and score(a', b) >= tau * score(a, b)
    For B-side expansion: pick all b' such that hit(b', a) exists and passes thresholds
                          and score(b', a) >= tau * score(a, b)
    """
    a, b, score_ab = anchor
    thr = tau * score_ab

    # A' : q in A, s == b
    df_a = df_all[(df_all["sseqid"] == b)]
    # thresholds
    df_a = df_a[
        (df_a["bitscore"] >= thr) &
        (df_a["cov_both"] >= min_cov) &
        (df_a["len_ratio"] >= len_ratio) &
        (df_a["length"] >= min_align) &
        (df_a["evalue"] <= evalue)
    ]
    A_members = set(df_a["qseqid"].tolist())
    if a not in A_members:
        A_members.add(a)

    # B' : q in B, s == a
    df_b = df_all[(df_all["sseqid"] == a)]
    df_b = df_b[
        (df_b["bitscore"] >= thr) &
        (df_b["cov_both"] >= min_cov) &
        (df_b["len_ratio"] >= len_ratio) &
        (df_b["length"] >= min_align) &
        (df_b["evalue"] <= evalue)
    ]
    B_members = set(df_b["qseqid"].tolist())
    if b not in B_members:
        B_members.add(b)

    return sorted(A_members), sorted(B_members)

def compute_pairwise_groups(df, ids_A, ids_B, tau, min_cov, len_ratio, min_align, evalue, pair_tag):
    # cross-subsets
    AB = subset_cross_species(df, ids_A, ids_B)
    BA = subset_cross_species(df, ids_B, ids_A)

    # apply basic filters up front to stabilize best hits
    filt = lambda x: x[
        (x["evalue"] <= evalue) &
        (x["cov_both"] >= min_cov) &
        (x["len_ratio"] >= len_ratio) &
        (x["length"] >= min_align)
    ].copy()
    ABf = filt(AB)
    BAf = filt(BA)

    bestAB = best_hits(ABf, metric="bitscore")
    bestBA = best_hits(BAf, metric="bitscore")
    anchors = rbh_pairs(bestAB, bestBA)

    groups = []
    for i,(a,b,score_ab) in enumerate(anchors, start=1):
        A_members, B_members = expand_inparalogs(
            pd.concat([ABf, BAf], ignore_index=True),
            (a,b,score_ab),
            tau=tau, min_cov=min_cov, len_ratio=len_ratio,
            min_align=min_align, evalue=evalue
        )
        groups.append({
            "pair": pair_tag,
            "group_id": f"{pair_tag}.g{i}",
            "anchor_a": a,
            "anchor_b": b,
            "anchor_score": score_ab,
            "A_members": A_members,
            "B_members": B_members,
        })
    return groups

def write_groups(groups, outdir, pair_tag):
    path = os.path.join(outdir, f"{pair_tag}.rbh_inparalog.tsv")
    rows = []
    for g in groups:
        rows.append({
            "pair": g["pair"],
            "group_id": g["group_id"],
            "anchor_a": g["anchor_a"],
            "anchor_b": g["anchor_b"],
            "anchor_score": g["anchor_score"],
            "A_count": len(g["A_members"]),
            "B_count": len(g["B_members"]),
            "A_members": ",".join(g["A_members"]),
            "B_members": ",".join(g["B_members"]),
        })
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)
    return path

def build_mapping_to_target(groups, A_label, T_label):
    """
    Build mapping dict: each A gene -> sorted list of T genes
    Here, A is the 'A_members' side and T is the 'B_members' side,
    assuming groups were computed for pair A:T.
    """
    amap = defaultdict(set)
    for g in groups:
        # A_members map to B_members
        for a in g["A_members"]:
            for t in g["B_members"]:
                amap[a].add(t)
    # sort
    return {a: sorted(list(v)) for a,v in amap.items()}

def main():
    ap = argparse.ArgumentParser(description="RBH + in-paralog expansion (ortholog co-groups) from DIAMOND/BLAST tabular hits.")
    ap.add_argument("--hits", required=True, help="All-vs-all hits (outfmt 6) with extra fields: qlen slen qstart qend sstart send")
    ap.add_argument("--species", action="append", nargs=2, metavar=("LABEL","ID_LIST"),
                    help="Species label and 1-column file of protein IDs. Use labels like A, B, T etc. Example: --species A arabidopsis.ids")
    ap.add_argument("--pairs", nargs="+", required=True,
                    help="Pairs to compute as LabelX:LabelY (direction-agnostic). Example: A:B A:T B:T")
    ap.add_argument("--tau", type=float, default=0.7, help="In-paralog expansion threshold as fraction of anchor bitscore (default 0.7)")
    ap.add_argument("--min_cov", type=float, default=0.6, help="Min coverage (min(query,subject)) per hit (default 0.6)")
    ap.add_argument("--len_ratio", type=float, default=0.7, help="Min min(qlen,slen)/max(qlen,slen) (default 0.7)")
    ap.add_argument("--evalue", type=float, default=1e-5, help="Max e-value (default 1e-5)")
    ap.add_argument("--min_align", type=int, default=50, help="Min alignment length (aa) (default 50)")
    ap.add_argument("--out", required=True, help="Output directory")
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)

    # species maps
    if not args.species or len(args.species) < 2:
        print("ERROR: Provide at least two --species LABEL ID_LIST pairs.", file=sys.stderr)
        sys.exit(1)

    label_to_ids = {}
    for label, idfile in args.species:
        label_to_ids[label] = read_id_list(idfile)

    # load table
    df = load_hits_table(args.hits)

    # compute groups for requested pairs
    pair_groups = {}
    for pair_tag in args.pairs:
        try:
            x,y = pair_tag.split(":")
        except ValueError:
            print(f"ERROR: pair must be LabelX:LabelY, got {pair_tag}", file=sys.stderr)
            sys.exit(1)
        if x not in label_to_ids or y not in label_to_ids:
            print(f"ERROR: unknown species label in {pair_tag}. Known: {list(label_to_ids.keys())}", file=sys.stderr)
            sys.exit(1)
        groups = compute_pairwise_groups(
            df,
            label_to_ids[x], label_to_ids[y],
            tau=args.tau, min_cov=args.min_cov, len_ratio=args.len_ratio,
            min_align=args.min_align, evalue=args.evalue,
            pair_tag=f"{x}:{y}"
        )
        pair_groups[pair_tag] = groups
        path = write_groups(groups, args.out, f"{x}:{y}")
        print(f"[OK] Wrote groups: {path}")

    # Optional: if we have A:T and/or B:T pairs, emit mapping tables
    # Detect target label 'T' if provided as species label; else skip
    labels = set(label_to_ids.keys())
    if "T" in labels:
        for src in labels - {"T"}:
            tag = f"{src}:T"
            if tag in pair_groups:
                mapping = build_mapping_to_target(pair_groups[tag], src, "T")
                outp = os.path.join(args.out, f"map_{src}_to_T.tsv")
                with open(outp, "w", encoding="utf-8") as f:
                    f.write("src_id\ttarget_ids\n")
                    for a, tgts in sorted(mapping.items()):
                        f.write(f"{a}\t{','.join(tgts)}\n")
                print(f"[OK] Wrote mapping {src}->T: {outp}")

if __name__ == "__main__":
    main()
