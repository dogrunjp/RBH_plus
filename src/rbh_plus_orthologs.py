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


def parse_hit_line(line):
    cols = line.rstrip("\n").split("\t")
    if len(cols) < len(COLS):
        return None
    qseqid, sseqid = cols[0], cols[1]
    try:
        pident = float(cols[2])
        length = float(cols[3])
        evalue = float(cols[4])
        bitscore = float(cols[5])
        qlen = float(cols[6])
        slen = float(cols[7])
    except ValueError:
        return None
    cov_q = length / qlen if qlen else float("nan")
    cov_s = length / slen if slen else float("nan")
    cov_both = min(cov_q, cov_s) if (cov_q == cov_q and cov_s == cov_s) else float("nan")
    len_ratio = min(qlen, slen) / max(qlen, slen) if max(qlen, slen) > 0 else 0.0
    return {
        "qseqid": qseqid,
        "sseqid": sseqid,
        "pident": pident,
        "length": length,
        "evalue": evalue,
        "bitscore": bitscore,
        "qlen": qlen,
        "slen": slen,
        "cov_both": cov_both,
        "len_ratio": len_ratio,
    }


def pass_filters(hit, min_cov, len_ratio, min_align, evalue):
    return (
        hit["evalue"] <= evalue
        and hit["cov_both"] >= min_cov
        and hit["len_ratio"] >= len_ratio
        and hit["length"] >= min_align
    )


def best_hits(path, ids_q, ids_s, min_cov, len_ratio, min_align, evalue):
    best = {}
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            hit = parse_hit_line(line)
            if not hit:
                continue
            if hit["qseqid"] not in ids_q or hit["sseqid"] not in ids_s:
                continue
            if not pass_filters(hit, min_cov, len_ratio, min_align, evalue):
                continue
            prev = best.get(hit["qseqid"])
            if prev is None or hit["bitscore"] > prev[1] or (
                hit["bitscore"] == prev[1] and hit["evalue"] < prev[2]
            ):
                best[hit["qseqid"]] = (hit["sseqid"], hit["bitscore"], hit["evalue"])
    return best


def index_filtered_hits(path, ids_from, ids_to, min_cov, len_ratio, min_align, evalue):
    index = defaultdict(list)
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            hit = parse_hit_line(line)
            if not hit:
                continue
            if hit["qseqid"] not in ids_from or hit["sseqid"] not in ids_to:
                continue
            if not pass_filters(hit, min_cov, len_ratio, min_align, evalue):
                continue
            index[hit["sseqid"]].append((hit["qseqid"], hit["bitscore"]))
    return index


def rbh_pairs(bestAB, bestBA):
    anchors = []
    for a, (b, score_ab, _) in bestAB.items():
        if b in bestBA:
            a2, score_ba, _ = bestBA[b]
            if a2 == a:
                anchors.append((a, b, float(score_ab)))
    return anchors


def expand_inparalogs(index_ab, index_ba, anchor, tau):
    """
    anchor: (a, b, score_ab)
    For A-side expansion: pick all a' such that hit(a', b) exists and passes thresholds
                          and score(a', b) >= tau * score(a, b)
    For B-side expansion: pick all b' such that hit(b', a) exists and passes thresholds
                          and score(b', a) >= tau * score(a, b)
    """
    a, b, score_ab = anchor
    thr = tau * score_ab

    A_members = {a}
    for qid, score in index_ab.get(b, []):
        if score >= thr:
            A_members.add(qid)

    B_members = {b}
    for qid, score in index_ba.get(a, []):
        if score >= thr:
            B_members.add(qid)

    return sorted(A_members), sorted(B_members)


def compute_pairwise_groups(path, ids_A, ids_B, tau, min_cov, len_ratio, min_align, evalue, pair_tag):
    bestAB = best_hits(path, ids_A, ids_B, min_cov, len_ratio, min_align, evalue)
    bestBA = best_hits(path, ids_B, ids_A, min_cov, len_ratio, min_align, evalue)
    anchors = rbh_pairs(bestAB, bestBA)

    index_ab = index_filtered_hits(path, ids_A, ids_B, min_cov, len_ratio, min_align, evalue)
    index_ba = index_filtered_hits(path, ids_B, ids_A, min_cov, len_ratio, min_align, evalue)

    groups = []
    for i, (a, b, score_ab) in enumerate(anchors, start=1):
        A_members, B_members = expand_inparalogs(
            index_ab, index_ba, (a, b, score_ab), tau=tau
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
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("pair\tgroup_id\tanchor_a\tanchor_b\tanchor_score\tA_count\tB_count\tA_members\tB_members\n")
        for row in rows:
            handle.write(
                f"{row['pair']}\t{row['group_id']}\t{row['anchor_a']}\t{row['anchor_b']}\t{row['anchor_score']}\t"
                f"{row['A_count']}\t{row['B_count']}\t{row['A_members']}\t{row['B_members']}\n"
            )
    return path


def build_mapping_to_target(groups, A_label, T_label):
    """
    Build mapping dict: each A gene -> sorted list of T genes
    Here, A is the 'A_members' side and T is the 'B_members' side,
    assuming groups were computed for pair A:T.
    """
    amap = defaultdict(set)
    for g in groups:
        for a in g["A_members"]:
            for t in g["B_members"]:
                amap[a].add(t)
    return {a: sorted(list(v)) for a, v in amap.items()}


def main():
    ap = argparse.ArgumentParser(description="RBH + in-paralog expansion (ortholog co-groups) from DIAMOND/BLAST tabular hits.")
    ap.add_argument("--hits", required=True, help="All-vs-all hits (outfmt 6) with extra fields: qlen slen qstart qend sstart send")
    ap.add_argument("--species", action="append", nargs=2, metavar=("LABEL", "ID_LIST"),
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

    if not args.species or len(args.species) < 2:
        print("ERROR: Provide at least two --species LABEL ID_LIST pairs.", file=sys.stderr)
        sys.exit(1)

    label_to_ids = {}
    for label, idfile in args.species:
        label_to_ids[label] = read_id_list(idfile)

    pair_groups = {}
    for pair_tag in args.pairs:
        try:
            x, y = pair_tag.split(":")
        except ValueError:
            print(f"ERROR: pair must be LabelX:LabelY, got {pair_tag}", file=sys.stderr)
            sys.exit(1)
        if x not in label_to_ids or y not in label_to_ids:
            print(f"ERROR: unknown species label in {pair_tag}. Known: {list(label_to_ids.keys())}", file=sys.stderr)
            sys.exit(1)
        groups = compute_pairwise_groups(
            args.hits,
            label_to_ids[x], label_to_ids[y],
            tau=args.tau, min_cov=args.min_cov, len_ratio=args.len_ratio,
            min_align=args.min_align, evalue=args.evalue,
            pair_tag=f"{x}:{y}"
        )
        pair_groups[pair_tag] = groups
        path = write_groups(groups, args.out, f"{x}:{y}")
        print(f"[OK] Wrote groups: {path}")

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
