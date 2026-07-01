#!/usr/bin/env python3
"""Convert an RBH_PLUS mapping table into PathLift's provided-table format."""

from __future__ import annotations

import argparse
import csv
import re
import sys
from dataclasses import dataclass
from pathlib import Path


SYMBOL_PATTERNS = (
    re.compile(r"(?:^|\s)gene_symbol:([^\s]+)"),
    re.compile(r"\[gene=([^\]]+)\]"),
    re.compile(r"(?:^|\s)GN=([^\s]+)"),
)


class ConversionError(ValueError):
    """Raised when an input cannot be converted safely."""


@dataclass(frozen=True)
class ConversionStats:
    mapping_rows: int
    source_ids_with_symbol: int
    source_ids_without_symbol: int
    target_links: int
    output_rows: int
    target_ids_without_mapping: int


def _symbol_from_header(header: str) -> str | None:
    for pattern in SYMBOL_PATTERNS:
        match = pattern.search(header)
        if match:
            symbol = match.group(1).strip()
            if symbol:
                return symbol
    return None


def load_source_symbols(path: str | Path) -> dict[str, str]:
    """Load protein ID -> gene symbol from Ensembl/NCBI/UniProt FASTA headers."""
    symbols: dict[str, str] = {}
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            header = line[1:].strip()
            source_id = header.split(maxsplit=1)[0] if header else ""
            symbol = _symbol_from_header(header)
            if source_id and symbol:
                previous = symbols.get(source_id)
                if previous is not None and previous != symbol:
                    raise ConversionError(
                        f"source FASTA contains conflicting symbols for {source_id}: "
                        f"{previous!r}, {symbol!r}"
                    )
                symbols[source_id] = symbol
    return symbols


def _versionless(identifier: str) -> str:
    base, separator, version = identifier.rpartition(".")
    return base if separator and version.isdigit() else identifier


def load_target_id_map(path: str | Path) -> tuple[dict[str, str], dict[str, str]]:
    """Load a TSV with protein_id and gene_id columns.

    Exact accessions are preferred. Versionless accessions are accepted only when
    they map unambiguously to one gene.
    """
    exact: dict[str, str] = {}
    versionless_values: dict[str, set[str]] = {}
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ConversionError(f"target ID map is empty: {path}")
        missing = {"protein_id", "gene_id"} - set(reader.fieldnames)
        if missing:
            raise ConversionError(
                "target ID map requires columns protein_id and gene_id; "
                f"missing: {', '.join(sorted(missing))}"
            )
        for line_number, row in enumerate(reader, start=2):
            protein_id = (row.get("protein_id") or "").strip()
            gene_id = (row.get("gene_id") or "").strip()
            if not protein_id or not gene_id:
                continue
            previous = exact.get(protein_id)
            if previous is not None and previous != gene_id:
                raise ConversionError(
                    f"target ID map line {line_number}: {protein_id} maps to both "
                    f"{previous!r} and {gene_id!r}"
                )
            exact[protein_id] = gene_id
            versionless_values.setdefault(_versionless(protein_id), set()).add(gene_id)

    versionless = {
        protein_id: next(iter(gene_ids))
        for protein_id, gene_ids in versionless_values.items()
        if len(gene_ids) == 1
    }
    return exact, versionless


def _mapped_target(
    target_id: str,
    exact: dict[str, str] | None,
    versionless: dict[str, str] | None,
) -> str | None:
    if exact is None or versionless is None:
        return target_id
    return exact.get(target_id) or versionless.get(_versionless(target_id))


def convert_mapping(
    mapping_path: str | Path,
    source_fasta: str | Path,
    output_path: str | Path,
    *,
    target_id_map: str | Path | None = None,
    missing_target: str = "error",
) -> ConversionStats:
    """Convert map_*_to_T.tsv into target_id/source_symbol/source_pid rows."""
    if missing_target not in {"error", "skip", "keep"}:
        raise ValueError("missing_target must be error, skip, or keep")

    source_symbols = load_source_symbols(source_fasta)
    exact = versionless = None
    if target_id_map is not None:
        exact, versionless = load_target_id_map(target_id_map)

    output_rows: set[tuple[str, str, str]] = set()
    mapping_rows = 0
    with_symbol = 0
    without_symbol = 0
    target_links = 0
    missing_targets: set[str] = set()

    with open(mapping_path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ConversionError(f"RBH mapping is empty: {mapping_path}")
        missing_columns = {"src_id", "target_ids"} - set(reader.fieldnames)
        if missing_columns:
            raise ConversionError(
                "RBH mapping requires columns src_id and target_ids; "
                f"missing: {', '.join(sorted(missing_columns))}"
            )

        for row in reader:
            source_id = (row.get("src_id") or "").strip()
            targets = [
                value.strip()
                for value in (row.get("target_ids") or "").split(",")
                if value.strip()
            ]
            if not source_id:
                continue
            mapping_rows += 1
            symbol = source_symbols.get(source_id)
            if not symbol:
                without_symbol += 1
                continue
            with_symbol += 1
            for target_id in targets:
                target_links += 1
                converted_target = _mapped_target(target_id, exact, versionless)
                if converted_target is None:
                    missing_targets.add(target_id)
                    if missing_target == "skip":
                        continue
                    if missing_target == "keep":
                        converted_target = target_id
                    else:
                        continue
                output_rows.add((converted_target, symbol, source_id))

    if missing_targets and missing_target == "error":
        sample = ", ".join(sorted(missing_targets)[:5])
        raise ConversionError(
            f"{len(missing_targets)} target protein IDs are absent from the target ID map "
            f"(examples: {sample}). Use --missing-target skip or keep to override."
        )

    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    with open(output, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("target_id", "source_symbol", "source_pid"))
        writer.writerows(sorted(output_rows, key=lambda row: (row[1], row[0], row[2])))

    return ConversionStats(
        mapping_rows=mapping_rows,
        source_ids_with_symbol=with_symbol,
        source_ids_without_symbol=without_symbol,
        target_links=target_links,
        output_rows=len(output_rows),
        target_ids_without_mapping=len(missing_targets),
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Convert an RBH_PLUS mapping table for PathLift provided_table input."
    )
    parser.add_argument("--mapping", required=True, help="RBH_PLUS map_*_to_T.tsv")
    parser.add_argument(
        "--source-fasta",
        required=True,
        help="Source protein FASTA containing gene_symbol:, [gene=...], or GN= metadata",
    )
    parser.add_argument("--out", required=True, help="Output PathLift-compatible TSV")
    parser.add_argument(
        "--target-id-map",
        help="Optional TSV with protein_id and gene_id columns; collapses target proteins to genes",
    )
    parser.add_argument(
        "--missing-target",
        choices=("error", "skip", "keep"),
        default="error",
        help="Handling of proteins absent from --target-id-map (default: error)",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        stats = convert_mapping(
            args.mapping,
            args.source_fasta,
            args.out,
            target_id_map=args.target_id_map,
            missing_target=args.missing_target,
        )
    except (ConversionError, OSError) as error:
        print(f"Error: {error}", file=sys.stderr)
        return 2

    print(f"[OK] PathLift table: {args.out}")
    print(
        "  RBH rows: "
        f"{stats.mapping_rows}; with symbol: {stats.source_ids_with_symbol}; "
        f"without symbol: {stats.source_ids_without_symbol}"
    )
    print(f"  target links: {stats.target_links}; output rows: {stats.output_rows}")
    if args.target_id_map is None:
        print(
            "  Note: target IDs remain protein accessions. Supply --target-id-map "
            "for gene-level PathLift output.",
            file=sys.stderr,
        )
    elif stats.target_ids_without_mapping:
        print(
            f"  target proteins without gene mapping: {stats.target_ids_without_mapping}",
            file=sys.stderr,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
