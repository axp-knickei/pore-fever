#!/usr/bin/env python3

from pathlib import Path
import re
import sys
from typing import Any

import pandas as pd


PROJECT = Path(__file__).resolve().parents[2]
SAMPLES_FILE = PROJECT / "config" / "samples.tsv"
SYLPH_FILE = PROJECT / "results" / "taxonomy" / "sylph_species_profile.tsv"
OUTDIR = PROJECT / "results" / "taxonomy"
REQUIRED_METADATA_COLUMNS = {"sample_id", "treatment", "timepoint", "replicate"}

SAMPLE_COLUMN_CANDIDATES = [
    "Sample_file",
    "sample_file",
    "sample",
    "Sample",
    "query",
    "Query",
]
TAXON_COLUMN_CANDIDATES = [
    "Taxonomy",
    "taxonomy",
    "Genome_name",
    "genome_name",
    "Genome_file",
    "genome_file",
    "Species",
    "species",
    "Name",
    "name",
]
ABUNDANCE_COLUMN_CANDIDATES = [
    "Taxonomic_abundance",
    "taxonomic_abundance",
    "Sequence_abundance",
    "sequence_abundance",
    "Relative_abundance",
    "relative_abundance",
    "Abundance",
    "abundance",
]


def require_file(path: Path) -> None:
    if not path.exists() or path.stat().st_size == 0:
        sys.exit(f"ERROR: Missing or empty file: {path}")


def normalize_sample_token(value: Any) -> str:
    text = str(value).strip()
    name = Path(text).name
    name = re.sub(r"\.(fastq|fq)(\.gz)?$", "", name, flags=re.IGNORECASE)
    name = re.sub(r"\.host_removed$", "", name)
    return name


def infer_sample_id(value: Any, valid_ids: set[str]) -> str | None:
    normalized = normalize_sample_token(value)
    if normalized in valid_ids:
        return normalized
    return None


def find_first_column(columns: list[str], candidates: list[str]) -> str | None:
    lowered = {col.lower(): col for col in columns}
    for candidate in candidates:
        if candidate.lower() in lowered:
            return lowered[candidate.lower()]
    return None


def validate_metadata(metadata: pd.DataFrame) -> None:
    missing = REQUIRED_METADATA_COLUMNS - set(metadata.columns)
    if missing:
        sys.exit(f"ERROR: samples.tsv missing columns: {sorted(missing)}")

    for column in REQUIRED_METADATA_COLUMNS:
        empty = metadata[column].astype(str).str.strip() == ""
        if empty.any():
            rows = ", ".join(str(idx + 2) for idx in metadata.index[empty].tolist())
            sys.exit(f"ERROR: samples.tsv has empty {column} values on line(s): {rows}")

    duplicated = metadata["sample_id"].duplicated(keep=False)
    if duplicated.any():
        sample_ids = sorted(metadata.loc[duplicated, "sample_id"].astype(str).unique())
        sys.exit(f"ERROR: Duplicate sample_id values in samples.tsv: {sample_ids}")


def write_diagnostic(
    path: Path,
    message: str,
    tax: pd.DataFrame,
    detected_columns: dict[str, str | None] | None = None,
    extra: str = "",
) -> None:
    lines = [message, "", "Observed columns:", repr(list(tax.columns)), ""]
    if detected_columns is not None:
        lines.extend(["Detected columns:", repr(detected_columns), ""])
    if extra:
        lines.extend([extra, ""])
    lines.extend(["First rows:", tax.head(10).to_string(index=False)])
    path.write_text("\n".join(lines) + "\n")


def detect_columns(tax: pd.DataFrame) -> dict[str, str | None]:
    return {
        "sample": find_first_column(list(tax.columns), SAMPLE_COLUMN_CANDIDATES),
        "taxon": find_first_column(list(tax.columns), TAXON_COLUMN_CANDIDATES),
        "abundance": find_first_column(list(tax.columns), ABUNDANCE_COLUMN_CANDIDATES),
    }


def summarize_sample_sums(tidy: pd.DataFrame) -> pd.DataFrame:
    return (
        tidy.groupby("sample_id", as_index=False)["abundance"]
        .sum()
        .rename(columns={"abundance": "total_abundance"})
        .sort_values("sample_id")
    )


def write_run_summary(
    path: Path,
    metadata: pd.DataFrame,
    tax: pd.DataFrame,
    tidy_raw_rows: int,
    tidy_rows: int,
    dropped_unmatched_sample: int,
    dropped_missing_taxon: int,
    dropped_invalid_abundance: int,
    duplicate_rows_aggregated: int,
    detected_columns: dict[str, str | None],
    sample_sums: pd.DataFrame,
    output_files: list[Path],
) -> None:
    matched_samples = sorted(tax.attrs.get("matched_samples", []))
    missing_samples = sorted(set(metadata["sample_id"].astype(str)) - set(matched_samples))

    lines = [
        "Taxonomy table preparation summary",
        f"Input metadata: {SAMPLES_FILE}",
        f"Input Sylph profile: {SYLPH_FILE}",
        "",
        "Detected columns:",
        f"  sample: {detected_columns['sample']}",
        f"  taxon: {detected_columns['taxon']}",
        f"  abundance: {detected_columns['abundance']}",
        "",
        "Row counts:",
        f"  metadata samples: {metadata['sample_id'].nunique()}",
        f"  raw sylph rows: {len(tax)}",
        f"  selected rows before drops: {tidy_raw_rows}",
        f"  final tidy rows: {tidy_rows}",
        f"  dropped unmatched sample rows: {dropped_unmatched_sample}",
        f"  dropped missing/empty taxon rows: {dropped_missing_taxon}",
        f"  dropped invalid abundance rows: {dropped_invalid_abundance}",
        f"  duplicate sample/taxon rows aggregated by summing: {duplicate_rows_aggregated}",
        "",
        "Samples:",
        f"  matched samples: {len(matched_samples)}",
        f"  metadata samples without taxa rows: {len(missing_samples)}",
    ]

    if missing_samples:
        lines.append(f"  missing sample IDs: {', '.join(missing_samples)}")

    lines.extend(
        [
            "",
            "Per-sample abundance totals:",
            sample_sums.to_string(index=False),
            "",
            "Output files:",
        ]
    )
    lines.extend(f"  {path}" for path in output_files)
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    require_file(SAMPLES_FILE)
    require_file(SYLPH_FILE)

    metadata = pd.read_csv(SAMPLES_FILE, sep="\t", comment="#", dtype=str, keep_default_na=False)
    validate_metadata(metadata)
    valid_ids = set(metadata["sample_id"].astype(str))

    tax = pd.read_csv(SYLPH_FILE, sep="\t", comment="#")
    if tax.empty:
        sys.exit(f"ERROR: Sylph profile is empty after parsing: {SYLPH_FILE}")

    detected_columns = detect_columns(tax)
    if not all(detected_columns.values()):
        note = OUTDIR / "taxonomy_table_preparation_DIAGNOSTIC.txt"
        write_diagnostic(
            note,
            "Sylph output was loaded, but the script could not infer sample, taxon, and abundance columns automatically.",
            tax,
            detected_columns,
            extra="Edit the candidate column lists in scripts/python/03_prepare_taxonomy_table.py after confirming the Sylph output format.",
        )
        sys.exit(f"ERROR: Could not infer Sylph columns. See {note}")

    sample_col = detected_columns["sample"]
    taxon_col = detected_columns["taxon"]
    abundance_col = detected_columns["abundance"]
    assert sample_col is not None and taxon_col is not None and abundance_col is not None

    tidy = tax[[sample_col, taxon_col, abundance_col]].copy()
    tidy.columns = ["sample_raw", "taxon", "abundance"]
    tidy_raw_rows = len(tidy)

    tidy["sample_id"] = tidy["sample_raw"].map(lambda value: infer_sample_id(value, valid_ids))
    raw_taxon = tidy["taxon"]
    tidy["taxon"] = tidy["taxon"].fillna("").astype(str).str.strip()
    tidy["abundance"] = pd.to_numeric(tidy["abundance"], errors="coerce")

    dropped_unmatched_sample = int(tidy["sample_id"].isna().sum())
    dropped_missing_taxon = int(raw_taxon.isna().sum() + (tidy["taxon"] == "").sum())
    dropped_invalid_abundance = int(tidy["abundance"].isna().sum())

    unmatched_examples = tidy.loc[tidy["sample_id"].isna(), "sample_raw"].astype(str).drop_duplicates().head(20)
    if dropped_unmatched_sample:
        note = OUTDIR / "taxonomy_table_preparation_DIAGNOSTIC.txt"
        write_diagnostic(
            note,
            "Some Sylph rows could not be matched exactly to sample_id values from config/samples.tsv.",
            tax,
            detected_columns,
            extra="Unmatched sample field examples:\n" + unmatched_examples.to_string(index=False),
        )
        sys.exit(f"ERROR: {dropped_unmatched_sample} rows had unmatched samples. See {note}")

    tidy = tidy.dropna(subset=["sample_id", "abundance"])
    tidy = tidy[tidy["taxon"] != ""]

    if tidy.empty:
        sys.exit("ERROR: No usable taxonomy rows remained after sample/taxon/abundance parsing.")

    if (tidy["abundance"] < 0).any():
        negative_count = int((tidy["abundance"] < 0).sum())
        sys.exit(f"ERROR: Abundance column contains {negative_count} negative values.")

    tax.attrs["matched_samples"] = set(tidy["sample_id"].unique())
    duplicate_mask = tidy.duplicated(subset=["sample_id", "taxon"], keep=False)
    duplicate_rows_aggregated = int(duplicate_mask.sum())
    if duplicate_rows_aggregated:
        print(
            f"WARNING: Aggregating {duplicate_rows_aggregated} duplicate sample/taxon rows by summing abundances.",
            file=sys.stderr,
        )

    sample_sums = summarize_sample_sums(tidy)
    if sample_sums["total_abundance"].isna().any():
        sys.exit("ERROR: Per-sample abundance totals contain NA values.")
    if (sample_sums["total_abundance"] <= 0).any():
        print("WARNING: At least one sample has non-positive total abundance after parsing.", file=sys.stderr)
    if (sample_sums["total_abundance"] > 1000).any():
        print(
            "WARNING: Some per-sample abundance totals exceed 1000. Confirm the selected abundance column and units.",
            file=sys.stderr,
        )

    matrix = (
        tidy.pivot_table(index="taxon", columns="sample_id", values="abundance", aggfunc="sum", fill_value=0)
        .reset_index()
        .rename_axis(None, axis=1)
    )

    matrix_file = OUTDIR / "taxa_abundance_matrix.tsv"
    tidy_file = OUTDIR / "taxa_abundance_long.tsv"
    metadata_file = OUTDIR / "taxa_abundance_long_with_metadata.tsv"
    sample_sums_file = OUTDIR / "taxa_abundance_sample_sums.tsv"
    summary_file = OUTDIR / "taxonomy_table_preparation_summary.txt"

    tidy_out = tidy[["sample_id", "taxon", "abundance"]].sort_values(["sample_id", "taxon"])
    metadata_out = tidy_out.merge(metadata, on="sample_id", how="left")

    matrix.to_csv(matrix_file, sep="\t", index=False)
    tidy_out.to_csv(tidy_file, sep="\t", index=False)
    metadata_out.to_csv(metadata_file, sep="\t", index=False)
    sample_sums.to_csv(sample_sums_file, sep="\t", index=False)

    output_files = [matrix_file, tidy_file, metadata_file, sample_sums_file, summary_file]
    write_run_summary(
        summary_file,
        metadata,
        tax,
        tidy_raw_rows,
        len(tidy_out),
        dropped_unmatched_sample,
        dropped_missing_taxon,
        dropped_invalid_abundance,
        duplicate_rows_aggregated,
        detected_columns,
        sample_sums,
        output_files,
    )

    print(f"Wrote abundance matrix: {matrix_file}")
    print(f"Wrote long table: {tidy_file}")
    print(f"Wrote metadata-joined table: {metadata_file}")
    print(f"Wrote sample abundance sums: {sample_sums_file}")
    print(f"Wrote run summary: {summary_file}")


if __name__ == "__main__":
    main()
