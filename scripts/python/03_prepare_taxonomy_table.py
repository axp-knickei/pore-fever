#!/usr/bin/env python3

from pathlib import Path
import re
import sys

import pandas as pd


PROJECT = Path.cwd()
SAMPLES_FILE = PROJECT / "config" / "samples.tsv"
SYLPH_FILE = PROJECT / "results" / "taxonomy" / "sylph_species_profile.tsv"
OUTDIR = PROJECT / "results" / "taxonomy"


def require_file(path: Path) -> None:
    if not path.exists() or path.stat().st_size == 0:
        sys.exit(f"ERROR: Missing or empty file: {path}")


def infer_sample_id(value: str, valid_ids: set[str]) -> str | None:
    text = str(value)
    stem = Path(text).name
    stem = re.sub(r"\.(fastq|fq)(\.gz)?$", "", stem)
    stem = stem.replace(".host_removed", "")
    if stem in valid_ids:
        return stem
    for sample_id in valid_ids:
        if sample_id in text:
            return sample_id
    return None


def find_first_column(columns: list[str], candidates: list[str]) -> str | None:
    lowered = {col.lower(): col for col in columns}
    for candidate in candidates:
        if candidate.lower() in lowered:
            return lowered[candidate.lower()]
    return None


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    require_file(SAMPLES_FILE)
    require_file(SYLPH_FILE)

    metadata = pd.read_csv(SAMPLES_FILE, sep="\t")
    required_cols = {"sample_id", "treatment", "timepoint", "replicate"}
    missing = required_cols - set(metadata.columns)
    if missing:
        sys.exit(f"ERROR: samples.tsv missing columns: {sorted(missing)}")

    valid_ids = set(metadata["sample_id"].astype(str))
    tax = pd.read_csv(SYLPH_FILE, sep="\t", comment="#")

    sample_col = find_first_column(
        list(tax.columns),
        ["Sample_file", "sample_file", "sample", "Sample", "query", "Query"],
    )
    taxon_col = find_first_column(
        list(tax.columns),
        ["Taxonomy", "taxonomic_abundance", "Taxonomic_abundance", "Genome_name", "Name"],
    )
    abundance_col = find_first_column(
        list(tax.columns),
        ["Sequence_abundance", "sequence_abundance", "abundance", "Abundance", "Adjusted_ANI_5-95_percentile"],
    )

    if not all([sample_col, taxon_col, abundance_col]):
        note = OUTDIR / "taxonomy_table_preparation_NOTE.txt"
        note.write_text(
            "Sylph output was loaded, but the script could not infer sample, taxon, "
            "and abundance columns automatically.\n\n"
            f"Observed columns:\n{list(tax.columns)}\n"
        )
        sys.exit(f"ERROR: Could not infer sylph columns. See {note}")

    tidy = tax[[sample_col, taxon_col, abundance_col]].copy()
    tidy.columns = ["sample_raw", "taxon", "abundance"]
    tidy["sample_id"] = tidy["sample_raw"].map(lambda value: infer_sample_id(value, valid_ids))
    tidy["abundance"] = pd.to_numeric(tidy["abundance"], errors="coerce")
    tidy = tidy.dropna(subset=["sample_id", "taxon", "abundance"])

    if tidy.empty:
        sys.exit("ERROR: No usable taxonomy rows remained after sample/taxon/abundance parsing.")

    matrix = (
        tidy.pivot_table(index="taxon", columns="sample_id", values="abundance", aggfunc="sum", fill_value=0)
        .reset_index()
        .rename_axis(None, axis=1)
    )

    matrix_file = OUTDIR / "taxa_abundance_matrix.tsv"
    tidy_file = OUTDIR / "taxa_abundance_long.tsv"
    metadata_file = OUTDIR / "taxa_abundance_long_with_metadata.tsv"

    matrix.to_csv(matrix_file, sep="\t", index=False)
    tidy[["sample_id", "taxon", "abundance"]].to_csv(tidy_file, sep="\t", index=False)
    tidy.merge(metadata, on="sample_id", how="left").to_csv(metadata_file, sep="\t", index=False)

    print(f"Wrote abundance matrix: {matrix_file}")
    print(f"Wrote long table: {tidy_file}")
    print(f"Wrote metadata-joined table: {metadata_file}")


if __name__ == "__main__":
    main()
