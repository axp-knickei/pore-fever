#!/usr/bin/env python3

from pathlib import Path
import gzip
import re
import sys
from typing import Any

import pandas as pd


PROJECT = Path(__file__).resolve().parents[2]
SAMPLES_FILE = PROJECT / "config" / "samples.tsv"
HOST_REMOVED_DIR = PROJECT / "data" / "host_removed_fastq"
HOST_STATS_DIR = PROJECT / "results" / "qc" / "host_removal"
OUTDIR = PROJECT / "results" / "qc"
REQUIRED_METADATA_COLUMNS = {"sample_id", "treatment", "timepoint", "replicate"}


def require_file(path: Path) -> None:
    if not path.exists() or path.stat().st_size == 0:
        sys.exit(f"ERROR: Missing or empty file: {path}")


def fail_fastq(path: Path, read_number: int, message: str) -> None:
    sys.exit(f"ERROR: Malformed FASTQ in {path}, read {read_number}: {message}")


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


def fastq_basic_stats(fastq_gz: Path, max_reads: int = 100000) -> dict[str, float | int]:
    """Count all reads and estimate length statistics from up to max_reads reads."""
    lengths: list[int] = []
    total_reads = 0

    with gzip.open(fastq_gz, "rt", errors="replace") as handle:
        while True:
            header = handle.readline()
            if not header:
                break

            total_reads += 1
            seq = handle.readline()
            plus = handle.readline()
            qual = handle.readline()

            if not seq or not plus or not qual:
                fail_fastq(fastq_gz, total_reads, "truncated record")
            if not header.startswith("@"):
                fail_fastq(fastq_gz, total_reads, "header line does not start with '@'")
            if not plus.startswith("+"):
                fail_fastq(fastq_gz, total_reads, "plus line does not start with '+'")

            seq = seq.rstrip("\n\r")
            qual = qual.rstrip("\n\r")
            if len(seq) != len(qual):
                fail_fastq(
                    fastq_gz,
                    total_reads,
                    f"sequence length ({len(seq)}) does not match quality length ({len(qual)})",
                )

            if total_reads <= max_reads:
                lengths.append(len(seq))

    if lengths:
        mean_len = sum(lengths) / len(lengths)
        min_len = min(lengths)
        max_len = max(lengths)
    else:
        mean_len = min_len = max_len = 0

    return {
        "total_reads": total_reads,
        "mean_length_sampled_reads": mean_len,
        "min_length_sampled_reads": min_len,
        "max_length_sampled_reads": max_len,
        "reads_used_for_length_stats": len(lengths),
    }


def first_two_counts(line: str) -> tuple[int, int] | None:
    match = re.match(r"^(\d+) \+ (\d+) ", line)
    if not match:
        return None
    return int(match.group(1)), int(match.group(2))


def first_percent(line: str) -> float | None:
    match = re.search(r"\((\d+(?:\.\d+)?)%", line)
    if not match:
        return None
    return float(match.group(1))


def parse_flagstat(flagstat_file: Path) -> dict[str, Any]:
    """Parse key samtools flagstat values from the host-removal BAM."""
    stats: dict[str, Any] = {
        "host_flagstat_file": str(flagstat_file),
        "host_flagstat_found": False,
        "host_total_alignments": pd.NA,
        "host_primary_alignments": pd.NA,
        "host_secondary_alignments": pd.NA,
        "host_supplementary_alignments": pd.NA,
        "mouse_mapped_alignments": pd.NA,
        "mouse_mapped_percent": pd.NA,
        "mouse_primary_mapped_reads": pd.NA,
        "mouse_primary_mapped_percent": pd.NA,
        "host_primary_unmapped_reads": pd.NA,
    }

    if not flagstat_file.exists() or flagstat_file.stat().st_size == 0:
        print(f"WARNING: Missing host-removal flagstat file: {flagstat_file}", file=sys.stderr)
        return stats

    stats["host_flagstat_found"] = True
    for line in flagstat_file.read_text().splitlines():
        counts = first_two_counts(line)
        if counts is None:
            continue
        passed, failed = counts
        total = passed + failed

        if " in total " in line:
            stats["host_total_alignments"] = total
        elif line.endswith(" primary"):
            stats["host_primary_alignments"] = total
        elif line.endswith(" secondary"):
            stats["host_secondary_alignments"] = total
        elif line.endswith(" supplementary"):
            stats["host_supplementary_alignments"] = total
        elif " primary mapped (" in line:
            stats["mouse_primary_mapped_reads"] = total
            stats["mouse_primary_mapped_percent"] = first_percent(line)
        elif " mapped (" in line:
            stats["mouse_mapped_alignments"] = total
            stats["mouse_mapped_percent"] = first_percent(line)

    primary = stats["host_primary_alignments"]
    primary_mapped = stats["mouse_primary_mapped_reads"]
    if pd.notna(primary) and pd.notna(primary_mapped):
        stats["host_primary_unmapped_reads"] = int(primary) - int(primary_mapped)

    return stats


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    require_file(SAMPLES_FILE)
    metadata = pd.read_csv(SAMPLES_FILE, sep="\t", comment="#", dtype=str, keep_default_na=False)
    validate_metadata(metadata)

    records = []
    for _, row in metadata.iterrows():
        sample_id = row["sample_id"]
        fq = HOST_REMOVED_DIR / f"{sample_id}.host_removed.fastq.gz"
        flagstat_file = HOST_STATS_DIR / f"{sample_id}.host_removal_stats.txt"
        require_file(fq)

        stats: dict[str, Any] = fastq_basic_stats(fq)
        stats.update(parse_flagstat(flagstat_file))
        stats.update(
            {
                "sample_id": sample_id,
                "treatment": row["treatment"],
                "timepoint": row["timepoint"],
                "replicate": row["replicate"],
                "fastq": str(fq),
            }
        )
        records.append(stats)

    qc = pd.DataFrame(records)
    qc = qc[
        [
            "sample_id",
            "treatment",
            "timepoint",
            "replicate",
            "total_reads",
            "mean_length_sampled_reads",
            "min_length_sampled_reads",
            "max_length_sampled_reads",
            "reads_used_for_length_stats",
            "host_total_alignments",
            "host_primary_alignments",
            "host_secondary_alignments",
            "host_supplementary_alignments",
            "mouse_mapped_alignments",
            "mouse_mapped_percent",
            "mouse_primary_mapped_reads",
            "mouse_primary_mapped_percent",
            "host_primary_unmapped_reads",
            "host_flagstat_found",
            "fastq",
            "host_flagstat_file",
        ]
    ]

    outfile = OUTDIR / "host_removed_fastq_qc_summary.tsv"
    qc.to_csv(outfile, sep="\t", index=False)

    print(f"QC summary written to: {outfile}")
    print(qc.head())


if __name__ == "__main__":
    main()
