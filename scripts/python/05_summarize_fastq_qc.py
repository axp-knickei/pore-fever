#!/usr/bin/env python3

from pathlib import Path
import gzip
import sys

import pandas as pd


PROJECT = Path.cwd()
SAMPLES_FILE = PROJECT / "config" / "samples.tsv"
HOST_REMOVED_DIR = PROJECT / "data" / "host_removed_fastq"
OUTDIR = PROJECT / "results" / "qc"


def require_file(path: Path) -> None:
    if not path.exists() or path.stat().st_size == 0:
        sys.exit(f"ERROR: Missing or empty file: {path}")


def fastq_basic_stats(fastq_gz: Path, max_reads: int = 100000) -> dict[str, float | int]:
    """Estimate length statistics from up to max_reads reads while counting all reads."""
    lengths: list[int] = []
    read_count = 0

    with gzip.open(fastq_gz, "rt", errors="replace") as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            seq = handle.readline().strip()
            plus = handle.readline()
            qual = handle.readline()

            if not seq or not plus or not qual:
                sys.exit(f"ERROR: Truncated FASTQ record in {fastq_gz}")

            read_count += 1
            if read_count <= max_reads:
                lengths.append(len(seq))

    if lengths:
        mean_len = sum(lengths) / len(lengths)
        min_len = min(lengths)
        max_len = max(lengths)
    else:
        mean_len = min_len = max_len = 0

    return {
        "estimated_total_reads": read_count,
        "mean_length_first_reads": mean_len,
        "min_length_first_reads": min_len,
        "max_length_first_reads": max_len,
        "reads_used_for_length_estimate": len(lengths),
    }


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    require_file(SAMPLES_FILE)
    metadata = pd.read_csv(SAMPLES_FILE, sep="\t")

    records = []
    for _, row in metadata.iterrows():
        sample_id = row["sample_id"]
        fq = HOST_REMOVED_DIR / f"{sample_id}.host_removed.fastq.gz"
        require_file(fq)

        stats = fastq_basic_stats(fq)
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
            "estimated_total_reads",
            "mean_length_first_reads",
            "min_length_first_reads",
            "max_length_first_reads",
            "reads_used_for_length_estimate",
            "fastq",
        ]
    ]

    outfile = OUTDIR / "host_removed_fastq_qc_summary.tsv"
    qc.to_csv(outfile, sep="\t", index=False)
    print(f"QC summary written to: {outfile}")
    print(qc.head())


if __name__ == "__main__":
    main()
