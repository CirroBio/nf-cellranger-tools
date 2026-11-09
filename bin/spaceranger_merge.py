#!/usr/bin/env python3

import gzip
import logging
import sys
from pathlib import Path
from typing import Dict
import pandas

logger = logging.getLogger(__name__)


def read_custom() -> Dict[str, str]:
    """Read the sequences from custom/"""

    return {
        header: seq
        for file in Path("custom").rglob("*")
        for header, seq in parse_fasta(file)
    }


def format_fasta(recs: Dict[str, str]):
    return "\n".join([
        f">{header}\n{seq}"
        for header, seq in recs.items()
    ])


def format_gtf(recs: Dict[str, str]):
    return "\n".join([
        format_gtf_line(header, seq)
        for header, seq in recs.items()
    ])


def format_gtf_line(header, seq):
    return "\t".join([
        header,
        "unknown",
        "exon",
        str(1),
        str(len(seq)),
        ".",
        "+",
        ".",
        f'gene_id "{header}"; transcript_id "{header}"; gene_name "{header}"; gene_biotype "protein_coding"'
    ])


def safe_open(file: Path, mode: str):
    if file.name.endswith(".gz"):
        return gzip.open(file, mode)
    else:
        return open(file, mode)


def parse_fasta(file: Path):

    if not file.exists():
        logger.info(f"File does not exist: {file}")
        return

    logger.info(f"Reading FASTA records from {file}")

    header = None
    seq = []

    with safe_open(file, "r") as handle:
        for line in handle:
            if line[0] == ">":
                if header is not None and len(seq) > 0:
                    yield header, "".join(seq)
                    logger.info(f"Read {header} ({len(''.join(seq)):,}bp)")
                header = line[1:].split(" ")[0].rstrip("\n")
                seq = []
            else:
                if len(line) > 1:
                    seq.append(line.rstrip("\n"))

    if header is not None and len(seq) > 0:
        yield header, "".join(seq)    
        logger.info(f"Read {header} ({len(''.join(seq)):,}bp)")


def append(file_in, file_out, footer):
    with open(file_in, "r") as handle_in:
        with open(file_out, "w") as handle_out:
            for line in handle_in:
                handle_out.write(line)
            handle_out.write(footer)


def merge_probes(file_in, file_out, recs: Dict[str, str], genome_name: str):
    """Add each custom sequence to the probe set as the probe for its own gene.

    The rows follow the probe set's header, as Space Ranger rejects a row with fewer
    columns. Each custom gene is a single exon, so its probe does not span a splice
    junction ("unspliced"). Space Ranger also requires every probe to be as long as
    the others, and #reference_genome to name the genome mkref builds (CG000621).
    """

    with open(file_in, "r") as handle_in:
        lines = [line.rstrip("\n") for line in handle_in if line.strip()]

    rows = [line for line in lines if not line.startswith("#")]
    columns = rows[0].split(",")
    probe_length = len(rows[1].split(",")[columns.index("probe_seq")])

    wrong_length = [header for header, seq in recs.items() if len(seq) != probe_length]
    if wrong_length:
        raise ValueError(
            f"Custom sequences must be {probe_length} bp, the length of every probe "
            f"in {file_in}: {', '.join(wrong_length)}"
        )

    with open(file_out, "w") as handle_out:
        for line in lines:
            if line.startswith("#reference_genome="):
                line = f"#reference_genome={genome_name}"
            handle_out.write(line + "\n")
        for header, seq in recs.items():
            values = dict(
                gene_id=header,
                probe_seq=seq,
                probe_id=header,
                included="TRUE",
                region="unspliced"
            )
            handle_out.write(",".join(values[column] for column in columns) + "\n")


def main():
    custom = read_custom()

    append(
        "input.fasta",
        "merged.fasta",
        format_fasta(custom)
    )

    append(
        "input.gtf",
        "merged.gtf",
        format_gtf(custom)
    )

    merge_probes(
        "input.csv",
        "probes.csv",
        custom,
        sys.argv[1]
    )


if __name__ == "__main__":
    main()
