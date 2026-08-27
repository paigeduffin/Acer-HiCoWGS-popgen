#!/usr/bin/env python3
"""Filter FASTQ reads using the quality criteria applied to Panama samples.

Reads are retained when at least 90% of bases have a Phred quality score
of 20 or greater. This reproduces the behavior of:

    fastq_quality_filter -q 20 -p 90

Input FASTQ files must be uncompressed. The script writes an uncompressed
FASTQ file containing only reads that pass the filter.

Usage:
    python3 custom_qual.filt_90perc.over.20.py INPUT.fastq OUTPUT.clean
"""

import argparse


def phred_quality_filter(
    input_file: str,
    output_file: str,
    quality_threshold: int = 20,
    min_percent_high_quality: float = 0.90,
) -> None:
    """Retain reads meeting the specified Phred-quality threshold."""
    with open(input_file, "r", encoding="utf-8") as infile, open(
        output_file, "w", encoding="utf-8"
    ) as outfile:
        record_number = 0

        while True:
            header = infile.readline()
            if not header:
                break

            sequence = infile.readline()
            plus = infile.readline()
            quality = infile.readline()
            record_number += 1

            if not sequence or not plus or not quality:
                raise ValueError(
                    f"Incomplete FASTQ record {record_number} in {input_file}"
                )

            header = header.rstrip("\r\n")
            sequence = sequence.rstrip("\r\n")
            plus = plus.rstrip("\r\n")
            quality = quality.rstrip("\r\n")

            if len(sequence) != len(quality):
                raise ValueError(
                    f"Sequence and quality lengths differ in FASTQ record "
                    f"{record_number} in {input_file}"
                )

            high_quality_bases = sum(
                (ord(character) - 33) >= quality_threshold
                for character in quality
            )
            high_quality_fraction = high_quality_bases / len(quality)

            if high_quality_fraction >= min_percent_high_quality:
                outfile.write(f"{header}\n{sequence}\n{plus}\n{quality}\n")


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Retain FASTQ reads for which at least 90% of bases have "
            "Phred quality scores of 20 or greater."
        )
    )
    parser.add_argument("input_file", help="Uncompressed input FASTQ file")
    parser.add_argument("output_file", help="Filtered output FASTQ file")
    args = parser.parse_args()

    phred_quality_filter(args.input_file, args.output_file)


if __name__ == "__main__":
    main()
