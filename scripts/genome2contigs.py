#!/usr/bin/env python3

"""
Script to cut a genome into multiple contigs.

The script expects at least 3 arguments, the path to a FASTA file with
the reference genome, the path to the output FASTA file and a list of
intervals to be removed from the reference. If N intervals were given
as parameters the script will write N+1 contigs to the output file.

Usage:
    python3 genomes2contigs.py reference.fasta output.fasta [x1-y1]+

Args:
    reference_path: the path to the reference genome
    xa-yb: intervals to be eliminated from the reference sequence,
        both values are inclusive
"""

import sys


def _cyclic_print(string, output_file, limit):
    start = 0

    while start < len(string):
        output_file.write(string[start:start + limit])
        output_file.write("\n")
        start = start + limit


def main(reference_path, output_path, cut_intervals):
    with open(reference_path, "r") as reference_file:
        name = reference_file.readline()[1:-1]
        reference = []

        for line in reference_file:
            if line:
                reference.append(line.strip())
            else:
                break

        reference = "".join(reference)

        keep_intervals = []

    last = 0

    for interval in cut_intervals:
        parts = interval.split("-")

        try:
            start, end = int(parts[0]), int(parts[1])
            length = end - start + 1

            if length > 1:
                keep_intervals.append((last, start - 1))
                last = end + 1
            else:
                raise ValueError("Illegal interval length!")
        except (IndexError, ValueError):
            print("Illegal remove interval format!")
            exit(1)

    if len(reference) - last > 0:
        keep_intervals.append((last, len(reference) - 1))

    with open(output_path, "w") as output_file:
        contig_id = 0

        for start, end in keep_intervals:
            output_file.write(">%s, contig %d\n" % (name, contig_id))
            _cyclic_print(reference[start:end + 1], output_file, 70)
            contig_id += 1


if __name__ == "__main__":
    if len(sys.argv) < 4:
        print(__doc__)
        exit(1)

    main(sys.argv[1], sys.argv[2], sys.argv[3:])
