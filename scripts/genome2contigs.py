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


def main(reference_path, output_path, remove_intervals):
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

    for interval in remove_intervals:
        parts = interval.split("-")

        try:
            x, y = int(parts[0]), int(parts[1])
            length = y - x + 1

            if length > 1:
                keep_intervals.append((last, x - 1))
                last = y + 1
        except (IndexError, ValueError):
            print("Illegal remove interval format!")
            exit(1)

    keep_intervals.append((last, len(reference) - 1))
    print(keep_intervals)
    with open(output_path, "w") as output_file:
        contig_id = 0

        for start, end in keep_intervals:
            output_file.write(">%s, contig %d\n" % (name, contig_id))
            output_file.write(reference[start:end + 1])
            output_file.write("\n")

            contig_id += 1


if __name__ == "__main__":
    if len(sys.argv) < 4:
        print(__doc__)
        exit(1)

    main(sys.argv[1], sys.argv[2], sys.argv[3:])
