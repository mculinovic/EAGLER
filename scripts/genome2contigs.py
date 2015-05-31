#!/usr/bin/env python3

"""
Cut a genome into multiple contigs. The script expects at least 3 arguments,
the path to a FASTA file with the reference genome, the path to the output
FASTA file and a list of intervals to be removed (or kept) from the reference.
If N intervals were given as parameters the script will write N+1 contigs to
the output file (or N contigs if --keep is used).
"""

from argparse import ArgumentParser

from shared.bio_utils import cyclic_print


def main(reference_path, output_path, cut_intervals, keep):
    with open(reference_path, "r") as reference_file:
        name = reference_file.readline()[1:-1]
        reference = []

        for line in reference_file:
            if line and line != "\n" and not line.startswith(">"):
                reference.append(line.strip())
            else:
                break

        reference = "".join(reference)

    name = name.split(maxsplit=1)[0]
    if not name.endswith("|"):
        name += "|"

    keep_intervals = []
    last = 0

    last_interval = cut_intervals[-1].split("-")
    if len(last_interval) == 2 and last_interval[1] == "END":
        cut_intervals[-1] = last_interval[0] + "-" + str(len(reference) - 1)

    for interval in cut_intervals:
        parts = interval.split("-")

        try:
            start, end = int(parts[0]), int(parts[1])
            cut_length = end - start + 1

            if cut_length > 1:
                if keep:
                    keep_intervals.append((start, end))
                else:
                    if start != 0:
                        keep_intervals.append((last, start - 1))
                    last = end + 1
            else:
                raise ValueError("Illegal interval length!")
        except (IndexError, ValueError):
            print("Illegal remove interval format!")
            exit(1)

    if not keep and len(reference) - last > 0:
        keep_intervals.append((last, len(reference) - 1))

    with open(output_path, "w") as output_file:
        contig_id = 0

        for start, end in keep_intervals:
            output_file.write(">%s%d|\n" % (name, contig_id))
            cyclic_print(reference[start:end + 1], output_file, 70)
            contig_id += 1


if __name__ == "__main__":
    parser = ArgumentParser(prog=__file__[:-3], description=__doc__)

    parser.add_argument(
        "-v", "--version",
        action="version",
        version="%(prog)s v0.3"
    )

    parser.add_argument(
        "-k", "--keep",
        action="store_true",
        default=False,
        help="keep the specified intervals of the genome instead of "
             "cutting them out"
    )

    parser.add_argument("reference", help="path to the reference genome")

    parser.add_argument("output", help="path to the output file with the "
                                       "extracted contigs")

    parser.add_argument("intervals", nargs="+", help="x-y formatted inclusive"
                        " intervals to be eliminated from the reference "
                        "sequence, y may be substituted with the keyword END "
                        "to represent the last base of a reference genome")

    args = parser.parse_args()

    main(args.reference, args.output, args.intervals, args.keep)
