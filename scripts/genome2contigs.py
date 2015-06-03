#!/usr/bin/env python3

"""
Cut a genome into multiple contigs. The script expects at least 3 arguments,
the path to a FASTA file with the reference genome, the path to the output
FASTA file and a list of intervals to be removed (or kept) from the reference.
If N intervals were given as parameters the script will write N+1 contigs to
the output file (or N contigs if --keep is used).
"""

from argparse import ArgumentParser

from shared.bio_structs import SequenceCollection


def main(reference_path, output_path, cut_intervals, keep):
    reference = SequenceCollection.load_from_fasta(reference_path)[0]

    reference.name = reference.name.split(maxsplit=1)[0]
    if not reference.name.endswith("|"):
                reference.name += "|"

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

    contigs = SequenceCollection()

    for contig_id in range(len(keep_intervals)):
        start, end = keep_intervals[contig_id]
        contiga_name = "%s%d|" % (reference.name, contig_id)
        contigs.append(contiga_name, reference[start:end + 1])

    contigs.dump_to_fasta(output_path)


if __name__ == "__main__":
    parser = ArgumentParser(prog=__file__[:-3], description=__doc__)

    parser.add_argument(
        "-v", "--version",
        action="version",
        version="%(prog)s v0.4"
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
