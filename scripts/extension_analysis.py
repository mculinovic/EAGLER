#!/usr/bin/env python3

"""
Runs various statistics on extension created by ONTscaffolder. The scripts
expects at 2 arguments. The first argument is FASTA file containing the
reference genome of the analysed organism. All further arguments will be
treated as SAM files containing alignments of extensions created by the
scaffolder. A table-like format will be outputted to stdout with each SAM file
in a separate table and each alignment in a separate row.
"""

from argparse import ArgumentParser

from shared.bio_structs import SequenceCollection, AlignmentRecord


def _run_statistics(reference, alignments):
    for alignment in alignments:
        print("%40s\tID\tI\tD\tI+D\tLEN" % (alignment.read_name))


def main(reference_path, sam_files):
    reference = SequenceCollection.load_from_fasta(reference_path)[0]
    print("Reference genome: %s" % reference.name)

    for sam_file in sam_files:
        alignments = AlignmentRecord.load_from_sam(sam_file)
        print("\n%40s\tID\tI\tD\tI+D\tLEN" % sam_file)

        _run_statistics(reference, alignments)


if __name__ == "__main__":
    parser = ArgumentParser(prog=__file__[:-3], description=__doc__)

    parser.add_argument(
        "-v", "--version",
        action="version",
        version="%(prog)s v0.1"
    )

    parser.add_argument("reference", help="path to the reference genome")

    parser.add_argument("sam_file", nargs="+", help="one or more SAM files "
                        "containing the alignments of the extension of the "
                        "draft genome on the reference genome")

    args = parser.parse_args()

    main(args.reference, args.sam_file)
