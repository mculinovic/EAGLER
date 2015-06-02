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


class AlignmentStats(object):
    def __init__(self, alignment, reference_sequence):
        self.alignment = alignment

        self.matches = 0
        self.mismatches = 0
        self.insertions = 0
        self.deletions = 0

        reference_position = alignment.reference_position
        read_position = 0

        for pair in alignment.cigar:
            if pair.symbol == "H":
                pass
            elif pair.symbol == "S":
                read_position += pair.count
            elif pair.symbol == "I":
                self.insertions += pair.count
                read_position += pair.count
            elif pair.symbol == "D":
                self.deletions += pair.count
                reference_position += pair.count
            elif pair.symbol == "M":
                for index in range(pair.count):
                    read_base = alignment.sequence[read_position]
                    reference_base = reference_sequence[reference_position]

                    if read_base == reference_base:
                        self.matches += 1
                    else:
                        self.mismatches += 1

                    read_position += 1
                    reference_position += 1

        self.indels = self.insertions + self.deletions
        self.length_error = abs(self.insertions - self.deletions)
        self.identity = 100.0 * self.matches / (self.matches + self.mismatches
                                                + self.deletions)

    def __str__(self):
        return "%40s\t%9.4f%%\t%10d\t%10d\t%10d\t%10d\t%10d\t%10d\t%10d" % (
            self.alignment.read_name,
            self.identity,
            self.matches,
            self.mismatches,
            self.insertions,
            self.deletions,
            self.indels,
            self.length_error,
            len(self.alignment.sequence)
        )


def main(reference_path, sam_files):
    reference = SequenceCollection.load_from_fasta(reference_path)[0]
    print("Reference genome: %s" % reference.name)

    for sam_file in sam_files:
        alignments = AlignmentRecord.load_from_sam(sam_file)
        print("\n%40s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s\t%10s" % (
            sam_file, "ID", "MTCH", "MISM", "I", "D", "I+D", "|I-D|", "LEN"))

        for alignment in alignments:
            stats = AlignmentStats(alignment, reference.sequence)
            print(str(stats))


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
