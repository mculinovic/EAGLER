#!/usr/bin/env python3

"""
Generates the reverse complement of the given sequences. The scripts
expects at least 2 arguments. The first argument is a FASTA file containing the
input sequences, while the second argument is the path of the output FASTA
file. If only 2 arguments are given, all sequences in the input file will be
processed. Any additional argument will be treated as an index in the input
file and the sequence at that index will be reversed. All sequences whose index
was not listed will be forwarded to the output in the original format.
"""

from argparse import ArgumentParser

from shared.bio_structs import SequenceCollection


def main(input_path, output_path, sequence_index):
    sequences = SequenceCollection.load_from_fasta(input_path)
    sequence_index = set([int(index) for index in sequence_index])

    if len(sequence_index) == 0:
        sequence_index = set(range(len(sequences)))

    for index in range(len(sequences)):
        if index in sequence_index:
            sequences[index].reverse_complement()
            sequences[index].name += "|RC"

    sequences.dump_to_fasta(output_path)


if __name__ == "__main__":
    parser = ArgumentParser(prog=__file__[:-3], description=__doc__)

    parser.add_argument(
        "-v", "--version",
        action="version",
        version="%(prog)s v0.1"
    )

    parser.add_argument("input_path", help="path to the input FASTA file")

    parser.add_argument("output_path", help="path to the output FASTA file")

    parser.add_argument("sequence_index", nargs="*", help="indexes of the "
                        "sequences to be reversed by the script ")

    args = parser.parse_args()

    main(args.input_path, args.output_path, args.sequence_index)
