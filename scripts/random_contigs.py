#!/usr/bin/env python3

import random
import os

from argparse import ArgumentParser

from shared.bio_structs import SequenceCollection

def main(reference_path, output_path, num_gaps, len_interval):
    reference = SequenceCollection.load_from_fasta(reference_path)[0]

    lengths = len_interval.split("-")
    min_size, max_size = int(lengths[0]), int(lengths[1])

    THRESHOLD = 100000

    if min_size >= max_size:
        print("Illegal length interval")
        exit(1)

    if max_size >= len(reference):
        print("Max gap length greater than reference length")
        exit(1)

    diff = int (len(reference) / (int(num_gaps) + 1))
    #print(diff)

    intervals = []
    for i in range(int(num_gaps)):
        start = random.randint(i * diff, (i + 1) * diff - THRESHOLD - max_size)
        end = random.randint(start + min_size, start + max_size)
        intervals.append((start, end))

    for i in range(len(intervals)):
        print("Creating gap at indices: " + str(intervals[i]))

    command = "python3 genome2contigs.py " + reference_path + " " + output_path
    for i in range(len(intervals)):
        command += " " + str(intervals[i][0]) + "-" + str(intervals[i][1])

    print("Executing script genome2contigs.py with command:")
    print(command)
    os.system(command)



if __name__ == "__main__":
    parser = ArgumentParser(prog=__file__[:-3], description=__doc__)

    parser.add_argument(
        "-v", "--version",
        action="version",
        version="%(prog)s v0.1"
    )

    parser.add_argument("reference", help="path to the reference genome")

    parser.add_argument("output", help="path to the output file with the "
                                   "extracted contigs")

    parser.add_argument("num_gaps", help="number of randomly created gaps in genome")

    parser.add_argument("len_interval", help="x-y formatted interval where x is minimum"
                        " gap length and y is maximum gap length. Each randomly created gap"
                        " has length within this interval")

    args = parser.parse_args()

    main(args.reference, args.output, args.num_gaps, args.len_interval)
