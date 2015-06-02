"""
Data structures commonly used in bioinformatics.
"""

import re

from shared.bio_utils import cyclic_print


class SequenceRead(object):
    def __init__(self, name, sequence):
        self.name = name
        self.sequence = sequence

    def __repr__(self):
        if len(self.sequence) <= 70:
            return ">%s\n%s" % (self.name, self.sequence)
        else:
            return ">%s\n%s..." % (self.name, self.sequence[:67])


class SequenceCollection(object):
    def __init__(self, entries):
        self.entries = entries

    def __len__(self):
        return len(self.entries)

    def __getitem__(self, index):
        return self.entries[index]

    def __iter__(self):
        return iter(self.entries)

    def __repr__(self):
        return "\n\n".join([repr(entry) for entry in self])

    def dump_to_fasta(self, path):
        with open(path, "w") as output_file:
            for index in range(len(self)):
                output_file.write(">")
                output_file.write(self[index].name)
                output_file.write("\n")
                cyclic_print(self[index].sequence, output_file, 70)
                output_file.write("\n" if index < len(self) - 1 else "")

    @staticmethod
    def load_from_fasta(path):
        entries = []
        name = None
        sequence = []

        with open(path, "r") as fasta_file:
            for line in fasta_file:
                if line == "\n":
                    continue
                elif line.startswith(">"):
                    if len(sequence) > 0:
                        entries.append(SequenceRead(name, "".join(sequence)))
                        sequence = []

                    name = line[1:-1].split(maxsplit=1)[0]
                else:
                    sequence.append(line.strip())

            if len(sequence) > 0:
                entries.append(SequenceRead(name, "".join(sequence)))

        return SequenceCollection(entries)


class CigarString(object):

    class CigarPair(object):
        def __init__(self, count, symbol):
            self.count = count
            self.symbol = symbol

        def __repr__(self):
            return "%d%s" % (self.count, self.symbol)

    def __init__(self, string):
        pattern = re.compile("([0-9]+)([MIDNSHPX=])")
        pairs = pattern.findall(string)
        self.pairs = [CigarString.CigarPair(int(c), s) for c, s in pairs]

    def __len__(self):
        return len(self.pairs)

    def __getitem__(self, index):
        return self.pairs[index]

    def __iter__(self):
        return iter(self.pairs)

    def __repr__(self):
        return "".join([repr(pair) for pair in self])


class AlignmentRecord(object):

    class Alignment(object):
        def __init__(self, string):
            parts = string.split("\t")

            self.read_name = parts[0]
            self.flag = int(parts[1])
            self.reference_name = parts[2]
            self.reference_position = int(parts[3]) - 1
            self.cigar = CigarString(parts[5])
            self.sequence = parts[9]

        def __repr__(self):
            buffer = ["%s: %s" % (k, v) for k, v in self.__dict__.items()]
            return "\n".join(buffer)

    def __init__(self, headers, alignments):
        self.headers = headers
        self.alignments = alignments

    def __len__(self):
        return len(self.alignments)

    def __getitem__(self, index):
        return self.alignments[index]

    def __iter__(self):
        return iter(self.alignments)

    def __repr__(self):
        buffer = []
        buffer.extend(self.headers)
        buffer.extend([repr(alignment) for alignment in self])
        return "\n".join(buffer)

    @staticmethod
    def load_from_sam(path):
        headers = []
        alignments = []

        with open(path, "r") as sam_file:
            for line in sam_file:
                if line == "\n":
                    continue
                elif line.startswith("@"):
                    headers.append(line.strip())
                else:
                    alignment = AlignmentRecord.Alignment(line.strip())
                    if alignment.flag & 0x4 == 0:
                        alignments.append(alignment)

        return AlignmentRecord(headers, alignments)
