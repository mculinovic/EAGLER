"""
Data structures commonly used in bioinformatics.
"""

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
