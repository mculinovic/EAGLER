"""
Various utility functions used with bioinformatic data.
"""


def complement_base(base):
    if base == "A":
        return "T"
    elif base == "T":
        return "A"
    elif base == "C":
        return "G"
    elif base == "G":
        return "C"
    else:
        raise ValueError("Illegal base")


def cyclic_print(string, output_file, limit):
    start = 0

    while start < len(string):
        output_file.write(string[start:start + limit])
        output_file.write("\n")
        start = start + limit
