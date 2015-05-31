"""
Various utility functions used with bioinformatic data.
"""


def cyclic_print(string, output_file, limit):
    start = 0

    while start < len(string):
        output_file.write(string[start:start + limit])
        output_file.write("\n")
        start = start + limit
