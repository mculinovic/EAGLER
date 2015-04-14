// Copyright @mculinovic
#ifndef UTILITY_H
#define UTILITY_H

#include <seqan/seq_io.h>
#include <iostream>
#include <exception>
#include <cstdio>

using std::remove;

using seqan::StringSet;
using seqan::CharString;
using seqan::Dna5String;
using seqan::SeqFileIn;
using seqan::SeqFileOut;

using seqan::open;
using seqan::readRecords;
using seqan::writeRecord;

using std::exception;

namespace utility {

// reads sequences data from fasta file and
// stores it in two sets: sequences ids and
// sequences
void read_fasta(StringSet<CharString>* pids,
                StringSet<Dna5String>* pseqs,
                char *ont_reads_filename) {
    auto& ids = *pids;
    auto& seqs = *pseqs;

    // opening input file
    SeqFileIn input_file;
    if (!open(input_file, ont_reads_filename)) {
        std::cerr << "Error: Could not open file "
                  << ont_reads_filename << "\n";
        exit(1);
    }

    // read all reads in file
    try {
        readRecords(ids, seqs, input_file);
    } catch(exception const& e) {
        std::cout << "Error: " << e.what() << std::endl;
        exit(1);
    }
}

// writes given sequence id and sequence to file with filename
// @filename
void write_fasta(const CharString &id, const Dna5String &seq,
                 const char*& filename) {
    // opening output file
    SeqFileOut out_file;
    if (!open(out_file, filename)) {
        std::cerr << "Error: Could not open file "
                  << filename << "\n";
        exit(1);
    }

    try {
        writeRecord(out_file, id, seq);
    } catch(exception const& e) {
        std::cout << "Error: " << e.what() << std::endl;
        exit(1);
    }
}


void delete_file(const char*& filename) {
    if (remove(filename)) {
        std::cerr << "Error deleting file" << std :: endl;
        exit(1);
    } else {
        std::cout << "Temp file successfully deleted" << std::endl;
    }
}

};  // namespace utility

#endif  // UTILITY_H
