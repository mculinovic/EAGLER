// Copyright @mculinovic
#include <seqan/seq_io.h>
#include <iostream>
#include <exception>

using seqan::StringSet;
using seqan::CharString;
using seqan::Dna5String;
using seqan::SeqFileIn;

using seqan::open;
using seqan::readRecords;

using std::exception;

void read_fasta(StringSet<CharString>* pids,
                      StringSet<Dna5String>* pseqs,
                      char *ont_reads_filename) {
    auto& ids = *pids;
    auto& seqs = *pseqs;

    // opening input file
    SeqFileIn reads_file;
    if (!open(reads_file, ont_reads_filename)) {
        std::cerr << "Error: Could not open file "
                  << ont_reads_filename << "\n";
        exit(1);
    }

    // read all reads in file
    try {
        readRecords(ids, seqs, reads_file);
    } catch(exception const& e) {
        std::cout << "Error: " << e.what() << std::endl;
        exit(1);
    }
}

