// Copyright @mculinovic
#include <seqan/seq_io.h>
#include <seqan/bam_io.h>
#include <iostream>
#include <exception>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <string>
#include <stdexcept>

#include "./utility.h"

using std::vector;
using std::remove;
using std::string;
using std::runtime_error;

using seqan::StringSet;
using seqan::CharString;
using seqan::Dna5String;
using seqan::SeqFileIn;
using seqan::SeqFileOut;
using seqan::BamHeader;
using seqan::BamAlignmentRecord;
using seqan::BamFileIn;

using seqan::open;
using seqan::readRecords;
using seqan::writeRecord;
using seqan::readHeader;
using seqan::readRecord;
using seqan::atEnd;

using std::exception;

namespace utility {


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


void write_fasta(const CharString &id, const Dna5String &seq,
                 const char* filename) {
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


void write_fasta(const StringSet<CharString>& ids,
                const StringSet<Dna5String>& seqs,
                const char *filename) {
    // opening output file
    SeqFileOut out_file;
    if (!open(out_file, filename)) {
        std::cerr << "Error: Could not open file "
                  << filename << "\n";
        exit(1);
    }

    for (uint32_t i = 0; i < length(ids); ++i) {
        try {
            writeRecord(out_file, ids[i], seqs[i]);
        } catch(exception const& e) {
            std::cout << "Error: " << e.what() << std::endl;
            exit(1);
        }
    }
}


void read_sam(BamHeader* pheader, vector<BamAlignmentRecord>* precords,
              const char* filename) {
    auto& header = *pheader;
    auto& records = *precords;

    BamFileIn input_file;
    if (!open(input_file, filename)) {
        std::cerr << "Error: Could not open " << filename << std::endl;
        exit(1);
    }

    try {
        // Copy header
        readHeader(header, input_file);

        // Copy records.
        BamAlignmentRecord record;
        while (!atEnd(input_file)) {
            readRecord(record, input_file);
            records.emplace_back(record);
        }
    } catch(exception const & e) {
        std::cout << "Error: " << e.what() << std::endl;
        exit(1);
    }
}


void delete_file(const char* filename) {
    if (remove(filename)) {
        std::cerr << "Error deleting file" << std :: endl;
        exit(1);
    } else {
        std::cout << "Temp file successfully deleted" << std::endl;
    }
}


void execute_command(const string& command) {
    int ret = system(command.c_str());
    if (ret != 0) {
        string desc = "Command \"";
        desc += command;
        desc += "\" failed with exit status ";
        char numstr[21];
        snprintf(numstr, 21, "%d", ret);
        desc += numstr;
        desc += "!";
        throw runtime_error(desc);
    }
}

}  // namespace utility