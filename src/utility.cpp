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


void map_alignments(const char *filename,
                    alignment_collection *pcollection) {
    auto& collection = *pcollection;

    BamHeader header;
    vector<BamAlignmentRecord> records;
    read_sam(&header, &records, filename);

    std::cout << "creating map" << std::endl;
    for (auto& record: records) {
        if (record.rID != BamAlignmentRecord::INVALID_REFID &&
            (record.flag & UNMAPPED) == 0) {
            collection[record.rID].emplace_back(record);
        }
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

// get array index for given base
int base_to_idx(char c) {
    switch (c) {
        case 'A': return 0;
        case 'T': return 1;
        case 'G': return 2;
        case 'C': return 3;
    }
    throw invalid_argument("Illegal base character.");
}


// get genome base from array index
char idx_to_base(int idx) {
    switch (idx) {
        case 0: return 'A';
        case 1: return 'T';
        case 2: return 'G';
        case 3: return 'C';
    }
    throw invalid_argument("Illegal base ID.");
}

}  // namespace utility
