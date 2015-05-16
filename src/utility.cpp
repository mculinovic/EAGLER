/**
 * @file utility.cpp
 * @author Marko Culinovic <marko.culinovic@fer.hr>
 * @author Luka Sterbic <luka.sterbic@fer.hr>
 * @brief Implementation of various utility functions
 * @details Implementation of utility functions for genomic data I/O, shell
 * commands execution and data conversion.
 */

#include <seqan/seq_io.h>
#include <seqan/bam_io.h>
#include <iostream>
#include <exception>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <string>
#include <stdexcept>
#include <thread>
#include <algorithm>
#include <cstdarg>

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

using std::max;
using std::exception;


namespace utility {


const unsigned int hardware_concurrency = max(
    1u, std::thread::hardware_concurrency());


char command_buffer[COMMAND_BUFFER_SIZE] = { 0 };


char error_buffer[ERROR_BUFFER_SIZE] = { 0 };


void read_fasta(StringSet<CharString>* pids, StringSet<Dna5String>* pseqs,
                char *ont_reads_filename) {
    auto& ids = *pids;
    auto& seqs = *pseqs;

    // opening input file
    SeqFileIn input_file;
    if (!open(input_file, ont_reads_filename)) {
        exit_with_message("Could not open file %s", ont_reads_filename);
    }

    // read all reads in file
    try {
        readRecords(ids, seqs, input_file);
    } catch(exception const& e) {
        exit_with_message(e.what());
    }
}


void write_fasta(const CharString &id, const Dna5String &seq,
                 const char* filename) {
    // opening output file
    SeqFileOut out_file;
    if (!open(out_file, filename)) {
        exit_with_message("Could not open file %s", filename);
    }

    try {
        writeRecord(out_file, id, seq);
    } catch(exception const& e) {
        exit_with_message(e.what());
    }
}


void write_fasta(const StringSet<CharString>& ids,
                const StringSet<Dna5String>& seqs,
                const char *filename) {
    // opening output file
    SeqFileOut out_file;
    if (!open(out_file, filename)) {
        exit_with_message("Could not open file %s", filename);
    }

    for (uint32_t i = 0; i < length(ids); ++i) {
        try {
            writeRecord(out_file, ids[i], seqs[i]);
        } catch(exception const& e) {
            exit_with_message(e.what());
        }
    }
}


void read_sam(BamHeader* pheader, vector<BamAlignmentRecord>* precords,
              const char* filename) {
    auto& header = *pheader;
    auto& records = *precords;

    BamFileIn input_file;
    if (!open(input_file, filename)) {
        exit_with_message("could not open file %s", filename);
    }

    try {
        // Copy header
        readHeader(header, input_file);

        // Copy records
        BamAlignmentRecord record;
        while (!atEnd(input_file)) {
            readRecord(record, input_file);
            records.emplace_back(record);
        }
    } catch(exception const& e) {
        exit_with_message(e.what());
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


void execute_command(const char *format, ...) {
    va_list args_list;
    va_start(args_list, format);

    vsnprintf(command_buffer, COMMAND_BUFFER_SIZE, format, args_list);

    va_end(args_list);

    int exit_value = system(command_buffer);

    if (exit_value != 0) {
        throw_exception<runtime_error>(
            "command \"%s\" failed with exit status %d",
            command_buffer,
            exit_value
        );
    }
}


int base_to_idx(char base) {
    switch (base) {
        case 'A': return 0;
        case 'T': return 1;
        case 'G': return 2;
        case 'C': return 3;
    }

    throw invalid_argument("Illegal base character.");
}


char idx_to_base(int idx) {
    switch (idx) {
        case 0: return 'A';
        case 1: return 'T';
        case 2: return 'G';
        case 3: return 'C';
    }

    throw invalid_argument("Illegal base ID.");
}

template<typename T>
void throw_exception(const char *format, ...) {
    va_list args_list;
    va_start(args_list, format);

    vsnprintf(error_buffer, ERROR_BUFFER_SIZE, format, args_list);
    va_end(args_list);

    throw T(error_buffer);
}


void exit_with_message(const char *format, ...) {
    va_list args_list;
    va_start(args_list, format);

    fprintf(stderr, "[ERROR] ");
    vfprintf(stderr, format, args_list);

    va_end(args_list);
    exit(1);
}


}  // namespace utility
