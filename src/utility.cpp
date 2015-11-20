/**
 * @file utility.cpp
 * @copyright Marko Culinovic <marko.culinovic@fer.hr>
 * @copyright Luka Sterbic <luka.sterbic@fer.hr>
 * @brief Implementation of various utility functions.
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
#include <algorithm>
#include <regex>
#include <cstdarg>

#include "utility.h"


using std::vector;
using std::remove;
using std::string;
using std::runtime_error;
using std::reverse;
using std::max;
using std::exception;
using std::regex;

using seqan::StringSet;
using seqan::CharString;
using seqan::Dna5String;
using seqan::SeqFileIn;
using seqan::SeqFileOut;
using seqan::BamHeader;
using seqan::BamAlignmentRecord;
using seqan::BamFileIn;
using seqan::String;
using seqan::CStyle;

using seqan::open;
using seqan::readRecords;
using seqan::writeRecord;
using seqan::readHeader;
using seqan::readRecord;
using seqan::atEnd;



namespace utility {


unsigned int hardware_concurrency = max(
    1u, std::thread::hardware_concurrency());


char command_buffer[COMMAND_BUFFER_SIZE] = { 0 };


char error_buffer[ERROR_BUFFER_SIZE] = { 0 };


char seq_id_buffer[SEQ_ID_BUFFER_SIZE] = { 0 };


unsigned int get_concurrency_level() {
    return hardware_concurrency;
}


void set_concurrency_level(int threads) {
    if (threads > 0) {
        hardware_concurrency = threads;
    } else {
        exit_with_message("Illegal concurrency level");
    }
}


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

    // attempt write
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

    // attempt write
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
        // copy header
        readHeader(header, input_file);

        // copy records
        BamAlignmentRecord record;
        while (!atEnd(input_file)) {
            readRecord(record, input_file);
            records.emplace_back(record);
        }
    } catch(exception const& e) {
        exit_with_message(e.what());
    }
}


void read_sam(BamFileIn* pinput_file, BamHeader* pheader,
              vector<BamAlignmentRecord>* precords) {
    auto& header = *pheader;
    auto& records = *precords;
    auto& input_file = *pinput_file;

    try {
        // copy header
        readHeader(header, input_file);

        // copy records
        BamAlignmentRecord record;
        while (!atEnd(input_file)) {
            readRecord(record, input_file);
            records.emplace_back(record);
        }
    } catch(exception const& e) {
        exit_with_message(e.what());
    }
}


void map_alignments(const char *filename, AlignmentCollection *pcollection,
                    const unordered_map<string, uint32_t>& contig_name_to_id) {
    auto& collection = *pcollection;

    BamFileIn input_file;
    if (!open(input_file, filename)) {
        exit_with_message("could not open file %s", filename);
    }

    BamHeader header;
    vector<BamAlignmentRecord> records;
    read_sam(&input_file, &header, &records);

    for (auto& record : records) {
        if (record.rID != BamAlignmentRecord::INVALID_REFID &&
            (record.flag & UNMAPPED) == 0) {
            // get contig name for the current read
            CharString contig_id = getContigName(record, input_file);
            string contig_name = CharString_to_string(contig_id);

            // fetch ID record of the contig in the draft genome file
            auto draft_contig_id = contig_name_to_id.find(contig_name);

            // associate the read to the correct contig
            collection[draft_contig_id->second].emplace_back(record);
        }
    }
}


void execute_command(const char *format, ...) {
    // add quatation marks around all string arguments
    regex argument_re("([^%])%th");
    string escaped_fmt = regex_replace(format, argument_re, "$1\"%s\"");

    va_list args_list;
    va_start(args_list, format);

    vsnprintf(command_buffer, COMMAND_BUFFER_SIZE, escaped_fmt.c_str(),
              args_list);

    va_end(args_list);

    DEBUG(command_buffer);
    int exit_value = system(command_buffer);

    if (exit_value != 0) {
        throw_exception<runtime_error>(
            "command \"%s\" failed with exit status %d",
            command_buffer, exit_value);
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


void exit_with_message(const char *format, ...) {
    va_list args_list;
    va_start(args_list, format);

    fprintf(stderr, "[ERROR] ");
    vfprintf(stderr, format, args_list);
    fprintf(stderr, "\n");

    va_end(args_list);
    exit(1);
}


string CharString_to_string(const CharString& str) {
    String<char, CStyle> tmp = str;
    string cppstr(tmp);
    return cppstr;
}


string Dna5String_to_string(const Dna5String& str) {
    String<char, CStyle> tmp = str;
    string cppstr(tmp);
    return cppstr;
}


int contributes_to_seq_len(char c) {
    switch (c) {
        case 'M': return 1;  // alignment match
        case 'I': return 1;  // insertion to reference
        case 'S': return 1;  // soft clipping
        case 'X': return 1;  // sequence mismatch
        case '=': return 1;  // sequence match
        default: return 0;
    }
}


int contributes_to_contig_len(char c) {
    switch (c) {
        case 'M': return 1;  // alignment match
        case 'D': return 1;  // deletion from reference
        case 'X': return 1;  // sequence mismatch
        case '=': return 1;  // sequence match
        default: return 0;
    }
}


string reverse_complement(const Dna5String& seq) {
    string tmp = Dna5String_to_string(seq);
    reverse(tmp.begin(), tmp.end());

    auto complement = [](char base) -> char {
        switch (base) {
            case 'A': return 'T';
            case 'T': return 'A';
            case 'C': return 'G';
            case 'G': return 'C';
            default: throw_exception<runtime_error>("Illegal base");
        }
        return 'N';  // unknown base
    };

    for (size_t i = 0; i < tmp.length(); ++i) {
        tmp[i] = complement(tmp[i]);
    }

    return tmp;
}


string create_seq_id(const char *format, ...) {
    va_list args_list;
    va_start(args_list, format);

    vsnprintf(seq_id_buffer, SEQ_ID_BUFFER_SIZE, format, args_list);

    va_end(args_list);

    return string(seq_id_buffer);
}


bool is_command_available(const char* command) {
    snprintf(command_buffer, COMMAND_BUFFER_SIZE, "type \"%s\" >/dev/null 2>&1",
             command);
    return system(command_buffer) ? false : true;
}

}  // namespace utility
