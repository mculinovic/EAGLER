// Copyright @mculinovic
#ifndef SCAFFOLDER_H
#define SCAFFOLDER_H

#include <seqan/sequence.h>
#include <algorithm>
#include <vector>
#include <string>
#include <iostream>
#include "./utility.h"

#define UNMAPPED 0x4
#define NUM_BASES 4

using std::vector;
using std::string;
using std::reverse;
using std::invalid_argument;

using seqan::CharString;
using seqan::Dna5String;
using seqan::StringSet;
using seqan::toCString;
using seqan::String;
using seqan::CStyle;
using seqan::length;


namespace scaffolder {

// method finds substrings of reads which extend contig on both ends
void find_possible_extensions(const vector<BamAlignmentRecord>& aln_records,
                              vector<string>* pleft_extensions,
                              vector<string>* pright_extensions,
                              uint64_t contig_len) {
    auto& left_extensions = *pleft_extensions;
    auto& right_extensions = *pright_extensions;

    // used to determine read length from cigar string
    auto contributes_to_seq_len = [](char c) -> int {
        switch (c) {
            case 'M': return 1;  // alignment match
            case 'I': return 1;  // insertion to reference
            case 'S': return 1;  // soft clipping
            case 'X': return 1;  // sequence mismatch
            case '=': return 1;  // sequence match
            default: return 0;
        }
    };

    // used to determine length of contig part to which
    // read is aligned
    auto contributes_to_contig_len = [](char c) -> int {
        switch (c) {
            case 'M': return 1;  // alignment match
            case 'D': return 1;  // deletion from reference
            case 'X': return 1;  // sequence mismatch
            case '=': return 1;  // sequence match
            default: return 0;
        }
    };

    for (auto const& record: aln_records) {
        // record is suitable for extending contig

        // if it is soft clipped and clipped part extends
        // left of contig start
        // example:
        // contig ->     ------------
        // read ->  ----------
        if ((record.flag & UNMAPPED) == 0 &&
            record.cigar[0].operation == 'S' &&
            record.cigar[0].count > record.beginPos) { // <-- ?????????????
            int len = record.cigar[0].count - record.beginPos;
            String<char, CStyle> tmp = record.seq;
            string seq(tmp);
            // std::cout << record.qName << " " << len << std::endl;
            string extension = seq.substr(0, len);
            // reverse it because when searching for next base
            // in contig extension on left side we're moving
            // in direction right to left: <--------
            reverse(extension.begin(), extension.end());
            left_extensions.emplace_back(extension);
        }

        // if it is soft clipped and clipped part extends
        // right of contig start
        // example:
        // contig ->  ------------
        // read ->            ----------
        if ((record.flag & UNMAPPED) == 0 &&
            record.cigar[length(record.cigar) - 1].operation == 'S') {
            // iterate over cigar string to get lengths of
            // read and contig parts used in alignment
            int used_read_size = 0;
            int used_contig_size = 0;
            for (auto const& e: record.cigar) {
                if (contributes_to_seq_len(e.operation)) {
                    used_read_size += e.count;
                }
                if (contributes_to_contig_len(e.operation)) {
                    used_contig_size += e.count;
                }
                // std::cout << e.operation << " " << used_contig_size << std::endl;
            }
            int right_clipping_len = record.cigar[length(record.cigar) - 1].count;
            used_read_size -= right_clipping_len;
            int len = right_clipping_len -
                      (contig_len - (record.beginPos + used_contig_size));

            // if read doesn't extend right of contig skip it
            if (len <= 0)
                continue;
            String<char, CStyle> tmp = record.seq;
            string seq(tmp);
            // std::cout << record.qName << " " << seq.length() << " " << used_read_size + (right_clipping_len - len) << " " << len << std::endl;
            string extension = seq.substr(used_read_size + (right_clipping_len - len), len);
            right_extensions.emplace_back(extension);
        }
    }
}


string get_extension(const vector<string> extensions, int k) {
    // get array index for given base
    auto base_to_idx = [](char c) -> int {
        switch (c) {
            case 'A': return 0;
            case 'T': return 1;
            case 'G': return 2;
            case 'C': return 3;
        }
        throw invalid_argument("Illegal base character.");
    };

    // get genome base from array index
    auto idx_to_base = [](int idx) -> char {
        switch (idx) {
            case 0: return 'A';
            case 1: return 'T';
            case 2: return 'G';
            case 3: return 'C';
        }
        throw invalid_argument("Illegal base ID.");
    };

    // calculate extension by majority vote
    string extension("");
    int coverage = extensions.size();
    if (coverage >= k) {
        unsigned int i = 0;
        do {
            coverage = 0;
            vector<int> base(4, 0);
            int max_idx = -1;
            int max = -1;
            for (auto const& read: extensions) {
                if (read.size() <= i)
                    continue;
                ++coverage;
                int idx = base_to_idx(read[i]);
                base[idx]++;
                if (base[idx] > max) {
                    max_idx = idx;
                    max = base[idx];
                }
            }
            if (coverage >= k) {
                extension.push_back(idx_to_base(max_idx));
            }
            ++i;
        } while (coverage >= k);
    }
    return extension;
}


// method tries to extend contig on both sides using read alignments
Dna5String extend_contig(const CharString& contig_id,
                   const Dna5String& contig_seq,
                   char *alignment_filename) {
    // read alignment data
    BamHeader header;
    vector<BamAlignmentRecord> aln_records;
    utility::read_sam(&header, &aln_records, alignment_filename);

    vector<string> left_extensions;
    vector<string> right_extensions;

    find_possible_extensions(aln_records,
                            &left_extensions,
                            &right_extensions,
                            length(contig_seq));

    int k = 5;  // minimum coverage for position
    string left_extension = get_extension(left_extensions, k);
    // reverse it because we want to have it in direction
    // left to right -------->
    reverse(left_extension.begin(), left_extension.end());
    string right_extension = get_extension(right_extensions, k);

    Dna5String extended_contig = left_extension;
    extended_contig += contig_seq;
    extended_contig += right_extension;
    return extended_contig;
}

}

#endif  // SCAFFOLDER_H
