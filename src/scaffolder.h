// Copyright @mculinovic
#ifndef SCAFFOLDER_H
#define SCAFFOLDER_H

#include <seqan/sequence.h>
#include <vector>
#include "./utility.h"

using std::vector;

using seqan::CharString;
using seqan::Dna5String;
using seqan::StringSet;


namespace scaffolder {

    // method tries to extend contig on both sides using read alignments
    void extend_contig(const CharString& contig_id,
                       const Dna5String& contig_seq,
                       const StringSet<CharString>& read_ids,
                       const StringSet<Dna5String>& read_seqs,
                       char *alignment_filename) {
        // read alignment data
        BamHeader header;
        vector<BamAlignmentRecord> aln_records;
        utility::read_sam(&header, &aln_records, alignment_filename);

        // test output
        for (int i = 0; i < aln_records.size(); ++i) {
            std::cout << "*** Proccessing new read ***" << std::endl;
            std::cout << aln_records[i].qName << std::endl;
            for (auto const& e: aln_records[i].cigar) {
                std::cout << e.count << e.operation;
            }
            std::cout << std::endl;
        }
    }
};

#endif  // SCAFFOLDER_H
