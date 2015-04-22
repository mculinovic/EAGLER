// Copyright @mculinovic
#ifndef SCAFFOLDER_H
#define SCAFFOLDER_H

#include <seqan/sequence.h>
#include <vector>
#include <string>

using std::vector;
using std::string;

using seqan::Dna5String;
using seqan::toCString;
using seqan::BamAlignmentRecord;


namespace scaffolder {

// method finds substrings of reads which extend contig on both ends
void find_possible_extensions(const vector<BamAlignmentRecord>& aln_records,
                              vector<string>* pleft_extensions,
                              vector<string>* pright_extensions,
                              uint64_t contig_len);


// method finds contig extension using majority vote on each
// position of possible extensions while coverage >= k
// @extensions - possible contig extensions strings
// @k - minimum coverage for position
string get_extension(const vector<string> extensions, int k);


// method tries to extend contig on both sides using read alignments
Dna5String extend_contig(const Dna5String& contig_seq,
                   const char *alignment_filename);


}  // namespace scaffolder

#endif  // SCAFFOLDER_H
