// Copyright @mculinovic
#ifndef SCAFFOLDER_H
#define SCAFFOLDER_H

#include <seqan/sequence.h>
#include <vector>
#include <string>
#include <utility>
#include <unordered_map>
#include <memory>

#include "./extension.h"

using std::vector;
using std::string;
using std::pair;
using std::unordered_map;
using std::shared_ptr;

using seqan::Dna5String;
using seqan::toCString;
using seqan::BamAlignmentRecord;


namespace scaffolder {

/**
 * method finds substrings of reads which extend contig on both ends
 * @param aln_records records from sam file
 * @param pleft_extensions pointer to possible left end extensions
 * @param pright_extensions pointer to possible right end extensions
 * @param contig_len length of contig
 */
void find_possible_extensions(const vector<BamAlignmentRecord>& aln_records,
                              vector<string>* pleft_extensions,
                              vector<string>* pright_extensions,
                              uint64_t contig_len);


/* method finds contig extension using majority vote on each
 * position of possible extensions while coverage >= k
 * @param extensions possible contig extensions strings*/
string get_extension_mv_simple(const vector<string>& extensions);


/* method finds contig extension using majority vote on each
 * position of possible extensions while coverage >= k and tries
 * to realign reads which base isn't elected as majority
 * @param extensions possible contig extensions strings*/
string get_extension_mv_realign(const vector<shared_ptr<Extension>>& extensions);


/**
 * method tries to extend contig on both sides with given alignment
 * records
 * @param  contig_seq the sequence of the contig to be extended
 * @param  aln_records alignment records from sam file
 * @param  read_name_to_id map read name to integer id
 * @return extended contig
 */
Dna5String extend_contig(Dna5String& contig_seq,
                         const vector<BamAlignmentRecord>& aln_records,
                         const unordered_map<string, uint32_t>& read_name_to_id,
                         const StringSet<CharString>& read_ids,
                         const StringSet<Dna5String>& read_seqs);


Dna5String extend_contig_poa(const Dna5String& contig_seq,
                    const vector<BamAlignmentRecord>& aln_records,
                    const unordered_map<string, uint32_t>& read_name_to_id);

}  // namespace scaffolder

#endif  // SCAFFOLDER_H
