/**
 * @file scaffolder.h
 * @copyright Marko Culinovic <marko.culinovic@gmail.com>
 * @copyright Luka Sterbic <luka.sterbic@gmail.com>
 * @brief Header file for scaffolder namespace.
 * @details Header file for scaffolder namespace. It provides
 * functions for different methods of contig extension.
 */

#ifndef SCAFFOLDER_H
#define SCAFFOLDER_H

#include <seqan/sequence.h>
#include <vector>
#include <string>
#include <utility>
#include <unordered_map>
#include <memory>

#include "extension.h"
#include "contig.h"

using std::vector;
using std::string;
using std::pair;
using std::unordered_map;
using std::shared_ptr;

using seqan::Dna5String;
using seqan::toCString;
using seqan::BamAlignmentRecord;


/**
 * @brief Scaffolder namespace provides
 * functions for different methods of contig extension.
 * @details Scaffolder namespace provides
 * functions for different methods of contig extension.
 * Simple majority vote method, local realignment with
 * global realignment method and POA consensus method.
 */
namespace scaffolder {


/**
 * @brief Methods sets maximum extension length.
 * @details Method limits length of extension to
 * value given as parameter for all methods.
 *
 * @param length Maximum extension length.
 */
void set_max_extension_len(int length);



/**
 * @brief Sets the inner margin.
 * @details Sets the value in BP for the margin of a read to be considered
 * directly in the extension base computation.
 *
 * @param margin the desired inner margin in base pairs
 */
void set_inner_margin(int margin);


/**
 * @brief Sets the outher margin.
 * @details Sets the value in BP for the margin of a read to be considered for
 * potential global realignment.
 *
 * @param margin the desired outer margin in base pairs
 */
void set_outer_margin(int margin);


/**
 * @brief Method sets the minimum coverage.
 * @details Sets the minimum numbers of reads that have to be present in the
 * mojority vote to output an extension base.
 *
 * @param coverage the desired coverage
 */
void set_min_coverage(int coverage);


/**
 * @brief Method finds substrings of reads which extend contig
 * on both ends.
 * @details Every record/read mapped to contig needs to be checked
 * if it is a possible extension. A record is suitable for extending
 * contig if it is soft clipped and clipped part extends left of contig
 * start, or if it is soft clipped and clipped part extends right
 * of contig start. Additional criteria is introduced to check if
 * alignment record, i.e. read, is feasible for contig extension and
 * these are inner and outer margin. Reads whose alignment starts
 * (left extension) or ends (right extension) within the inner margin is
 * immediately suitable for extension. Reads whose alignment
 * starts or ends within the outer margin are called dropped reads and
 * are later used in global realignment method.
 *
 * @param aln_records Records from SAM file
 * @param pleft_extensions Pointer to possible left end extensions
 * @param pright_extensions Pointer to possible right end extensions
 * @param contig_len Length of contig
 */
void find_possible_extensions(const vector<BamAlignmentRecord>& aln_records,
                              vector<string>* pleft_extensions,
                              vector<string>* pright_extensions,
                              uint64_t contig_len);


/**
 * @brief Method finds contig extension using majority vote on each
 * position of possible extensions while coverage >= k
 * @deprecated
 *
 * @param extensions Possible contig extensions strings
 * @return Contig extension.
 */
string get_extension_mv_simple(const vector<string>& extensions);


/**
 * @brief Method finds contig extension using majority vote on each
 * position of possible extensions while coverage >= k and tries
 * to realign reads which base isn't elected as majority
 * @details Variation of simple extension majority vote method is to
 * locally realign reads which bases at current extension position are
 * not the same as base calculated with majority vote. Local realignment
 * "looks" one move ahead, and checks if one of the alignment operations
 * (match, mismatch, deletion, insertion) could correct read alignment to
 * (extended) contig. "Looking ahead" is done by calculating the next
 * base by majority vote, but only reads withcorrect base at current
 * position are considered eligible for counting.
 *
 * @param extensions Possible contig extensions.
 * @return Resulting contig extension.
 */
string get_extension_mv_realign(const vector<shared_ptr<Extension>>&
                                extensions);


/**
 * @brief Method tries to extend contig using global realignment
 * method on both sides with given alignment records.
 * @details This process is iterative - in each step, after contig
 * extension is found using local realignment method, dropped reads are
 * globally realigned and new possible extensions of already extended
 * contig are found. If the coverage of left and right possible
 * extensions are both below the minimum coverage, contig cannot be
 * extended anymore and process is stopped.
 *
 * @param contig_seq the Sequence of the contig to be extended
 * @param aln_records Alignment records from SAM file
 * @param read_name_to_id Mapping from read name to integer ID.
 * @param read_ids Reads names / string IDs.
 * @param read_seqs Reads sequnces.
 *
 * @return Contig extended on both sides
 */
Contig* extend_contig(Dna5String& contig_seq,
                      const vector<BamAlignmentRecord>& aln_records,
                      const unordered_map<string, uint32_t>& read_name_to_id,
                      const StringSet<CharString>& read_ids,
                      const StringSet<Dna5String>& read_seqs);


/**
 * @brief Method tries to extend contig using POA consensus
 * method on both sides with given alignment records.
 * @details Using POA generated consesus to extend contigs is actually
 * a very simple idea. Local realignment method is used to get possible
 * extension reads, substrings of defined length are created for each
 * extension and these substrings are passed to the POA algorithm for
 * generating consensus of these sequences. Afterwards,
 * consensus sequences are appended to contig as left and right
 * extensions.
 *
 * @param contig_seq the Sequence of the contig to be extended
 * @param aln_records Alignment records from SAM file
 * @param read_name_to_id Mapping from read name to integer ID.
 *
 * @return Contig extended on both sides.
 */
Contig* extend_contig_poa(const Dna5String& contig_seq,
                          const vector<BamAlignmentRecord>& aln_records,
                          const unordered_map<string, uint32_t>&
                          read_name_to_id);


}  // namespace scaffolder

#endif  // SCAFFOLDER_H
