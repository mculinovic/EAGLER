// Copyright @mculinovic
#ifndef UTILITY_H
#define UTILITY_H

#include <seqan/bam_io.h>
#include <vector>
#include <string>
#include <unordered_map>

using std::vector;
using std::string;
using std::invalid_argument;
using std::unordered_map;

using seqan::StringSet;
using seqan::CharString;
using seqan::Dna5String;
using seqan::BamHeader;
using seqan::BamAlignmentRecord;

#define UNMAPPED 0x4
typedef unordered_map<int, vector<BamAlignmentRecord>> alignment_collection;


namespace utility {

/**
 * reads sequences data from fasta file and
 * stores it in two sets: sequences ids and
 * sequences
 */
void read_fasta(StringSet<CharString>* pids,
                StringSet<Dna5String>* pseqs,
                char *ont_reads_filename);

// writes given sequence id and sequence to file with filename
// @filename
void write_fasta(const CharString &id, const Dna5String &seq,
                 const char* filename);


// writes set of sequence ids and sequences to file
// with filename @filename
void write_fasta(const StringSet<CharString>& ids,
                const StringSet<Dna5String>& seqs,
                const char *filename);


// reads alignment data from sam file and stores it in containers given
// as function arguments
void read_sam(BamHeader* pheader, vector<BamAlignmentRecord>* precords,
              const char* filename);

/**
 * Reads sam file and maps aligned reads to contigs
 * @param filename sam file with alignments
 * @param collection hash map for storage
 */
void map_alignments(const char *filename,
                    alignment_collection *collection);


// wrapper for system() call
void execute_command(const string& command);


// get array index for given base
int base_to_idx(char c);


// get genome base from array index
char idx_to_base(int idx);

}  // namespace utility

#endif  // UTILITY_H
