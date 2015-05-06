// Copyright @mculinovic
#ifndef BWA_H
#define BWA_H


#include <seqan/sequence.h>

using seqan::CharString;
using seqan::Dna5String;


namespace aligner {

extern const char *alignment_filename;
extern const char *tmp_reference_filename;
extern const char *contig_tmp_filename;

/**
 * creates bwa index for temporary contig using command "bwa index"
 * @param filename fasta file for creating index
 */
void bwa_index(const char* filename);


/**
 * alignes reads to contig using "bwa mem" command
 * @param reference_file fasta file with reference sequence
 * @param reads_file fasta file with reads
 */
void bwa_mem(const char *reference_file, const char *reads_file);


/**
 * align reads from file to contig using bwa
 * alignemnt info is stored in aln.sam file
 * @param id contig id
 * @param contig contig sequence
 * @param reads_filename name of .fasta file with reads
 */
void align(const CharString &id, const Dna5String &contig,
           const char *reads_filename);

}  // namespace aligner

#endif  // BWA_H
