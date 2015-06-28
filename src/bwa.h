/**
 * @file bwa.h
 * @author Marko Culinovic <marko.culinovic@gmail.com>
 * @author Luka Sterbic <luka.sterbic@gmail.com>
 * @brief Header file for aligner namespace.
 * @details Header file for aligner namespace. It consists of various
 * bwa tool wrapper functions which make system calls to execute
 * these bwa commands.
 */
#ifndef BWA_H
#define BWA_H


#include <seqan/sequence.h>


using seqan::CharString;
using seqan::Dna5String;

/**
 * Aligner namespace. It consists of various
 * bwa tool wrapper functions which make system calls to execute
 * these bwa commands.
 */
namespace aligner {

/** 
 * Temporary alignment filename. Used for storing results of
 * bwa mem command.
 */
extern const char *tmp_alignment_filename;


/**
 * Temporary reference filename. Used when copying draft genome file 
 * to temporary folder to avoid data folder polution when executing
 * bwa index and bwa mem commands.
 */
extern const char *tmp_reference_filename;

/**
 * Temporary contig filename. Used when contig sequence needs to
 * be written in temporary file for executing bwa index and
 * bwa mem commands.
 */
extern const char *tmp_contig_filename;


/**
 * @brief Bwa index command wrapper.
 * @details Method makes system call to execute bwa index command.
 * Command creates index for sequences stored in file.
 * 
 * @param filename FASTA file with sequences for creating index.
 */
void bwa_index(const char* filename);


/**
 * @brief Bwa mem command wrapper.
 * @details Method makes system call to execute bwa mem command.
 * Aligns reads from read file to reference in reference file (
 * index has to be created for this reference before). Alignment
 * results are stored in tmp_alignment_filename and bwa mem command
 * is executed in verbose mode.
 * 
 * @param reference_file FASTA file with reference sequence(s).
 * @param reads_file FASTA file with reads.
 */
void bwa_mem(const char *reference_file, const char *reads_file);


/**
 * @brief Bwa mem command wrapper.
 * @details Method makes system call to execute bwa mem command.
 * Aligns reads from read file to reference in reference file (
 * index has to be created for this reference before).
 * 
 * @param reference_file FASTA file with reference sequence(s).
 * @param reads_file FASTA file with reads.
 * @param sam_file SAM file for storing results of bwa mem command.
 * @param only_primary Flag for executing bwa mem in verbose mode.
 * Verbose mode is turned on if flag is set to false.
 */
void bwa_mem(const char *reference_file, const char *reads_file,
    const char *sam_file, bool only_primary);


/**
 * @brief Bwa mem command wrapper.
 * @details Method makes system call to execute bwa mem command.
 * Aligns reads from read file to reference in reference file (
 * index has to be created for this reference before). Bwa mem
 * command is executed in verbose mode.
 * 
 * @param reference_file FASTA file with reference sequence(s).
 * @param reads_file FASTA file with reads.
 * @param sam_file SAM file for storing results of bwa mem command.
 */
void bwa_mem(const char *reference_file, const char *reads_file,
    const char *sam_file);


/**
 * @brief Align reads from file to contig using bwa.
 * @details Method writes contig sequence to tmp_contig_filename,
 * creates index using bwa index command and align reads
 * to contig using bwa mem command.
 * @deprecated
 * 
 * @param id Contig id.
 * @param contig Contig sequence.
 * @param reads_filename FASTA file with reads.
 */
void align(const CharString &id, const Dna5String &contig,
           const char *reads_filename);


}  // namespace aligner

#endif  // BWA_H
