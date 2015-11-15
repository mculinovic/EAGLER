/**
 * @file bwa.cpp
 * @copyright Marko Culinovic <marko.culinovic@gmail.com>
 * @copyright Luka Sterbic <luka.sterbic@gmail.com>
 * @brief Implementation file for aligner namespace.
 * @details Implementation file for aligner namespace. It consists of various
 * bwa tool wrapper functions which make system calls to execute
 * these bwa commands.
 */
#include <seqan/sequence.h>
#include <cstdlib>
#include <string>

#include "./bwa.h"
#include "./utility.h"

using std::string;


void BwaAligner::index(const char *filename) {
    utility::execute_command("bwa index %th 2> /dev/null", filename);
}


void BwaAligner::align(const char *reference_file, const char *reads_file,
    const char *sam_file, bool only_primary) {
    utility::execute_command(
        "bwa mem -t %d -x pacbio %s %th %th > %th 2> /dev/null",
        utility::get_concurrency_level(),
        only_primary ? "" : "-Y",
        reference_file,
        reads_file,
        sam_file);
}


void BwaAligner::align(const char *reference_file, const char *reads_file,
    const char *sam_file) {
    align(reference_file, reads_file, sam_file, false);
}


void BwaAligner::align(const char *reference_file, const char *reads_file) {
    align(reference_file, reads_file, Aligner::get_tmp_alignment_filename(),
          false);
}


void BwaAligner::align(const CharString &id, const Dna5String &contig,
           const char *reads_filename) {
    // write contig to temporary .fasta file
    utility::write_fasta(id, contig, Aligner::get_tmp_contig_filename());

    // create index for contig
    index(Aligner::get_tmp_contig_filename());

    // align reads to conting
    align(Aligner::get_tmp_contig_filename(), reads_filename);
}
