// Copyright @mculinovic
#include <seqan/sequence.h>
#include <cstdlib>
#include <string>

#include "./bwa.h"
#include "./utility.h"

using std::string;

namespace aligner {


const char *tmp_alignment_filename = "./tmp/aln.sam";
const char *tmp_reference_filename = "./tmp/reference.fasta";
const char *tmp_contig_filename = "./tmp/contig_tmp.fasta";


void bwa_index(const char *filename) {
    utility::execute_command("bwa index %s 2> /dev/null", filename);
}


void bwa_mem(const char *reference_file, const char *reads_file, const char *sam_file) {
    utility::execute_command(
        "bwa mem -t %d -x pacbio -Y %s %s > %s 2> /dev/null",
        utility::hardware_concurrency,
        reference_file,
        reads_file,
        sam_file
    );
}


void bwa_mem(const char *reference_file, const char *reads_file) {
    bwa_mem(reference_file, reads_file, tmp_alignment_filename);
}


void align(const CharString &id, const Dna5String &contig,
           const char *reads_filename) {
    // write contig to temporary .fasta file
    utility::write_fasta(id, contig, tmp_contig_filename);

    // create index for contig
    bwa_index(tmp_contig_filename);

    // align reads to conting
    bwa_mem(tmp_contig_filename, reads_filename);
}


}  // namespace aligner
