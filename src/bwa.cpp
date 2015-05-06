// Copyright @mculinovic
#include <seqan/sequence.h>
#include <cstdlib>
#include <string>

#include "./bwa.h"
#include "./utility.h"

using std::string;

namespace aligner {

const char *alignment_filename = "./tmp/aln.sam";
const char *tmp_reference_filename = "./tmp/reference.fasta";
const char *contig_tmp_filename = "./tmp/contig_tmp.fasta";


void bwa_index(const char *filename) {
    string command("bwa index ");
    command += filename;
    utility::execute_command(command);
}


void bwa_mem(const char *reference_file, const char *reads_file) {
    // string command("bwa mem -x ont2d ");
    string command("bwa mem -t 4 -x pacbio -Y ");
    command += reference_file;
    command += " ";
    command += reads_file;
    command += " > ";
    command += alignment_filename;
    command += "2> /dev/null";

    utility::execute_command(command);
}


void align(const CharString &id, const Dna5String &contig,
           const char *reads_filename) {
    // write contig to temporary .fasta file
    utility::write_fasta(id, contig, contig_tmp_filename);

    // create index for contig
    bwa_index(contig_tmp_filename);

    // align reads to conting
    bwa_mem(contig_tmp_filename, reads_filename);
}


}  // namespace aligner
