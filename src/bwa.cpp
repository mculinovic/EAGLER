// Copyright @mculinovic
#include <seqan/sequence.h>
#include <cstdlib>
#include <string>

#include "./bwa.h"
#include "./utility.h"

using std::string;


const char* contig_tmp_filename = "./tmp/contig_tmp.fasta";
const char* alignemnt_filename = "./tmp/aln.sam";

namespace aligner {


void bwa_index() {
    std::cout << "Creating bwa index for contig" << std::endl;
    string command("bwa index ");
    command += contig_tmp_filename;
    utility::execute_command(command);
    std::cout << "Bwa index created" << std::endl;
}


void bwa_mem(char* ont_reads_filename) {
    std::cout << "Aligning reads to contig" << std::endl;
    string command("bwa mem -x ont2d ");
    command += contig_tmp_filename;
    command += " ";
    command += ont_reads_filename;
    command += " > ";
    command += alignemnt_filename;
    utility::execute_command(command);
    std::cout << "Alignment finished" << std::endl;
}


void align(const CharString &id, const Dna5String &contig,
           char* ont_reads_filename) {
    // write contig to temporary .fasta file
    utility::write_fasta(id, contig, contig_tmp_filename);

    // create index for contig
    bwa_index();

    // align reads to conting
    bwa_mem(ont_reads_filename);

    // delete temporary file
    utility::delete_file(contig_tmp_filename);
}


}  // namespace aligner
