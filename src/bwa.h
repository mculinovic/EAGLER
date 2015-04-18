// Copyright @mculinovic
#ifndef BWA_H
#define BWA_H


#include <seqan/sequence.h>
#include <cstdlib>
#include <string>

#include "./utility.h"

using std::string;


const char* contig_tmp_filename = "../data/tmp/contig_tmp.fasta";
const char* alignemnt_filename = "../data/tmp/aln.sam";

namespace aligner {

// creates bwa index for temporary contig using command "bwa index"
void bwa_index() {
    std::cout << "Creating bwa index for contig" << std::endl;
    string command("bwa index ");
    command += contig_tmp_filename;
    system(command.c_str());
    std::cout << "Bwa index created" << std::endl;
}


// alignes reads to contig using "bwa mem" command
void bwa_mem(char* ont_reads_filename) {
    std::cout << "Aligning reads to conting" << std::endl;
    string command("bwa mem -x ont2d ");
    command += contig_tmp_filename;
    command += " ";
    command += ont_reads_filename;
    command += " > ";
    command += alignemnt_filename;
    system(command.c_str());
    std::cout << "Alignment finished" << std::endl;
}


// align reads from file to contig using bwa
// alignemnt info is stored in aln.sam file
// arguments:
// -> @id - contig id
// -> @contig - contig sequence
// -> @ont_reads_filename - name of .fasta file with ONT reads
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

} // namespace aligner

#endif  // BWA_H
