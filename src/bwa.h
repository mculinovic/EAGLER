// Copyright @mculinovic
#ifndef BWA_H
#define BWA_H


#include <seqan/sequence.h>

using seqan::CharString;
using seqan::Dna5String;

namespace aligner {

// creates bwa index for temporary contig using command "bwa index"
void bwa_index();


// alignes reads to contig using "bwa mem" command
void bwa_mem(char* ont_reads_filename);


// align reads from file to contig using bwa
// alignemnt info is stored in aln.sam file
// arguments:
// -> @id - contig id
// -> @contig - contig sequence
// -> @ont_reads_filename - name of .fasta file with ONT reads
void align(const CharString &id, const Dna5String &contig,
           char* ont_reads_filename);

}  // namespace aligner

#endif  // BWA_H
