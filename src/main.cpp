// Copyright @mculinovic
#include <seqan/sequence.h>
#include <parsero/parsero.h>
#include <iostream>
#include "./utility.h"
#include "./bwa.h"
#include "./scaffolder.h"


using seqan::StringSet;
using seqan::CharString;
using seqan::Dna5String;
using seqan::appendValue;


char *reads_filename = nullptr;
char *draft_genome_filename = nullptr;
char *result_filename = nullptr;

// using parsero library for command line settings
void setup_cmd_interface(int argc, char **argv) {
    // argument - oxford nanopore reads in fasta format
    parsero::add_argument("ont_reads.fasta",
        [] (char *filename) { reads_filename = filename; });
    // argument - draft genome in fasta format
    parsero::add_argument("draft_genome.fasta",
        [] (char *filename) { draft_genome_filename = filename; });
    // argument - output file in fasta format
    parsero::add_argument("output_file.fasta",
        [] (char *filename) { result_filename = filename; });

    parsero::parse(argc, argv);
}


int main(int argc, char **argv) {
    setup_cmd_interface(argc, argv);

    if (reads_filename == nullptr || draft_genome_filename == nullptr) {
        parsero::help(argv[0]);
        exit(1);
    }

    // read oxford nanopore long reads from file
    StringSet<CharString> read_ids;
    StringSet<Dna5String> read_seqs;
    utility::read_fasta(&read_ids, &read_seqs, reads_filename);

    // read contigs from draft genome file
    StringSet<CharString> contig_ids;
    StringSet<Dna5String> contig_seqs;
    utility::read_fasta(&contig_ids, &contig_seqs, draft_genome_filename);
    /*
    // copy file to temporary folder to avoid data folder polution
    utility::write_fasta(contig_ids, contig_seqs,
                         aligner::tmp_reference_filename);

    // align reads to draft genome
    aligner::bwa_index(aligner::tmp_reference_filename);
    aligner::bwa_mem(aligner::tmp_reference_filename, reads_filename);
    */

    std::cout << "map - start" << std::endl;
    alignment_collection contig_alignments;
    utility::map_alignments(aligner::alignment_filename, &contig_alignments);
    std::cout << "map - end" << std::endl;

    StringSet<Dna5String> result_contig_seqs;
    int contigs_size = length(contig_ids);
    for (int i = 0; i < contigs_size; ++i) {
        // for every contig do following

        // first method when aligning all reads to contig
/*        // 1. align reads to it using bwa
        aligner::align(contig_ids[i], contig_seqs[i], reads_filename);

        // 2. try to extend it using alignments
        // TODO(mculinovic, lukasterbic): extension algorithm
        Dna5String contig = scaffolder::extend_contig(
                                contig_seqs[i],
                                aligner::alignment_filename);*/

        std::cout << i << " " << contig_alignments[i].size() << std::endl;
        Dna5String contig = scaffolder::extend_contig(contig_seqs[i],
                                                contig_alignments[i]);
        appendValue(result_contig_seqs, contig);
    }

    utility::write_fasta(contig_ids, result_contig_seqs, result_filename);

    return 0;
}
