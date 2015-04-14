// Copyright @mculinovic
#include <seqan/sequence.h>
#include <parsero/parsero.h>
#include <iostream>
#include "./utility.h"
#include "./bwa.h"


using seqan::StringSet;
using seqan::CharString;
using seqan::Dna5String;


char *ont_reads_filename = nullptr;
char *draft_genome_filename = nullptr;

// using parsero library for command line settings
void setup_cmd_interface(int argc, char **argv) {
    // argument - oxford nanopore reads in fasta format
    parsero::add_argument("ont_reads.fasta",
        [] (char *filename) { ont_reads_filename = filename; });
    // argument - draft genome in fasta format
    parsero::add_argument("draft_genome.fasta",
        [] (char * filename) { draft_genome_filename = filename; });

    parsero::parse(argc, argv);
}


int main(int argc, char **argv) {
    setup_cmd_interface(argc, argv);

    if (ont_reads_filename == nullptr || draft_genome_filename == nullptr) {
        parsero::help(argv[0]);
        exit(1);
    }

    // read oxford nanopore long reads from file
    StringSet<CharString> read_ids;
    StringSet<Dna5String> read_seqs;
    utility::read_fasta(&read_ids, &read_seqs, ont_reads_filename);

    // read contigs from draft genome file
    StringSet<CharString> contig_ids;
    StringSet<Dna5String> contig_seqs;
    utility::read_fasta(&contig_ids, &contig_seqs, draft_genome_filename);

    int contigs_size = length(contig_ids);
    for (int i = 0; i < contigs_size; ++i) {
        // for every contig do following

        // 1. align reads to it using bwa
        aligner::align(contig_ids[i], contig_seqs[i], ont_reads_filename);

        // 2. try to extend it using alignments
        // TODO(mculinovic, lukasterbic): extension algorithm
    }

    return 0;
}
