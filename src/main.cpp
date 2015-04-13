// Copyright @mculinovic
#include <seqan/sequence.h>
#include <parsero/parsero.h>
#include <iostream>
#include "./fasta_IO.h"


using seqan::StringSet;
using seqan::CharString;
using seqan::Dna5String;


char *tmp_contig_filename = "../data/contig_tmp.fasta";
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

    StringSet<CharString> read_ids;
    StringSet<Dna5String> read_seqs;
    read_fasta(&read_ids, &read_seqs, ont_reads_filename);

    StringSet<CharString> contig_ids;
    StringSet<Dna5String> contig_seqs;
    read_fasta(&contig_ids, &contig_seqs, draft_genome_filename);


    for (int i = 0; i < length(read_ids); ++i) {
        std::cout << read_ids[i] << '\t' << read_seqs[i] << '\n';
    }

    return 0;
}
