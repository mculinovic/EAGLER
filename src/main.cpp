/**
 * @file main.cpp
 * @author Marko Culinovic <marko.culinovic@gmail.com>
 * @author Luka Sterbic <luka.sterbic@gmail.com>
 * @brief Entry point of program and main program of project.
 */
#include <seqan/sequence.h>
#include <parsero/parsero.h>
#include <iostream>
#include <unordered_map>
#include <string>
#include <vector>
#include <utility>

#include "./utility.h"
#include "./aligner.h"
#include "./scaffolder.h"
#include "./contig.h"
#include "./connector.h"

using std::cout;
using std::endl;
using std::unordered_map;
using std::string;
using std::pair;

using seqan::StringSet;
using seqan::CharString;
using seqan::Dna5String;
using seqan::appendValue;
using seqan::String;
using seqan::CStyle;


char *reads_filename = nullptr;
char *draft_genome_filename = nullptr;
char *result_filename = nullptr;
char *extensions_filename = nullptr;
bool use_POA_consensus = false;
bool use_graphmap_aligner = false;

// using parsero library for command line settings
void setup_cmd_interface(int argc, char **argv) {
    // set footer
    string footer;
    footer += "Copyright (C) by Marko Culinovic, Luka Sterbic and Mile Sikic";
    footer += "\nEAGLER is licensed under the GNU General Public License.";
    parsero::set_footer(footer);

    // option - enable poa
    parsero::add_option("p", "use POA consensus algorithm [flag]",
        [] (char *) { use_POA_consensus = true; });
    // option - enable poa
    parsero::add_option("g", "use GraphMap aligner [flag]",
        [] (char *) { use_graphmap_aligner = true; });
    // option - set number of threads
    parsero::add_option("t:", "number of parallel threads [int]",
        [] (char *option) { utility::set_concurrency_level(atoi(option)); });
    // option - set extension size
    parsero::add_option("s:", "the maximum extension size in base pairs [int]",
        [] (char *option) { scaffolder::set_max_extension_len(atoi(option)); });

    // argument - oxford nanopore reads in fasta format
    parsero::add_argument("ont_reads.fasta",
        [] (char *filename) { reads_filename = filename; });
    // argument - draft genome in fasta format
    parsero::add_argument("draft_genome.fasta",
        [] (char *filename) { draft_genome_filename = filename; });
    // argument - output file in fasta format
    parsero::add_argument("output_file.fasta",
        [] (char *filename) { result_filename = filename; });
    // argument - exntensions output file in fasta format
    parsero::add_argument("output_extensions.fasta",
        [] (char *filename) { extensions_filename = filename; });

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

    unordered_map<string, uint32_t> read_name_to_id;
    for (uint32_t id = 0; id < length(read_ids); ++id) {
        String<char, CStyle> tmp = read_ids[id];
        string name(tmp);
        read_name_to_id[name] = id;
    }


    // read contigs from draft genome file
    StringSet<CharString> contig_ids;
    StringSet<Dna5String> contig_seqs;
    utility::read_fasta(&contig_ids, &contig_seqs, draft_genome_filename);


    // copy file to temporary folder to avoid data folder polution
    utility::write_fasta(contig_ids, contig_seqs,
                         Aligner::get_tmp_reference_filename());


    // initialize Aligner
    Aligner::init(use_graphmap_aligner);

    // create index for all contigs in draft genome
    cout << "[Aligner] creating index..." << endl;

    Aligner::get_instance().index(draft_genome_filename);

    // align all reads to the draft genome
    cout << "[Aligner] aligning reads to draft genome using ";
    cout << utility::get_concurrency_level() << " threads..." << endl;

    Aligner::get_instance().align(draft_genome_filename, reads_filename);

    cout << "[INFO] creating alignments map..." << endl;
    AlignmentCollection contig_alns;
    utility::map_alignments(Aligner::get_tmp_alignment_filename(), &contig_alns);

    StringSet<Dna5String> result_contig_seqs;
    StringSet<Dna5String> extensions;
    StringSet<CharString> ext_ids;
    vector< Contig* > contigs;
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

        std::cout << i << " " << contig_alns[i].size() << std::endl;
        Dna5String contig_seq;
        Contig *contig = nullptr;
        if (use_POA_consensus) {
            contig = scaffolder::extend_contig_poa(contig_seqs[i],
                                                   contig_alns[i],
                                                   read_name_to_id);
        } else {
            cout << "### len before: " << length(contig_seqs[i]) << endl;
            contig = scaffolder::extend_contig(contig_seqs[i],
                                               contig_alns[i],
                                               read_name_to_id,
                                               read_ids,
                                               read_seqs);
            cout << "### len after: " << length(contig->seq()) << endl;
        }
        contigs.emplace_back(contig);
        contig->set_id(contig_ids[i]);
        appendValue(result_contig_seqs, contig->seq());

        // extensions
        appendValue(ext_ids, contig->left_id());
        appendValue(extensions, contig->ext_left());

        appendValue(ext_ids, contig->right_id());
        appendValue(extensions, contig->ext_right());
    }

    utility::write_fasta(ext_ids, extensions, extensions_filename);
    utility::write_fasta(contig_ids, result_contig_seqs, result_filename);

    Connector connector(contigs);
    connector.connect_contigs();

    connector.dump_scaffolds();

    // free memory
    for (auto contig : contigs) {
        delete contig;
    }

    return 0;
}
