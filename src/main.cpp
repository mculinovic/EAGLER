/**
 * @file main.cpp
 * @copyright Marko Culinovic <marko.culinovic@gmail.com>
 * @copyright Luka Sterbic <luka.sterbic@gmail.com>
 * @brief Entry point of program and main program of project.
 */
#include <seqan/sequence.h>
#include <parsero/parsero.h>
#include <iostream>
#include <unordered_map>
#include <string>
#include <vector>
#include <utility>

#include "aligners/aligner.h"
#include "utility.h"
#include "scaffolder.h"
#include "contig.h"
#include "connector.h"

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
read_type::ReadType use_tech_type = read_type::PacBio;


// using parsero library for command line settings
void setup_cmd_interface(int argc, char **argv) {
    // set header
    string header;

    header += "EAGLER is a scaffolding tool for long reads. The scaffolder ";
    header += "takes as input a draft genome created by any NGS assembler and ";
    header += "a set of long reads. The long reads are used to extend the ";
    header += "contigs present in the NGS draft and possibly join overlapping ";
    header += "contigs. EAGLER supports both PacBio and Oxford Nanopore reads.";

    parsero::set_header(header);

    // set footer
    string footer;

    footer += "Copyright (C) by Marko Culinovic, Luka Sterbic and Mile Sikic";
    footer += "\nEAGLER is licensed under the GNU General Public License.";

    parsero::set_footer(footer);

    // option - enable poa, hack to avoid unused variable warning
    parsero::add_option("p", "use POA consensus algorithm [flag]",
        [] (char *option) { use_POA_consensus = true || option; });
    // option - enable graphmap aligner, hack to avoid unused variable warning
    parsero::add_option("g", "use GraphMap aligner [flag]",
        [] (char *option) { use_graphmap_aligner = true || option; });
    // option - set read type
    parsero::add_option("x:",
        "input reads type, by default set to PacBio [pacbio, ont]",
        [] (char *option) {
            use_tech_type = read_type::string_to_read_type(option); });
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

    cout << "[INPUT] Reading draft genome: " << draft_genome_filename
        << endl;

    // read contigs from draft genome file
    StringSet<CharString> contig_ids;
    StringSet<Dna5String> contig_seqs;
    utility::read_fasta(&contig_ids, &contig_seqs, draft_genome_filename);

    // create map<contig_str_name, contig_int_id>
    unordered_map<string, uint32_t> contig_name_to_id;
    for (uint32_t id = 0; id < length(contig_ids); ++id) {
        string contig_name = utility::CharString_to_string(contig_ids[id]);
        contig_name_to_id[contig_name] = id;
    }

    cout << "[INPUT] Reading long reads: " << reads_filename << endl;

    // read long reads from file
    StringSet<CharString> read_ids;
    StringSet<Dna5String> read_seqs;
    utility::read_fasta(&read_ids, &read_seqs, reads_filename);

    // create map<read_str_name, read_int_id>
    unordered_map<string, uint32_t> read_name_to_id;
    for (uint32_t id = 0; id < length(read_ids); ++id) {
        string read_name = utility::CharString_to_string(read_ids[id]);
        read_name_to_id[read_name] = id;
    }

    // copy file to temporary folder to avoid data folder polution
    utility::write_fasta(contig_ids, contig_seqs,
                         Aligner::get_tmp_reference_filename());


    // initialize Aligner
    Aligner::init(use_graphmap_aligner, use_tech_type);
    const char *aligner_name = Aligner::get_instance().get_name().c_str();

    if (!utility::is_command_available(aligner_name)) {
        utility::exit_with_message("The %s aligner has not been detected!",
                                   aligner_name);
    }

    cout << "[ALIGNER] Initializing "<< aligner_name << " aligner..." << endl;

    // create index for all contigs in draft genome
    cout << "[ALIGNER] Creating index..." << endl;

    Aligner::get_instance().index(draft_genome_filename);

    // align all reads to the draft genome
    cout << "[ALIGNER] Aligning reads to draft genome using ";
    cout << utility::get_concurrency_level() << " threads..." << endl;

    Aligner::get_instance().align(draft_genome_filename, reads_filename);

    cout << "[ALIGNER] Creating alignments map..." << endl;
    AlignmentCollection contig_alns;
    utility::map_alignments(Aligner::get_tmp_alignment_filename(), &contig_alns,
                            contig_name_to_id);

    StringSet<Dna5String> result_contig_seqs;
    StringSet<Dna5String> extensions;
    StringSet<CharString> ext_ids;

    vector< Contig* > contigs;
    int contigs_size = length(contig_ids);

    // attempt to extend each contig
    for (int i = 0; i < contigs_size; ++i) {
        Dna5String contig_seq;
        Contig *contig = nullptr;

        cout << "[EXTENDER] Starting extension procedure for contig [" << i + 1
            << "/" << contigs_size << "]: " << contig_ids[i] << endl;

        if (use_POA_consensus) {
            contig = scaffolder::extend_contig_poa(contig_seqs[i],
                                                   contig_alns[i],
                                                   read_name_to_id);
        } else {
            contig = scaffolder::extend_contig(contig_seqs[i], contig_alns[i],
                                               read_name_to_id, read_ids,
                                               read_seqs);
        }

        cout << "\tLeft extension: " << contig->total_ext_left() << " BP"
            << endl;
        cout << "\tRight extension: " << contig->total_ext_right() << " BP"
            << endl;
        cout << "\tExtended conitg length: " << contig->total_len() << " BP"
            << endl;

        // store extended contig
        contigs.emplace_back(contig);
        contig->set_id(contig_ids[i]);
        appendValue(result_contig_seqs, contig->seq());

        // store extensions
        appendValue(ext_ids, contig->left_id());
        appendValue(extensions, contig->ext_left());

        appendValue(ext_ids, contig->right_id());
        appendValue(extensions, contig->ext_right());
    }

    cout << "[CONNECTOR] Attempting to connect extended contigs..." << endl;

    // attempt to cennect extended contigs
    Connector connector(contigs);
    connector.connect_contigs();

    // write all output files
    cout << "[OUTPUT] Writing extended contigs to file: " << result_filename
        << endl;
    utility::write_fasta(contig_ids, result_contig_seqs, result_filename);

    cout << "[OUTPUT] Writing extensions to file: " << extensions_filename
        << endl;
    utility::write_fasta(ext_ids, extensions, extensions_filename);

    cout << "[OUTPUT] Writing scaffolds to file..." << endl;
    connector.dump_scaffolds();

    // cleanup contigs
    for (auto contig : contigs) {
        delete contig;
    }

    return 0;
}
