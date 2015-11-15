/**
 * @file graphmap.cpp
 * @copyright Marko Culinovic <marko.culinovic@gmail.com>
 * @copyright Luka Sterbic <luka.sterbic@gmail.com>
 */

#include "graphmap.h"
#include "utility.h"


void GraphMapAligner::index(const char* filename) {
    utility::execute_command("graphmap -v 0 -I -r %th", filename);
}


void GraphMapAligner::align(const char* reference_file,
                            const char* reads_file) {
    align(reference_file, reads_file, get_tmp_alignment_filename(), false);
}


void GraphMapAligner::align(const char* reference_file,
                            const char* reads_file,
                            const char* sam_file,
                            bool only_primary) {
    utility::execute_command(
        "graphmap -v 0 -t %d %s -a anchor -r %th -d %th -o %th",
        utility::get_concurrency_level(),
        only_primary ? "" : "-Z",
        reference_file,
        reads_file,
        sam_file);
}


void GraphMapAligner::align(const char* reference_file,
                            const char* reads_file,
                            const char* sam_file) {
    align(reference_file, reads_file, sam_file, false);
}


void GraphMapAligner::align(const CharString& id,
                            const Dna5String& contig,
                            const char* reads_filename) {
    // write contig to temporary .fasta file
    utility::write_fasta(id, contig, get_tmp_contig_filename());

    // create index for contig
    index(get_tmp_contig_filename());

    // align reads to conting
    align(get_tmp_contig_filename(), reads_filename);
}
