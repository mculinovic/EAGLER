#include "graphmap.h"
#include "utility.h"

void GraphMapAligner::index(const char* filename) {
	// utility::execute_command("graphmap -v 0 -I -r %s", filename);
}

void GraphMapAligner::align(const char* reference_file,
                            const char* reads_file) {
    utility::execute_command(
        "graphmap -v 0 -t %d -r %s -d %s -o %s",
        utility::get_concurrency_level(),
        reference_file,
        reads_file,
        get_tmp_alignment_filename()
    );
}

void GraphMapAligner::align(const char* reference_file,
                            const char* reads_file,
                            const char* sam_file,
                            bool only_primary) {
    // TODO handle only primary
    utility::execute_command(
        "graphmap -v 0 -t %d %s -r %s -d %s -o %s",
        only_primary ? "" : "-Z",
        utility::get_concurrency_level(),
        reference_file,
        reads_file,
        sam_file
    );
}

void GraphMapAligner::align(const char* reference_file,
                            const char* reads_file,
                            const char* sam_file) {
    // TODO handle only primary
    utility::execute_command(
        "graphmap -v 0 -t %d -r %s -d %s -o %s",
        utility::get_concurrency_level(),
        reference_file,
        reads_file,
        sam_file
    );
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
