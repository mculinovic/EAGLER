#include "aligner.h"
#include "bwa.h"
#include "graphmap.h"

const char *Aligner::tmp_alignment_filename = "./tmp/aln.sam";
const char *Aligner::tmp_reference_filename = "./tmp/reference.fasta";
const char *Aligner::tmp_contig_filename = "./tmp/contig_tmp.fasta";

Aligner *Aligner::instance = nullptr;

void Aligner::init(bool use_graphmap_aligner) {
    if (use_graphmap_aligner) {
        instance = new GraphMapAligner();
    } else {
        instance = new BwaAligner();
    }
}

Aligner& Aligner::get_instance() {
    if (instance == nullptr) {
        // throw exception
    }
    return *instance;
}

const char *Aligner::get_tmp_alignment_filename() {
    return tmp_alignment_filename;
}

const char *Aligner::get_tmp_reference_filename() {
    return tmp_reference_filename;
}

const char *Aligner::get_tmp_contig_filename() {
    return tmp_contig_filename;
}
