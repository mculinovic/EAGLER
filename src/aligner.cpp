#include "aligner.h"
#include "bwa.h"

const char *Aligner::tmp_alignment_filename = "./tmp/aln.sam";
const char *Aligner::tmp_reference_filename = "./tmp/reference.fasta";
const char *Aligner::tmp_contig_filename = "./tmp/contig_tmp.fasta";
Aligner *Aligner::sInstance;

void Aligner::init(bool use_graphmap_aligner) {
    if (use_graphmap_aligner) {
        //
    } else {
        sInstance = new BwaAligner();
    }
}

Aligner& Aligner::get_instance() {
    if (sInstance == nullptr) {
        // throw exception
    }
    return *sInstance;
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
