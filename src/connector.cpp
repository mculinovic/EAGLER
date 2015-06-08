#include <vector>
#include "./connector.h"
#include "./utility.h"
#include "./bwa.h"

using std::vector;

const char* Connector::reference_file = "./tmp/connector.fasta";
const char* Connector::anchors_file = "./tmp/anchors.fasta";

Connector::Connector(const vector<Contig*>& contigs):
                    contigs_(contigs) {}


void Connector::connect_contigs() {
    Contig::dump_anchors(contigs_);

    curr = *contigs_[0];
    used_ids_.insert(utility::CharString_to_string(curr.id()));
    do {
        connect_next();
    } while (used_ids_.size() < contigs_.size());
}

void Connector::connect_next() {
    utility::write_fasta(curr.id(), curr.seq(), reference_file);
    aligner::bwa_index(reference_file);
    aligner::bwa_mem(reference_file, anchors_file);
    // TODO(mculinovic, lsterbic) read .sam file and extend contigs
}
