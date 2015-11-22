/**
 * @file connector.cpp
 * @copyright Marko Culinovic <marko.culinovic@gmail.com>
 * @copyright Luka Sterbic <luka.sterbic@gmail.com>
 * @brief Implementation file for Connector class
 * @details Implementation file for Connector class. Class provides
 * functionality for connecting extended contigs into scaffolds if they mutually
 * overlap.
 */
#include <seqan/bam_io.h>
#include <seqan/sequence.h>
#include <vector>
#include <algorithm>
#include <string>
#include <stdexcept>

#include "aligners/aligner.h"
#include "connector.h"
#include "utility.h"


using std::runtime_error;
using std::vector;
using std::max;
using std::min;
using std::string;
using std::cout;
using std::endl;

using seqan::BamHeader;
using seqan::BamAlignmentRecord;
using seqan::length;


const char* Connector::reference_file = "./tmp/connector.fasta";
const char* Connector::anchors_file = "./tmp/anchors.fasta";
const char* Connector::aln_file = "./tmp/connector_aln.sam";


Connector::Connector(const vector<Contig*>& contigs):
                    contigs_(contigs) {
    for (auto contig : contigs_) {
        string id = utility::CharString_to_string(contig->id());
        unused_contigs[id] = contig;
    }
}


Connector::~Connector() {
    curr = nullptr;

    for (auto scaffold : scaffolds) {
        delete scaffold;
    }
}


void Connector::connect_contigs() {
    cout << "\tWriting contig anchors to file..." << endl;
    Contig::dump_anchors(contigs_);

    curr = create_scaffold();
    scaffolds.emplace_back(curr);

    bool found = false;
    while (true) {
        found = connect_next();

        if (!found) {
            curr = create_scaffold();
            if (curr != nullptr) {
                scaffolds.emplace_back(curr);
            } else {
                break;
            }
        }
    }

    cout << "\tCorrecting circular genome scaffolds..." << endl;

    for (uint32_t i = 0; i < scaffolds.size(); i++) {
        cout << "\t\tExamining scaffold [" << i + 1 << "/" << scaffolds.size()
            << "]... ";

        bool did_correct = correct_circular_scaffold(scaffolds[i]);

        cout << (did_correct ? "CORRECTED" : "UNTOUCHED") << endl;
    }
}


bool Connector::connect_next() {
    Contig *curr_contig = curr->last_contig();
    string curr_contig_id = utility::CharString_to_string(curr_contig->id());

    DEBUG("Current contig: " << curr_contig->id() << endl)

    utility::write_fasta(curr_contig->id(), curr_contig->seq(), reference_file);
    Aligner::get_instance().index(reference_file);
    Aligner::get_instance().align(reference_file, anchors_file, aln_file, true);

    BamHeader header;
    vector<BamAlignmentRecord> records;
    utility::read_sam(&header, &records, aln_file);

    for (auto const& record : records) {
        DEBUG("Examining record for anchor: " << record.qName)

        if ((record.flag & UNMAPPED) || (record.flag & SECONDARY_LINE)) {
            continue;
        }

        if (used_ids_.count(utility::CharString_to_string(record.qName)) > 0) {
            continue;
        }

        string anchor_id = utility::CharString_to_string(record.qName);
        string next_id = anchor_id.substr(0, anchor_id.length() - 1);

        // if next contig is the same as current
        // do not extend with itself
        if (next_id == curr_contig_id) {
            DEBUG("Id's are same [nxt, cur]: " << next_id << curr_contig_id)
            continue;
        }

        // if next contig is already in current scaffold
        if (curr->contains(next_id)) {
            break;
        }

        // if contig inside other scaffold
        bool merge_scaffold = false;
        bool is_first = false;
        if (contig_to_scaffold.find(next_id) != contig_to_scaffold.end()) {
            Scaffold *next_scaffold = contig_to_scaffold[next_id];

            is_first = *(next_scaffold->first_contig()) != *(curr_contig);

            if (!is_first) {
                continue;
            }

            merge_scaffold = true;
        }

        // ovo je mozda problematicno - provjeriti!
        if (!should_connect(curr_contig, record)) {
            continue;
        }

        DEBUG("Attempting merge for anchor: " << record.qName)

        Contig *next = find_contig(next_id);

        cout << "\t\tConnecting contig: " << next->id() << endl;

        int merge_start = max(curr_contig->right_ext_pos(), record.beginPos);
        bool reverse_complement = record.flag & COMPLEMENT;

        if (next == nullptr) {
            utility::throw_exception<runtime_error>("Contig invalid id");
        }

        if (reverse_complement) {
            next->reverse_complement();
        }

        int right_ext_len = curr_contig->total_len() - merge_start;
        int next_start = min(right_ext_len, next->total_ext_left());
        int merge_end = next_start + record.beginPos;

        int merge_len = merge_end - merge_start;

        curr->add_contig(next, merge_start + merge_len / 2,
                         next_start - merge_len / 2);

        contig_to_scaffold[next_id] = curr;

        used_ids_.insert(anchor_id);
        used_ids_.insert(utility::CharString_to_string(
            curr_contig->right_id()));

        if (merge_scaffold) {
            Scaffold *next_scaffold = contig_to_scaffold[next_id];
            curr->merge(next_scaffold);

            // refresh contig to scaffold mapping
            for (auto contig : next_scaffold->get_contigs()) {
                string tmp_id = utility::CharString_to_string(contig->id());
                contig_to_scaffold[tmp_id] = curr;
            }

            // remove scaffold from list
            for (size_t i = 0; i < scaffolds.size(); ++i) {
                if (scaffolds[i] == next_scaffold) {
                    scaffolds.erase(scaffolds.begin() + i);
                    break;
                }
            }
            delete next_scaffold;
        }

        if (!merge_scaffold) {
            unused_contigs.erase(next_id);
        }

        return true;
    }

    return false;
}


bool Connector::should_connect(Contig *contig,
                               const BamAlignmentRecord& record) {
    // iterate over cigar string to get lengths of
    // read and contig parts used in alignment
    int cigar_len = length(record.cigar);
    int used_read_size = 0;
    int used_contig_size = 0;

    // if the anchor is not soft clipped skip it
    if (record.cigar[cigar_len - 1].operation != 'S') {
        return false;
    }

    for (auto const& e : record.cigar) {
        if (utility::contributes_to_seq_len(e.operation)) {
            used_read_size += e.count;
        }
        if (utility::contributes_to_contig_len(e.operation)) {
            used_contig_size += e.count;
        }
    }

    int right_clipping_len = record.cigar[cigar_len - 1].count;
    used_read_size -= right_clipping_len;
    int len = right_clipping_len -
              (contig->total_len() - (record.beginPos + used_contig_size));

    // if read doesn't extend right of contig skip it
    return len > ANCHOR_THRESHOLD * ANCHOR_LEN;
}


Contig* Connector::find_contig(const string& id) {
    CharString contig_id = id;

    for (auto contig : contigs_) {
        if (contig->id() == contig_id) {
            return contig;
        }
    }

    return nullptr;
}


Scaffold* Connector::create_scaffold() {
    if (unused_contigs.empty()) {
        return nullptr;
    }

    auto it = unused_contigs.begin();
    Scaffold *scaffold = new Scaffold(it->second);

    cout << "\tCreated scaffold with base contig: " << it->second->id() << endl;

    contig_to_scaffold[it->first] = scaffold;

    unused_contigs.erase(it);
    return scaffold;
}


void Connector::dump_scaffolds() {
    StringSet<Dna5String> scaffold_seqs;
    StringSet<CharString> ids;

    int i = 0;
    for (auto scaffold : scaffolds) {
        appendValue(ids, utility::create_seq_id("scaffold|%d", i));
        appendValue(scaffold_seqs, scaffold->get_combined_sequence());
        ++i;
    }

    utility::write_fasta(ids, scaffold_seqs, "./tmp/genome.fasta");
}


bool Connector::correct_circular_scaffold(Scaffold *scaffold) {
    Contig *last_contig = scaffold->last_contig();
    string contig_id = utility::CharString_to_string(last_contig->id());

    utility::write_fasta(last_contig->id(), last_contig->seq(), reference_file);

    Aligner::get_instance().index(reference_file);
    Aligner::get_instance().align(reference_file, anchors_file,
                                  aln_file, false);

    Contig *first_contig = scaffold->first_contig();
    string left_id = utility::CharString_to_string(first_contig->left_id());

    BamHeader header;
    vector<BamAlignmentRecord> records;
    utility::read_sam(&header, &records, aln_file);

    for (auto const& record : records) {
        if ((record.flag & UNMAPPED) || (record.flag & SECONDARY_LINE)) {
            continue;
        }

        string anchor_id = utility::CharString_to_string(record.qName);

        if (anchor_id == left_id && should_connect(last_contig, record)) {
            int trim_right_idx = max(last_contig->right_ext_pos(),
                                     record.beginPos);

            int trim_left_idx = 0;
            if (record.beginPos < last_contig->right_ext_pos()) {
                trim_left_idx = last_contig->right_ext_pos() - record.beginPos;
            }

            scaffold->trim(trim_left_idx, trim_right_idx);
            return true;
        }
    }

    return false;
}
