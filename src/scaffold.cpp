#include <string>

#include "./scaffold.h"
#include "./utility.h"

Scaffold::Scaffold(Contig *first_contig) {
    contigs.emplace_back(first_contig);
    contig_ids.insert(utility::CharString_to_string(first_contig->id()));
    contributions.emplace_back(0, first_contig->total_len());
}

void Scaffold::add_contig(Contig *contig, int last_end, int this_start) {
    contigs.emplace_back(contig);
    contig_ids.insert(utility::CharString_to_string(contig->id()));
    contributions[contributions.size() - 1].second = last_end;
    contributions.emplace_back(this_start, contig->total_len());
}

Dna5String Scaffold::get_combined_sequence() {
    string seq;
    for (size_t i = 0; i < contigs.size(); ++i) {
        auto &contig = contigs[i];
        string contig_seq = utility::Dna5String_to_string(contig->seq());
        seq += contig_seq.substr(contributions[i].first,
                         contributions[i].second - contributions[i].first);
    }
    Dna5String dna_string = seq;
    return dna_string;
}


void Scaffold::merge(Scaffold *scaffold) {
    int size = num_contigs();

    contributions[size - 1].second = scaffold->contributions[0].second;
    contributions.insert(contributions.end(),
                         scaffold->contributions.begin() + 1,
                         scaffold->contributions.end());

    contigs.insert(contigs.end(),
                   scaffold->contigs.begin() + 1,
                   scaffold->contigs.end());

    contig_ids.insert(scaffold->contig_ids.begin(),
                      scaffold->contig_ids.end());
}


void Scaffold::trim(int left_start_pos, int right_end_pos) {
    contributions[0].first = left_start_pos;
    contributions[contributions.size() - 1].second = right_end_pos;
}
