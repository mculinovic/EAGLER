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

void Scaffold::circular_genome_trim() {
    int n = num_contigs();
    Contig *first = contigs[0];
    Contig *last = contigs[n - 1];
    contributions[0].first = first->total_ext_left();
    contributions[n - 1].second = last->total_len() - last->total_ext_right();
}
