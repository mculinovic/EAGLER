#ifndef SCAFFOLD_H
#define SCAFFOLD_H

#include <vector>
#include <unordered_set>
#include <utility>

#include "./contig.h"

using std::vector;
using std::unordered_set;
using std::pair;

class Scaffold {
 public:
    Scaffold(Contig *first_contig);
    void add_contig(Contig *contig, int last_end, int this_start);
    Contig* last_contig() { return contigs[contigs.size() - 1]; }
    bool contains(const string& id) { return contig_ids.count(id) > 0; }
    Dna5String get_combined_sequence();
    void trim_ends();
    int num_contigs() { return contigs.size(); }

 private:
    vector<Contig *> contigs;
    vector<pair<int, int>> contributions;
    unordered_set<string> contig_ids;
};


#endif  // SCAFFOLD_H
