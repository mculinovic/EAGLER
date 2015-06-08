#ifndef CONNECTOR_H
#define CONNECTOR_H

#include <vector>
#include <unordered_set>
#include "./contig.h"

using std::vector;
using std::unordered_set;

class Connector {
 public:
    Connector(const vector<Contig*>& contigs);
    void connect_contigs();

 private:
    void connect_next();

    static const char* reference_file;
    static const char* anchors_file;

    const vector<Contig*>& contigs_;
    unordered_set<string> used_ids_;
    Contig curr;
};


#endif  // CONNECTOR_H
