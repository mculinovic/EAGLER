#ifndef CONNECTOR_H
#define CONNECTOR_H

#include <seqan/bam_io.h>
#include <vector>
#include <unordered_set>
#include <string>
#include "./contig.h"
#include "./scaffold.h"

using std::vector;
using std::string;
using std::unordered_set;
using seqan::BamAlignmentRecord;

class Connector {
 public:
    Connector(const vector<Contig*>& contigs);
    ~Connector();
    void connect_contigs();
    const vector<Scaffold*>& get_scaffolds() { return scaffolds; }

 private:
    bool connect_next();
    bool should_connect(Contig *contig, const BamAlignmentRecord& record);
    Contig* find_contig(const string& id);

    static const char* reference_file;
    static const char* anchors_file;

    const vector<Contig*>& contigs_;
    unordered_set<string> used_ids_;
    Scaffold* curr;
    vector<Scaffold*> scaffolds;
};


#endif  // CONNECTOR_H
