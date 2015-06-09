#ifndef CONTIG_H
#define CONTIG_H

#include <seqan/sequence.h>
#include <string>
#include <vector>
#include "./utility.h"

using std::string;
using std::vector;

using seqan::toCString;
using seqan::String;
using seqan::CStyle;
using seqan::length;
using seqan::Dna5String;
using seqan::CharString;
using seqan::StringSet;


#define ANCHOR_LEN 10000


class Contig {
 public:
    Contig();
    Contig(Dna5String& seq, int total_ext_left, int total_ext_right);
    Contig(const Dna5String& contig_seq, string& left_extension, string &right_extension);

    Dna5String& seq() { return seq_; }
    CharString& id() { return id_; }
    int total_len() { return length(seq_); }

    int total_ext_left() { return total_ext_left_; }
    int total_ext_right() { return total_ext_right_; }
    int right_ext_pos() { return total_len() - total_ext_right(); }

    string& ext_left() { return ext_left_; }
    string& ext_right() { return ext_right_; }

    CharString& left_id() { return left_id_; }
    CharString& right_id() { return right_id_; }

    bool is_reversed() { return need_reversal_; }
    bool is_complemented() { return need_complement_; }

    void set_id(const CharString& id);
    void reverse() { need_reversal_ = !need_reversal_; }
    void complement() { need_reversal_ = !need_reversal_; }

    static void dump_anchors(const vector< Contig *>& contigs) {
        StringSet<Dna5String> anchors;
        StringSet<CharString> ids;

        for (auto contig : contigs) {
            appendValue(ids, contig->left_id());
            appendValue(anchors, contig->anchor_left());
            appendValue(ids, contig->right_id());
            appendValue(anchors, contig->anchor_right());
        }

        utility::write_fasta(ids, anchors, "./tmp/anchors.fasta");
    }

 private:
    Dna5String anchor_left();
    Dna5String anchor_right();

    CharString id_;
    Dna5String seq_;

    CharString left_id_;
    int total_ext_left_;
    string ext_left_;

    CharString right_id_;
    int total_ext_right_;
    string ext_right_;

    bool need_reversal_;
    bool need_complement_;
};


#endif  // CONTIG_H
