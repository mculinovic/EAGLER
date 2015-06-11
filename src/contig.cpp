#include <string>
#include <algorithm>
#include <utility>
#include <cstring>
#include "./contig.h"

using std::string;
using std::swap;


Contig::Contig() {}


Contig::Contig(Dna5String& seq,
               int total_ext_left,
               int total_ext_right): seq_(seq) {
    total_ext_left_ = total_ext_left;
    total_ext_right_ = total_ext_right;

    // creating string extensions
    String<char, CStyle> tmp = seq_;
    string contig_str(tmp);
    ext_left_ = contig_str.substr(0, total_ext_left_);
    ext_right_ = contig_str.substr(contig_str.length() - total_ext_right_,
                                   total_ext_right_);
}

Contig::Contig(const Dna5String& contig_seq,
               string& left_extension,
               string &right_extension) {
    seq_ = left_extension;
    seq_ += contig_seq;
    seq_ += right_extension;
    ext_left_ = left_extension;
    ext_right_ = right_extension;
    total_ext_left_ = ext_left_.length();
    total_ext_right_ = ext_right_.length();
}


void Contig::set_id(const CharString& id) {
    id_ = id;

    // set left extension id
    String<char, CStyle> ltmp = id_;
    string lid(ltmp);
    lid += "L";
    left_id_ = lid;

    // set right extension id
    String<char, CStyle> rtmp = id_;
    string rid(rtmp);
    rid += "R";
    right_id_ = rid;
}

Dna5String Contig::anchor_left() {
    String<char, CStyle> tmp = seq_;
    string seq(tmp);
    Dna5String anchor = seq.substr(0, total_ext_left() + ANCHOR_LEN);
    return anchor;
}


Dna5String Contig::anchor_right() {
    String<char, CStyle> tmp = seq_;
    string seq(tmp);
    Dna5String anchor = seq.substr(seq.length() - total_ext_right() - ANCHOR_LEN);
    return anchor;
}

void Contig::reverse_complement() {
    string contig_str = utility::reverse_complement(seq_);
    seq_ = contig_str;

    swap(total_ext_left_, total_ext_right_);
    swap(left_id_, right_id_);

    ext_left_ = contig_str.substr(0, total_ext_left_);
    ext_right_ = contig_str.substr(contig_str.length() - total_ext_right_,
                                   total_ext_right_);
}
