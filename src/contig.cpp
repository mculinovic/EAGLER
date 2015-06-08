#include <string>
#include "./contig.h"

using std::string;

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

