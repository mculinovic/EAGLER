#include <seqan/sequence.h>
#include <string>

using std::string;

using seqan::toCString;
using seqan::String;
using seqan::CStyle;
using seqan::length;
using seqan::Dna5String;

class Contig {
 public:
    Contig(Dna5String& seq, int total_ext_left, int total_ext_right);
    Contig(const Dna5String& contig_seq, string& left_extension, string &right_extension);

    Dna5String& seq() { return seq_; }
    int total_ext_left() { return total_ext_left_; }
    int total_ext_right() { return total_ext_right_; }
    string& ext_left() { return ext_left_; }
    string& ext_right() { return ext_right_; }

 private:
    Dna5String seq_;
    int total_ext_left_;
    int total_ext_right_;
    string ext_left_;
    string ext_right_;
};
