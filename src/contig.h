/**
 * @file contig.h
 * @author Marko Culinovic <marko.culinovic@gmail.com>
 * @author Luka Sterbic <luka.sterbic@gmail.com>
 * @brief Header file for Contig class
 * @details Header file for Contig class. Contig class is wrapper for
 * original contigs. It provides functionality for accessing extensions
 * of contig found in extension process. It also provides functionality
 * for accessing contig anchors - subsequences of contig used in contigs
 * merge process.
 */
#ifndef CONTIG_H
#define CONTIG_H

#include <seqan/sequence.h>
#include <string>
#include <vector>

#include "utility.h"


using std::string;
using std::vector;

using seqan::toCString;
using seqan::String;
using seqan::CStyle;
using seqan::length;
using seqan::Dna5String;
using seqan::CharString;
using seqan::StringSet;

// length of subsequence in original contig sequence to be
// inserted in anchor
#define ANCHOR_LEN 10000

/**
 * @brief Contig class represent contig sequences.
 * @details Contig class is wrapper for
 * original contigs. It provides functionality for accessing extensions
 * of contig found in extension process. It also provides functionality
 * for accessing contig anchors - subsequences of contig used in contigs
 * merge process.
 */
class Contig {
 public:
    /**
     * @brief Contig class default constructor.
     */
    Contig();


    /**
     * @brief Contig class constructor
     * @details Constructor used for creating Contig object
     * from already extended contig and lengths of right and
     * left extensions.
     *
     * @param seq Extended contig sequence.
     * @param total_ext_left Length of left extension.
     * @param total_ext_right Length of right extension.
     */
    Contig(Dna5String& seq, int total_ext_left, int total_ext_right);


    /**
     * @brief Contig class constructor
     * @details Constructor used for creating Contig object
     * from original contig sequence and sequences of left
     * and right extensions.
     *
     * @param contig_seq Original contig sequence.
     * @param left_extension Left extension sequence.
     * @param right_extension Right extension sequence.
     */
    Contig(const Dna5String& contig_seq, string& left_extension, string &right_extension);


    /**
     * @brief Getter for extended contig sequence.
     * @return Contig extended sequence.
     */
    Dna5String& seq() { return seq_; }


    /**
     * @brief Getter for contig id.
     * @return Contig id.
     */
    CharString& id() { return id_; }


    /**
     * @brief Getter for total length of extended contig sequence.
     * @return Length of extended contig sequence.
     */
    int total_len() { return length(seq_); }


    /**
     * @brief Getter for length of left extension.
     * @return Length of left extension
     */
    int total_ext_left() { return total_ext_left_; }


    /**
     * @brief Getter for length of right extension.
     * @return Length of right extension.
     */
    int total_ext_right() { return total_ext_right_; }


    /**
     * @brief Getter for right extension start position in
     * extended contig sequence.
     * @return Start index of right extension in extended
     * contig sequence.
     */
    int right_ext_pos() { return total_len() - total_ext_right(); }


    /**
     * @brief Getter for left extension.
     * @return Contig left extension.
     */
    string& ext_left() { return ext_left_; }


    /**
     * @brief Getter for right extension.
     * @return Contig right extension.
     */
    string& ext_right() { return ext_right_; }


    /**
     * @brief Getter for left extension id.
     * @return Left extension id.
     */
    CharString& left_id() { return left_id_; }


    /**
     * @brief Getter for right extension id.
     * @return Right extension id.
     */
    CharString& right_id() { return right_id_; }


    /**
     * @brief Setter for contig id.
     *
     * @param id Contig id.
     */
    void set_id(const CharString& id);


    /**
     * @brief Method reverse complements this contig.
     * @details Extended sequence is reversed and complemented.
     * Left and right extensions data is updated from this
     * reverse complemented sequence.
     */
    void reverse_complement();


    /**
     * @brief Method creates left and right anchor from contigs
     * and dumps them to file.
     *
     * @param contigs contigs used for creation of anchors
     * @param anchors_file output file to dump the anchors
     */
    static void dump_anchors(const vector< Contig *>& contigs,
                             const char *anchors_file) {
        StringSet<Dna5String> anchors;
        StringSet<CharString> ids;

        for (auto contig : contigs) {
            appendValue(ids, contig->left_id());
            appendValue(anchors, contig->anchor_left());
            appendValue(ids, contig->right_id());
            appendValue(anchors, contig->anchor_right());
        }

        utility::write_fasta(ids, anchors, anchors_file);
    }


    /**
     * @brief Operator == overload
     * @return True if contigs are equal, false otherwise.
     */
    bool operator==(const Contig& other) const {
        return strcmp(toCString(id_), toCString(other.id_)) == 0;
    }


    /**
     * @brief Operator != overload
     * @return True if contigs are different, false otherwise.
     */
    bool operator!=(const Contig& other) const {
        return !(*this == other);
    }

 private:
    /**
     * @brief Creates left anchor of contig.
     * @details Left anchor is created so that subsequence[0, ANCHOR_LEN]
     * of original contig sequence is appended on left extension.
     * @return Left contig anchor.
     */
    Dna5String anchor_left();


    /**
     * @brief Creates right anchor of contig.
     * @details Right extension is created so that right extension
     * is appended on subsequence[LEN - ANCHOR_LEN, LEN] of original
     * contig sequence.
     * @return Right contig anchor.
     */
    Dna5String anchor_right();

    // Contig id
    CharString id_;
    // Extended contig sequence
    Dna5String seq_;

    // Left extension id
    CharString left_id_;
    // Left extension length
    int total_ext_left_;
    // Left extension sequence
    string ext_left_;

    // Right extension id
    CharString right_id_;
    // Right extension length
    int total_ext_right_;
    // Right extension sequence.
    string ext_right_;
};


#endif  // CONTIG_H
