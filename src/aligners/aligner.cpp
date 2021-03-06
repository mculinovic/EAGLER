/**
 * @file aligner.cpp
 * @copyright Marko Culinovic <marko.culinovic@gmail.com>
 * @copyright Luka Sterbic <luka.sterbic@gmail.com>
 * @brief Implementation file for the abstract Aligner class.
 */

#include "aligners/bwa.h"
#include "aligners/graphmap.h"
#include "aligner.h"
#include "utility.h"


using std::runtime_error;


const char *Aligner::tmp_alignment_filename = "./tmp/aln.sam";
const char *Aligner::tmp_reference_filename = "./tmp/reference.fasta";
const char *Aligner::tmp_contig_filename = "./tmp/contig_tmp.fasta";

Aligner *Aligner::instance = nullptr;


read_type::ReadType read_type::string_to_read_type(const char *tech_type) {
    string s_tech_type(tech_type);

    if (s_tech_type == "pacbio") {
        return PacBio;
    } else if (s_tech_type == "ont") {
        return ONT;
    } else {
        utility::exit_with_message("Unknown read type.");
        // silence compiler warning for no return value
        return PacBio;
    }
}


void Aligner::init(bool use_graphmap_aligner, read_type::ReadType tech_type) {
    if (instance != nullptr) {
        utility::throw_exception<runtime_error>(
            "The init method should not be called more than once.");
    }

    if (use_graphmap_aligner) {
        instance = new GraphMapAligner(tech_type);
    } else {
        instance = new BwaAligner(tech_type);
    }
}


Aligner& Aligner::get_instance() {
    if (instance == nullptr) {
        utility::throw_exception<runtime_error>(
            "The init method should be called before get_instance.");
    }

    return *instance;
}


const char *Aligner::get_tmp_alignment_filename() {
    return tmp_alignment_filename;
}


const char *Aligner::get_tmp_reference_filename() {
    return tmp_reference_filename;
}


const char *Aligner::get_tmp_contig_filename() {
    return tmp_contig_filename;
}


const std::string& Aligner::get_name() const {
    return this->name;
}
