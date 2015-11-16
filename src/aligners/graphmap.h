/**
 * @file graphmap.cpp
 * @copyright Marko Culinovic <marko.culinovic@fer.hr>
 * @copyright Luka Sterbic <luka.sterbic@fer.hr>
 * @brief Declaration of the GraphMap class.
 */

#ifndef ALIGNER_GRAPHMAP_H
#define ALIGNER_GRAPHMAP_H

#include <seqan/sequence.h>

#include "aligner.h"

using seqan::CharString;
using seqan::Dna5String;


class GraphMapAligner : public Aligner {
 public:
    explicit GraphMapAligner(read_type::ReadType tech_type)
        : Aligner("graphmap", tech_type) {}
    virtual ~GraphMapAligner() = default;
    virtual void index(const char* filename);
    virtual void align(const char* reference_file,
                       const char* reads_file);
    virtual void align(const char* reference_file,
                       const char* reads_file,
                       const char* sam_file,
                       bool only_primary);
    virtual void align(const char* reference_file,
                       const char* reads_file,
                       const char* sam_file);
    virtual void align(const CharString& id,
                       const Dna5String& contig,
                       const char* reads_filename);
};

#endif  // ALIGNER_GRAPHMAP_H
