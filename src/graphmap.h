#ifndef GRAPHMAP_H
#define GRAPHMAP_H

#include <seqan/sequence.h>

#include "aligner.h"

using seqan::CharString;
using seqan::Dna5String;


class GraphMapAligner: public Aligner {
public:
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

#endif // GRAPHMAP_H
