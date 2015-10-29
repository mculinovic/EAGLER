#ifndef ALIGNER_H
#define ALIGNER_H

#include <seqan/sequence.h>


using seqan::CharString;
using seqan::Dna5String;


class Aligner {
private:
    static const char *tmp_alignment_filename;
    static const char *tmp_reference_filename;
    static const char *tmp_contig_filename;
public:
	virtual ~Aligner() = default;
	virtual void index(const char* filename) = 0;
	virtual void align(const char* reference_file,
	                   const char* reads_file) = 0;
	virtual void align(const char* reference_file,
	                   const char* reads_file,
	                   const char* sam_file,
	                   bool only_primary) = 0;
	virtual void align(const char* reference_file,
	                   const char* reads_file,
	                   const char* sam_file) = 0;
	virtual void align(const CharString& id,
	                   const Dna5String& contig,
	                   const char* reads_filename) = 0;

    static const char *get_tmp_alignment_filename();
    static const char *get_tmp_reference_filename();
    static const char *get_tmp_contig_filename();
};

#endif // ALIGNER_H
