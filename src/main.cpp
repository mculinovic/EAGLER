// Copyright @mculinovic
#include <iostream>
#include <seqan/bam_io.h>

using seqan::BamFileIn;
using seqan::BamHeader;
using seqan::readHeader;
using seqan::BamAlignmentRecord;
using seqan::atEnd;
using seqan::FormattedFileContext;

using seqan::contigNames;
using seqan::contigLengths;

int main() {
    BamFileIn bamFileIn("../data/E-Coli/e-coli.sam");

    BamHeader header;
    readHeader(header, bamFileIn);

    typedef FormattedFileContext<BamFileIn, void>::Type TBamContext;
    TBamContext const & bamContext = context(bamFileIn);

    for (unsigned i = 0; i < length(contigNames(bamContext)); ++i)
        std::cout << contigNames(bamContext)[i] << '\t'
                  << contigLengths(bamContext)[i] << '\n';

    BamAlignmentRecord record;
    while (!atEnd(bamFileIn)) {
        readRecord(record, bamFileIn);
        std::cout << record.qName << std::endl;
        for (auto& e: record.cigar) {
            std::cout << e.count << e.operation << " ";
        }
        std::cout << std::endl;
        std::cout << record.seq << std::endl;
    }
    return 0;
}
