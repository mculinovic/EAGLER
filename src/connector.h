/**
 * @file connector.h
 * @copyright Marko Culinovic <marko.culinovic@gmail.com>
 * @copyright Luka Sterbic <luka.sterbic@gmail.com>
 * @brief Header file for Connector class
 * @details Header file for Connector class. Class provides functionality
 * for connecting extended contigs into scaffolds if they mutually overlap.
 */
#ifndef CONNECTOR_H
#define CONNECTOR_H

#include <seqan/bam_io.h>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <string>
#include "./contig.h"
#include "./scaffold.h"

using std::vector;
using std::string;
using std::unordered_set;
using std::unordered_map;
using seqan::BamAlignmentRecord;


/**
 * @brief SAM file flag denoting a non-primary alignment
 */
#define SECONDARY_LINE 0x900

/**
 * @brief SAM file flag denoting an unmapped query
 */
#define UNMAPPED 0x4

/**
 * @brief [brief description] SAM file flag denoting a query from the complement
 * strand of the reference genome
 */
#define COMPLEMENT 0x10

/**
 * @brief Minimum percentage of anchor sequence that must extend the end of the
 * current scaffold to extend it
 */
#define ANCHOR_THRESHOLD 0.66


/**
 * @brief Connector class
 * @details Class provides functionality
 * for connecting extended contigs into scaffolds if they mutually overlap.
 */
class Connector {
 public:
    /**
     * @brief Connector class constructor
     * @details Constructor initializes Connector class with
     * contigs which are afterwards merged into scaffolds.
     *
     * @param contigs Vector of contigs.
     */
    explicit Connector(const vector<Contig*>& contigs);


    /**
     * @brief Connector class destructor
     */
    ~Connector();


    /**
     * @brief Method connects contigs passed as input
     * into scaffolds if this is possible.
     * @details Anchors are created from each contigs.
     * For each contig, overlap between anchors and this contig
     * are found (if any). This contig is afterwards merged
     * with contig that anchor belongs to. Procedure is repeated
     * until all contigs are processed.
     */
    void connect_contigs();


    /**
     * @brief Getter for scaffolds.
     * @return Vector of scaffolds.
     */
    const vector<Scaffold*>& get_scaffolds() { return scaffolds; }


    /**
     * @brief Output scaffolds to file.
     * @details Scaffolds are outputed in FASTA format file
     * genome.fasta
     *
     * @param output_file path to the output file
     */
    void dump_scaffolds(const char *output_file);


 private:
    /**
     * Contig as reference filename during extension
     * process for bwa tool.
     */
    static const char* tmp_reference_file;

    /**
     * Filename for storing anchor sequences.
     */
    static const char* tmp_anchors_file;


    /**
     * Filename for storing alignment results
     * obtained by bwa tool.
     */
    static const char* tmp_alignment_file;


    /**
     * Vector of contigs
     */
    const vector<Contig*>& contigs_;

    /**
     * @brief IDs set of used contigs in process.
     */
    unordered_set<string> used_ids_;

    /**
     * @brief Map from ID to Contig. Unused contigs
     * are only stored in this map.
     */
    unordered_map<string, Contig*> unused_contigs;


    /**
     * @brief Current scaffold to merge next contigs into.
     */
    Scaffold* curr;

    /**
     * @brief Vector of scaffolds. Empty before connect_contigs
     * function is called.
     */
    vector<Scaffold*> scaffolds;

    /**
     * @brief Map from ContigID to Scaffold. So that information
     * about which contig belongs to which scaffold is memorized.
     */
    unordered_map<string, Scaffold*> contig_to_scaffold;

    /**
     * @brief Creates new Scaffold from next unused Contig
     * @details New Scaffold object is created if there is
     * at least one unused contig in previously created scaffolds.
     * @return New scaffold object or nullptr if all contigs are used.
     */
    Scaffold* create_scaffold();


    /**
     * @brief Method finds Contig with id given as parameter.
     *
     * @param id Contig id.
     * @return Contig object if exists, nullptr otherwise.
     */
    Contig* find_contig(const string& id);


    /**
     * @brief Method merges existing scaffold with next contig.
     * @details Merthod checks if there exists any overlap between
     * contig anchors and last contig in scaffold. If overlap is found
     * contig that anchor belongs to is merged into scaffold.
     * @return True if some contig is found and merged into scaffold,
     * false otherwise.
     */
    bool connect_next();


    /**
     * @brief Method checks if contig should be connected with
     * contig represented by record in alignment file.
     * @details Check is done by iterating over cigar string
     * and calculating lengths of used part of sequences in
     * reference and in record. If record extends right of
     * contig end, overlap is found, and contigs should be
     * merged.
     *
     * @param contig Contig used as reference in bwa.
     * @param record Record representing anchor of another contig.
     *
     * @return True if contigs should be connected, false otherwise.
     */
    bool should_connect(Contig *contig, const BamAlignmentRecord& record);


    /**
     * @brief Method trims scaffold at its ends if scaffold is circular.
     *
     * @param scaffold Scaffold to check for cicularity and correct
     * if neccessary.
     * @return true if the scaffold has been corrected, false otherwise
     */
    bool correct_circular_scaffold(Scaffold *scaffold);
};


#endif  // CONNECTOR_H
