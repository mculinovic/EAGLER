/**
 * @file scaffold.h
 * @author Marko Culinovic <marko.culinovic@gmail.com>
 * @author Luka Sterbic <luka.sterbic@gmail.com>
 * @brief Header file for Scaffold class
 * @details Header file for Scaffold class. It is used as representation
 * of scaffolds. Scaffold consists of multiple Contigs and
 * memorizes contributions of each contig into scaffold sequence.
 */
#ifndef SCAFFOLD_H
#define SCAFFOLD_H

#include <vector>
#include <unordered_set>
#include <utility>
#include <string>

#include "./contig.h"


using std::vector;
using std::unordered_set;
using std::pair;
using std::string;

/**
 * @brief Scaffold class is representation of scaffold structure.
 * @details Scaffold class is used as representation
 * of scaffolds. Scaffold consists of multiple Contigs and
 * memorizes contributions of each contig into scaffold sequence.
 */
class Scaffold {
 public:
    /**
     * @brief Scaffold constructor.
     * @details Constructor creates scaffold from initial
     * contig.
     * 
     * @param first_contig First contig scaffold consists of.
     */
    explicit Scaffold(Contig *first_contig);


    /**
     * @brief Method adds next contig to scaffold structure.
     * @details Contig is added as last contig in a scaffold.
     * Because of overlaps between previously last contig and
     * the one being added contributions of both contigs to
     * scaffold sequence have to be updated.
     * 
     * @param contig Contig to be added.
     * @param last_end End contribution index of previously
     * last contig.
     * @param this_start Start contribution index of contig
     * being added.
     */
    void add_contig(Contig *contig, int last_end, int this_start);


    /**
     * @brief Checks if scaffold contains contig with given Id.
     * 
     * @param id Contig Id.
     * @return True if contains, false otherwise.
     */
    bool contains(const string& id) { return contig_ids.count(id) > 0; }


    /**
     * @brief Method creates unified scaffold sequence.
     * @details One scaffold sequence is generated from
     * multiple contig sequences. Depending on contributions
     * of each contig, subsequences from contig sequences
     * are extracted and connected with each other. 
     * @return Scaffold sequence.
     */
    Dna5String get_combined_sequence();


    /**
     * @brief Getter for first contig in scaffold.
     * @return First contig in scaffold.
     */
    Contig* first_contig() { return contigs[0]; }


    /**
     * @brief Getter for last contig in scaffold.
     * @return Last contig in scaffold.
     */
    Contig* last_contig() { return contigs[contigs.size() - 1]; }
    

    /**
     * @brief Getter for all contigs in scaffold.
     * @return Contigs scaffold consists of.
     */
    const vector<Contig*>& get_contigs() { return contigs; }


    /**
     * @brief Getter for number of contigs in scaffold.
     * @return Number of contigs in scaffold.
     */
    int num_contigs() { return contigs.size(); }


    /**
     * @brief Merge this scaffold with scaffold given as parameter.
     * @details All contigs from scaffold given as parameter
     * are inserted at the end of vector of contigs of this
     * scaffold. Also, contigs contributions are updated.
     * 
     * @param scaffold Scaffold to merge into this one.
     */
    void merge(Scaffold *scaffold);


    /**
     * @brief Trim scaffold if it is circular.
     * @details Only contributions of first and last contig
     * in scaffold have to be modified accordingly to
     * parameters passed.
     * 
     * @param left_start_pos Start index in sequence of
     * first contig.
     * @param right_end_pos End index in sequence of
     * last contig.
     */
    void trim(int left_start_pos, int right_end_pos);

 private:
    // contigs scaffold consists of
    vector<Contig *> contigs;


    // Start and end indices of each contig sequences.
    // Only subsequences in this interval contribute
    // to scaffold sequence.
    vector<pair<int, int>> contributions;


    // Ids of contigs scaffold consits of.
    unordered_set<string> contig_ids;
};


#endif  // SCAFFOLD_H
