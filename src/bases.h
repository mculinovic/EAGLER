/**
 * @file bases.h
 * @copyright Coypright 2015 Marko Culinovic, Luka Sterbic
 * @author Marko Culinovic <marko.culinovic@gmail.com>
 * @author Luka Sterbic <luka.sterbic@gmail.com>
 * @brief Header file for BasesCounter class and various functions which
 * calculate statistics over a sequence of bases.
 * @details Header file for BasesCounter class and supporting functions. It is
 * used for calculating various statistics at specific positions in the contig
 * extension process.
 */
#ifndef BASES_H
#define BASES_H

#include <string>
#include <memory>
#include <vector>
#include <utility>
#include <functional>

#include "./extension.h"


using std::string;
using std::vector;
using std::pair;
using std::shared_ptr;


/**
 * @brief The number of different bases in a DNA sequence
 */
#define NUM_BASES 4


/**
 * @brief Type used to represent a boolean function over a chara cter
 */
typedef std::function<bool(char)> bool_predicate;


/**
 * @brief Namespace for BasesCounter class and supporting functions
 */
namespace bases {

/**
 * @brief Used for calculating various statistics at specific positions in the
 * contig extension process.
 * @details The BasesCounter class is used for storing the frequency of
 * appearances of each base at a specific position. The class also calculates
 * the coverage and the most frequent base at the given position.
 */
class BasesCounter {
 public:
    /**
     * @brief array used to store the number of appearances of each base
     */
    uint32_t count[NUM_BASES];

    /**
     * @brief number of bases at this position, sum(count)
     */
    uint32_t coverage;

    /**
     * @brief index in the count array of base with most appearances
     */
    uint32_t max_idx;

    /**
     * @brief BasesCounter class default constructor
     */
    BasesCounter();

    /**
     * @brief Proccesing a single base
     * @details Converts the given base to it's corresponding index in the count
     * array and increments the count at that index by one.
     *
     * @param base character base in a DNA sequence
     */
    void digest_base(char base);

    /**
     * @brief Refresh all member variables
     * @details The count array is used to compute the coverage and the most
     * frequent base.
     */
    void refresh_stats();
};


/**
 * @brief Count bases at a specific position in the given extensions.
 * @details Creates a BasesCounter object and digests one base from each
 * extending sequence.
 *
 * @param extensions Extension sequences, each extension stores the index of the
 * first unprocessed index in it's sequence
 * @param is_read_eligible lambda function called over the active base of each
 * extension, should return true if the base should be processed
 * @param offset the offset from current index in all Extension objects of the
 * base to be digested
 *
 * @return BasesCounter object with summarized data from one base from each
 * extension
 */
BasesCounter count_bases(const vector<shared_ptr<Extension>>& extensions,
                         bool_predicate is_read_eligible,
                         int offset);


/**
 * @brief Count bases at specfic position in extensions.
 * @details Wrapper method for other count_bases method. All extensions
 * are eligible for usage.
 *
 * @param extensions Extension sequences, each extension stores the index of the
 * first unprocessed index in it's sequence
 *
 * @return BasesCounter object with summarized data from one base from each
 * extension
 */
BasesCounter count_bases(const vector<shared_ptr<Extension>>& extensions);


}  // namespace bases


#endif  // BASES_H
