/**
 * @file bases.h
 * @author Marko Culinovic <marko.culinovic@gmail.com>
 * @author Luka Sterbic <luka.sterbic@gmail.com>
 * @brief Header file for BasesCounter class, and various functions which
 * calculate data for object of this class.
 * @details Header file for BasesCounter class and supporting functions. It is 
 * used for calculating various data at specific position in contig extension process
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


// number of different bases in sequence
#define NUM_BASES 4


typedef std::function<bool(char)> bool_predicate;


/**
 * namespace for BasesCounter class and supporting functions
 */
namespace bases {

/**
 * @brief BasesCounter class. It is used for calculating
 * various data at specific position in contig extension process.
 * @details BasesCounter class is used for storing information
 * about appearances of each base at specific position, coverage
 * at this position and which base appeared maximum number of
 * times at this position.
 */
class BasesCounter {
 public:
    // array used for countig appearances of each base
    uint32_t count[NUM_BASES];
    // number of bases at this position - sum(count)
    uint32_t coverage;
    // index in count of base with maximum appearances
    uint32_t max_idx;

    /**
     * @brief BasesCounter class default constructor
     */
    BasesCounter();

    /**
     * @brief Method for proccesing base.
     * @details Method converts base character to
     * index in count array and increases its
     * counter for one.
     * 
     * @param base base at specific position in sequence
     */
    void digest_base(char base);

    /**
     * @brief Method updates data in class
     * @details Array count is used for calculating other
     * data (coverage and maximum index) which thiss class wrapps.
     */
    void refresh_stats();
};


/**
 * @brief Count bases at specific position in extensions.
 * @details Firstly, BasesCounter object is created. Base at
 * specific position in each extension is digested by this object.
 * Data is summarized in BasesCounter object.
 * 
 * @param extensions Extension sequences. Every extension has information
 * about which index is current/active in its sequence.
 * @param is_read_eligible Lambda function which returns bool value which
 * reoresents if this extension is still eligible for usage.
 * @param offset offset from current/active index in extension sequence
 * @return BasesCounter object with summarized data from extensions at
 * specific position
 */
BasesCounter count_bases(const vector<shared_ptr<Extension>>& extensions,
                    bool_predicate is_read_eligible,
                    int offset);


/**
 * @brief Count bases at specfic position in extensions.
 * @details Wrapper method for other count_bases method. All extensions
 * are eligible for usage.
 * 
 * @param extensions Extension sequences.
 * @return BasesCounter object with summarized data from extensions at
 * specific position
 */
BasesCounter count_bases(const vector<shared_ptr<Extension>>& extensions);



}  // namespace bases


#endif
