/**
 * @file extension.h
 * @author Marko Culinovic <marko.culinovic@gmail.com>
 * @author Luka Sterbic <luka.sterbic@gmail.com>
 * @brief Header file for Extension class
 * @details Header file for Extension class. It is used as representation
 * of possible extension reads of contig. It provides functionality
 * for local realignment method used in contig extension process.
 */

#ifndef EXTENSION_H
#define EXTENSION_H

#include <string>


using std::string;


/**
 * Enum Operation is used for
 * handling moves in local alignment.
 */
enum Operation {
    match = 1,
    mismatch = 1,
    insertion_1 = 2,
    deletion_1 = 0
};


/**
 * @brief Extension class.
 * @details Extension class is used as representation
 * of possible extension reads of contig. It provides functionality
 * for local realignment method used in process.
 */
class Extension {
 public:
    /**
     * @brief Extension class constructor.
     *
     * @param read_id Id of read that is possible extension.
     * @param seq Subsequence of read sequence that is possible extension.
     * @param drop Bool value that representes if this read is dropped.
     */
    Extension(uint32_t read_id, const string& seq, bool drop);


    /**
     * @brief Getter for read Id.
     * @return Read Id.
     */
    uint32_t read_id() { return read_id_; }

    /**
     * @brief Getter for possible extension sequence.
     * @return Extension sequence.
     */
    const string& seq() { return seq_; }


    /**
     * @brief Getter for current position in extension
     * durring extension process.
     * @return Index of current position in extension sequence.
     */
    uint32_t curr_pos() { return curr_pos_; }


    /**
     * @brief Local realignment operation executor
     * @details Depending on operation current position
     * in extension sequence is moved ahead by 1 or 2 or
     * remains unchanged.
     *
     * @param op Alignment operation.
     */
    void do_operation(const Operation& op);

    // Flag represents if read is dropped.
    bool is_droped;

 private:
    // Read Id
    uint32_t read_id_;

    // Sequence that is possible extension.
    string seq_;

    // Current position in sequence during
    // contig extension process.
    uint32_t curr_pos_;
};

#endif  // EXTENSION_H
