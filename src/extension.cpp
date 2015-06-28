/**
 * @file extension.cpp
 * @author Marko Culinovic <marko.culinovic@gmail.com>
 * @author Luka Sterbic <luka.sterbic@gmail.com>
 * @brief Implementation file for Extension class
 * @details Implementation file for Extension class. It is used as representation
 * of possible extension reads of contig. It provides functionality
 * for local realignment method used in contig extension process.
 */
#include "./extension.h"

#include <string>

using std::string;

Extension::Extension(uint32_t read_id, const string& seq,
                     bool drop): is_droped(drop),
                                 read_id_(read_id),
                                 seq_(seq),
                                 curr_pos_(0) {}


void Extension::do_operation(const Operation& op) {
    curr_pos_ += op;
}
