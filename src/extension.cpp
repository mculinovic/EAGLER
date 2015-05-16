// Copyright
#include "./extension.h"

#include <string>

using std::string;

Extension::Extension(uint32_t id, const string& seq,
                     bool drop): is_droped(drop),
                                 id_(id),
                                 seq_(seq),
                                 curr_pos_(0) {}


void Extension::do_operation(const Operation& op) {
    curr_pos_ += op;
}
