// Copyright
#ifndef EXTENSION_H
#define EXTENSION_H

#include <string>

using std::string;


enum Operation {
    match = 1,
    mismatch = 1,
    insertion_1 = 2,
    deletion_1 = 0
};

class Extension {
public:
    Extension(uint32_t read_id, const string& seq, bool drop);
    uint32_t read_id() { return read_id_; }
    const string& seq() { return seq_; }
    uint32_t curr_pos() { return curr_pos_; }
    void do_operation(const Operation& op);

    bool is_droped;

private:
    uint32_t read_id_;
    string seq_;
    uint32_t curr_pos_;
};

#endif  // EXTENSION_H
