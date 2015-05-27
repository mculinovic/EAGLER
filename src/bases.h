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


#define NUM_BASES 4


typedef std::function<bool(char)> bool_predicate;

namespace bases {


class BasesCounter {
 public:
    uint32_t count[NUM_BASES];
    uint32_t coverage;
    uint32_t max_idx;

    BasesCounter();
    void digest_base(char base);
    void refresh_stats();
};


BasesCounter count_bases(const vector<shared_ptr<Extension>>& extensions,
                    bool_predicate is_read_eligible,
                    int offset);


BasesCounter count_bases(const vector<shared_ptr<Extension>>& extensions);



}  // namespace bases


#endif
