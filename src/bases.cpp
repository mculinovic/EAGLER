/**
 * @file bases.cpp
 * @copyright Coypright 2015 Marko Culinovic, Luka Sterbic
 * @author Marko Culinovic <marko.culinovic@gmail.com>
 * @author Luka Sterbic <luka.sterbic@gmail.com>
 * @brief Implementation file for the BasesCounter class and associated factory
 * functions.
 * @details Implementation file for the BasesCounter class and supporting
 * associated factory functions. It is used for calculating various statistics
 * at specific positions in the contig extension process.
 */
#include <vector>
#include <cstring>

#include "bases.h"
#include "utility.h"


namespace bases {


BasesCounter::BasesCounter() {
    std::memset(count, 0, NUM_BASES * sizeof(uint32_t));
    coverage = 0;
    max_idx = 0;
}


void BasesCounter::digest_base(char base) {
    count[utility::base_to_idx(base)]++;
}


void BasesCounter::refresh_stats() {
    coverage = count[0];
    max_idx = 0;

    for (uint32_t index = 1; index < NUM_BASES; ++index) {
        coverage += count[index];

        if (count[index] > count[max_idx]) {
            max_idx = index;
        }
    }
}


BasesCounter count_bases(const vector<shared_ptr<Extension>>& extensions,
                        bool_predicate is_read_eligible,
                        int offset) {
    BasesCounter counter;

    for (size_t j = 0; j < extensions.size(); ++j) {
        auto& extension = extensions[j];

        if (!extension->is_droped) {
            auto& seq = extension->seq();
            uint32_t position = extension->curr_pos();

            if (position + offset < seq.length() &&
                is_read_eligible(seq[position])) {
                counter.digest_base(seq[position + offset]);
            }
        }
    }

    counter.refresh_stats();
    return counter;
}


BasesCounter count_bases(const vector<shared_ptr<Extension>>& extensions) {
    auto is_read_eligible = [](char c) -> bool { (void) c; return true; };
    return count_bases(extensions, is_read_eligible, 0);
}


}  // namespace bases
