/**
 * @file bases.cpp
 * @author Marko Culinovic <marko.culinovic@gmail.com>
 * @author Luka Sterbic <luka.sterbic@gmail.com>
 * @brief Implementation file for BasesCounter class, and various functions which
 * calculate data for object of this class.
 * @details Implementation file for BasesCounter class and supporting functions. It is 
 * used for calculating various data at specific position in contig extension process
 */
#include <vector>
#include <cstring>

#include "./bases.h"
#include "./utility.h"


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
