#include <vector>

#include "./bases.h"
#include "./utility.h"

namespace bases {


vector<int> count_bases(const vector<shared_ptr<Extension>>&           extensions,
                        bool_predicate is_read_eligible,
                        int offset) {
	vector<int> bases(4, 0);

    for (size_t j = 0; j < extensions.size(); ++j) {
        auto& extension = extensions[j];

        if (!extension->is_droped) {
	        auto& seq = extension->seq();
	        uint32_t position = extension->curr_pos();

	        if (position + offset < seq.length() &&
	            is_read_eligible(seq[position])) {
	            int idx = utility::base_to_idx(seq[position + offset]);
	            bases[idx]++;
	        }
        }
    }

    return bases;
}


vector<int> count_bases(const vector<shared_ptr<Extension>>& extensions) {
    auto is_read_eligible = [](char c) -> bool { (void) c; return true; };
    return count_bases(extensions, is_read_eligible, 0);
}


vector<int> count_bases(const vector<string>& extensions,
                        const vector<uint32_t>& read_positions,
                        bool_predicate is_read_eligible,
                        int offset) {
    vector<int> bases(4, 0);
    for (size_t j = 0; j < extensions.size(); ++j) {
        auto& read = extensions[j];

        if (read_positions[j] + offset < read.length() &&
            is_read_eligible(read[read_positions[j]])) {
            int idx = utility::base_to_idx(read[read_positions[j] + offset]);
            bases[idx]++;
        }
    }
    return bases;
}


vector<int> count_bases(const vector<string>& extensions,
                        const vector<uint32_t>& read_positions) {
    auto is_read_eligible = [](char c) -> bool { (void) c; return true; };
    return count_bases(extensions, read_positions, is_read_eligible, 0);
}


vector<int> count_bases(const vector<string>& extensions, int pos) {
    vector<uint32_t> read_positions(extensions.size(), pos);
    return count_bases(extensions, read_positions);
}


pair<int, int> get_bases_stats(const vector<int>& bases) {
    int coverage = 0;
    int max_idx = 0;

    for (int i = 0; i < NUM_BASES; ++i) {
        coverage += bases[i];
        if (bases[i] > bases[max_idx]) {
            max_idx = i;
        }
    }

    return std::make_pair(coverage, max_idx);
}


}
