#ifndef BASES_H
#define BASES_H

#include <string>
#include <memory>
#include <vector>

#include "./extension.h"

using std::string;
using std::vector;
using std::pair;
using std::shared_ptr;


#define NUM_BASES 4


typedef std::function<bool(char)> bool_predicate;


namespace bases {


	vector<int> count_bases(vector<shared_ptr<Extension>>& extensions,
                        bool_predicate is_read_eligible,
                        int offset);


	vector<int> count_bases(vector<shared_ptr<Extension>>& extensions);


	vector<int> count_bases(const vector<string>& extensions,
                        const vector<uint32_t>& read_positions,
                        bool_predicate is_read_eligible,
                        int offset);


	vector<int> count_bases(const vector<string>& extensions,
	                        const vector<uint32_t>& read_positions);


	vector<int> count_bases(const vector<string>& extensions, int pos);


	pair<int, int> get_bases_stats(const vector<int>& bases);

}


#endif
