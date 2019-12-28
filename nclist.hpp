//============================================================================
// Name        : nclist.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Nested Containment Lists
//============================================================================

#include <algorithm>
#include <iostream>
#include <iterator>
#include <stack>
#include <string>
#include <vector>

typedef unsigned int uint; //could or should be bigger?

template <class CInterval>
class NCList {

	//tagged with sublist info (in H, or just make another struct or something to directly refer to indices?) and sublist id,
	//welp those are stuck on forever
	//could use unsigned ints, but then I need to use 0? 0-based coords also exist, so nah
	struct TaggedInterval : public CInterval{ //could edit to use composition rather than inheritance, use pointers or something
		//use size_type instead of int??????
		int H_idx; //-1 per paper, tells us where children live in h_list!
		int sub_id; //where (parent) sublist is in H, top-level is 0
		TaggedInterval(): H_idx(-1), sub_id(-1) {}
		TaggedInterval(const CInterval& c): CInterval(c), H_idx(-1), sub_id(-1) {}
		bool contains(const CInterval& i) const {
			return this->getStart()<=i.getStart() && i.getEnd()<this->getEnd();
		}
		static bool cmp(const TaggedInterval& a, const TaggedInterval& b) {
			return a.sub_id<b.sub_id;
		}
		bool operator<(const CInterval& o) const { //I should be nice and implement the other operators... this does the OL order per paper
			return (this->getStart()==o.getStart())?(this->getEnd()>o.getEnd()):(this->getStart()<o.getStart());
		}
	};

	struct H_list_el {
		int L_start_idx;
		int num_intervals;
		H_list_el(int idx=0, int n=0): L_start_idx(idx), num_intervals(n){}
		void set(int a, int b) {
			L_start_idx=a;
			num_intervals=b;
		}
	};

	//could be protected
	//stores vectors (would be std::array or std::tuple) of {first_index, length} of sublists
	//h_list[0] refers to top-level
	std::vector<H_list_el> h_list;
	std::vector<TaggedInterval> l_list;

	public:
	//should this be a constructor? could be a static method like that German version
	NCList(const std::vector<CInterval>& in){
		//copy stuff to l_list, there should be a way to accomplish this "tagging" without wrapping and copying :/
		//maybe with that composition not inheritance thing
		l_list.reserve(in.size());
		for (typename std::vector<CInterval>::const_iterator it=in.begin(); it!=in.end(); ++it)
			l_list.push_back(*it);

		h_list.push_back(H_list_el()); //blank to offset
		std::sort(l_list.begin(), l_list.end());
		std::stack<TaggedInterval> parents; //could be CInterval, but is that any more efficient?
		std::stack<int> h_indices; //holds the h_index of the parent
		h_indices.push(0);

		for (typename std::vector<TaggedInterval>::iterator it = l_list.begin(), dupe;
				it!=l_list.end();
				++it){
			while(parents.size()!=0 && !parents.top().contains(*it)) {
				parents.pop();
				h_indices.pop();
			}
			it->sub_id=h_indices.top();
			if (it!=l_list.end()-1 && it->contains(*++(dupe=it))) {
				parents.push(*it);
				h_indices.push(h_list.size());
				h_list.push_back(H_list_el());
				it->H_idx=h_indices.top();
			}
		}

		std::stable_sort(l_list.begin(), l_list.end(), TaggedInterval::cmp); //dunno if stable is necessary
		int curr_idx=0, l_start_idx=0, num_intervals=0;
		for (int i = 0; i<(int)l_list.size(); ++i){
			if (l_list[i].sub_id == curr_idx) {
				++num_intervals;
			} else {
				++curr_idx;
				l_start_idx = i;
				num_intervals = 1;
			}

			h_list[curr_idx].set(l_start_idx, num_intervals); // can do this part in-place probably, like no need for temp variables
		}
	}

	//so sketchy
	static bool compEnd(const CInterval& a, uint b) {
		return a.getEnd()<b;
	}

	void overlaps(const CInterval& query, std::vector<CInterval>& res) const {
		overlaps(query.getStart(), query.getEnd(), res);
	}

	void overlaps(uint start, uint end, std::vector<CInterval>& res, int h_idx=0) const {
		if (h_idx == -1) return;
		int subend = h_list[h_idx].L_start_idx + h_list[h_idx].num_intervals;
		typename std::vector<TaggedInterval>::const_iterator it = std::lower_bound(l_list.begin()+h_list[h_idx].L_start_idx, l_list.begin()+subend, start, compEnd);
		while (it-l_list.begin()<subend && it->getStart() <= end){ //i.e., while overlap
			res.push_back(*it);
			overlaps(start, end, res, it->H_idx);
			++it;
		}
	}

	void printNC() const {
		std::cout << "H: [";
		for(typename std::vector<H_list_el>::const_iterator it=h_list.begin(); it!=h_list.end(); ++it) {
			std::cout << "(" << it-> L_start_idx << ", " << it->num_intervals << ")" << (it!=h_list.end()-1?", ":"");
		}
		std::cout << "]\n";

		std::cout << "L: [";
		for(typename std::vector<TaggedInterval>::const_iterator it = l_list.begin(); it!=l_list.end(); ++it) {
			std::cout<<"("<< it->getStart()<< ", "<< it->getEnd() << ", "<< it->H_idx<<", "<<it->sub_id << ")"<<(it!=l_list.end()-1?", ":"");
		}
		std::cout << "]\n";
	}
};
