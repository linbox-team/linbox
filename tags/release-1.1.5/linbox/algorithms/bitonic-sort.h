/* -*-mode:C++ -*- */
/* Written by Zhendong Wan */
/* Implement bitonic sorting network */


#ifndef BITONIC_SORT_H__
#define BITONIC_SORT_H__
#include <algorithm>

namespace LinBox{


	/* end - begin must be a power of 2*/
	template <class Iterator, class Comparator>
		void bitonicSort(Iterator begin, Iterator end, const Comparator& comparator = Comparator());

	/* end - begin must be a power of 2*/
	template <class Iterator, class Comparator>
		void bitonicMerge(Iterator begin, Iterator end, const Comparator& comparator = Comparator());

	///
	template<class Iterator, class Comparator>
		void bitonicSort(Iterator begin, Iterator end, const Comparator& comparator){

		if (end - begin >= 2) {
			Iterator mid = begin + (end - begin) / 2;
			
			// Sort the first half
			bitonicSort(begin, mid, comparator);
			
			// Sort the second half
			bitonicSort(mid, end, comparator);
			
			// reverse the order of second half
			std::reverse(mid, end);
			
			// Bitonic merge two halves
			bitonicMerge(begin, end, comparator);
		}
		
	}
	
	template<class Iterator, class Comparator>
		void bitonicMerge(Iterator begin, Iterator end, const Comparator& comparator){
		
		if (end - begin >= 2) {
			
			Iterator mid = begin + (end - begin) / 2;
			
			Iterator p1 = begin;
			
			Iterator p2 = mid;
			
			// Compare network
			for (p1 = begin, p2 = mid; p2 != end; ++ p1, ++ p2)
				comparator(*p1, *p2);

			// Bitonic Merge first half
			bitonicMerge(begin, mid, comparator);

			// Bitonic Merge second half
			bitonicMerge(mid, end, comparator);
		}
	}
	
}

#endif
