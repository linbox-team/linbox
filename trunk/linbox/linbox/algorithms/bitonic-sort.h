/** -*-mode:C++ -*- */
/* Written by Zhendong Wan */
/* Implement bitonic sorting network */


#ifndef BITONIC_SORT_H__
#define BITONIC_SORT_H__
#include <algorithm>

namespace LinBox{


	/* end - begin must be a power of 2*/
	template <class Iterator, class Compare>
		void bitonicSort(Iterator begin, Iterator end, const Compare& compare = Compare());

	/* end - begin must be a power of 2*/
	template <class Iterator, class Compare>
		void bitonicMerge(Iterator begin, Iterator end, const Compare& compare = Compare());


	template<class Iterator, class Compare>
		void bitonicSort(Iterator begin, Iterator end, const Compare& compare){

		if (end - begin >= 2) {
			Iterator mid = begin + (end - begin) / 2;
			
			// Sort the first half
			bitonicSort(begin, mid, compare);
			
			// Sort the second half
			bitonicSort(mid, end, compare);
			
			// reverse the order of second half
			std::reverse(mid, end);
			
			// Bitonic merge two halves
			bitonicMerge(begin, end, compare);
		}
		
	}
	
	template<class Iterator, class Compare>
		void bitonicMerge(Iterator begin, Iterator end, const Compare& compare){
		
		if (end - begin >= 2) {
			
			Iterator mid = begin + (end - begin) / 2;
			
			Iterator p1 = begin;
			
			Iterator p2 = mid;
			
			// Compare network
			for (p1 = begin, p2 = mid; p2 != end; ++ p1, ++ p2)
				compare(*p1, *p2);

			// Bitonic Merge first half
			bitonicMerge(begin, mid, compare);

			// Bitonic Merge second half
			bitonicMerge(mid, end, compare);
		}
	}
	
}

#endif
