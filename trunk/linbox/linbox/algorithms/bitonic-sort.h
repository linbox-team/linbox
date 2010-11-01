/* Copyright (C) LinBox
 * Written by Zhendong Wan 
 *
 *
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#ifndef __LINBOX_bitonic_sort_H
#define __LINBOX_bitonic_sort_H

/*! @file algorithms/bitonic-sort.h
 * Implement bitonic sorting network 
 */

#include <algorithm>

namespace LinBox
{


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

#endif //__LINBOX_bitonic_sort_H
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
