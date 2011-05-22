/* -*- mode:C++ -*- */

/* File: random-matrix-traits.h
 *  Author: Zhendong Wan
 */
#ifndef __RANDOM_MATRIX_TRAITS_H__
#define __RANDOM_MATRIX_TRAITS_H__

namespace LinBox {

	template<class Matrix>
	class RandomMatrixTraits{
		public:
		typedef Matrix value_type;
	};
}

#endif
