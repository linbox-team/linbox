/** -*- mode:C++ -*- */

/** File: Compose-traits.h
 *  Author: Zhendong Wan
 *  the ComposeTrait class defines the return type of composition of three matrices
 */

#ifndef __LINBOX_COMPOSE_TRAITS_H__
#define __LINBOX_COMPOSE_TRAITS_H__

#include <linbox/blackbox/dense.h>
//#include <linbox/blackbox/sparse.h>
#include <linbox/blackbox/compose.h>

namespace LinBox{

	template<class IMatrix>
		class ComposeTraits {
		
		public:
		typedef Compose<IMatrix, IMatrix> value_type;
	};
		
	
	template<class Field>
		class ComposeTraits<DenseMatrix<Field> > {
		public:
		
		// define the return value type
		typedef DenseMatrix<Field> value_type;
	};	
}


#endif
