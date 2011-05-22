#ifndef __DENSITY_H__
#define __DENSITY_H__
#include <linbox/vector/vector-traits.h>

namespace LinBox {

	/** \brief Estimate nonzero entries in a vector, used in parallel elimination */
	template<class Vector>
	inline long density(const Vector& v) {
		
		return density(v, VectorTraits<Vector>::VectorCategory());
	}

	template<class Vector, class VectorCategory>
	inline long density(const Vector&, VectorCategory);

	template<class Vector>
	inline long density(const Vector& v, VectorCategories::DenseVectorTag) {

		return v.size();
	}


	template<class Vector>
	inline long density(const Vector& v, VectorCategories::SparseSequenceVectorTag) {

		return v.size();
	}


	template<class Vector>
	inline long density(const Vector& v, VectorCategories::SparseAssociativeVectorTag) {

		return v.size();
	}


	template<class Vector>
        inline long density(const Vector& v, VectorCategories::SparseParallelVectorTag) { 

		return v.first.size();
	}
}

#endif
