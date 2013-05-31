#ifndef __LINBOX_algorithm_iml_H
#define __LINBOX_algorithm_iml_H

// #include "linbox/matrix/blas-matrix.h"
// #include "linbox/integer.h"
// #include "linbox/algorithms/linbox-tags.h"
// #include "integer-matrix-tools.h"

namespace LinBox { namespace iml {

	void
	revseq(std::vector<size_t> & At,
	       const size_t r, const std::vector<size_t> & A)
	{

		size_t m = A.size();
		At.resize(m,0);
		for (size_t i = 0; i < m; i++) {
			At[i] = i;
		}
		for (size_t i = 1; i < r+1; i++) {
			if (A[i] != i)
			{
				size_t t = At[i-1];
				At[i-1] = At[A[i]-1];
				At[A[i]-1] = t;
			}
		}

		return ;
	}


} // iml
} // LinBox

#include "IML/row-echelon-transform.h"
#include "IML/RNS.h"
// #include "IML/p-adic-lift.h"
// #include "IML/nullspace.h"

#endif // __LINBOX_algorithm_iml_iml_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
