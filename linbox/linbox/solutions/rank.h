/* linbox/solutions/rank.h
 * Copyright(C) LinBox
 * ------------------------------------
 *
 * ========LICENCE========
 * This file is part of the library LinBox.
 *
 * LinBox is free software: you can redistribute it and/or modify
 * it under the terms of the  GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 *.
 */

#ifndef __LINBOX_rank_H
#define __LINBOX_rank_H

//#include "linbox/linbox-config.h"
#include "linbox/field/modular.h"
#include "linbox/randiter/random-prime.h"
#include "linbox/algorithms/matrix-hom.h"
#include "linbox/matrix/sparse-matrix.h"
#include "linbox/blackbox/diagonal.h"
#include "linbox/blackbox/diagonal-gf2.h"
#include "linbox/blackbox/compose.h"
#include "linbox/blackbox/permutation.h"
#include "linbox/blackbox/transpose.h"
#include "linbox/algorithms/blackbox-container-symmetrize.h"
#include "linbox/algorithms/blackbox-container-symmetric.h"
#include "linbox/algorithms/blackbox-container.h"
#include "linbox/algorithms/massey-domain.h"
#include "linbox/algorithms/gauss.h"
#include "linbox/algorithms/gauss-gf2.h"
#include "linbox/matrix/matrix-domain.h"
#include "linbox/algorithms/whisart_trace.h"
#include "linbox/matrix/dense-matrix.h"
#include "linbox/matrix/MatrixDomain/blas-matrix-domain.h"
#include "linbox/switch/cekstv.h"
#include "linbox/blackbox/butterfly.h"


#include "linbox/vector/vector-traits.h"
#include "linbox/solutions/trace.h"
#include "linbox/solutions/methods.h"


#include "linbox/util/debug.h"

// Namespace in which all LinBox library code resides
namespace LinBox
{


	/**
	 * Compute the rank of a linear transform A over a field by selected method.
	 * \ingroup solutions
	 * For very large and/or very sparse matrices the Wiedemann method will be faster
	 * (and it is memory efficient).
	 * For some sparse matrices SparseElimination may outperform Wiedemann.
	 * For small or dense matrices BlasElimination will be faster.
	 * \param[out] r  output rank of A.
	 * \param[in]  A linear transform, member of any blackbox class.
	 * \param[in]  M may be a \p Method::Hybrid (the default), a \p Method::Wiedemann, a  \p Method::BlasElimination, or a \p Method::SparseElimination..
	 * \param      tag UNDOC
	 * \return a reference to r.
	 */
	template <class Blackbox, class Method, class DomainCategory>
	inline unsigned long &rank (unsigned long                   &r,
				    const Blackbox                  &A,
				    const DomainCategory          &tag,
				    const Method                   &M);

	template <class Blackbox, class Method, class DomainCategory>
	inline unsigned long &rankin (unsigned long                   &r,
				      Blackbox                        &A,
				      const DomainCategory          &tag,
				      const Method                   &M);


	/**
	 * Compute the rank of a linear transform A over a field.
	 * \ingroup solutions
	 * The default method is Wiedemann(), using diagonal preconditioning and
	 * the minpoly.  For small or dense matrices BlasElimination will be faster.
	 * \param      A linear transform, member of any blackbox class.
	 * \param[out] r rank of \p A
	 * \return     \p r rank of \p A.
	  */
	template <class Blackbox>
	inline unsigned long &rank (unsigned long                   &r,
				    const Blackbox                  &A)
	{
		return rank(r, A, typename FieldTraits<typename Blackbox::Field>::categoryTag(), Method::Hybrid());
	}

	/**
	 * Compute the rank of a linear transform A over a field.
	 * \ingroup solutions
	 *
	 * The default method is \p Wiedemann(), using diagonal preconditioning and
	 * the minpoly.  For small or dense matrices \p BlasElimination will be faster.
	 * \return \p r rank of \p A.
	 * \param A linear transform, member of any blackbox class.
	 * @param[out] r rank of \p A
	 * @param M method (see ???)
	 */
	template <class Blackbox, class Method>
	inline unsigned long &rank (unsigned long                   &r,
				    const Blackbox                  &A,
				    const Method    &M)
	{
		return rank(r, A, typename FieldTraits<typename Blackbox::Field>::categoryTag(), M);
	}

	/** Rank of \p A.
	 * \p A may be modified
	 * @param A matrix
	 * @param r rank
	*/
	template <class Blackbox>
	inline unsigned long &rankin (unsigned long                   &r,
				      Blackbox &A)
	{
		//! @bug there is no Elimination() method there.
		return rankin(r, A, typename FieldTraits<typename Blackbox::Field>::categoryTag(), Method::SparseElimination());
	}



} // LinBox

#include "linbox/solutions/rank.inl"

#endif // __LINBOX_rank_H


// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
