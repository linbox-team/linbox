/* linbox/algorithms/cra-domain.h
 * Copyright (C) 1999-2010 The LinBox group
 *
 * Selector for ChineseRemainder
 * Parallel versions are transparent to the user
 * Time-stamp: <28 Oct 22 18:12:06 Jean-Guillaume.Dumas@imag.fr>
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 */
#ifndef __LINBOX_cra_domain_H
#define __LINBOX_cra_domain_H

/*! @file algorithms/cra-domain.h
 * @brief Wrapper around PAR/SEQ version of ChineseRemainder.
 * @ingroup algorithms
 * @ingroup CRA
 *
 * If __LINBOX_USE_OPENMP is defined, then we use ChineseRemainderParallel, else
 * we fall back to ChineseRemainderSequential
 */

#include "linbox/integer.h"
#include "linbox/field/rebind.h"
#include "linbox/vector/vector.h"

namespace LinBox
{
	/*! @brief Return type for CRA iteration.
	 *
	 * The function object passed to the ChineseRemainder operators
	 * should return one of these values on each iteration.
	 */
	enum struct IterationResult {
		CONTINUE, ///< successful iteration; add to the result and keep going
		SKIP, ///< "bad prime"; ignore this result and keep going
		RESTART ///< all previous iterations were bad; start over with this residue
	};

	/** \brief Type information for the residue in a CRA iteration.
	 *
	 * Here ResultType should be some kind of vector type whose default constructor
	 * results in zero-dimensional vectors where no init'ing is required.
	 */
	template <typename ResultType, typename Function>
	struct CRAResidue {
		template <typename Domain>
		using ResidueType = typename Rebind<ResultType,Domain>::other;

		template <typename Domain>
		static ResidueType<Domain> create(const Domain& d) {
			return ResidueType<Domain>(d);
		}
	};

	/** \brief Type information for the residue in a CRA iteration.
	 *
	 * This is the specialization for scalar types (namely Integer) where
	 * the residue type (such as a Modular element) must be init'ed.
	 */
	template <typename Function>
	struct CRAResidue<Integer, Function> {
		template <typename Domain>
		using ResidueType = typename Domain::Element;

		template <typename Domain>
		static ResidueType<Domain> create(const Domain& d) {
			ResidueType<Domain> r;
			d.init(r);
			return r;
		}
	};

	/** \brief Type information for the residue in a CRA iteration.
	 *
	 * This is the specialization for a vector of scalar types (namely Integer)
	 */
	template <typename Function>
	struct CRAResidue<std::vector<Integer>, Function> {
		template <typename Domain>
		using ResidueType = DenseVector<Domain>;

		template <typename Domain>
		static ResidueType<Domain> create(const Domain& d) {
            return ResidueType<Domain>(d);
		}
	};
}

#ifdef __LINBOX_USE_OPENMP

#include "linbox/algorithms/cra-domain-parallel.h"
namespace LinBox
{
	/*! @brief Wrapper around PAR/SEQ version of ChineseRemainderXXX<CRABase>.
	 * \ingroup CRA
	 *
	 * If __LINBOX_USE_OPENMP is defined, the we use ChineseRemainderParallel, else
	 * we fall back to ChineseRemainderSequential
	 *
	 * This is the Parallel version
	 */

	template<class CRABase>
        using ChineseRemainder = ChineseRemainderParallel<CRABase>;
}

#else

#include "linbox/algorithms/cra-domain-sequential.h"
namespace LinBox
{
	/*! @brief Wrapper around PAR/SEQ version of ChineseRemainderXXX<CRABase>.
	 * \ingroup CRA
	 *
	 * If __LINBOX_USE_OPENMP is defined, the we use ChineseRemainderParallel, else
	 * we fall back to ChineseRemainderSequential
	 *
	 * This is the SEQ version
	 */
	template<class CRABase>
        using ChineseRemainder = ChineseRemainderSequential<CRABase>;
}

#endif


#endif

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
