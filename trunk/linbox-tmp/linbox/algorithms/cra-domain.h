/* linbox/algorithms/cra-domain.h
 * Copyright (C) 1999-2010 The LinBox group
 *
 * Selector for ChineseRemainder
 * Parallel versions are transparent to the user
 * Time-stamp: <30 Mar 10 15:11:42 Jean-Guillaume.Dumas@imag.fr>
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
 * @brief Wrapper around OMP/SEQ version of ChineseRemainder.
 * @ingroup algorithms
 * @ingroup CRA
 *
 * If LINBOX_USES_OPENMP is defined, the we use ChineseRemainderOMP, else
 * we fall back to ChineseRemainderSeq
 */

#include "linbox/integer.h"
#ifdef LINBOX_USES_OPENMP

#include "linbox/algorithms/cra-domain-omp.h"
namespace LinBox
{
	/*! @brief Wrapper around OMP/SEQ version of ChineseRemainderXXX<CRABase>.
	 * \ingroup CRA
	 *
	 * If LINBOX_USES_OPENMP is defined, the we use ChineseRemainderOMP, else
	 * we fall back to ChineseRemainderSeq
	 *
	 * This is the OMP version
	 */

	template<class CRABase>
	struct ChineseRemainder : public ChineseRemainderOMP<CRABase> {
		typedef typename CRABase::Domain	Domain;
		typedef typename CRABase::DomainElement	DomainElement;

		template<class Param>
		ChineseRemainder(const Param& b) :
			ChineseRemainderOMP<CRABase>(b)
		{}

		ChineseRemainder(const CRABase& b) :
			ChineseRemainderOMP<CRABase>(b)
		{}
	};
}

#else

#include "linbox/algorithms/cra-domain-seq.h"
namespace LinBox
{
	/*! @brief Wrapper around OMP/SEQ version of ChineseRemainderXXX<CRABase>.
	 * \ingroup CRA
	 *
	 * If LINBOX_USES_OPENMP is defined, the we use ChineseRemainderOMP, else
	 * we fall back to ChineseRemainderSeq
	 *
	 * This is the SEQ version
	 */

	template<class CRABase>
	struct ChineseRemainder : public ChineseRemainderSeq<CRABase> {
		typedef typename CRABase::Domain	Domain;
		typedef typename CRABase::DomainElement	DomainElement;

		template<class Param>
		ChineseRemainder(const Param& b) :
			ChineseRemainderSeq<CRABase>(b)
		{}

		ChineseRemainder(const CRABase& b) :
			ChineseRemainderSeq<CRABase>(b)
		{}
	};
}

#endif


#endif

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

