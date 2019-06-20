/* Copyright (C) 2007 LinBox
 * Written by JG Dumas
 *
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

/*! @file algorithms/cra-builder-early-multip.h
 * @ingroup algorithms
 * @brief NO DOC
 */

#ifndef __LINBOX_cra_early_multip_H
#define __LINBOX_cra_early_multip_H

#include "linbox/util/timer.h"
#include <stdlib.h>
#include "linbox/integer.h"
#include "linbox/solutions/methods.h"
#include <vector>
#include <utility>

#include "linbox/algorithms/cra-builder-single.h"
#include "linbox/algorithms/cra-builder-full-multip.h"


namespace LinBox
{

	/*!  @brief NO DOC
	 * @ingroup CRA
	 *
	 */

	template<class Domain_Type>
	struct CRABuilderEarlyMultip : public CRABuilderEarlySingle<Domain_Type>, public CRABuilderFullMultip<Domain_Type> {
		typedef Domain_Type			Domain;
		typedef typename Domain::Element DomainElement;
		typedef CRABuilderEarlyMultip<Domain> 		Self_t;

	protected:
		// Random coefficients for a linear combination
		// of the elements to be reconstructed
		std::vector< size_t >      	randv;

		Integer& result(Integer &d) { std::cout << "should not be called" << std::endl; return d ;} ; // DON'T TOUCH
	public:

		CRABuilderEarlyMultip(const size_t EARLY=LINBOX_DEFAULT_EARLY_TERMINATION_THRESHOLD) :
			CRABuilderEarlySingle<Domain>(EARLY), CRABuilderFullMultip<Domain>()
		{}

		Integer& getModulus(Integer& m)
		{
			CRABuilderEarlySingle<Domain>::getModulus(m);
			return m;
		}
		Integer& getResidue(Integer& m)
		{
			CRABuilderEarlySingle<Domain>::getResidue(m);
			return m;
		}

		template<template<class T> class Vect>
		Vect<Integer>& getResidue(Vect<Integer>& m)
		{
			CRABuilderFullMultip<Domain>::getResidue(m);
			return m;
		}

		//! Init
		template<template<class T> class Vect>
		void initialize (const Integer& D, const Vect<Integer>& e)
		{
			srand48(BaseTimer::seed());
			randv. resize ( e.size() );
			for ( std::vector<size_t>::iterator int_p = randv. begin(); int_p != randv. end(); ++ int_p)
				*int_p = ((size_t)lrand48()) % 20000;
			Integer z;
			dot(z, D, e, randv);
			CRABuilderEarlySingle<Domain>::initialize(D, z);
			CRABuilderFullMultip<Domain>::initialize(D, e);
		}

		template<class Vect>
                //		template<template <class> class Alloc, template<class, class> class Vect>
		void initialize (const Domain& D, const Vect& e)
		{
			// Random coefficients for a linear combination
			// of the elements to be reconstructed
			srand48(BaseTimer::seed());
			randv. resize ( e.size() );
			for ( std::vector<size_t>::iterator int_p = randv. begin();
			      int_p != randv. end(); ++ int_p)
				*int_p = ((size_t)lrand48()) % 20000;
			DomainElement z;
			// Could be much faster
			// - do not compute twice the product of moduli
			// - reconstruct one element of e until Early Termination,
			//   then only, try a random linear combination.
			CRABuilderEarlySingle<Domain>::initialize (D, dot(z, D, e, randv));
			CRABuilderFullMultip<Domain>::initialize (D, e);
		}

		template<class OKDomain>
		void initialize (const Domain& D, const BlasVector<OKDomain>& e)
		{
			// Random coefficients for a linear combination
			// of the elements to be reconstructed
			srand48(BaseTimer::seed());
			randv. resize ( e.size() );
			for ( std::vector<size_t>::iterator int_p = randv. begin();
			      int_p != randv. end(); ++ int_p)
				*int_p = ((size_t)lrand48()) % 20000;
			DomainElement z;
			// Could be much faster
			// - do not compute twice the product of moduli
			// - reconstruct one element of e until Early Termination,
			//   then only, try a random linear combination.
			CRABuilderEarlySingle<Domain>::initialize(D,dot(z, D, e, randv) );
			CRABuilderFullMultip<Domain>::initialize(D, e);
		}

		//! Progress
		template<template<class T> class Vect>
		void progress (const Integer& D, const Vect<Integer>& e)
		{

			Integer z;
			CRABuilderEarlySingle<Domain>::progress(D, dot(z, D, e, randv));
			CRABuilderFullMultip<Domain>::progress(D, e);
		}

#if 1
		template<class Vect>
                //		template<template <class> class Alloc, template<class, class> class Vect>
		void progress (const Domain& D, const Vect& e)
		{
			// DomainElement z;
			/*!@todo Could be much faster
			  - do not compute twice the product of moduli
			  - reconstruct one element of e until Early Termination,
			  then only, try a random linear combination.
			*/
                        DomainElement z;
                        CRABuilderEarlySingle<Domain>::progress(D, dot(z, D, e, randv));
                        CRABuilderFullMultip<Domain>::progress(D, e);
		}
#endif

		template<class OKDomain>
		void progress (const Domain& D, const BlasVector<OKDomain>& e)
		{
			DomainElement z;
			/*!@todo Could be much faster
			  - do not compute twice the product of moduli
			  - reconstruct one element of e until Early Termination,
			  then only, try a random linear combination.
			*/
			CRABuilderEarlySingle<Domain>::progress(D, dot(z, D, e, randv));
			CRABuilderFullMultip<Domain>::progress(D, e);
		}

		//! Result
//		template<template <class> class Alloc, template<class, class> class Vect>
		template<class Vect>
                Vect& result(Vect& d)
		{
			return CRABuilderFullMultip<Domain>::result(d);
		}

		BlasVector<Givaro::ZRing<Integer> >& result(BlasVector<Givaro::ZRing<Integer> >& d)
		{
			return CRABuilderFullMultip<Domain>::result(d);
		}

		//! terminate
		bool terminated()
		{
			return CRABuilderEarlySingle<Domain>::terminated();
		}

		bool noncoprime(const Integer& i) const
		{
			return CRABuilderEarlySingle<Domain>::noncoprime(i);
		}

		bool changeVector()
		{
			for ( std::vector<size_t>::iterator int_p = randv. begin();int_p != randv. end(); ++ int_p)
				*int_p = ((size_t)lrand48()) % 20000;

			std::vector<Integer> e(randv.size());
			/* clear CRAEarlySingle; */
			CRABuilderEarlySingle<Domain>::occurency_ = 0;
			CRABuilderEarlySingle<Domain>::nextM_ = 1UL;
			CRABuilderEarlySingle<Domain>::primeProd_ = 1UL;
			CRABuilderEarlySingle<Domain>::residue_ = 0;

			/* Computation of residue_ */
            for (auto it = CRABuilderFullMultip<Domain>::shelves_begin();
                 it != CRABuilderFullMultip<Domain>::shelves_end();
                 ++it)
            {
                if (it->occupied) {
					Integer D = it->mod();
					std::vector<Integer> e_v(randv.size());
					e_v = it->residue;
					Integer z;
					dot(z,D, e_v, randv);
					Integer prev_residue_ = CRABuilderEarlySingle<Domain>::residue_;
					CRABuilderEarlySingle<Domain>::progress(D,z);
					if (prev_residue_ == CRABuilderEarlySingle<Domain>::residue_ )
						CRABuilderEarlySingle<Domain>::occurency_ += it->count;
					if ( CRABuilderEarlySingle<Domain>::terminated() ) {
						return true;
					}
                }
			}
			return false;
		}

	protected:

		/*! @bug why a dot product here ?
		 */
		template <class Vect1, class Vect2>
		Integer& dot (Integer& z, const Integer& D, const Vect1& v1, const Vect2& v2)
		{
			z = 0;
			typename Vect1::const_iterator v1_p;
			typename Vect2::const_iterator v2_p;
			for (v1_p  = v1. begin(), v2_p = v2. begin(); v1_p != v1. end(); ++ v1_p, ++ v2_p) {
				z = (z + (*v1_p)*(*v2_p))%D;
			}
			return z;
		}

		/*! @bug why a dot product here ?
		 */
		template <class Vect1, class Vect2>
		DomainElement& dot (DomainElement& z, const Domain& D,
				    const Vect1& v1,
				    const Vect2& v2)
		{

			D.assign(z,D.zero); DomainElement tmp;
			typename Vect1::const_iterator v1_p;
			typename Vect2::const_iterator v2_p;
			for (v1_p  = v1. begin(), v2_p = v2. begin();
			     v1_p != v1. end();
			     ++ v1_p, ++ v2_p)
				D.axpyin(z, (*v1_p), D.init(tmp, (*v2_p)));

			//             commentator().report(Commentator::LEVEL_ALWAYS, INTERNAL_DESCRIPTION) << "v: " << v2 << std::endl;
			//             commentator().report(Commentator::LEVEL_ALWAYS, INTERNAL_DESCRIPTION) << "z: " << z << std::endl;
			return z;
		}

		template <class Vect2, class OKDomain>
		DomainElement& dot (DomainElement& z, const Domain& D,
				    const BlasVector<OKDomain>& v1,
				    const Vect2& v2)
		{

			D.assign(z,D.zero); DomainElement tmp;
			typename BlasVector<Domain>::const_iterator v1_p;
			typename Vect2::const_iterator v2_p;
			for (v1_p  = v1. begin(), v2_p = v2. begin();
			     v1_p != v1. end();
			     ++ v1_p, ++ v2_p)
				D.axpyin(z, (*v1_p), D.init(tmp, (*v2_p)));

			//             commentator().report(Commentator::LEVEL_ALWAYS, INTERNAL_DESCRIPTION) << "v: " << v2 << std::endl;
			//             commentator().report(Commentator::LEVEL_ALWAYS, INTERNAL_DESCRIPTION) << "z: " << z << std::endl;
			return z;
		}

	};
}
#endif //__LINBOX_cra_early_multip_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
