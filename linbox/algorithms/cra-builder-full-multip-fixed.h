/* Copyright (C) 2010 LinBox
 * Written by <Jean-Guillaume.Dumas@imag.fr>
 *
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

/*! @file algorithms/cra-builder-full-multip-fixed.h
 * @ingroup algorithms
 * @brief CRA for multi-residues.
 *
 * An upper bound is given on the size of the data to reconstruct.
 */

#ifndef __LINBOX_cra_full_multip_fixed_H
#define __LINBOX_cra_full_multip_fixed_H

#include "linbox/util/timer.h"
#include <stdlib.h>
#include "linbox/integer.h"
#include "linbox/solutions/methods.h"
#include <vector>
#include <utility>

#include "linbox/algorithms/lazy-product.h"
#include "linbox/algorithms/cra-builder-full-multip.h"


namespace LinBox
{

	/*! @ingroup CRA
	 * @brief Chinese Remaindering Algorithm for multiple residues.
	 * An upper bound is given on the size of the data to reconstruct.
	 */
	template<class Domain_Type>
	struct CRABuilderFullMultipFixed : CRABuilderFullMultip<Domain_Type> {
		typedef Domain_Type			Domain;
		typedef typename Domain::Element 	DomainElement;
		typedef CRABuilderFullMultipFixed<Domain> 	Self_t;

		/*! Constructor.
		 * @param p is a pair such that
		 * - \c p.first is the size of a residue (ie. it would be 1 for \"FullSingle\")
		 * - \c p.second is the theoretical upperbound (log2)
		 *   on the size of the integers to reconstruct.
		 * .
		 */
		CRABuilderFullMultipFixed(const std::pair<size_t,double>& p ) :
			CRABuilderFullMultip<Domain>(p.second, p.first)
        { }

		/*! Intialize to the first residue/prime.
		 * @param D domain
		 * @param e iterator on the first residue
		 * @pre any CRA should first call \c initialize before \c progress
		 */
		template<class Iterator>
		inline void initialize (const Domain& D, Iterator& e)
        { this->initialize_iter(D, e, this->dimension_); }

		/*! Add a new residue (ie take into account a new prime).
		 * @param D domain
		 * @param e iterator for the new residue, for instance, a
		 * <code>std::vector<T>::iterator</code>.
		 */
		template<class Iterator>
		inline void progress (const Domain& D, Iterator& e)
		{ this->progress_iter(D, e, this->dimension_); }

		/*! Compute the result.
		 * moves low occupied shelves up.
		 * @param[out] d an iterator for the result.
		 */
		template<class Iterator>
		void result (Iterator &d) const
        { this->result_iter(d); }
    };
}

namespace LinBox
{
	/*! NO DOC..
	 * @ingroup CRA
	 * Version of LinBox::CRABuilderFullMultip for matrices.
	 */
	template<class Domain_Type>
	struct CRABuilderFullMultipMatrix : CRABuilderFullMultip<Domain_Type> {
		typedef Domain_Type			Domain;
		typedef typename Domain::Element 	DomainElement;
		typedef CRABuilderFullMultipMatrix<Domain> 	Self_t;

		/*! Constructor.
		 * @param p is a pair such that
		 * - \c p.first is the size of a residue, it would be 1 for \"FullSingle\"
		 * - \c p.second is the theoretical upperbound (natural
		 *   logarithm) on the size of the integers to reconstruct.
		 * .
		 */
		CRABuilderFullMultipMatrix(const std::pair<size_t,double>& p ) :
			CRABuilderFullMultip<Domain>(p.second, p.first)
        { }

		/*! Intialize to the first residue/prime.
		 * @param D domain
		 * @param e
		 * @pre any CRA should first call \c initialize before \c progress
		 */
		template<class Matrix>
		void initialize (const Domain& D, Matrix& e)
        { this->initialize_iter(D, e.Begin(), this->dimension_); }

		/*! Add a new residue (ie take into account a new prime).
		 * @param D domain
		 * @param e
		 */
		template<class Matrix>
		void progress (const Domain& D, Matrix& e)
		{ this->progress_iter(D, e.Begin(), this->dimension_); }

		/*! Compute the result.
		 * moves low occupied shelves up.
		 * @param[out] d
		 */
		template<class Matrix>
		Matrix& result (Matrix &d)
        {
            this->result_iter(d.Begin());
            return d;
        }
	};

}

#endif //__LINBOX_cra_full_multip_fixed_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
