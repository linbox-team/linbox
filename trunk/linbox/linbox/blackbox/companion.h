/* linbox/blackbox/companion.h
 * Copyright(c) LinBox
 *
 * Written by David Saunders <saunders@cis.udel.edu>
 *  ========LICENCE========
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
 */

#ifndef __LINBOX_companion_H
#define __LINBOX_companion_H

#include "linbox/blackbox/blackbox-interface.h"
#include "linbox/blackbox/triplesbb.h"
#include <vector>

namespace LinBox
{

	/** \ingroup blackbox
	  \brief %Companion matrix of a monic polynomial.
	  */
	template<class MatDom>
	struct Companion: public TriplesBB<MatDom> {
		typedef MatDom MatrixDomain;
		typedef typename MatrixDomain::Field Field;

		/// This is the n by n companion matrix of a given polynomial of degree n.
		template<class Polynomial>
		Companion(const Field& F = Field(), const Polynomial& P = Polynomial(1)) :
			TriplesBB<MatrixDomain>(F, P.size()-1, P.size()-1)
		{
			size_t n = P.size() - 1;
			for (size_t i = 1; i < n; ++i)
				this->setEntry(i, i-1, F.one);
			for (size_t i = 0; i < n; ++i) {
				typename Field::Element x;
				F.init(x, 0);
				F.neg(x, P[i]);
				this->setEntry(i, n-1, x);
			}
		}// Companion cstor




		/**
		 * \brief This constructs a random companion matrix.

		 Builds n by n matrix from degree n monic poly with other coefficients random.
		 */
		Companion(const Field& F, size_t n,
			  typename Field::RandIter r ) :
			TriplesBB<MatrixDomain>(F, n, n)
		{
			std::vector<typename Field::Element> p(n+1);
			for (typename std::vector<typename Field::Element>::iterator i = p.begin(); i != p.end(); ++i)
				r.random(*i); // we'll pretend p[n] == 1, ok?

			for (size_t i = 1; i < n; ++i) setEntry(i, i-1, F.one);
			for (size_t i = 0; i < n; ++i)
			{	typename Field::Element x;
				F.init(x, 0);
				F.neg(x, p[i]);
				setEntry(i, n-1, x);
			}

		}

		Companion(const Field& F, size_t n) :
			TriplesBB<MatrixDomain>(F,n,n)
		{
			typename Field::RandIter r(F);
			std::vector<typename Field::Element> p(n+1);
			for (typename std::vector<typename Field::Element>::iterator i = p.begin(); i != p.end(); ++i)
				r.random(*i); // we'll pretend p[n] == 1, ok?

			// const size_t indexbase = 1;
			for (size_t i = 1; i < n; ++i)
				this->setEntry(i, i-1, F.one);
			for (size_t i = 0; i < n; ++i)
			{	typename Field::Element x;
				F.init(x, 0);
				F.neg(x, p[i]);
				this->setEntry(i, n-1, x);
			}

		}

		template<typename _Tp1>
		struct rebind {
			typedef Companion<_Tp1> other;
		};




		// companion would be faster if built direct, using one axpy per entry: y_i = x_i-1 + p_i*x_n

	}; //Companion class

} //LinBox

#endif //__LINBOX_companion_H



// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
