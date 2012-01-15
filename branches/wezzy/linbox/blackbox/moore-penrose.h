/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* linbox/blackbox/moore-penrose.h
 * Copyright (C) 2001 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
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
 */

#ifndef __LINBOX_moore_penrose_H
#define __LINBOX_moore_penrose_H

#include "linbox/blackbox/blackbox-interface.h"
#include "linbox/blackbox/submatrix.h"
#include "linbox/blackbox/inverse.h"
#include "linbox/blackbox/transpose.h"
#include "linbox/blackbox/compose.h"
#include "linbox/vector/vector-traits.h"
#include "linbox/util/debug.h"
#include "linbox/util/error.h"

// Namespace in which all LinBox library code resides
namespace LinBox
{

	/** \brief Generalized inverse of a blackbox.  Efficiency concerns when many applications are used.
	 *
	 \ingroup blackbox
	 * Given an arbitrary matrix in black box representation, this black box
	 * represents the Moore-Penrose inverse of the matrix.
	 *
	 * This implementation assumes that A already has a nonsingular
	 * principal r x r minor. It is the caller's responsibility to ensure
	 * that that condition holds.
	 *
	 * Given MoorePenrose M(A, r), and vector b, we have that M.apply(u, b) provides
	 * the least norm, least squares solution x = u to Ax = b.
	 *
	 * TODO: remove the requirement that lpm is nonsingular.  Specialize for dense matrices.
	 */
	template <class Blackbox>
	class MoorePenrose : public BlackboxInterface {
	public:

		typedef typename Blackbox::Field Field;
		typedef typename Blackbox::Element Element;
		template<typename _Tp1>
		struct rebind {
		       	typedef MoorePenrose<typename Blackbox::template rebind<_Tp1>::other> other;
		};


		/** Constructor from field and dense vector of field elements.
		 * -param BB   Black box from which to extract the submatrix
		 * -param row  First row of the submatrix to extract (1.._BB->rowdim ())
		 * -param col  First column of the submatrix to extract (1.._BB->coldim ())
		 * -param rowdim Row dimension
		 * -param coldim Column dimension
		 *  @param A
		 *  @param rank
		 */
		MoorePenrose (const Blackbox *A, size_t rank) :
			_matA (A), _rank (rank)
		{
			_matB1     = new Submatrix<Blackbox> (_matA, 0, 0, rank, rank);
			_matF      = new Submatrix<Blackbox> (_matA, 0, 0, _matA->rowdim (), rank);
			_matGG     = new Submatrix<Blackbox> (_matA, 0, 0, rank, _matA->coldim ());
			_matFT     = new Transpose<Submatrix<Blackbox> > (_matF);
			_matGT     = new Transpose<Submatrix<Blackbox> > (_matGG);
			_matFTF    = new Compose<Transpose<Submatrix<Blackbox> >,Submatrix<Blackbox> > (_matFT, _matF);
			_matGGT    = new Compose<Submatrix<Blackbox>, Transpose<Submatrix<Blackbox> > > (_matGG, _matGT);
			_matFTFinv = new Inverse<Compose<Transpose<Submatrix<Blackbox> >,Submatrix<Blackbox> > > ( _matFTF);
			_matGGTinv = new Inverse<Compose<Submatrix<Blackbox>, Transpose<Submatrix<Blackbox> > > > ( _matGGT);
		}

		/** Copy constructor
		*/
		MoorePenrose (const MoorePenrose &A) :
			_matA      ( A._matA),
			_matB1     ( A._matB1),
			_matF      ( A._matF),
			_matGG     ( A._matGG),
			_matFT     ( A._matFT),
			_matGT     ( A._matGT),
			_matFTF    ( A._matFTF),
			_matGGT    ( A._matGGT),
			_matFTFinv ( A._matFTFinv),
			_matGGTinv ( A._matGGTinv),
			_rank      ( A._rank)
		{}

		/** Destructor
		*/
		~MoorePenrose ()
		{
			delete _matGGTinv;
			delete _matFTFinv;
			delete _matGGT;
			delete _matFTF;
			delete _matGT;
			delete _matFT;
			delete _matGG;
			delete _matF;
			delete _matB1;
		}



		/** Application of BlackBox matrix.
		 * <code>y= A*x</code>.
		 * Requires one vector conforming to the \ref LinBox
		 * vector @link Archetypes archetype@endlink.
		 * Required by abstract base class.
		 * @return reference to vector y containing output.
		 * @param  x constant reference to vector to contain input
		 * @param y
		 */
		template <class OutVector, class InVector>
		OutVector& apply (OutVector &y, const InVector& x) const
		{
			InVector _z1 (_rank);
			InVector _z2 (_rank);

			_matF->applyTranspose (_z1, x);
			_matFTFinv->apply (_z2, _z1);
			_matB1->apply (_z1, _z2);
			_matGGTinv->apply (_z2, _z1);
			_matGG->applyTranspose (y, _z2);

			return y;
		}

		/** Application of BlackBox matrix transpose.
		 * <code>y= transpose(A)*x</code>.
		 * Requires one vector conforming to the \ref LinBox
		 * vector @link Archetypes archetype@endlink.
		 * Required by abstract base class.
		 * @return reference to vector y containing output.
		 * @param  x constant reference to vector to contain input
		 * @param y
		 */
		template <class OutVector, class InVector>
		OutVector& applyTranspose (OutVector &y, const InVector& x) const
		{
			InVector _z1 (_rank);
			InVector _z2 (_rank);

			_matGG->apply (_z1, x);
			_matGGTinv->applyTranspose (_z2, _z1);
			_matB1->applyTranspose (_z1, _z2);
			_matFTFinv->applyTranspose (_z2, _z1);
			_matF->apply (y, _z2);

			return y;
		}

		/** Retreive _row dimensions of BlackBox matrix.
		 * This may be needed for applying preconditioners.
		 * Required by abstract base class.
		 * @return integer number of _rows of black box matrix.
		 */
		size_t rowdim (void) const
		{
			return _matA->coldim ();
		}

		/** Retreive _column dimensions of BlackBox matrix.
		 * Required by abstract base class.
		 * @return integer number of _columns of black box matrix.
		 */
		size_t coldim (void) const
		{
			return _matA->rowdim ();
		}

		const Field& field() { return _matA -> field(); }

	private:

		const Blackbox       * _matA;
		Submatrix<Blackbox>  * _matB1;
		Submatrix<Blackbox>  * _matF;
		Submatrix<Blackbox>  * _matGG;
		Transpose<Submatrix<Blackbox> >  *_matFT;
		Transpose<Submatrix<Blackbox> >  *_matGT;
		Compose<Transpose<Submatrix<Blackbox> >,Submatrix<Blackbox> >   * _matFTF;
		Compose<Submatrix<Blackbox>, Transpose<Submatrix<Blackbox> > >  * _matGGT;
		Inverse<Compose<Transpose<Submatrix<Blackbox> >,Submatrix<Blackbox> > >    * _matFTFinv;
		Inverse<Compose<Submatrix<Blackbox>, Transpose<Submatrix<Blackbox> >  > >  * _matGGTinv;

		size_t     _rank;
	};

} // namespace LinBox

#endif // __LINBOX_moore_penrose_H

