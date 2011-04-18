/* linbox/blackbox/moore-penrose.h
 * Copyright (C) 2001 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#ifndef __LINBOX_moore_penrose_H
#define __LINBOX_moore_penrose_H

#include <linbox/blackbox/blackbox-interface.h>
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
	class MoorePenrose : public BlackboxInterface
	{
	    public:

		typedef typename Blackbox::Field Field;
		typedef typename Blackbox::Element Element;
            template<typename _Tp1>
            struct rebind
            { typedef MoorePenrose<typename Blackbox::template rebind<_Tp1>::other> other; };


		/** Constructor from field and dense vector of field elements.
		 * @param BB   Black box from which to extract the submatrix
		 * @param row  First row of the submatrix to extract (1.._BB->rowdim ())
		 * @param col  First column of the submatrix to extract (1.._BB->coldim ())
		 * @param rowdim Row dimension
		 * @param coldim Column dimension
		 */
		MoorePenrose (const Blackbox *A, size_t rank)
			: _A (A), _rank (rank)
		{
			_B1 = new Submatrix<Blackbox> (_A, 0, 0, rank, rank);
			_F = new Submatrix<Blackbox> (_A, 0, 0, _A->rowdim (), rank);
			_GG = new Submatrix<Blackbox> (_A, 0, 0, rank, _A->coldim ());
			_FT = new Transpose<Submatrix<Blackbox> > (_F);
			_GT = new Transpose<Submatrix<Blackbox> > (_GG);
			_FTF = new Compose<Transpose<Submatrix<Blackbox> >,Submatrix<Blackbox> > (_FT, _F);
			_GGT = new Compose<Submatrix<Blackbox>, Transpose<Submatrix<Blackbox> > > (_GG, _GT);
			_FTFinv = new Inverse<Compose<Transpose<Submatrix<Blackbox> >,Submatrix<Blackbox> > > ( _FTF);
			_GGTinv = new Inverse<Compose<Submatrix<Blackbox>, Transpose<Submatrix<Blackbox> > > > ( _GGT);
		}

		/** Copy constructor
		 */
		MoorePenrose (const MoorePenrose &A)
			: _A (A._A),
			_B1 (A._B1),
			_F (A._F),
			_GG (A._GG),
			_FT (A._FT),
			_GT (A._GT),
			_FTF (A._FTF),
			_GGT (A._GGT),
			_FTFinv (A._FTFinv),
			_GGTinv (A._GGTinv),
			_rank (A._rank)
			{}

		/** Destructor
		 */
		~MoorePenrose ()
		{
			delete _GGTinv;
			delete _FTFinv;
			delete _GGT;
			delete _FTF;
			delete _GT;
			delete _FT;
			delete _GG;
			delete _F;
			delete _B1;
		}

	

		/** Application of BlackBox matrix.
		 * y= A*x.
		 * Requires one vector conforming to the \ref{LinBox}
		 * vector {@link Archetypes archetype}.
		 * Required by abstract base class.
		 * @return reference to vector y containing output.
		 * @param  x constant reference to vector to contain input
		 */
		template <class OutVector, class InVector>
	        OutVector& apply (OutVector &y, const InVector& x) const
		{
			InVector _z1 (_rank);
			InVector _z2 (_rank);

			_F->applyTranspose (_z1, x);
			_FTFinv->apply (_z2, _z1);
			_B1->apply (_z1, _z2);
			_GGTinv->apply (_z2, _z1);
			_GG->applyTranspose (y, _z2);

			return y;
		}

		/** Application of BlackBox matrix transpose.
		 * y= transpose(A)*x.
		 * Requires one vector conforming to the \ref{LinBox}
		 * vector {@link Archetypes archetype}.
		 * Required by abstract base class.
		 * @return reference to vector y containing output.
		 * @param  x constant reference to vector to contain input
		 */
		template <class OutVector, class InVector>
		OutVector& applyTranspose (OutVector &y, const InVector& x) const
		{
			InVector _z1 (_rank);
			InVector _z2 (_rank);

			_GG->apply (_z1, x);
			_GGTinv->applyTranspose (_z2, _z1);
			_B1->applyTranspose (_z1, _z2);
			_FTFinv->applyTranspose (_z2, _z1);
			_F->apply (y, _z2);

			return y;
		}

		/** Retreive _row dimensions of BlackBox matrix.
		 * This may be needed for applying preconditioners.
		 * Required by abstract base class.
		 * @return integer number of _rows of black box matrix.
		 */
		size_t rowdim (void) const
			{ return _A->coldim (); }
    
		/** Retreive _column dimensions of BlackBox matrix.
		 * Required by abstract base class.
		 * @return integer number of _columns of black box matrix.
		 */
		size_t coldim (void) const
			{ return _A->rowdim (); }

		const Field& field() { return _A -> field(); }

	    private:

		const Blackbox  *_A;
		Submatrix<Blackbox>  *_B1;
		Submatrix<Blackbox>  *_F;
		Submatrix<Blackbox>  *_GG;
		Transpose<Submatrix<Blackbox> >  *_FT;
		Transpose<Submatrix<Blackbox> >  *_GT;
		Compose<Transpose<Submatrix<Blackbox> >,Submatrix<Blackbox> >  *_FTF;
		Compose<Submatrix<Blackbox>, Transpose<Submatrix<Blackbox> > >  *_GGT;
		Inverse<Compose<Transpose<Submatrix<Blackbox> >,Submatrix<Blackbox> > >  *_FTFinv;
		Inverse<Compose<Submatrix<Blackbox>, Transpose<Submatrix<Blackbox> >  > >  *_GGTinv;

		size_t     _rank;
	}; 

} // namespace LinBox

#endif // __LINBOX_moore_penrose_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
