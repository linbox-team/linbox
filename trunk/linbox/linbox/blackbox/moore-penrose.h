/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

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

#ifndef __MOORE_PENROSE_H
#define __MOORE_PENROSE_H

#include "linbox/blackbox/archetype.h"
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

	/** Blackbox MoorePenrose.
	 *
	 * Given an arbitrary matrix in black box representation, this black box
	 * represents the Moore-Penrose inverse of the matrix.
	 *
	 * This implementation assumes that A already has a nonsingular
	 * principal r x r minor. It is the caller's responsibility to ensure
	 * that that condition holds.
	 */
	template <class Field, class Vector>
	class MoorePenrose : public BlackboxArchetype<Vector>
	{
		public:

		typedef BlackboxArchetype<Vector> Blackbox;

		/** Constructor from field and dense vector of field elements.
		 * @param BB   Black box from which to extract the submatrix
		 * @param row  First row of the submatrix to extract (1.._BB->rowdim ())
		 * @param col  First column of the submatrix to extract (1.._BB->coldim ())
		 * @param rowdim Row dimension
		 * @param coldim Column dimension
		 */
		MoorePenrose (Field &F, const Blackbox *A, size_t rank)
			: _A (A->clone ()), _rank (rank)
		{
			_B1 = new Submatrix<Vector> (_A, 0, 0, rank, rank);
			_F = new Submatrix<Vector> (_A, 0, 0, _A->rowdim (), rank);
			_G = new Submatrix<Vector> (_A, 0, 0, rank, _A->coldim ());
			_FT = new Transpose<Vector> (_F);
			_GT = new Transpose<Vector> (_G);
			_FTF = new Compose<Vector> (_FT, _F);
			_GGT = new Compose<Vector> (_G, _GT);
			_FTFinv = new Inverse<Field, Vector> (F, _FTF);
			_GGTinv = new Inverse<Field, Vector> (F, _GGT);
		}

		/** Copy constructor
		 */
		MoorePenrose (const MoorePenrose &A)
			: _A (A._A->clone ()),
			_B1 (A._B1->clone ()),
			_F (A._F->clone ()),
			_G (A._G->clone ()),
			_FT (A._FT->clone ()),
			_GT (A._GT->clone ()),
			_FTF (A._FTF->clone ()),
			_GGT (A._GGT->clone ()),
			_FTFinv (A._FTFinv->clone ()),
			_GGTinv (A._GGTinv->clone ()),
			_rank (A._rank)
			{}

		/** Destructor
		 */
		virtual ~MoorePenrose ()
		{
			delete _GGTinv;
			delete _FTFinv;
			delete _GGT;
			delete _FTF;
			delete _GT;
			delete _FT;
			delete _G;
			delete _F;
			delete _A;
			delete _B1;
		}

		/** Virtual constructor.
		 * Required because constructors cannot be virtual.
		 * Make a copy of the BlackboxArchetype object.
		 * Required by abstract base class.
		 * @return pointer to new blackbox object
		 */
		Blackbox *clone () const
			{ return new MoorePenrose (*this); }

		/** Application of BlackBox matrix.
		 * y= A*x.
		 * Requires one vector conforming to the \Ref{LinBox}
		 * vector {@link Archetypes archetype}.
		 * Required by abstract base class.
		 * @return reference to vector y containing output.
		 * @param  x constant reference to vector to contain input
		 */
	        Vector& apply (Vector &y, const Vector& x) const
		{
			Vector _z1 (_rank);
			Vector _z2 (_rank);

			_F->applyTranspose (_z1, x);
			_FTFinv->apply (_z2, _z1);
			_B1->apply (_z1, _z2);
			_GGTinv->apply (_z2, _z1);
			_G->applyTranspose (y, _z2);

			return y;
		}

		/** Application of BlackBox matrix transpose.
		 * y= transpose(A)*x.
		 * Requires one vector conforming to the \Ref{LinBox}
		 * vector {@link Archetypes archetype}.
		 * Required by abstract base class.
		 * @return reference to vector y containing output.
		 * @param  x constant reference to vector to contain input
		 */
		Vector& applyTranspose (Vector &y, const Vector& x) const
		{
			Vector _z1 (_rank);
			Vector _z2 (_rank);

			_G->apply (_z1, x);
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

	    private:

		Blackbox  *_A;
		Blackbox  *_B1;
		Blackbox  *_F;
		Blackbox  *_G;
		Blackbox  *_FT;
		Blackbox  *_GT;
		Blackbox  *_FTF;
		Blackbox  *_GGT;
		Blackbox  *_FTFinv;
		Blackbox  *_GGTinv;

		size_t     _rank;
	}; // template <Vector> class MoorePenrose

} // namespace LinBox

#endif // __MOORE_PENROSE_H
