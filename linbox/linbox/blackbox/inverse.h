/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/blackbox/inverse.h
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#ifndef __INVERSE_H
#define __INVERSE_H

#include "linbox/blackbox/archetype.h"
#include "linbox/blackbox/transpose.h"
#include "linbox/field/vector-domain.h"
#include "linbox/solutions/minpoly.h"
#include "linbox/vector/vector-traits.h"

// Namespace in which all LinBox library code resides
namespace LinBox
{

	/** Blackbox Inverse.  This represents the inverse of a nonsingular
	 * matrix.
	 *
	 * The matrix itself is not stored in memory.  Rather, its apply
	 * methods use a vector of {@link Fields field} elements, which are 
	 * used to "multiply" the matrix to a vector.
	 * 
	 * This class has three template parameters.  The first is the field in 
	 * which the arithmetic is to be done.  The second is the type of 
	 * \Ref{LinBox} vector to which to apply the matrix.  The 
	 * third is chosen be default to be the \Ref{LinBox} vector trait
	 * of the vector.  This class is then specialized for dense and sparse 
	 * vectors.
	 *
	 * @param Field \Ref{LinBox} field
	 * @param Vector \Ref{LinBox} dense or sparse vector of field elements
	 * @param Trait  Marker whether to use dense or sparse LinBox vector 
	 *               implementation.  This is chosen by a default parameter 
	 *               and partial template specialization.  */
	template <class Field, class Vector>
	class Inverse : public BlackboxArchetype<Vector>
	{
	    public:

		typedef BlackboxArchetype<Vector> Blackbox;
		typedef typename Field::Element   Element;
		typedef std::vector<Element>      Polynomial;

		/** Constructor from field and dense vector of field elements.
		 * @param __BB   Black box of which to get the inverse
		 */
		Inverse (const Field &F, const Blackbox *BB)
		    : _F (F), _VD (F), _BB (BB->clone ())
		{
			linbox_check (BB->rowdim () == BB->coldim ());

			_minpoly.clear ();
			_transposeMinpoly.clear ();

			VectorWrapper::ensureDim (_z, _BB->coldim ());
		}

		/** Copy constructor, so that we don't have to recompute the
		 * minimal polynomial every time this black box is used inside
		 * another black box
		 */
		Inverse (const Inverse &BB)
		    : _F (BB._F), _VD (BB._F), _BB (BB._BB->clone ()), _minpoly (BB._minpoly)
		{
			_z.resize (_BB->coldim ());
		}

		/** Virtual constructor.
		 * Required because constructors cannot be virtual.
		 * Make a copy of the BlackboxArchetype object.
		 * Required by abstract base class.
		 * @return pointer to new blackbox object
		 */
	        Blackbox *clone () const
	        {
			return new Inverse (*this);
		}

		/** Application of BlackBox matrix.
		 * y= A*x.
		 * Requires one vector conforming to the \Ref{LinBox}
		 * vector {@link Archetypes archetype}.
		 * Required by abstract base class.
		 * @return reference to vector y containing output.
		 * @param  y reference to vector into which to store the result
		 * @param  x constant reference to vector to contain input
		 */
	        Vector& apply (Vector &y, const Vector& x) const
	        {
			int i;

			if (_minpoly.empty ()) {
				Polynomial _mp1;
				Element a0;

				minpoly (_mp1, *_BB, _F);

				_minpoly.resize (_mp1.size () - 1);

				_F.inv (a0, _mp1[0]);
				_F.negin (a0);

				for (i = 1; i < (int) _mp1.size (); i++)
					_F.mul (_minpoly[i-1], _mp1[i], a0);
			}

			int n = _minpoly.size () - 1;

			_VD.mul (y, x, _minpoly[n]);

			for (i = n - 1; i >= 0; i--) {
				_BB->apply (_z, y);
				_VD.axpy (y, _minpoly[i], x, _z);
			}

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
			int i;

			if (_transposeMinpoly.empty ()) {
				Polynomial _mp1;
				Element a0;

				Transpose<Vector> BBT (_BB);

				minpoly (_mp1, BBT, _F);

				_transposeMinpoly.resize (_mp1.size () - 1);

				_F.inv (a0, _mp1[0]);
				_F.negin (a0);

				for (i = 1; i < (int) _mp1.size (); i++)
					_F.mul (_transposeMinpoly[i-1], _mp1[i], a0);
			}

			int n = _transposeMinpoly.size () - 1;

			_VD.mul (y, x, _transposeMinpoly[n]);

			for (i = n - 1; i >= 0; i--) {
				_BB->applyTranspose (_z, y);
				_VD.axpy (y, _transposeMinpoly[i], x, _z);
			}

			return y;
		}

		/** Retreive row dimensions of BlackBox matrix.
		 * This may be needed for applying preconditioners.
		 * Required by abstract base class.
		 * @return integer number of rows of black box matrix.
		 */
		size_t rowdim (void) const
		{
			return _BB->rowdim ();
		}
    
		/** Retreive column dimensions of BlackBox matrix.
		 * Required by abstract base class.
		 * @return integer number of columns of black box matrix.
		 */
		size_t coldim(void) const
		{
			return _BB->coldim ();
		}

	    private:

		const Field               &_F;
		const VectorDomain<Field>  _VD;
		Blackbox                  *_BB;

		mutable Polynomial         _minpoly;
		mutable Polynomial         _transposeMinpoly;

		// Temporary for reducing necessary memory allocation
		mutable Vector             _z;

	}; // template <Field, Vector> class Inverse

} // namespace LinBox

#endif // __INVERSE_H
