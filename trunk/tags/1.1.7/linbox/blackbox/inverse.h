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

#ifndef __LINBOX_inverse_H
#define __LINBOX_inverse_H

#include <linbox/blackbox/blackbox-interface.h>
#include "linbox/blackbox/transpose.h"
#include "linbox/vector/vector-domain.h"
#include "linbox/solutions/minpoly.h"
#include "linbox/vector/vector-traits.h"

// Namespace in which all LinBox library code resides
namespace LinBox
{

	/** \brief A Blackbox for the inverse.  Not efficient if many applications are used.
	 * \ingroup blackbox
	 *
	 * The matrix itself is not stored in memory.  Rather, its apply
	 * methods use a vector of {@link Fields field} elements, which are 
	 * used to "multiply" the matrix to a vector.
	 * 
	 * This class has three template parameters.  The first is the field in 
	 * which the arithmetic is to be done.  The second is the type of 
	 * \ref{LinBox} vector to which to apply the matrix.  The 
	 * third is chosen be default to be the \ref{LinBox} vector trait
	 * of the vector.  This class is then specialized for dense and sparse 
	 * vectors.
	 *
	 * @param Field \ref{LinBox} field
	 * @param Vector \ref{LinBox} dense or sparse vector of field elements
	 * @param Trait  Marker whether to use dense or sparse LinBox vector 
	 *               implementation.  This is chosen by a default parameter 
	 *               and partial template specialization.  */
	template <class Blackbox>
	class Inverse : public BlackboxInterface
	{
	    public:

		typedef typename Blackbox::Field Field;
		typedef typename Field::Element   Element;
		typedef std::vector<Element>      Polynomial;

		/** Constructor from field and dense vector of field elements.
		 * @param __BB   Black box of which to get the inverse
		 */
		Inverse (const Blackbox *BB)
		    :  _VD (BB->field()), _BB (BB)
		{
			linbox_check ((BB->rowdim ()) == (BB->coldim ()));

			_minpoly.clear ();
			_transposeMinpoly.clear ();

			VectorWrapper::ensureDim (_z, _BB->coldim ());
		}

		/** Copy constructor, so that we don't have to recompute the
		 * minimal polynomial every time this black box is used inside
		 * another black box
		 */
		Inverse (const Inverse &BB)
		    : _VD (BB->field()), _BB (BB._BB), _minpoly (BB._minpoly)
		{
			_z.resize (_BB->coldim ());
		}


		/** Application of BlackBox matrix.
		 * y= A*x.
		 * Requires one vector conforming to the \ref{LinBox}
		 * vector {@link Archetypes archetype}.
		 * Required by abstract base class.
		 * @return reference to vector y containing output.
		 * @param  y reference to vector into which to store the result
		 * @param  x constant reference to vector to contain input
		 */
		template<class OutVector, class InVector>
	        OutVector& apply (OutVector &y, const InVector& x) const
	        {
			int i;

			if (_minpoly.empty ()) {
				Polynomial _mp1;
				Element a0;

				minpoly (_mp1, *_BB);

				_minpoly.resize (_mp1.size () - 1);

				field().inv (a0, _mp1[0]);
				field().negin (a0);

				for (i = 1; i < (int) _mp1.size (); i++)
					field().mul (_minpoly[i-1], _mp1[i], a0);
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
		 * Requires one vector conforming to the \ref{LinBox}
		 * vector {@link Archetypes archetype}.
		 * Required by abstract base class.
		 * @return reference to vector y containing output.
		 * @param  x constant reference to vector to contain input
		 */
		template<class OutVector, class InVector>
		OutVector& applyTranspose (OutVector &y, const InVector& x) const
		{
			int i;

			if (_transposeMinpoly.empty ()) {
				Polynomial _mp1;
				Element a0;

				Transpose<Blackbox> BBT (_BB);

				minpoly (_mp1, BBT);

				_transposeMinpoly.resize (_mp1.size () - 1);

				field().inv (a0, _mp1[0]);
				field().negin (a0);

				for (i = 1; i < (int) _mp1.size (); i++)
					field().mul (_transposeMinpoly[i-1], _mp1[i], a0);
			}

			int n = _transposeMinpoly.size () - 1;

			_VD.mul (y, x, _transposeMinpoly[n]);

			for (i = n - 1; i >= 0; i--) {
				_BB->applyTranspose (_z, y);
				_VD.axpy (y, _transposeMinpoly[i], x, _z);
			}

			return y;
		}

            template<typename _Tp1>
            struct rebind
            { typedef Inverse<typename Blackbox::template rebind<_Tp1>::other> other; };


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

		const Field& field() const { return _BB->field();}
	    private:

		const VectorDomain<Field>  _VD;
		const Blackbox             *_BB;

		mutable Polynomial         _minpoly;
		mutable Polynomial         _transposeMinpoly;

		// Temporary for reducing necessary memory allocation
		mutable Polynomial         _z;

	};

} // namespace LinBox

#endif // __LINBOX_inverse_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
