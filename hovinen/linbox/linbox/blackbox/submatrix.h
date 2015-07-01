/* -*- mode: c; style: linux -*- */

/* linbox/src/library/objects/blackbox/submatrix.h
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

#ifndef __SUBMATRIX_H
#define __SUBMATRIX_H

#include "linbox/blackbox/archetype.h"
#include "linbox/vector/vector-traits.h"

// Namespace in which all LinBox library code resides
namespace LinBox
{

	/** Blackbox Submatrix.  This black box allows the extraction of a
	 * leading principal minor of an existing matrix in a black box fashion.
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
	class Submatrix : public Blackbox_archetype<Vector>
	{
	    public:

		typedef Blackbox_archetype<Vector> Blackbox;

		/** Constructor from field and dense vector of field elements.
		 * @param __BB   Black box from which to extract the submatrix
		 * @param __row  First row of the submatrix to extract (1.._BB->rowdim ())
		 * @param __col  First column of the submatrix to extract (1.._BB->coldim ())
		 */
		Submatrix (const Blackbox *BB,
			   const size_t &row,
			   const size_t &col)
		    : _BB(BB->clone ()), _row(row), _col(col)
		{
			_z.resize (_BB->_rowdim ());
			_zt.resize (_BB->_coldim ());
		}

		/** Virtual constructor.
		 * Required because constructors cannot be virtual.
		 * Make a copy of the Blackbox_archetype object.
		 * Required by abstract base class.
		 * @return pointer to new blackbox object
		 */
	        Blackbox *clone () const
	        {
			return new Submatrix (_BB, _row, _col);
		}

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
			Vector xtmp;

			copy (x.begin (), x.end (), _z.begin ());  // Copying. Yuck.
			_BB->apply (y, _z);
			y.resize (_col);
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
			Vector xtmp;

			copy (x.begin (), x.end (), _zt.begin ());  // Copying. Yuck.
			_BB->applyTranspose (y, _zt);
			y.resize (_col);
			return y;
		}

		/** Retreive _row dimensions of BlackBox matrix.
		 * This may be needed for applying preconditioners.
		 * Required by abstract base class.
		 * @return integer number of _rows of black box matrix.
		 */
		size_t _rowdim (void) const
		{
			return _row;
		}
    
		/** Retreive _column dimensions of BlackBox matrix.
		 * Required by abstract base class.
		 * @return integer number of _columns of black box matrix.
		 */
		size_t _coldim(void) const
		{
			return _col;
		}

	    private:

		Blackbox *_BB;
		size_t    _row;
		size_t    _col;

	        // Temporaries for reducing the amount of memory allocation we do
	        Vector    _z;
	        Vector    _zt;

	}; // template <Field, Vector> class Submatrix

} // namespace LinBox

#endif // __SUBMATRIX_H
