/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/blackbox/archetype.h
 * Copyright (C) 1999-2001 William J Turner,
 *               2001 Bradford Hovinen
 *
 * Written by William J Turner <wjturner@math.ncsu.edu>,
 *            Bradford Hovinen <hovinen@cis.udel.edu>
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

#ifndef __BLACKBOX_ARCHETYPE_H
#define __BLACKBOX_ARCHETYPE_H

#include <cstddef>

#include "linbox/util/error.h"
#include "linbox-config.h"

#ifdef __LINBOX_XMLENABLED

#include "linbox/util/xml/linbox-writer.h"
#include <iostream>

using std::ostream;

#endif

namespace LinBox
{

	/*-  put this info in archetypes dxx file -bds
	 * Found in file \URL{linbox/blackbox/archetype.h}.
	 * Base class from which derived concrete blackbox classes.
	 * Unlike the LinBox field common object interface,
	 * the common object interface for LinBox BlackBoxes does not require
	 * the copy constructor.  All object management routines are given 
	 * through virtual clone and killclone methods.  This allows the base
	 * object to be the archetype, which is not possible for LinBox fields. 
	 *
	 * In general, there are three uses of archetypes:
	 * \begin{enumerate}
	 * \item To define the interface, i.e., document what an
	 *       explicitly designed field class must have.  This is
	 *       useful, for instance, in building a wrapper or adaptor
	 *       to an existing library.
	 * \item To distribute compiled code and ease the testing of
	 *       library components.
	 * \item To control code bloat.
	 * \end{enumerate}
	 * Because of their use of virtual member funtions, these archetypes can be 
	 * inefficient.
	 *
	 * @param Vector \Ref{LinBox} dense or sparse vector of field elements
	 */

	/*- 
	@memo BlackBox base class and archetype 
	@doc 
	This archetype is an abstract base class for blackbox matrix classes.
	The key member functions are {\tt apply, applyTranspose, rodwim, coldim}.
	They are pure virtual, and hence are implemented in each child class.  
	
	Concrete classes inheriting from the archetype
	use a variety of representation schemes for matrices internally. 
	All provide the blackbox interface described here and can be used 
	interchangably in blackbox algorithms.
	Some also implement a dense matrix or sparse matrix interface to support elimination
	techniques.  Each has unique constructor(s) reflecting it's specific scheme for representing
	a linear operator.

	Algorithms written with a Blackbox template parameter 
	may be compiled against any of these classes specifically or may be separately compiled
	against the archetype.  Algorithms may also be written with a BlackboxArchetype parameter
	and then called with an instance of a concrete blackbox class. 
	In contrast with the situation for \Ref{Field}s there is 
	negligible performance cost for separate compilation here.
	
	{\bf Template Parameter:} Vector - A type meeting the LinBox \Ref{VectorArchetype} interface.
	Vectors of this type are the normal arguments to {\tt apply} and {\tt applyTranspose}.
	
	@see \Ref{../archetypes} for general discussion of LinBox archetypes.
	*/
	template <class Vector>
	class BlackboxArchetype 
	{
	public:

		/// Deallocates the memory used for the matrix representation.
		virtual ~BlackboxArchetype (void) {}

		/*- Serves in place of copy constructor.
		 * Required because constructors cannot be virtual.
		 * Make a copy of the BlackboxArchetype object.
		 * @return pointer to new blackbox object
		 */
		virtual BlackboxArchetype* clone () const = 0;
		// Should we have clone conform more to copy construction?
		// clone(A) = make this a copy of A. -bds

	// apply variants //

		/*- Matrix vector product. y := Ax.
		The vector x must be of size A.coldim(), where A is this blackbox.
		On entry to apply, the vector y must be of size A.rowdim().
		Neither vector has it's size or capacity modified by apply.  Apply is not
		responsible for the validity of the sizes, which may or may not be checked.
		The two vectors may not overlap in memory.
		@param y it's entries are set and a reference to it is also returned to allow for 
		use in nested expressions.
		@param x it's entries are the input data.
		*/
		virtual Vector &apply (Vector &y, const Vector &x) const = 0;

		/*- y := Ax, using a handle for ...
		The handle serves as "protection from the future".  The idea is that the handle
		could allow the blackbox to operate more as a pure container, with the field
		(or other functionality such as dot product) provided through the handle.

		However, there are no known current uses (2003 june).  
		*/
		virtual Vector &apply (Vector &y, const Vector &x, void *handle) const 
		{
			if (handle == 0)
				return apply (y, x);
			else
				throw LinboxError ("no handle handled in this blackbox");
		}

		/*- In-place application of BlackBox matrix.  x := Ax.
		This matrix must be square and x a vector of matching size.

		@param x is modified.  On exit the entries are those of Ax', where x' is the value of x on entry. and a reference to x is returned for possible use in nested expressions.
		*/
		/* Delete - this should only exist for special classes such as Diagonal
		virtual Vector& applyIn (Vector &x) const 
		{
			Vector y (x);
			return apply (x, y);
		}
		*/

		/*- Return new vector Ax.
		 * Requires one vector conforming to the \Ref{LinBox}
		 * vector {@link Archetypes archetype}.
		 * Purely virtual.
		 * @return reference to vector y containing output.
		 * @param  x constant reference to vector to contain input
		 */
		/* Delete -- not part of current interface spec, and not used anywhere so far as I know.  -bds
		virtual Vector& apply (const Vector &x) const 
		{
			Vector *y = new Vector(rowdim());
			return apply (*y, x);
		}
		*/

	// applyTranspose variants //

		/*- Application of BlackBox matrix transpose. y := xA.
		 * (Or taking the column vector view (y = A^T x.)
		 * Requires two vectors conforming to the \Ref{LinBox}
		 * vector {@link Archetypes archetype}.
		 * Virtual.
		 * @return reference to vector containing output.
		 * @param  y reference to vector to contain output
		 * @param  x constant reference to vector to contain input
		 */
		/*- Matrix transpose times vector product. y := xA.
		(or in column vector view: y := A^T x.)
		The vector x must be of size A.rowdim(), where A is this blackbox.
		On entry to apply, the vector y must be of size A.coldim().
		Neither vector has it's size or capacity modified by applyTranspose.  ApplyTranspose is not
		responsible for the validity of the sizes, which may or may not be checked.
		The two vectors may not overlap in memory.
		@param y it's entries are set and a reference to it is also returned to allow for 
		use in nested expressions.
		@param x it's entries are the input data.
		*/
		virtual Vector &applyTranspose (Vector &y, const Vector &x) const = 0;

		/*- y := xA, using a handle for ...
		The handle serves as "protection from the future".  The idea is that the handle
		could allow the blackbox to operate more as a pure container, with the field
		(or other functionality such as dot product) provided through the handle.

		However, there are no known current uses (2003 june).  
		*/
		virtual Vector &applyTranspose (Vector &y, const Vector &x, void *handle) const 
		{
			if (handle == 0)
				return applyTranspose (y, x);
			else
				throw LinboxError ("no handle handled in this blackbox");
		}



		/*- Application of BlackBox matrix Transpose. new y := xA.
		 * return transpose (A)*x.
		 * Requires one vector conforming to the \Ref{LinBox}
		 * vector {@link Archetypes archetype}.
		 * Purely virtual.
		 * @return reference to vector y containing output.
		 * @param  x constant reference to vector to contain input
		 */
		/* Delete -- not part of current interface spec, and not used anywhere so far as I know.  -bds
		virtual Vector& applyTranspose (const Vector &x) const 
		{
			Vector *y = new Vector;
			return applyTranspose (*y, x);
		}
		*/

		/*- In-place application of BlackBox matrix tranpose.  x := A^T x.
		This matrix must be square and x a vector of length equal to it's dimension.

		Depreciated:  There are no known uses of in-place apply so far and no blackbox class
		 provides a more efficient impl than the archetype.

		 * @param x is mutated to the output value and a reference to it is returned.
		 */
		 /* now have inplace forms only in special cases like diagonal
		virtual Vector& applyTransposeIn (Vector &x) const 
		{
			Vector y (x);
			return applyTranspose (x, y);
		}
		*/

		/*- Returns the number of rows of the matrix.
		This may be zero or greater.  Currently matrix size beyond size_t is not supported.
		*/
		virtual size_t rowdim (void) const = 0;

		/*- Returns the number of columns of the matrix.
		This may be zero or greater.  Currently matrix size beyond size_t is not supported.
		 */
		virtual size_t coldim (void) const = 0;

#ifdef __LINBOX_XMLENABLED
		virtual ostream &write(ostream &) const = 0;
		virtual bool toTag(Writer &W) const = 0;
#endif


		/*
	    protected:

		/- Default constructor.
		 * Developer: The default constructor is required by derived classes, but because 
		 * this class is abstract - it contains purely virtual functions - it 
		 * should never be called without a derived class.
		 -/
		BlackboxArchetype (void) {}
		*/

	}; // BlackBox Archetype

} // namespace LinBox

#endif // __BLACKBOX_ARCHETYPE_H










