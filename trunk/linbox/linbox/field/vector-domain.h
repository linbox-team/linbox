/* -*- mode: c; style: linux -*- */

/* linbox/field/vector-domain.h
 * Copyright (C) 2001-2002 Bradford Hovinen
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

#ifndef __FIELD_VECTOR_DOMAIN_H
#define __FIELD_VECTOR_DOMAIN_H

#include <iostream>

#include "linbox/field/archetype.h"
#include "linbox/vector/vector-traits.h"

#include "linbox/util/debug.h"

#define VectorDomainType(tag1, tag2) \
	VectorDomain<Field, Vector1, Vector2, VectorCategories::tag1##VectorTag, VectorCategories::tag2##VectorTag>

#define VectorDomainSimpleType(tag) \
	VectorDomain<Field, Vector1, Vector2, VectorCategories::tag##VectorTag>

namespace LinBox
{
	/** Vector Domain.
	 * Archetype for the vector domain \Ref{LinBox}.
	 *
	 * This is a generic wrapper around classes matching the
	 * \Ref{FieldArchetype} interface. It implements vector-vector
	 * operations such as axpy, mul, and dotprod. It also contains an
	 * interface to the underlying field whereby calls simply pass
	 * through. Template specializations permit optimizations to be done on
	 * these operations based on the characteristics of the field.
	 *
	 * This class is usable by itself. Simply supply any preexisting field
	 * as a template parameter and it will work as intended, though its
	 * operation may not be fully optimized.
	 */
	template <class Field, class Vector1, class Vector2,
		  class Trait1 = VectorTraits<Vector1>::VectorCategory,
		  class Trait2 = VectorTraits<Vector2>::VectorCategory>
	class VectorDomain
	{
		public:
    
		typedef typename Field::element         element;

		/** Copy constructor.
		 * Constructs VectorDomain_archetype object by copying the domain.
		 * This is required to allow matrix domain objects to be passed
		 * by value into functions.
		 * @param  MD VectorDomain_archetype object.
		 */
		VectorDomain (const VectorDomain &MD) 
			: _F (MD._F)
			{}
    
		/** Assignment operator.
		 * Assigns VectorDomain object MD to field.
		 * @param  MD VectorDomain object.
		 */
		VectorDomain &operator = (const VectorDomain &MD)
			{ _F = MD._F; return *this; }
    
		/** Retrieve the underlying field
		 * Return a reference to the field that this matrix domain
		 * object uses
		 * @return reference to field
		 */

		const Field &field () const
			{ return _F; }
    
		//@} Object Management

		/** Vector input/output operations
		 * These routines are useful for reading and writing vectors to
		 * and from file streams. They are analagous to field read and
		 * write operations.
		 */

		//@{

		/** Print vector of field elements.
		 * This function assumes the field element has already been 
		 * constructed and initialized.
		 * In this implementation, this means for the _elem_ptr for x 
		 * exists and does not point to null.
		 * @return output stream to which field element is written.
		 * @param  os  output stream to which field element is written.
		 * @param  x   field element.
		 */
		ostream &write (ostream &os, const Vector1 &x) const;

		/** Read vector of field elements.
		 * This function assumes the field element has already been 
		 * constructed and initialized.
		 * In this implementation, this means for the _elem_ptr for x 
		 * exists and does not point to null.
		 * @return input stream from which field element is read.
		 * @param  is  input stream from which field element is read.
		 * @param  x   field element.
		 */
		istream &read (istream &is, Vector1 &x) const;
    
		//@} Input/Output Operations

		/** @name Vector arithmetic operations
		 * These routes are analogs of field arithmetic operations, but
		 * they take vectors of elements as input. Vector-vector dot
		 * product and vector-vector axpy are supported here.
		 */

		//@{

		/** Vector-vector dot product
		 * @param res element into which to store result
		 * @param v1 Input vector
		 * @param v2 Input vector
		 */
		element &dotprod (element &res, const Vector1 &v1, const Vector2 &v2) const;

		/** Scalar-vector multiplication
		 * res <- a * x
		 * @param res Vector into which to store result
		 * @param x Input vector x
		 * @param a Input element a
		 */
		Vector1 &mul (Vector1 &res, const Vector1 &x, const element &a) const;

		/** In-place scalar-vector multiplication
		 * a <- a * x
		 * @param res Vector into which to store result
		 * @param x Input vector x
		 * @param a Input element a
		 */
		Vector1 &mulin (Vector1 &x, const element &a) const;

		/** Vector axpy
		 * res <- y + a*x
		 * @param res Vector into which to store result
		 * @param y Input vector y
		 * @param a Scalar element a
		 * @param x Input vector x
		 */
		Vector1 &axpy (Vector1 &res, const Vector1 &y, const element &a, const Vector1 &x) const;

		/** Vector in-place axpy
		 * y <- y + a*x
		 * @param y Input vector y; result is stored here
		 * @param a Scalar element a
		 * @param x Input vector x
		 */
		Vector1 &axpyin (Vector1 &y, const element &a, const Vector1 &x) const;
    
		//@} Common Object Interface
    
		/** @name Implementation-Specific Methods.
		 * These methods are not required of all \Ref{LinBox Fields}
		 * and are included only for this implementation of the archetype.
		 */
		//@{

		/** Construct from a field
		 * @param F Field from which to construct
		 */
		VectorDomain (const Field &F)
			: _F (F)
		{}

		//@} Implementation-Specific Methods
    
	    private:

		Field _F;

	}; // class VectorDomain

	template <class Field, class Vector1, class Vector2>
	class VectorDomainType(Dense, Dense) 
	{
	    public:
    
		typedef typename Field::element         element;

		VectorDomain (const VectorDomain &MD) : _F (MD._F)                  {}
		VectorDomain &operator = (const VectorDomain &MD)                   { _F = MD._F; return *this; }
  		Field &field () const                                               { return _F; }
		ostream &write (ostream &os, Vector1 &x) const;
		ostream &read (ostream &os, Vector1 &x) const;
		element &dotprod (element &res, const Vector1 &v1, const Vector2 &v2) const;
		Vector1 &mul (Vector1 &res, const Vector1 &x, const element &a) const;
		Vector1 &mulin (Vector1 &x, const element &a) const;
		Vector1 &axpy (Vector1 &res, const Vector1 &y, const element &a, const Vector1 &x) const;
		Vector1 &axpyin (Vector1 &y, const element &a, const Vector1 &x) const;

		VectorDomain (const Field &F) : _F (F)                              {}

	    private:

		Field _F;
	};

	template <class Field, class Vector1, class Vector2>
	class VectorDomainType(SparseSequence, Dense) 
	{
	    public:
    
		typedef typename Field::element         element;

		VectorDomain (const VectorDomain &MD) : _F (MD._F)                  {}
		VectorDomain &operator = (const VectorDomain &MD)                   { _F = MD._F; return *this; }
  		Field &field () const                                               { return _F; }
		ostream &write (ostream &os, Vector1 &x) const;
		ostream &read (ostream &os, Vector1 &x) const;
		element &dotprod (element &res, const Vector1 &v1, const Vector2 &v2) const;
		Vector1 &mul (Vector1 &res, const Vector1 &x, const element &a) const;
		Vector1 &mulin (Vector1 &x, const element &a) const;
		Vector1 &axpy (Vector1 &res, const Vector1 &y, const element &a, const Vector1 &x) const;
		Vector1 &axpyin (Vector1 &y, const element &a, const Vector1 &x) const;

		VectorDomain (const Field &F) : _F (F)                              {}

	    private:

		Field _F;
	};

} // namespace LinBox

#include "linbox/field/vector-domain.C"

#endif // __FIELD_MATRIX_DOMAIN_H
