/* -*- mode: C++; style: linux -*- */

/* linbox/field/vector-domain.h
 * Copyright (C) 2001-2002 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
 * Modified by Zhendong Wan <wan@mail.eecis.udel.edu>
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

#ifndef __DENSE_VECTOR_DOMAIN_H
#define __DENSE_VECTOR_DOMAIN_H

#include <iostream>

namespace LinBox
{
	/**Dense Vector Domain.
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
        template<class Field>
	class DenseVectorDomain
	{
		public:
    
		typedef typename Field::Element         Element;
		
		/** Construct from a field
		 * @param F Field from which to construct
		 */
		explicit DenseVectorDomain(const Field &F =Field())
			: _F (F)
		  {}

		/** Copy constructor.
		 * Constructs VectorDomain_archetype object by copying the domain.
		 * This is required to allow matrix domain objects to be passed
		 * by value into functions.
		 * @param  MD VectorDomain_archetype object.
		 */
		DenseVectorDomain (const DenseVectorDomain &MD) 
			: _F (MD._F)
		  {}
    
		/** Assignment operator.
		 * Assigns VectorDomain object MD to field.
		 * @param  MD VectorDomain object.
		 */
		DenseVectorDomain &operator = (const DenseVectorDomain &MD)
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
		template<class Vect>
		ostream &write (ostream &os, const Vect &x) const;

		/** Read vector of field elements.
		 * This function assumes the field element has already been 
		 * constructed and initialized.
		 * In this implementation, this means for the _elem_ptr for x 
		 * exists and does not point to null.
		 * @return input stream from which field element is read.
		 * @param  is  input stream from which field element is read.
		 * @param  x   field element.
		 */
		template<class Vector1>
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
		template<class Vector1, class Vector2>
		Element& dotproduct (Element &res, const Vector1 &v1, const Vector2 &v2) const;

		/** Scalar-vector multiplication
		 * res <- a * x
		 * @param res Vector into which to store result
		 * @param x Input vector x
		 * @param a Input element a
		 */
		template<class Vector1, class Vector2>
		Vector1 &mul (Vector1 &res, const Vector2 &x, const Element &a) const;

		/** In-place scalar-vector multiplication
		 * a <- a * x
		 * @param res Vector into which to store result
		 * @param x Input vector x
		 * @param a Input element a
		 */
		template<class Vector1>
		Vector1 &mulin (Vector1 &x, const Element &a) const;

		/** Vector axpy
		 * res <- y + a*x
		 * @param res Vector into which to store result
		 * @param y Input vector y
		 * @param a Scalar element a
		 * @param x Input vector x
		 */
		template<class Vector1, class Vector2, class Vector3>
		Vector1& axpy (Vector1 &res, const Element &a, const Vector2 &x, const Vector3 &y) const;

		/** Vector in-place axpy
		 * y <- y + a*x
		 * @param y Input vector y; result is stored here
		 * @param a Scalar element a
		 * @param x Input vector x
		 */
		template<class Vector1, class Vector2>
		Vector1 &axpyin (Vector1 &y, const Element &a, const Vector2 &x) const;
    
		//@} Common Object Interface
    
		/** @name Implementation-Specific Methods.
		 * These methods are not required of all \Ref{LinBox Fields}
		 * and are included only for this implementation of the archetype.
		 */
		//@{


		//@} Implementation-Specific Methods
    
	    private:

		Field _F;

	}; // class VectorDomain

} // namespace LinBox
#include "dense-vector-domain.C"

#endif // __FIELD_MATRIX_DOMAIN_H
