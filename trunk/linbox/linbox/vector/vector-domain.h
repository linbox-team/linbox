/* -*- mode: c; style: linux -*- */

/* linbox/field/vector-domain.h
 * Copyright (C) 2001-2002 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * ------------------------------------
 * 2002-04-10 Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * There are now two vector domain types: VectorDomainBase and
 * VectorDomain. VectorDomainBase, which in principle can be used independently,
 * contains all functions that require only one vector type (such as axpy, mul,
 * read, and write). VectorDomain inherits VectorDomainBase and implements
 * dotprod, which requires two vector types.
 * 
 * ------------------------------------
 * Modified by Dmitriy Morozov <linbox@foxcub.org>. May 27, 2002.
 *
 * Added the modifications for categories and vector traits that were designed
 * at the Rootbeer meeting. Added parametrization of VectorTags by VectorTraits.
 * 
 * ------------------------------------
 * 2002-06-04 Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * Reverted change of 2002-04-10, reintegrating VectorDomain and
 * VectorDomainBase. Now using template specialization on the functions, now
 * that I know how to do it.
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __FIELD_VECTOR_DOMAIN_H
#define __FIELD_VECTOR_DOMAIN_H

#include <iostream>

#include "linbox/field/archetype.h"
#include "linbox/vector/vector-traits.h"

#include "linbox/util/debug.h"

namespace LinBox
{
	/** Vector Domain
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
	template <class Field>
	class VectorDomain
	{
		public:
    
		typedef typename Field::Element         Element;

		/** Copy constructor.
		 * Constructs VectorDomain_archetype object by copying the domain.
		 * This is required to allow matrix domain objects to be passed
		 * by value into functions.
		 * @param  MD VectorDomain_archetype object.
		 */
		VectorDomain (const VectorDomain &VD)
			: _F (VD._F) 
			{}
    
		/** Assignment operator.
		 * Assigns VectorDomain object MD to field.
		 * @param  MD VectorDomain object.
		 */
		VectorDomain &operator = (const VectorDomain &VD)
			{ _F = VD._F; }
    
		/** Retrieve the underlying field
		 * Return a reference to the field that this matrix domain
		 * object uses
		 * @return reference to field
		 */

		const Field &field () const;
    
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
		template <class Vector>
		inline ostream &write (ostream &os, const Vector &x) const
			{ writeSpecialized (os, x, VectorTraits<Vector>::VectorCategory ()); }

		/** Read vector of field elements.
		 * This function assumes the field element has already been 
		 * constructed and initialized.
		 * In this implementation, this means for the _elem_ptr for x 
		 * exists and does not point to null.
		 * @return input stream from which field element is read.
		 * @param  is  input stream from which field element is read.
		 * @param  x   field element.
		 */
		template <class Vector>
		inline istream &read (istream &is, Vector &x) const
    			{ readSpecialized (is, x, VectorTraits<Vector>::VectorCategory ()); }

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
		template <class Vector1, class Vector2>
		inline Element &dot (Element &res, const Vector1 &v1, const Vector2 &v2) const
			{ dotSpecialized (res, v1, v2,
					  VectorTraits<Vector1>::VectorCategory (),
					  VectorTraits<Vector2>::VectorCategory ()); }

		/* Alias for the above, to avoid source incompatibility */
		template <class Vector1, class Vector2>
		inline Element &dotprod (Element &res, const Vector1 &v1, const Vector2 &v2) const
			{ dot (res, v1, v2); }

		/** Scalar-vector multiplication
		 * res <- a * x
		 * @param res Vector into which to store result
		 * @param x Input vector x
		 * @param a Input element a
		 */
		template <class Vector>
		inline Vector &mul (Vector &res, const Vector &x, const Element &a) const
			{ mulSpecialized (res, x, a, VectorTraits<Vector>::VectorCategory ()); }

		/** In-place scalar-vector multiplication
		 * a <- a * x
		 * @param res Vector into which to store result
		 * @param x Input vector x
		 * @param a Input element a
		 */
		template <class Vector>
		inline Vector &mulin (Vector &x, const Element &a) const
			{ mulinSpecialized (x, a, VectorTraits<Vector>::VectorCategory ()); }

		/** Vector axpy
		 * res <- y + a*x
		 * @param res Vector into which to store result
		 * @param y Input vector y
		 * @param a Scalar element a
		 * @param x Input vector x
		 */
		template <class Vector>
		inline Vector &axpy (Vector &res, const Vector &y, const Element &a, const Vector &x) const
			{ axpySpecialized (res, y, a, x, VectorTraits<Vector>::VectorCategory ()); }

		/** Vector in-place axpy
		 * y <- y + a*x
		 * @param y Input vector y; result is stored here
		 * @param a Scalar element a
		 * @param x Input vector x
		 */
		template <class Vector>
		inline Vector &axpyin (Vector &y, const Element &a, const Vector &x) const
			{ axpyinSpecialized (y, a, x, VectorTraits<Vector>::VectorCategory ()); }
    
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
    
	    protected:

		Field _F;

		// Specialized function implementations
		template <class Vector, class Trait>
		ostream &writeSpecialized (ostream &os, const Vector &x,
					   VectorCategories::DenseVectorTag<Trait> tag) const;
		template <class Vector, class Trait>
		ostream &writeSpecialized (ostream &os, const Vector &x,
					   VectorCategories::SparseSequenceVectorTag<Trait> tag) const;
		template <class Vector, class Trait>
		ostream &writeSpecialized (ostream &os, const Vector &x,
					   VectorCategories::SparseAssociativeVectorTag<Trait> tag) const;

		template <class Vector, class Trait>
		istream &readSpecialized (istream &is, const Vector &x,
					  VectorCategories::DenseVectorTag<Trait> tag) const;
		template <class Vector, class Trait>
		istream &readSpecialized (istream &is, const Vector &x,
					  VectorCategories::SparseSequenceVectorTag<Trait> tag) const;
		template <class Vector, class Trait>
		istream &readSpecialized (istream &is, const Vector &x,
					  VectorCategories::SparseAssociativeVectorTag<Trait> tag) const;

		template <class Vector1, class Trait1, class Vector2, class Trait2>
		Element &dotSpecialized (Element &res, const Vector1 &v1, const Vector2 &v2,
					 VectorCategories::DenseVectorTag<Trait1> tag1,
					 VectorCategories::DenseVectorTag<Trait2> tag2) const;
		template <class Vector1, class Trait1, class Vector2, class Trait2>
		Element &dotSpecialized (Element &res, const Vector1 &v1, const Vector2 &v2,
					 VectorCategories::SparseSequenceVectorTag<Trait1> tag1,
					 VectorCategories::DenseVectorTag<Trait2> tag2) const;
		template <class Vector1, class Trait1, class Vector2, class Trait2>
		Element &dotSpecialized (Element &res, const Vector1 &v1, const Vector2 &v2,
					 VectorCategories::SparseAssociativeVectorTag<Trait1> tag1,
					 VectorCategories::DenseVectorTag<Trait2> tag2) const;

		template <class Vector1, class Trait1, class Vector2, class Trait2>
		Element &dotSpecialized (Element &res, const Vector1 &v1, const Vector2 &v2,
					 VectorCategories::DenseVectorTag<Trait1> tag1,
					 VectorCategories::SparseSequenceVectorTag<Trait2> tag2) const
			{ return dot (res, v2, v1); }
		template <class Vector1, class Trait1, class Vector2, class Trait2>
		Element &dotSpecialized (Element &res, const Vector1 &v1, const Vector2 &v2,
					 VectorCategories::SparseSequenceVectorTag<Trait1> tag1,
					 VectorCategories::SparseSequenceVectorTag<Trait2> tag2) const;
		template <class Vector1, class Trait1, class Vector2, class Trait2>
		Element &dotSpecialized (Element &res, const Vector1 &v1, const Vector2 &v2,
					 VectorCategories::SparseAssociativeVectorTag<Trait1> tag1,
					 VectorCategories::SparseSequenceVectorTag<Trait2> tag2) const;

		template <class Vector1, class Trait1, class Vector2, class Trait2>
		Element &dotSpecialized (Element &res, const Vector1 &v1, const Vector2 &v2,
					 VectorCategories::DenseVectorTag<Trait1> tag1,
					 VectorCategories::SparseAssociativeVectorTag<Trait2> tag2) const
			{ return dot (res, v2, v1); }
		template <class Vector1, class Trait1, class Vector2, class Trait2>
		inline Element &dotSpecialized (Element &res, const Vector1 &v1, const Vector2 &v2,
						VectorCategories::SparseSequenceVectorTag<Trait1> tag1,
						VectorCategories::SparseAssociativeVectorTag<Trait2> tag2) const
			{ return dot (res, v2, v1); }
		template <class Vector1, class Trait1, class Vector2, class Trait2>
		Element &dotSpecialized (Element &res, const Vector1 &v1, const Vector2 &v2,
					 VectorCategories::SparseAssociativeVectorTag<Trait1> tag1,
					 VectorCategories::SparseAssociativeVectorTag<Trait2> tag2) const;

		template <class Vector, class Trait>
		Vector &mulSpecialized (Vector &res, const Vector &x, const Element &a,
					VectorCategories::DenseVectorTag<Trait> tag) const;
		template <class Vector, class Trait>
		Vector &mulSpecialized (Vector &res, const Vector &x, const Element &a,
					VectorCategories::SparseSequenceVectorTag<Trait> tag) const;
		template <class Vector, class Trait>
		Vector &mulSpecialized (Vector &res, const Vector &x, const Element &a,
					VectorCategories::SparseAssociativeVectorTag<Trait> tag) const;

		template <class Vector, class Trait>
		Vector &mulinSpecialized (Vector &x, const Element &a,
					  VectorCategories::DenseVectorTag<Trait> tag) const;
		template <class Vector, class Trait>
		Vector &mulinSpecialized (Vector &x, const Element &a,
					  VectorCategories::SparseSequenceVectorTag<Trait> tag) const;
		template <class Vector, class Trait>
		Vector &mulinSpecialized (Vector &x, const Element &a,
					  VectorCategories::SparseAssociativeVectorTag<Trait> tag) const;

		template <class Vector, class Trait>
		Vector &axpySpecialized (Vector &res, const Vector &y, const Element &a, const Vector &x,
					 VectorCategories::DenseVectorTag<Trait> tag) const;
		template <class Vector, class Trait>
		Vector &axpySpecialized (Vector &res, const Vector &y, const Element &a, const Vector &x,
					 VectorCategories::SparseSequenceVectorTag<Trait> tag) const;
		template <class Vector, class Trait>
		Vector &axpySpecialized (Vector &res, const Vector &y, const Element &a, const Vector &x,
					 VectorCategories::SparseAssociativeVectorTag<Trait> tag) const;

		template <class Vector, class Trait>
		Vector &axpyinSpecialized (Vector &y, const Element &a, const Vector &x,
					   VectorCategories::DenseVectorTag<Trait> tag) const;
		template <class Vector, class Trait>
		Vector &axpyinSpecialized (Vector &y, const Element &a, const Vector &x,
					   VectorCategories::SparseSequenceVectorTag<Trait> tag) const;
		template <class Vector, class Trait>
		Vector &axpyinSpecialized (Vector &y, const Element &a, const Vector &x,
					   VectorCategories::SparseAssociativeVectorTag<Trait> tag) const;
	}; // class VectorDomain

} // namespace LinBox

#include "linbox/field/vector-domain.C"

#endif // __FIELD_MATRIX_DOMAIN_H
