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

#define VectorDomainBaseType(tag) \
	VectorDomainBase<Field, Vector, VectorCategories::tag##VectorTag>

#define VectorDomainType(tag1, tag2) \
	VectorDomain<Field, Vector1, Vector2, VectorCategories::tag1##VectorTag, VectorCategories::tag2##VectorTag>

namespace LinBox
{
	/** Vector Domain base
	 * Vector domain functions that are parameterized by only one vector
	 * type.
	 */
	template <class Field, class Vector, class Trait = VectorTraits<Vector>::VectorCategory>
	class VectorDomainBase
	{
		public:
    
		typedef typename Field::Element         Element;

		/** Copy constructor.
		 * Constructs VectorDomain_archetype object by copying the domain.
		 * This is required to allow matrix domain objects to be passed
		 * by value into functions.
		 * @param  MD VectorDomain_archetype object.
		 */
		VectorDomainBase (const VectorDomainBase &VD);
    
		/** Assignment operator.
		 * Assigns VectorDomain object MD to field.
		 * @param  MD VectorDomain object.
		 */
		VectorDomainBase &operator = (const VectorDomainBase &VD);
    
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
		ostream &write (ostream &os, const Vector &x) const;

		/** Read vector of field elements.
		 * This function assumes the field element has already been 
		 * constructed and initialized.
		 * In this implementation, this means for the _elem_ptr for x 
		 * exists and does not point to null.
		 * @return input stream from which field element is read.
		 * @param  is  input stream from which field element is read.
		 * @param  x   field element.
		 */
		istream &read (istream &is, Vector &x) const;
    
		//@} Input/Output Operations

		/** @name Vector arithmetic operations
		 * These routes are analogs of field arithmetic operations, but
		 * they take vectors of elements as input. Vector-vector dot
		 * product and vector-vector axpy are supported here.
		 */

		//@{

		/** Scalar-vector multiplication
		 * res <- a * x
		 * @param res Vector into which to store result
		 * @param x Input vector x
		 * @param a Input element a
		 */
		Vector &mul (Vector &res, const Vector &x, const Element &a) const;

		/** In-place scalar-vector multiplication
		 * a <- a * x
		 * @param res Vector into which to store result
		 * @param x Input vector x
		 * @param a Input element a
		 */
		Vector &mulin (Vector &x, const Element &a) const;

		/** Vector axpy
		 * res <- y + a*x
		 * @param res Vector into which to store result
		 * @param y Input vector y
		 * @param a Scalar element a
		 * @param x Input vector x
		 */
		Vector &axpy (Vector &res, const Vector &y, const Element &a, const Vector &x) const;

		/** Vector in-place axpy
		 * y <- y + a*x
		 * @param y Input vector y; result is stored here
		 * @param a Scalar element a
		 * @param x Input vector x
		 */
		Vector &axpyin (Vector &y, const Element &a, const Vector &x) const;
    
		//@} Common Object Interface
    
		/** @name Implementation-Specific Methods.
		 * These methods are not required of all \Ref{LinBox Fields}
		 * and are included only for this implementation of the archetype.
		 */
		//@{

		/** Construct from a field
		 * @param F Field from which to construct
		 */
		VectorDomainBase (const Field &F);

		//@} Implementation-Specific Methods
    
	    protected:

		Field _F;

	}; // class VectorDomainBase

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
	template <class Field,
		  class Vector1, class Vector2,
		  class Trait1 = VectorTraits<Vector1>::VectorCategory,
		  class Trait2 = VectorTraits<Vector2>::VectorCategory>
	class VectorDomain : public VectorDomainBase<Field, Vector1, Trait1>
	{
	    public:

		/** Copy constructor.
		 * Constructs VectorDomain_archetype object by copying the domain.
		 * This is required to allow matrix domain objects to be passed
		 * by value into functions.
		 * @param  MD VectorDomain_archetype object.
		 */
		VectorDomain (const VectorDomain<Field, Vector1, Vector2, Trait1, Trait2> &VD);
    
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
		Element &dotprod (Element &res, const Vector1 &v1, const Vector2 &v2) const;

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
			: VectorDomainBase<Field, Vector1, Trait1> (F)
		{}

		//@} Implementation-Specific Methods
    
	}; // class VectorDomain

	template <class Field, class Vector>
	class VectorDomainBaseType(Dense) 
	{
	    public:
    
		typedef typename Field::Element         Element;

		VectorDomainBase (const VectorDomainBase &VD) : _F (VD._F) {}
		VectorDomainBase &operator = (const VectorDomainBase &VD) { _F = VD._F; return *this; }
  		Field &field () const { return _F; }
		ostream &write (ostream &os, Vector &x) const;
		ostream &read (ostream &os, Vector &x) const;
		Vector &mul (Vector &res, const Vector &x, const Element &a) const;
		Vector &mulin (Vector &x, const Element &a) const;
		Vector &axpy (Vector &res, const Vector &y, const Element &a, const Vector &x) const;
		Vector &axpyin (Vector &y, const Element &a, const Vector &x) const;

		VectorDomainBase (const Field &F) : _F (F) {}

	    protected:

		Field _F;
	};

	template <class Field, class Vector>
	class VectorDomainBaseType(SparseSequence)
	{
	    public:
    
		typedef typename Field::Element         Element;

		VectorDomainBase (const VectorDomainBase &MD) : _F (MD._F) {}
		VectorDomainBase &operator = (const VectorDomainBase &MD)
			{ _F = MD._F; return *this; }
  		Field &field () const { return _F; }
		ostream &write (ostream &os, Vector &x) const;
		ostream &read (ostream &os, Vector &x) const;
		Vector &mul (Vector &res, const Vector &x, const Element &a) const;
		Vector &mulin (Vector &x, const Element &a) const;
		Vector &axpy (Vector &res, const Vector &y, const Element &a, const Vector &x) const;
		Vector &axpyin (Vector &y, const Element &a, const Vector &x) const;

		VectorDomainBase (const Field &F) : _F (F) {}

	    protected:

		Field _F;
	};

	template <class Field, class Vector1, class Vector2>
	class VectorDomainType(Dense, Dense) : public VectorDomainBase<Field, Vector1, VectorCategories::DenseVectorTag>
	{
	    public:
		VectorDomain (const VectorDomainType(Dense, Dense) &VD)
			: VectorDomainBase<Field, Vector1, VectorCategories::DenseVectorTag> (VD) {}
		Element &dotprod (Element &res, const Vector1 &v1, const Vector2 &v2) const;
		VectorDomain (const Field &F) : VectorDomainBase<Field, Vector1, VectorCategories::DenseVectorTag> (F) {}
	};

	template <class Field, class Vector1, class Vector2>
	class VectorDomainType(SparseSequence, Dense)
		: public VectorDomainBase<Field, Vector1, VectorCategories::SparseSequenceVectorTag>
	{
	    public:
		VectorDomain (const VectorDomainType(SparseSequence, Dense) &VD)
			: VectorDomainBase<Field, Vector1, VectorCategories::SparseSequenceVectorTag> (VD) {}
		Element &dotprod (Element &res, const Vector1 &v1, const Vector2 &v2) const;
		VectorDomain (const Field &F) : VectorDomainBase<Field, Vector1, VectorCategories::SparseSequenceVectorTag> (F) {}
	};

	template <class Field, class Vector1, class Vector2>
	class VectorDomainType(SparseSequence, SparseSequence)
		: public VectorDomainBase<Field, Vector1, VectorCategories::SparseSequenceVectorTag>
	{
	    public:
		VectorDomain (const VectorDomainType(SparseSequence, Dense) &VD)
			: VectorDomainBase<Field, Vector1, VectorCategories::SparseSequenceVectorTag> (VD) {}
		Element &dotprod (Element &res, const Vector1 &v1, const Vector2 &v2) const;
		VectorDomain (const Field &F) : VectorDomainBase<Field, Vector1, VectorCategories::SparseSequenceVectorTag> (F) {}
	};

} // namespace LinBox

#include "linbox/field/vector-domain.C"

#endif // __FIELD_MATRIX_DOMAIN_H
