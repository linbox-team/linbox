/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/vector/vector-domain.h
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
 * 2002-06-21 Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * Added methods add, addin, sub, subin, areEqual, isZero, and copy.
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
#include "linbox/util/field-axpy.h"
#include "linbox/util/debug.h"

namespace LinBox
{
	/** @name Vector domain base
	 *
	 * @doc
	 * This class provides a virtual base for the VectorDomain and the
	 * DotProductDomain. Its purpose is to provide the field of computation
	 * in a single location.
	 */
	template <class Field>
	class VectorDomainBase 
	{
	public:
		VectorDomainBase (const Field &F)
			: _F (F), accu(F)
		{}

	protected:
		Field _F;
		mutable FieldAXPY<Field> accu;
	};

	/** @name Dot product domain
	 * @memo Performance-critical dot products
	 *
	 * @doc
	 * This class contains all of the "high-performance" dot product types
	 * for a vector domain, i.e. dense/dense and dense/sparse parallel. It
	 * allows these dot products to be specialized with very highly tuned
	 * and tightly optimized field-specific implementations, possibly in
	 * assembly language. The other dot products are not considered as
	 * critical for performance, so their implementations are in the derived
	 * class VectorDomain.
	 */
	template <class Field>
	class DotProductDomain : public virtual VectorDomainBase<Field>
	{
	public:

		typedef typename Field::Element Element;

		DotProductDomain (const Field &F)
			: VectorDomainBase<Field> (F)
		{}

	protected:
		template <class Vector1, class Vector2>
		inline Element &dotSpecializedDD (Element &res, const Vector1 &v1, const Vector2 &v2) const;

		template <class Vector1, class Vector2>
		inline Element &dotSpecializedDSP (Element &res, const Vector1 &v1, const Vector2 &v2) const;

	};


	/** @name Vector Domain
	 * @memo Vector arithmetic
	 *
	 * @doc
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
// JGD 01.10.2003 : Why inherit twice from VectorDomainBase<Field> ???
	template <class Field>
	class VectorDomain : public virtual DotProductDomain<Field>, public virtual VectorDomainBase<Field>
	{
      	public:
    
		typedef typename Field::Element         Element;

		/** Copy constructor.
		 * Constructs VectorDomain object by copying the domain.
		 * This is required to allow matrix domain objects to be passed
		 * by value into functions.
		 * @param  MD VectorDomain object.
		 */
		VectorDomain (const VectorDomain &VD)
			: DotProductDomain<Field> (VD._F), VectorDomainBase<Field> (VD._F)
		{}
    
		/** Assignment operator.
		 * Assigns VectorDomain object MD to field.
		 * @param  MD VectorDomain object.
		 */
		VectorDomain &operator = (const VectorDomain &VD)
		{ _F = VD._F; accu = VD.accu;}

		/** Retrieve the underlying field
		 * Return a reference to the field that this matrix domain
		 * object uses
		 * @return reference to field
		 */

		const Field &field () const
		{ return _F; }
    
		/** Vector input/output operations
		 * These routines are useful for reading and writing vectors to
		 * and from file streams. They are analagous to field read and
		 * write operations.
		 */

		//@{

		/** Print vector of field elements.
		 * This function assumes the field element has already been 
		 * constructed and initialized.
		 * @return output stream to which field element is written.
		 * @param  os  output stream to which field element is written.
		 * @param  x   field element.
		 */
		template <class Vector>
		inline std::ostream &write (std::ostream &os, const Vector &x) const
		{ return writeSpecialized (os, x, VectorTraits<Vector>::VectorCategory ()); }

		/** Read vector of field elements.
		 * This function assumes the field element has already been 
		 * constructed and initialized.
		 * @return input stream from which field element is read.
		 * @param  is  input stream from which field element is read.
		 * @param  x   field element.
		 */
		template <class Vector>
		inline std::istream &read (std::istream &is, Vector &x) const
		{ return readSpecialized (is, x, VectorTraits<Vector>::VectorCategory ()); }

		//@} Input/Output Operations

		/** @name Vector arithmetic operations
		 * These routes are analogs of field arithmetic operations, but
		 * they take vectors of elements as input. Vector-vector dot
		 * product and vector-vector axpy are supported here.
		 */

		//@{

		/** Vector copy
		 * Copy a vector to another vector, possibly converting to a
		 * different representation
		 * @param res Output vector
		 * @param v Input vector
		 * @returns reference to output vector
		 */
		template <class Vector1, class Vector2>
		inline Vector1 &copy (Vector1 &res, const Vector2 &v) const
		{ return copySpecialized (res, v,
					  VectorTraits<Vector1>::VectorCategory (),
					  VectorTraits<Vector2>::VectorCategory ()); }

		/** Vector copy
		 * Copy a vector to a portion of another vector, possibly
		 * converting to a different representation
		 * @param res Output vector
		 * @param i Index to which to copy
		 * @param len Length (in indices) of output vector to
		 *            invalidate. If len == 0, then copy the whole vector v
		 * @param v Input vector
		 * @returns reference to output vector
		 */
		template <class Vector1, class Vector2>
		inline Vector1 &copy (Vector1 &res, const Vector2 &v, size_t i, size_t len = 0) const
		{ return copySpecialized (res, v, i, len,
					  VectorTraits<Vector1>::VectorCategory ()); }

		/** Vector equality
		 * @param v1 Input vector
		 * @param v2 Input vector
		 * @returns true if and only if the vectors v1 and v2 are equal
		 */
		template <class Vector1, class Vector2>
		inline bool areEqual (const Vector1 &v1, const Vector2 &v2) const
		{ return areEqualSpecialized (v1, v2,
					      VectorTraits<Vector1>::VectorCategory (),
					      VectorTraits<Vector2>::VectorCategory ()); }

		/** Vector equality with zero
		 * @param v Input vector
		 * @returns true if and only if the vector v is zero
		 */
		template <class Vector>
		inline bool isZero (const Vector &v) const
		{ return isZeroSpecialized (v, VectorTraits<Vector>::VectorCategory ()); }

		/** Vector-vector dot product
		 * @param res element into which to store result
		 * @param v1 Input vector
		 * @param v2 Input vector
		 */
		template <class Vector1, class Vector2>
		inline Element &dot (Element &res, const Vector1 &v1, const Vector2 &v2) const
		{ return dotSpecialized (res, v1, v2,
					 VectorTraits<Vector1>::VectorCategory (),
					 VectorTraits<Vector2>::VectorCategory ()); }

		/* Alias for the above, to avoid source incompatibility */
		template <class Vector1, class Vector2>
		inline Element &dotprod (Element &res, const Vector1 &v1, const Vector2 &v2) const
		{ return dot (res, v1, v2); }

		/** Vector add
		 * res <- y + x
		 * @param res Vector into which to store result
		 * @param y Input vector y
		 * @param x Input vector x
		 */
		template <class Vector1, class Vector2, class Vector3>
		inline Vector1 &add (Vector1 &res, const Vector2 &y, const Vector3 &x) const
		{ return addSpecialized (res, y, x,
					 VectorTraits<Vector1>::VectorCategory (),
					 VectorTraits<Vector2>::VectorCategory (),
					 VectorTraits<Vector3>::VectorCategory ()); }

		/** Vector in-place add
		 * y <- y + x
		 * @param y Input vector y; result is stored here
		 * @param x Input vector x
		 */
		template <class Vector1, class Vector2>
		inline Vector1 &addin (Vector1 &y, const Vector2 &x) const
		{ return addinSpecialized (y, x,
					   VectorTraits<Vector1>::VectorCategory (),
					   VectorTraits<Vector2>::VectorCategory ()); }

		/** Vector subtract
		 * res <- y - x
		 * @param res Vector into which to store result
		 * @param y Input vector y
		 * @param x Input vector x
		 */
		template <class Vector1, class Vector2, class Vector3>
		inline Vector1 &sub (Vector1 &res, const Vector2 &y, const Vector3 &x) const
		{ return subSpecialized (res, y, x,
					 VectorTraits<Vector1>::VectorCategory (),
					 VectorTraits<Vector2>::VectorCategory (),
					 VectorTraits<Vector3>::VectorCategory ()); }

		/** Vector in-place subtract
		 * y <- y - x
		 * @param y Input vector y; result is stored here
		 * @param x Input vector x
		 */
		template <class Vector1, class Vector2>
		inline Vector1 &subin (Vector1 &y, const Vector2 &x) const
		{ return subinSpecialized (y, x,
					   VectorTraits<Vector1>::VectorCategory (),
					   VectorTraits<Vector2>::VectorCategory ()); }

		/** Vector negate
		 * res <- -x
		 * @param res Vector into which to store result
		 * @param x Input vector x
		 */
		template <class Vector1, class Vector2>
		inline Vector1 &neg (Vector1 &res, const Vector2 &x) const
		{ return negSpecialized (res, x,
					 VectorTraits<Vector1>::VectorCategory (),
					 VectorTraits<Vector2>::VectorCategory ()); }

		/** Vector in-place negate
		 * y <- -y
		 * @param y Input vector y; result is stored here
		 */
		template <class Vector>
		inline Vector &negin (Vector &y) const
		{ return neginSpecialized (y, VectorTraits<Vector>::VectorCategory ()); }

		/** Scalar-vector multiplication
		 * res <- a * x
		 * @param res Vector into which to store result
		 * @param x Input vector x
		 * @param a Input element a
		 */
		template <class Vector1, class Vector2>
		inline Vector1 &mul (Vector1 &res, const Vector2 &x, const Element &a) const
		{ return mulSpecialized (res, x, a, VectorTraits<Vector1>::VectorCategory ()); }

		/** In-place scalar-vector multiplication
		 * a <- a * x
		 * @param res Vector into which to store result
		 * @param x Input vector x
		 * @param a Input element a
		 */
		template <class Vector>
		inline Vector &mulin (Vector &x, const Element &a) const
		{ return mulinSpecialized (x, a, VectorTraits<Vector>::VectorCategory ()); }

		/** Vector axpy
		 * res <- y + a*x
		 * @param res Vector into which to store result
		 * @param a Scalar element a
		 * @param x Input vector x
		 * @param y Input vector y
		 */
		template <class Vector1, class Vector2, class Vector3>
		inline Vector1 &axpy (Vector1 &res, const Element &a, const Vector2 &x, const Vector3 &y) const
		{ return axpySpecialized (res, y, a, x, VectorTraits<Vector1>::VectorCategory ()); }

		/** Vector in-place axpy
		 * y <- y + a*x
		 * @param y Input vector y; result is stored here
		 * @param a Scalar element a
		 * @param x Input vector x
		 */
		template <class Vector1, class Vector2>
		inline Vector1 &axpyin (Vector1 &y, const Element &a, const Vector2 &x) const
		{ return axpyinSpecialized (y, a, x,
					    VectorTraits<Vector1>::VectorCategory (),
					    VectorTraits<Vector2>::VectorCategory ()); }

		//@} Vector arithmetic operations

		/** @name Implementation-Specific Methods.
		 * These methods are not required of all \Ref{LinBox Fields}
		 * and are included only for this implementation of the archetype.
		 */
		//@{

		/** Construct from a field
		 * @param F Field from which to construct
		 */
		VectorDomain (const Field &F)
			: 
                        VectorDomainBase<Field> (F), DotProductDomain<Field> (F)
		{}

		//@} Implementation-Specific Methods
    
	protected:

		// Specialized function implementations
		template <class Vector, class Trait>
		std::ostream &writeSpecialized (std::ostream &os, const Vector &x,
						VectorCategories::DenseVectorTag<Trait> tag) const;
		template <class Vector, class Trait>
		std::ostream &writeSpecialized (std::ostream &os, const Vector &x,
						VectorCategories::SparseSequenceVectorTag<Trait> tag) const;
		template <class Vector, class Trait>
		std::ostream &writeSpecialized (std::ostream &os, const Vector &x,
						VectorCategories::SparseAssociativeVectorTag<Trait> tag) const;
		template <class Vector, class Trait>
		std::ostream &writeSpecialized (std::ostream &os, const Vector &x,
						VectorCategories::SparseParallelVectorTag<Trait> tag) const;

		template <class Vector, class Trait>
		std::istream &readSpecialized (std::istream &is, const Vector &x,
					       VectorCategories::DenseVectorTag<Trait> tag) const;
		template <class Vector, class Trait>
		std::istream &readSpecialized (std::istream &is, const Vector &x,
					       VectorCategories::SparseSequenceVectorTag<Trait> tag) const;
		template <class Vector, class Trait>
		std::istream &readSpecialized (std::istream &is, const Vector &x,
					       VectorCategories::SparseAssociativeVectorTag<Trait> tag) const;
		template <class Vector, class Trait>
		std::istream &readSpecialized (std::istream &is, const Vector &x,
					       VectorCategories::SparseParallelVectorTag<Trait> tag) const;

		template <class Vector1, class Trait1, class Vector2, class Trait2>
		bool areEqualSpecialized (const Vector1 &v1, const Vector2 &v2,
					  VectorCategories::DenseVectorTag<Trait1> tag1,
					  VectorCategories::DenseVectorTag<Trait2> tag2) const;
		template <class Vector1, class Trait1, class Vector2, class Trait2>
		bool areEqualSpecialized (const Vector1 &v1, const Vector2 &v2,
					  VectorCategories::SparseSequenceVectorTag<Trait1> tag1,
					  VectorCategories::DenseVectorTag<Trait2> tag2) const;
		template <class Vector1, class Trait1, class Vector2, class Trait2>
		bool areEqualSpecialized (const Vector1 &v1, const Vector2 &v2,
					  VectorCategories::SparseAssociativeVectorTag<Trait1> tag1,
					  VectorCategories::DenseVectorTag<Trait2> tag2) const;
		template <class Vector1, class Trait1, class Vector2, class Trait2>
		bool areEqualSpecialized (const Vector1 &v1, const Vector2 &v2,
					  VectorCategories::SparseParallelVectorTag<Trait1> tag1,
					  VectorCategories::DenseVectorTag<Trait2> tag2) const;

		template <class Vector1, class Trait1, class Vector2, class Trait2>
		bool areEqualSpecialized (const Vector1 &v1, const Vector2 &v2,
					  VectorCategories::DenseVectorTag<Trait1> tag1,
					  VectorCategories::SparseSequenceVectorTag<Trait2> tag2) const
		{ return areEqual (v2, v1); }
		template <class Vector1, class Trait1, class Vector2, class Trait2>
		bool areEqualSpecialized (const Vector1 &v1, const Vector2 &v2,
					  VectorCategories::SparseSequenceVectorTag<Trait1> tag1,
					  VectorCategories::SparseSequenceVectorTag<Trait2> tag2) const;
		template <class Vector1, class Trait1, class Vector2, class Trait2>
		bool areEqualSpecialized (const Vector1 &v1, const Vector2 &v2,
					  VectorCategories::SparseAssociativeVectorTag<Trait1> tag1,
					  VectorCategories::SparseSequenceVectorTag<Trait2> tag2) const;
		template <class Vector1, class Trait1, class Vector2, class Trait2>
		bool areEqualSpecialized (const Vector1 &v1, const Vector2 &v2,
					  VectorCategories::SparseParallelVectorTag<Trait1> tag1,
					  VectorCategories::SparseSequenceVectorTag<Trait2> tag2) const;

		template <class Vector1, class Trait1, class Vector2, class Trait2>
		bool areEqualSpecialized (const Vector1 &v1, const Vector2 &v2,
					  VectorCategories::DenseVectorTag<Trait1> tag1,
					  VectorCategories::SparseAssociativeVectorTag<Trait2> tag2) const
		{ return areEqual (v2, v1); }
		template <class Vector1, class Trait1, class Vector2, class Trait2>
		inline bool areEqualSpecialized (const Vector1 &v1, const Vector2 &v2,
						 VectorCategories::SparseSequenceVectorTag<Trait1> tag1,
						 VectorCategories::SparseAssociativeVectorTag<Trait2> tag2) const
		{ return areEqual (v2, v1); }
		template <class Vector1, class Trait1, class Vector2, class Trait2>
		bool areEqualSpecialized (const Vector1 &v1, const Vector2 &v2,
					  VectorCategories::SparseAssociativeVectorTag<Trait1> tag1,
					  VectorCategories::SparseAssociativeVectorTag<Trait2> tag2) const;
		template <class Vector1, class Trait1, class Vector2, class Trait2>
		bool areEqualSpecialized (const Vector1 &v1, const Vector2 &v2,
					  VectorCategories::SparseParallelVectorTag<Trait1> tag1,
					  VectorCategories::SparseAssociativeVectorTag<Trait2> tag2) const;

		template <class Vector1, class Trait1, class Vector2, class Trait2>
		bool areEqualSpecialized (const Vector1 &v1, const Vector2 &v2,
					  VectorCategories::DenseVectorTag<Trait1> tag1,
					  VectorCategories::SparseParallelVectorTag<Trait2> tag2) const
		{ return areEqual (v2, v1); }
		template <class Vector1, class Trait1, class Vector2, class Trait2>
		inline bool areEqualSpecialized (const Vector1 &v1, const Vector2 &v2,
						 VectorCategories::SparseSequenceVectorTag<Trait1> tag1,
						 VectorCategories::SparseParallelVectorTag<Trait2> tag2) const
		{ return areEqual (v2, v1); }
		template <class Vector1, class Trait1, class Vector2, class Trait2>
		inline bool areEqualSpecialized (const Vector1 &v1, const Vector2 &v2,
						 VectorCategories::SparseAssociativeVectorTag<Trait1> tag1,
						 VectorCategories::SparseParallelVectorTag<Trait2> tag2) const
		{ return areEqual (v2, v1); }
		template <class Vector1, class Trait1, class Vector2, class Trait2>
		bool areEqualSpecialized (const Vector1 &v1, const Vector2 &v2,
					  VectorCategories::SparseParallelVectorTag<Trait1> tag1,
					  VectorCategories::SparseParallelVectorTag<Trait2> tag2) const;

		template <class Vector, class Trait>
		bool isZeroSpecialized (const Vector &v, VectorCategories::DenseVectorTag<Trait> tag) const;
		template <class Vector, class Trait>
		bool isZeroSpecialized (const Vector &v, VectorCategories::SparseSequenceVectorTag<Trait> tag) const;
		template <class Vector, class Trait>
		bool isZeroSpecialized (const Vector &v, VectorCategories::SparseAssociativeVectorTag<Trait> tag) const;
		template <class Vector, class Trait>
		bool isZeroSpecialized (const Vector &v, VectorCategories::SparseParallelVectorTag<Trait> tag) const;

		template <class Vector1, class Trait1, class Vector2, class Trait2>
		inline Vector1 &copySpecialized (Vector1 &res, const Vector2 &v,
						 VectorCategories::DenseVectorTag<Trait1> tag1,
						 VectorCategories::DenseVectorTag<Trait2> tag2) const
		{ std::copy (v.begin (), v.end (), res.begin ()); return res; }
		template <class Vector1, class Trait1, class Vector2, class Trait2>
		Vector1 &copySpecialized (Vector1 &res, const Vector2 &v,
					  VectorCategories::SparseSequenceVectorTag<Trait1> tag1,
					  VectorCategories::DenseVectorTag<Trait2> tag2) const;
		template <class Vector1, class Trait1, class Vector2, class Trait2>
		Vector1 &copySpecialized (Vector1 &res, const Vector2 &v,
					  VectorCategories::SparseAssociativeVectorTag<Trait1> tag1,
					  VectorCategories::DenseVectorTag<Trait2> tag2) const;
		template <class Vector1, class Trait1, class Vector2, class Trait2>
		Vector1 &copySpecialized (Vector1 &res, const Vector2 &v,
					  VectorCategories::SparseParallelVectorTag<Trait1> tag1,
					  VectorCategories::DenseVectorTag<Trait2> tag2) const;

		template <class Vector1, class Trait1, class Vector2, class Trait2>
		Vector1 &copySpecialized (Vector1 &res, const Vector2 &v,
					  VectorCategories::DenseVectorTag<Trait1> tag1,
					  VectorCategories::SparseSequenceVectorTag<Trait2> tag2) const;
		template <class Vector1, class Trait1, class Vector2, class Trait2>
		inline Vector1 &copySpecialized (Vector1 &res, const Vector2 &v,
						 VectorCategories::SparseSequenceVectorTag<Trait1> tag1,
						 VectorCategories::SparseSequenceVectorTag<Trait2> tag2) const
		{ res.resize (v.size ()); std::copy (v.begin (), v.end (), res.begin ()); return res; }
		template <class Vector1, class Trait1, class Vector2, class Trait2>
		Vector1 &copySpecialized (Vector1 &res, const Vector2 &v,
					  VectorCategories::SparseAssociativeVectorTag<Trait1> tag1,
					  VectorCategories::SparseSequenceVectorTag<Trait2> tag2) const;
		template <class Vector1, class Trait1, class Vector2, class Trait2>
		Vector1 &copySpecialized (Vector1 &res, const Vector2 &v,
					  VectorCategories::SparseParallelVectorTag<Trait1> tag1,
					  VectorCategories::SparseSequenceVectorTag<Trait2> tag2) const;

		template <class Vector1, class Trait1, class Vector2, class Trait2>
		Vector1 &copySpecialized (Vector1 &res, const Vector2 &v,
					  VectorCategories::DenseVectorTag<Trait1> tag1,
					  VectorCategories::SparseAssociativeVectorTag<Trait2> tag2) const;
		template <class Vector1, class Trait1, class Vector2, class Trait2>
		Vector1 &copySpecialized (Vector1 &res, const Vector2 &v,
					  VectorCategories::SparseSequenceVectorTag<Trait1> tag1,
					  VectorCategories::SparseAssociativeVectorTag<Trait2> tag2) const;
		template <class Vector1, class Trait1, class Vector2, class Trait2>
		Vector1 &copySpecialized (Vector1 &res, const Vector2 &v,
					  VectorCategories::SparseAssociativeVectorTag<Trait1> tag1,
					  VectorCategories::SparseAssociativeVectorTag<Trait2> tag2) const;
		template <class Vector1, class Trait1, class Vector2, class Trait2>
		Vector1 &copySpecialized (Vector1 &res, const Vector2 &v,
					  VectorCategories::SparseParallelVectorTag<Trait1> tag1,
					  VectorCategories::SparseAssociativeVectorTag<Trait2> tag2) const;

		template <class Vector1, class Trait1, class Vector2, class Trait2>
		Vector1 &copySpecialized (Vector1 &res, const Vector2 &v,
					  VectorCategories::DenseVectorTag<Trait1> tag1,
					  VectorCategories::SparseParallelVectorTag<Trait2> tag2) const;
		template <class Vector1, class Trait1, class Vector2, class Trait2>
		Vector1 &copySpecialized (Vector1 &res, const Vector2 &v,
					  VectorCategories::SparseSequenceVectorTag<Trait1> tag1,
					  VectorCategories::SparseParallelVectorTag<Trait2> tag2) const;
		template <class Vector1, class Trait1, class Vector2, class Trait2>
		Vector1 &copySpecialized (Vector1 &res, const Vector2 &v,
					  VectorCategories::SparseAssociativeVectorTag<Trait1> tag1,
					  VectorCategories::SparseParallelVectorTag<Trait2> tag2) const;
		template <class Vector1, class Trait1, class Vector2, class Trait2>
		Vector1 &copySpecialized (Vector1 &res, const Vector2 &v,
					  VectorCategories::SparseParallelVectorTag<Trait1> tag1,
					  VectorCategories::SparseParallelVectorTag<Trait2> tag2) const;

		template <class Vector1, class Trait, class Vector2>
		Vector1 &copySpecialized (Vector1 &res, const Vector2 &v, size_t i, size_t len,
					  VectorCategories::DenseVectorTag<Trait> tag) const;
		template <class Vector1, class Trait, class Vector2>
		Vector1 &copySpecialized (Vector1 &res, const Vector2 &v, size_t i, size_t len,
					  VectorCategories::SparseSequenceVectorTag<Trait> tag) const;
		template <class Vector1, class Trait, class Vector2>
		Vector1 &copySpecialized (Vector1 &res, const Vector2 &v, size_t i, size_t len,
					  VectorCategories::SparseAssociativeVectorTag<Trait> tag) const;
		template <class Vector1, class Trait, class Vector2>
		Vector1 &copySpecialized (Vector1 &res, const Vector2 &v, size_t i, size_t len,
					  VectorCategories::SparseParallelVectorTag<Trait> tag) const;

		// These versions are optimized for the case where one is
		// copying between vectors of the same type. It avoids
		// additional memory allocation and copying.
		template <class Vector, class Trait>
		Vector &copySpecialized (Vector &res, const Vector &v, size_t i, size_t len,
					 VectorCategories::DenseVectorTag<Trait> tag) const;
		template <class Vector, class Trait>
		Vector &copySpecialized (Vector &res, const Vector &v, size_t i, size_t len,
					 VectorCategories::SparseSequenceVectorTag<Trait> tag) const;
		template <class Vector, class Trait>
		Vector &copySpecialized (Vector &res, const Vector &v, size_t i, size_t len,
					 VectorCategories::SparseAssociativeVectorTag<Trait> tag) const;
		template <class Vector, class Trait>
		Vector &copySpecialized (Vector &res, const Vector &v, size_t i, size_t len,
					 VectorCategories::SparseParallelVectorTag<Trait> tag) const;

		template <class Vector1, class Trait1, class Vector2, class Trait2>
		inline Element &dotSpecialized (Element &res, const Vector1 &v1, const Vector2 &v2,
						VectorCategories::DenseVectorTag<Trait1> tag1,
						VectorCategories::DenseVectorTag<Trait2> tag2) const
		{ return DotProductDomain<Field>::dotSpecializedDD (res, v1, v2); }
		template <class Vector1, class Trait1, class Vector2, class Trait2>
		inline Element &dotSpecialized (Element &res, const Vector1 &v1, const Vector2 &v2,
						VectorCategories::SparseSequenceVectorTag<Trait1> tag1,
						VectorCategories::DenseVectorTag<Trait2> tag2) const;
		template <class Vector1, class Trait1, class Vector2, class Trait2>
		inline Element &dotSpecialized (Element &res, const Vector1 &v1, const Vector2 &v2,
						VectorCategories::SparseAssociativeVectorTag<Trait1> tag1,
						VectorCategories::DenseVectorTag<Trait2> tag2) const;
		template <class Vector1, class Trait1, class Vector2, class Trait2>
		inline Element &dotSpecialized (Element &res, const Vector1 &v1, const Vector2 &v2,
						VectorCategories::SparseParallelVectorTag<Trait1> tag1,
						VectorCategories::DenseVectorTag<Trait2> tag2) const
		{ return DotProductDomain<Field>::dotSpecializedDSP (res, v1, v2); }

		template <class Vector1, class Trait1, class Vector2, class Trait2>
		inline Element &dotSpecialized (Element &res, const Vector1 &v1, const Vector2 &v2,
						VectorCategories::DenseVectorTag<Trait1> tag1,
						VectorCategories::SparseSequenceVectorTag<Trait2> tag2) const
		{ return dot (res, v2, v1); }
		template <class Vector1, class Trait1, class Vector2, class Trait2>
		inline Element &dotSpecialized (Element &res, const Vector1 &v1, const Vector2 &v2,
						VectorCategories::SparseSequenceVectorTag<Trait1> tag1,
						VectorCategories::SparseSequenceVectorTag<Trait2> tag2) const;
		template <class Vector1, class Trait1, class Vector2, class Trait2>
		inline Element &dotSpecialized (Element &res, const Vector1 &v1, const Vector2 &v2,
						VectorCategories::SparseAssociativeVectorTag<Trait1> tag1,
						VectorCategories::SparseSequenceVectorTag<Trait2> tag2) const;
		template <class Vector1, class Trait1, class Vector2, class Trait2>
		inline Element &dotSpecialized (Element &res, const Vector1 &v1, const Vector2 &v2,
						VectorCategories::SparseParallelVectorTag<Trait1> tag1,
						VectorCategories::SparseSequenceVectorTag<Trait2> tag2) const;

		template <class Vector1, class Trait1, class Vector2, class Trait2>
		inline Element &dotSpecialized (Element &res, const Vector1 &v1, const Vector2 &v2,
						VectorCategories::DenseVectorTag<Trait1> tag1,
						VectorCategories::SparseAssociativeVectorTag<Trait2> tag2) const
		{ return dot (res, v2, v1); }
		template <class Vector1, class Trait1, class Vector2, class Trait2>
		inline Element &dotSpecialized (Element &res, const Vector1 &v1, const Vector2 &v2,
						VectorCategories::SparseSequenceVectorTag<Trait1> tag1,
						VectorCategories::SparseAssociativeVectorTag<Trait2> tag2) const
		{ return dot (res, v2, v1); }
		template <class Vector1, class Trait1, class Vector2, class Trait2>
		inline Element &dotSpecialized (Element &res, const Vector1 &v1, const Vector2 &v2,
						VectorCategories::SparseAssociativeVectorTag<Trait1> tag1,
						VectorCategories::SparseAssociativeVectorTag<Trait2> tag2) const;
		template <class Vector1, class Trait1, class Vector2, class Trait2>
		inline Element &dotSpecialized (Element &res, const Vector1 &v1, const Vector2 &v2,
						VectorCategories::SparseParallelVectorTag<Trait1> tag1,
						VectorCategories::SparseAssociativeVectorTag<Trait2> tag2) const;

		template <class Vector1, class Trait1, class Vector2, class Trait2>
		inline Element &dotSpecialized (Element &res, const Vector1 &v1, const Vector2 &v2,
						VectorCategories::DenseVectorTag<Trait1> tag1,
						VectorCategories::SparseParallelVectorTag<Trait2> tag2) const
		{ return dot (res, v2, v1); }
		template <class Vector1, class Trait1, class Vector2, class Trait2>
		inline Element &dotSpecialized (Element &res, const Vector1 &v1, const Vector2 &v2,
						VectorCategories::SparseSequenceVectorTag<Trait1> tag1,
						VectorCategories::SparseParallelVectorTag<Trait2> tag2) const
		{ return dot (res, v2, v1); }
		template <class Vector1, class Trait1, class Vector2, class Trait2>
		inline Element &dotSpecialized (Element &res, const Vector1 &v1, const Vector2 &v2,
						VectorCategories::SparseAssociativeVectorTag<Trait1> tag1,
						VectorCategories::SparseParallelVectorTag<Trait2> tag2) const
		{ return dot (res, v2, v1); }
		template <class Vector1, class Trait1, class Vector2, class Trait2>
		inline Element &dotSpecialized (Element &res, const Vector1 &v1, const Vector2 &v2,
						VectorCategories::SparseParallelVectorTag<Trait1> tag1,
						VectorCategories::SparseParallelVectorTag<Trait2> tag2) const;

		template <class Vector1, class Trait1, class Vector2, class Trait2, class Vector3, class Trait3>
		Vector1 &addSpecialized (Vector1 &res, const Vector2 &y, const Vector3 &x,
					 VectorCategories::DenseVectorTag<Trait1>,
					 VectorCategories::DenseVectorTag<Trait2>,
					 VectorCategories::DenseVectorTag<Trait3>) const;
		template <class Vector1, class Trait1, class Vector2, class Trait2, class Vector3, class Trait3>
		Vector1 &addSpecialized (Vector1 &res, const Vector2 &y, const Vector3 &x,
					 VectorCategories::SparseSequenceVectorTag<Trait1>,
					 VectorCategories::SparseSequenceVectorTag<Trait2>,
					 VectorCategories::SparseSequenceVectorTag<Trait3>) const;
		template <class Vector1, class Trait1, class Vector2, class Trait2, class Vector3, class Trait3>
		Vector1 &addSpecialized (Vector1 &res, const Vector2 &y, const Vector3 &x,
					 VectorCategories::SparseAssociativeVectorTag<Trait1>,
					 VectorCategories::SparseAssociativeVectorTag<Trait2>,
					 VectorCategories::SparseAssociativeVectorTag<Trait3>) const;
		template <class Vector1, class Trait1, class Vector2, class Trait2, class Vector3, class Trait3>
		Vector1 &addSpecialized (Vector1 &res, const Vector2 &y, const Vector3 &x,
					 VectorCategories::SparseParallelVectorTag<Trait1>,
					 VectorCategories::SparseParallelVectorTag<Trait2>,
					 VectorCategories::SparseParallelVectorTag<Trait3>) const;

		template <class Vector1, class Trait1, class Vector2, class Trait2>
		Vector1 &addinSpecialized (Vector1 &y, const Vector2 &x,
					   VectorCategories::DenseVectorTag<Trait1>,
					   VectorCategories::DenseVectorTag<Trait2>) const;
		template <class Vector1, class Trait1, class Vector2, class Trait2>
		Vector1 &addinSpecialized (Vector1 &y, const Vector2 &x,
					   VectorCategories::SparseSequenceVectorTag<Trait1>,
					   VectorCategories::SparseSequenceVectorTag<Trait2>) const;
		template <class Vector1, class Trait1, class Vector2, class Trait2>
		Vector1 &addinSpecialized (Vector1 &y, const Vector2 &x,
					   VectorCategories::SparseAssociativeVectorTag<Trait1>,
					   VectorCategories::SparseAssociativeVectorTag<Trait2>) const;
		template <class Vector1, class Trait1, class Vector2, class Trait2>
		Vector1 &addinSpecialized (Vector1 &y, const Vector2 &x,
					   VectorCategories::SparseParallelVectorTag<Trait1>,
					   VectorCategories::SparseParallelVectorTag<Trait2>) const;

		template <class Vector1, class Trait1, class Vector2, class Trait2, class Vector3, class Trait3>
		Vector1 &subSpecialized (Vector1 &res, const Vector2 &y, const Vector3 &x,
					 VectorCategories::DenseVectorTag<Trait1>,
					 VectorCategories::DenseVectorTag<Trait2>,
					 VectorCategories::DenseVectorTag<Trait3>) const;
		template <class Vector1, class Trait1, class Vector2, class Trait2, class Vector3, class Trait3>
		Vector1 &subSpecialized (Vector1 &res, const Vector2 &y, const Vector3 &x,
					 VectorCategories::SparseSequenceVectorTag<Trait1>,
					 VectorCategories::SparseSequenceVectorTag<Trait2>,
					 VectorCategories::SparseSequenceVectorTag<Trait3>) const;
		template <class Vector1, class Trait1, class Vector2, class Trait2, class Vector3, class Trait3>
		Vector1 &subSpecialized (Vector1 &res, const Vector2 &y, const Vector3 &x,
					 VectorCategories::SparseAssociativeVectorTag<Trait1>,
					 VectorCategories::SparseAssociativeVectorTag<Trait2>,
					 VectorCategories::SparseAssociativeVectorTag<Trait3>) const;
		template <class Vector1, class Trait1, class Vector2, class Trait2, class Vector3, class Trait3>
		Vector1 &subSpecialized (Vector1 &res, const Vector2 &y, const Vector3 &x,
					 VectorCategories::SparseParallelVectorTag<Trait1>,
					 VectorCategories::SparseParallelVectorTag<Trait2>,
					 VectorCategories::SparseParallelVectorTag<Trait3>) const;

		template <class Vector1, class Trait1, class Vector2, class Trait2>
		Vector1 &subinSpecialized (Vector1 &y, const Vector2 &x,
					   VectorCategories::DenseVectorTag<Trait1>,
					   VectorCategories::DenseVectorTag<Trait2>) const;
		template <class Vector1, class Trait1, class Vector2, class Trait2>
		Vector1 &subinSpecialized (Vector1 &y, const Vector2 &x,
					   VectorCategories::SparseSequenceVectorTag<Trait1>,
					   VectorCategories::SparseSequenceVectorTag<Trait2>) const;
		template <class Vector1, class Trait1, class Vector2, class Trait2>
		Vector1 &subinSpecialized (Vector1 &y, const Vector2 &x,
					   VectorCategories::SparseAssociativeVectorTag<Trait1>,
					   VectorCategories::SparseAssociativeVectorTag<Trait2>) const;
		template <class Vector1, class Trait1, class Vector2, class Trait2>
		Vector1 &subinSpecialized (Vector1 &y, const Vector2 &x,
					   VectorCategories::SparseParallelVectorTag<Trait1>,
					   VectorCategories::SparseParallelVectorTag<Trait2>) const;

		template <class Vector1, class Trait1, class Vector2, class Trait2>
		Vector1 &negSpecialized (Vector1 &res, const Vector2 &x,
					 VectorCategories::DenseVectorTag<Trait1>,
					 VectorCategories::DenseVectorTag<Trait2>) const;
		template <class Vector1, class Trait1, class Vector2, class Trait2>
		Vector1 &negSpecialized (Vector1 &res, const Vector2 &x,
					 VectorCategories::SparseSequenceVectorTag<Trait1>,
					 VectorCategories::SparseSequenceVectorTag<Trait2>) const;
		template <class Vector1, class Trait1, class Vector2, class Trait2>
		Vector1 &negSpecialized (Vector1 &res, const Vector2 &x,
					 VectorCategories::SparseAssociativeVectorTag<Trait1>,
					 VectorCategories::SparseAssociativeVectorTag<Trait2>) const;
		template <class Vector1, class Trait1, class Vector2, class Trait2>
		Vector1 &negSpecialized (Vector1 &res, const Vector2 &x,
					 VectorCategories::SparseParallelVectorTag<Trait1>,
					 VectorCategories::SparseParallelVectorTag<Trait2>) const;

		template <class Vector, class Trait>
		Vector &neginSpecialized (Vector &y,
					  VectorCategories::DenseVectorTag<Trait> tag) const;
		template <class Vector, class Trait>
		Vector &neginSpecialized (Vector &y,
					  VectorCategories::SparseSequenceVectorTag<Trait> tag) const;
		template <class Vector, class Trait>
		Vector &neginSpecialized (Vector &y,
					  VectorCategories::SparseAssociativeVectorTag<Trait> tag) const;
		template <class Vector, class Trait>
		Vector &neginSpecialized (Vector &y,
					  VectorCategories::SparseParallelVectorTag<Trait> tag) const;

		template <class Vector1, class Vector2, class Trait>
		Vector1 &mulSpecialized (Vector1 &res, const Vector2 &x, const Element &a,
					 VectorCategories::DenseVectorTag<Trait> tag) const;
		template <class Vector1, class Vector2, class Trait>
		Vector1 &mulSpecialized (Vector1 &res, const Vector2 &x, const Element &a,
					 VectorCategories::SparseSequenceVectorTag<Trait> tag) const;
		template <class Vector1, class Vector2, class Trait>
		Vector1 &mulSpecialized (Vector1 &res, const Vector2 &x, const Element &a,
					 VectorCategories::SparseAssociativeVectorTag<Trait> tag) const;
		template <class Vector1, class Vector2, class Trait>
		Vector1 &mulSpecialized (Vector1 &res, const Vector2 &x, const Element &a,
					 VectorCategories::SparseParallelVectorTag<Trait> tag) const;

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
		Vector &mulinSpecialized (Vector &x, const Element &a,
					  VectorCategories::SparseParallelVectorTag<Trait> tag) const;

		template <class Vector1, class Vector2, class Vector3, class Trait>
		Vector1 &axpySpecialized (Vector1 &res, const Vector2 &y,
					  const Element &a, const Vector3 &x,
					  VectorCategories::DenseVectorTag<Trait> tag) const;
		template <class Vector1, class Vector2, class Vector3, class Trait>
		Vector1 &axpySpecialized (Vector1 &res, const Vector2 &y,
					  const Element &a, const Vector3 &x,
					  VectorCategories::SparseSequenceVectorTag<Trait> tag) const;
		template <class Vector1, class Vector2, class Vector3, class Trait>
		Vector1 &axpySpecialized (Vector1 &res, const Vector2 &y,
					  const Element &a, const Vector3 &x,
					  VectorCategories::SparseAssociativeVectorTag<Trait> tag) const;
		template <class Vector1, class Vector2, class Vector3, class Trait>
		Vector1 &axpySpecialized (Vector1 &res, const Vector2 &y,
					  const Element &a, const Vector3 &x,
					  VectorCategories::SparseParallelVectorTag<Trait> tag) const;

		template <class Vector1, class Trait1, class Vector2, class Trait2>
		Vector1 &axpyinSpecialized (Vector1 &y, const Element &a, const Vector2 &x,
					    VectorCategories::DenseVectorTag<Trait1> tag1,
					    VectorCategories::DenseVectorTag<Trait2> tag2) const;
		template <class Vector1, class Trait1, class Vector2, class Trait2>
		Vector1 &axpyinSpecialized (Vector1 &y, const Element &a, const Vector2 &x,
					    VectorCategories::DenseVectorTag<Trait1> tag1,
					    VectorCategories::SparseSequenceVectorTag<Trait2> tag2) const;
		template <class Vector1, class Trait1, class Vector2, class Trait2>
		Vector1 &axpyinSpecialized (Vector1 &y, const Element &a, const Vector2 &x,
					    VectorCategories::DenseVectorTag<Trait1> tag1,
					    VectorCategories::SparseAssociativeVectorTag<Trait2> tag2) const;
		template <class Vector1, class Trait1, class Vector2, class Trait2>
		Vector1 &axpyinSpecialized (Vector1 &y, const Element &a, const Vector2 &x,
					    VectorCategories::DenseVectorTag<Trait1> tag1,
					    VectorCategories::SparseParallelVectorTag<Trait2> tag2) const;
		template <class Vector1, class Trait1, class Vector2, class Trait2>
		Vector1 &axpyinSpecialized (Vector1 &y, const Element &a, const Vector2 &x,
					    VectorCategories::SparseSequenceVectorTag<Trait1> tag1,
					    VectorCategories::SparseSequenceVectorTag<Trait2> tag2) const;
		template <class Vector1, class Trait1, class Vector2, class Trait2>
		Vector1 &axpyinSpecialized (Vector1 &y, const Element &a, const Vector2 &x,
					    VectorCategories::SparseAssociativeVectorTag<Trait1> tag1,
					    VectorCategories::SparseAssociativeVectorTag<Trait2> tag2) const;
		template <class Vector1, class Trait1, class Vector2, class Trait2>
		Vector1 &axpyinSpecialized (Vector1 &y, const Element &a, const Vector2 &x,
					    VectorCategories::SparseParallelVectorTag<Trait1> tag1,
					    VectorCategories::SparseParallelVectorTag<Trait2> tag2) const;

		// Specializations for the case where the two vectors are of
		// different representations. This is provided for the benefit
		// of MatrixDomain

		template <class Vector1, class Trait1, class Vector2, class Trait2, class Vector3, class Trait3>
		inline Vector1 &addSpecialized (Vector1 &res, const Vector2 &y, const Vector3 &x,
						VectorCategories::DenseVectorTag<Trait1>,
						VectorCategories::DenseVectorTag<Trait2>,
						VectorCategories::GenericVectorTag<Trait3>) const
		{
			typename LinBox::Vector<Field>::Dense v (res.size ());

			copy (v, x);
			add (res, y, v);

			return res;
		}

		template <class Vector1, class Trait1, class Vector2, class Trait2, class Vector3, class Trait3>
		inline Vector1 &addSpecialized (Vector1 &res, const Vector2 &y, const Vector3 &x,
						VectorCategories::DenseVectorTag<Trait1>,
						VectorCategories::GenericVectorTag<Trait2>,
						VectorCategories::DenseVectorTag<Trait3>) const
		{
			typename LinBox::Vector<Field>::Dense v (res.size ());

			copy (v, y);
			add (res, v, x);

			return res;
		}

		template <class Vector1, class Trait1, class Vector2, class Trait2, class Vector3, class Trait3>
		inline Vector1 &addSpecialized (Vector1 &res, const Vector2 &y, const Vector3 &x,
						VectorCategories::DenseVectorTag<Trait1>,
						VectorCategories::GenericVectorTag<Trait2>,
						VectorCategories::GenericVectorTag<Trait3>) const
		{
			typename LinBox::Vector<Field>::Dense v (res.size ());
			typename LinBox::Vector<Field>::Dense w (res.size ());

			copy (v, x);
			copy (w, y);
			add (res, w, v);

			return res;
		}

		template <class Vector1, class Trait1, class Vector2, class Trait2, class Vector3, class Trait3>
		inline Vector1 &addSpecialized (Vector1 &res, const Vector2 &y, const Vector3 &x,
						VectorCategories::GenericVectorTag<Trait1>,
						VectorCategories::GenericVectorTag<Trait2>,
						VectorCategories::GenericVectorTag<Trait3>) const
		{
			typename LinBox::Vector<Field>::Sparse v;
			typename LinBox::Vector<Field>::Sparse w;
			typename LinBox::Vector<Field>::Sparse u;

			copy (v, x);
			copy (w, y);
			add (u, w, v);
			copy (res, u);

			return u;
		}

		template <class Vector1, class Trait1, class Vector2, class Trait2>
		inline Vector1 &addinSpecialized (Vector1 &y, const Vector2 &x,
						  VectorCategories::DenseVectorTag<Trait1>,
						  VectorCategories::GenericVectorTag<Trait2>) const
		{
			typename LinBox::Vector<Field>::Dense v (y.size ());

			copy (v, x);
			addin (y, v);

			return y;
		}

		template <class Vector1, class Trait1, class Vector2, class Trait2>
		inline Vector1 &addinSpecialized (Vector1 &y, const Vector2 &x,
						  VectorCategories::GenericVectorTag<Trait1>,
						  VectorCategories::GenericVectorTag<Trait2>) const
		{
			typename LinBox::Vector<Field>::Sparse v;
			typename LinBox::Vector<Field>::Sparse w;

			copy (v, x);
			addin (w, v);
			copy (y, w);

			return y;
		}

		template <class Vector1, class Trait1, class Vector2, class Trait2, class Vector3, class Trait3>
		inline Vector1 &subSpecialized (Vector1 &res, const Vector2 &y, const Vector3 &x,
						VectorCategories::DenseVectorTag<Trait1>,
						VectorCategories::DenseVectorTag<Trait2>,
						VectorCategories::GenericVectorTag<Trait3>) const
		{
			typename LinBox::Vector<Field>::Dense v (res.size ());

			copy (v, x);
			sub (res, y, v);

			return res;
		}

		template <class Vector1, class Trait1, class Vector2, class Trait2, class Vector3, class Trait3>
		inline Vector1 &subSpecialized (Vector1 &res, const Vector2 &y, const Vector3 &x,
						VectorCategories::DenseVectorTag<Trait1>,
						VectorCategories::GenericVectorTag<Trait2>,
						VectorCategories::DenseVectorTag<Trait3>) const
		{
			typename LinBox::Vector<Field>::Dense v (res.size ());

			copy (v, y);
			sub (res, v, x);

			return res;
		}

		template <class Vector1, class Trait1, class Vector2, class Trait2, class Vector3, class Trait3>
		inline Vector1 &subSpecialized (Vector1 &res, const Vector2 &y, const Vector3 &x,
						VectorCategories::DenseVectorTag<Trait1>,
						VectorCategories::GenericVectorTag<Trait2>,
						VectorCategories::GenericVectorTag<Trait3>) const
		{
			typename LinBox::Vector<Field>::Dense v (res.size ());
			typename LinBox::Vector<Field>::Dense w (res.size ());

			copy (v, x);
			copy (w, y);
			sub (res, w, v);

			return res;
		}

		template <class Vector1, class Trait1, class Vector2, class Trait2, class Vector3, class Trait3>
		inline Vector1 &subSpecialized (Vector1 &res, const Vector2 &y, const Vector3 &x,
						VectorCategories::GenericVectorTag<Trait1>,
						VectorCategories::GenericVectorTag<Trait2>,
						VectorCategories::GenericVectorTag<Trait3>) const
		{
			typename LinBox::Vector<Field>::Sparse v;
			typename LinBox::Vector<Field>::Sparse w;
			typename LinBox::Vector<Field>::Sparse u;

			copy (v, x);
			copy (w, y);
			sub (u, w, v);
			copy (res, u);

			return u;
		}

		template <class Vector1, class Trait1, class Vector2, class Trait2>
		inline Vector1 &subinSpecialized (Vector1 &y, const Vector2 &x,
						  VectorCategories::DenseVectorTag<Trait1>,
						  VectorCategories::GenericVectorTag<Trait2>) const
		{
			typename LinBox::Vector<Field>::Dense v (y.size ());

			copy (v, x);
			subin (y, v);

			return y;
		}

		template <class Vector1, class Trait1, class Vector2, class Trait2>
		inline Vector1 &subinSpecialized (Vector1 &y, const Vector2 &x,
						  VectorCategories::GenericVectorTag<Trait1>,
						  VectorCategories::GenericVectorTag<Trait2>) const
		{
			typename LinBox::Vector<Field>::Sparse v;
			typename LinBox::Vector<Field>::Sparse w;

			copy (v, x);
			subin (w, v);
			copy (y, w);

			return y;
		}

		template <class Vector1, class Trait1, class Vector2, class Trait2>
		inline Vector1 &negSpecialized (Vector1 &y, const Vector2 &x,
						VectorCategories::DenseVectorTag<Trait1> tag1,
						VectorCategories::GenericVectorTag<Trait2> tag2) const
		{
			typename LinBox::Vector<Field>::Dense v (y.size ());

			copy (v, x);
			neg (y, v);

			return y;
		}

		template <class Vector1, class Trait1, class Vector2, class Trait2>
		inline Vector1 &negSpecialized (Vector1 &y, const Vector2 &x,
						VectorCategories::GenericVectorTag<Trait1>,
						VectorCategories::GenericVectorTag<Trait2>) const
		{
			typename LinBox::Vector<Field>::Sparse v;
			typename LinBox::Vector<Field>::Sparse w;

			copy (v, x);
			neg (w, v);
			copy (y, w);

			return y;
		}

		template <class Vector1, class Trait1, class Vector2, class Trait2>
		inline Vector1 &axpyinSpecialized (Vector1 &y, const Element &a, const Vector2 &x,
						   VectorCategories::GenericVectorTag<Trait1>,
						   VectorCategories::GenericVectorTag<Trait2>) const
		{
			typename LinBox::Vector<Field>::Sparse v;
			typename LinBox::Vector<Field>::Sparse w;

			copy (v, x);
			axpyin (w, a, v);
			copy (y, w);

			return y;
		}
	
	}; // class VectorDomain

} // namespace LinBox

#include "linbox/vector/vector-domain.inl"

#endif // __FIELD_MATRIX_DOMAIN_H
