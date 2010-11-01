/* tests/test-sparse0.C	(Formerly test-sparse-matrix.C)
 * Copyright (C) 2002 William J. Turner
 *
 * Written by William J. Turner <wjturner@acm.org>
 *
 * ------------------------------------
 * Modified by Dmitriy Morozov <linbox@foxcub.org>. May 28, 2002.
 *
 * Added parametrization of VectorCategory tags by VectorTraits. See
 * vector-traits.h for more details.
 *  
 * ------------------------------------
 *
 * See COPYING for lincense information.
 */

#ifndef __LINBOX_vector_random_H
#define __LINBOX_vector_random_H

#include <utility>
#include "linbox/integer.h"
#include "linbox/vector/vector-traits.h"

// Namespace in which all LinBox library code lives
namespace LinBox
{
	/** Random vector generator
	 * This templated function takes a field and a random field element
	 * generator and returns a vector of random field elements.
	 * The vector is dense in the field elements, even if the vector is
	 * a sparse LinBox vector.
	 * The funtion is templatized by the field and the vector types being used.
	 * This function calls another function by the same name with an additional
	 * parameter of the vector category of the vector it is called with.
	 * This mechanism is used because functions cannot have partial template
	 * specializations like classes can.
	 * This new, extended function can be specialized for specific fields
	 * and vectors to allow for better performance.
	 * @return	v	vector of random field elements
	 * @param	F	Field in which arithmetic is done
	 * @param	n	integer number of elements in vector
	 * @param	r	Random field element generator
	 \ingroup vector
	 */
	template <class Field, class Vector>
		inline
		Vector
		randomVector(Field& F, size_t n, typename Field::RandIter& r)
		{ 
			return 
				randomVector<Field, Vector>(
						F, 
						n, 
						r, 
						VectorTraits<Vector>::VectorCategory()
						); 
		}

	/* Random dense vector generator
	 * This templated function takes a field and a random field element
	 * generator and returns a vector of random field elements.
	 * The funtion is templatized by the field and the vector types being used.
	 * This function can be specialized for specific fields
	 * and vectors to allow for better performance.
	 * @return	v	vector of random field elements
	 * @param	F	Field in which arithmetic is done
	 * @param	n	integer number of elements in vector
	 * @param	r	Random field element generator
	 * @param tag	category of vector obtained from vector trait
	 */
	template <class Field, class Vector, class VectorTrait>
		inline
		Vector
		randomVector(
				Field& F, 
				size_t n,
				typename Field::RandIter& r,
				VectorCategories::DenseVectorTag<VectorTrait> tag
				)
		{
#ifdef TRACE
			cout << "Called dense random vector" << endl;
#endif // TRACE

			Vector v(n);
			typename Vector::iterator iter;

			for (iter = v.begin(); iter != v.end(); iter++)
				r.random(*iter);

			return v;
  
		} // randomVector(F, r, VectorCategories::DenseVectorTag& tag)

	/* Random sparse sequence vector generator
	 * This templated function takes a field and a random field element
	 * generator and returns a vector of random field elements.
	 * The funtion is templatized by the field and the vector types being used.
	 * This function can be specialized for specific fields
	 * and vectors to allow for better performance.
	 * @return	v	vector of random field elements
	 * @param	F	Field in which arithmetic is done
	 * @param	n	integer number of elements in vector
	 * @param	r	Random field element generator
	 * @param tag	category of vector obtained from vector trait
	 */
	template <class Field, class Vector, class VectorTrait>
		inline
		Vector
		randomVector(
				Field& F, 
				size_t n,
				typename Field::RandIter& r,
				VectorCategories::SparseSequenceVectorTag<VectorTrait> tag
				)
		{
#ifdef TRACE
			cout << "Called sparse sequence random vector" << endl
				<< "   return vector v = " << endl;
#endif // TRACE

			Vector v(n);
			typename Vector::iterator iter;
			size_t i = 0;
			for (iter = v.begin(); iter != v.end(); iter++, i++)
			{
				iter->first = i;
				r.random(iter->second);
#ifdef TRACE
				cout << "      v[" << iter->first << "] = ";
				F.write(cout, iter->second);
				cout << endl;
#endif // TRACE
			}

			return v;
  
		} // randomVector(F, r, VectorCategories::SparseSequenceVectorTag& tag)

	/* Random sparse associative vector generator
	 * This templated function takes a field and a random field element
	 * generator and returns a vector of random field elements.
	 * The funtion is templatized by the field and the vector types being used.
	 * This function can be specialized for specific fields
	 * and vectors to allow for better performance.
	 * @return	v	vector of random field elements
	 * @param	F	Field in which arithmetic is done
	 * @param	n	integer number of elements in vector
	 * @param	r	Random field element generator
	 * @param tag	category of vector obtained from vector trait
	 */
	template <class Field, class Vector, class VectorTrait>
		inline
		Vector
		randomVector(
				Field& F, 
				size_t n,
				typename Field::RandIter& r,
				VectorCategories::SparseAssociativeVectorTag<VectorTrait> tag
				)
		{
#ifdef TRACE
			cout << "Called sparse associative random vector" << endl
				<< "   return vector v = " << endl;
#endif // TRACE

			Vector v;
			typename Field::Element temp;
			
			for (size_t i = 0; i < size_t(n); i++)
			{
				r.random(temp);
				v.insert(make_pair(i, temp));
#ifdef TRACE
				cout << "      v[" << i << "] = ";
				F.write(cout, v[i]);
				cout << endl;
#endif // TRACE
			}

			return v;
  
		} // randomVector(F, r, VectorCategories::SparseAssociativeVectorTag& tag)

}

#endif // __LINBOX_vector_random_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
