/* -*- mode: c; style: linux -*- */

/* linbox/vector/vector-traits.h
 * Copyright (C) 1999-2001 William J Turner,
 *               2001 Bradford Hovinen
 *
 * Written by William J Turner <wjturner@math.ncsu.edu>,
 *            Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * ------------------------------------
 * Modified by Dmitriy Morozov <linbox@foxcub.org>. May 27, 2002.
 *
 * Implemented the Rootbeer meeting changes. Made VectorTag parametrized
 * and added typedef of those parameters to Traits. So now VectorCategory
 * for each vector has a "reference" back to the VectorTrait of each specific
 * vector (list of pairs, deque of pairs, etc.) through a typedef Trait. This
 * allows for generic manipulation of all vectors and placing the 
 * vector-implementation dependent code into VectorTraits only - as is done now
 * with the function sort.
 *
 * ------------------------------------ 
 */

#ifndef __VECTOR_TRAITS_H
#define __VECTOR_TRAITS_H

#include <vector>	// STL vectors
#include <list>		// STL lists
#include <deque>	// STL deques
#include <utility>      // STL pairs
#include <functional>   // STL functions
#include <map>          // STL maps

/** @name Vector traits.
 * Vector traits are use to allow template specialization to choose different
 * code for dense and sparse vectors.
 */
//@{

// Namespace in which all LinBox library code resides
namespace LinBox
{
	/** List of vector categories.
	 * This structure contains three structures: one relating to dense vectors,
	 * one relating to sparse vectors implemented as sequences of pairs, and 
	 * one relating to sparse vectors implemented as associative containers.
	 * These types allow us to use template specialization to use different 
	 * code for different types of vectors.
	 */
	struct VectorCategories
	{
		template <class T> struct DenseVectorTag { typedef T Traits; };
		template <class T> struct SparseSequenceVectorTag { typedef T Traits; };
		template <class T> struct SparseAssociativeVectorTag { typedef T Traits; };

	}; // struct VectorCategories

	// Helper structure used for various STL's sorts (std::list::sort and std::stable_sort) 
	// for comparison of two pairs of elements (by their first elements)
	template<class Element>
	struct SparseSequenceVectorPairLessThan: 
		public binary_function<const std::pair<size_t, Element>&, const std::pair<size_t, Element>&, bool >
	{
		bool operator() (const std::pair<size_t, Element>& p1, const std::pair<size_t, Element>& p2)
		{
			return p1.first < p2.first;
		}
	}; //struct SparseSequenceVectorPairLessThan


	/** Vector traits template structure.
	 * By default, it tries to take all information from the vector class,
	 * but it cannot usually do this.  For example, the vector_category
	 * type is not defined in STL types, so this must be done through 
	 * template specialization.
	 * @param Vector \Ref{LinBox} dense or sparse vector.
	 */
	template <class Vector> struct VectorTraits
	{
		typedef typename Vector::VectorCategory VectorCategory;
		typedef Vector VectorType;

		// These are defined for all STL vectors and sequence containers.

	};

	// Specialization for STL vectors
	template <class Element>
	struct VectorTraits< std::vector<Element> >
	{ 
		typedef std::vector<Element> VectorType;
		typedef typename VectorCategories::DenseVectorTag<VectorTraits<VectorType> > VectorCategory; 
	};

	// Specialization for STL vectors of pairs of size_t and elements
	template <class Element> 
	struct VectorTraits< std::vector< std::pair<size_t, Element> > >
	{ 
		typedef std::vector< std::pair<size_t, Element> > VectorType;
		typedef typename VectorCategories::SparseSequenceVectorTag<VectorTraits<VectorType> > VectorCategory; 

		static void sort (VectorType& v) { std::stable_sort(v.begin(), v.end(), SparseSequenceVectorPairLessThan<Element>()); }
	};

	// Specialization for STL lists of pairs of size_t and elements
	template <class Element> 
	struct VectorTraits< std::list< std::pair<size_t, Element> > >
	{ 
		typedef std::list< std::pair<size_t, Element> > VectorType;
		typedef typename VectorCategories::SparseSequenceVectorTag<VectorTraits<VectorType> > VectorCategory; 

		static void sort (VectorType& v) { v.sort(SparseSequenceVectorPairLessThan<Element>()); }
	};

	// Specialization for STL singly linked lists of pairs of size_t and elements
	template <class Element> 
	struct VectorTraits< std::deque< std::pair<size_t, Element> > >
	{ 
		typedef std::deque< std::pair<size_t, Element> > VectorType;
		typedef typename VectorCategories::SparseSequenceVectorTag<VectorTraits<VectorType> > VectorCategory; 

		static void sort (VectorType& v) { std::stable_sort(v.begin, v.end(), SparseSequenceVectorPairLessThan<Element>()); }
	};
  
	// Specialization for STL maps of size_t and elements
	template <class Element> 
	struct VectorTraits< std::map<size_t, Element> >
	{ 
		typedef std::map<size_t, Element> VectorType;
		typedef typename VectorCategories::SparseAssociativeVectorTag<VectorTraits<VectorType> > VectorCategory; 
	};

	// Namespace containing some useful generic functions

	namespace VectorWrapper 
	{
		template <class Field, class Vector, class Trait>
		inline typename Field::Element &refSpecialized
			(Vector &v, size_t i, VectorCategories::DenseVectorTag<Trait> tag)
			{ return v[i]; }

		template <class Field, class Vector, class Trait>
		inline typename Field::Element &refSpecialized
			(Vector &v, size_t i, VectorCategories::SparseSequenceVectorTag<Trait> tag)
		{
			typename Vector::iterator j;
			typename Field::Element zero;

			for (j = v.begin (); j != v.end () && (*j).first < i; j++);

			if (j == v.end () || (*j).first > i) {
				F.init (zero, 0);
				v.insert (j, pair<size_t, typename Field::Element> (i, zero));
				--j;
			}

			return (*j).second;
		}

		template <class Field, class Vector, class Trait>
		inline typename Field::Element &refSpecialized
			(Vector &v, size_t i, VectorCategories::SparseAssociativeVectorTag<Trait> tag)
			{ return v[i]; }

		template <class Field, class Vector>
		inline typename Field::Element &ref (Vector &v, size_t i) 
			{ return refSpecialized<Field, Vector> (v, i, VectorTraits<Vector>::VectorCategory()); }

		template <class Field, class Vector, class Trait>
		inline const typename Field::Element &constRefSpecialized
			(Vector &v, size_t i, VectorCategories::DenseVectorTag<Trait> tag)
			{ return v[i]; }

		template <class Field, class Vector, class Trait>
		inline const typename Field::Element &constRefSpecialized
			(Vector &v, size_t i, VectorCategories::SparseSequenceVectorTag<Trait> tag)
		{
			typename Vector::iterator j;
			static typename Field::Element zero;

			for (j = v.begin (); j != v.end () && (*j).first < i; j++);

			if (j == v.end () || (*j).first > i) {
				F.init (zero, 0);
				return zero;
			}
			else
				return (*j).second;
		}

		template <class Field, class Vector, class Trait>
		inline const typename Field::Element &constRefSpecialized
			(Vector &v, size_t i, VectorCategories::SparseAssociativeVectorTag<Trait> tag)
			{ return v[i]; }

		template <class Field, class Vector>
		inline const typename Field::Element &constRef (Vector &v, size_t i) 
			{ return constRefSpecialized<Field, Vector> (v, i, VectorTraits<Vector>::VectorCategory()); }

		template <class Vector, class Trait>
		inline void ensureDimSpecialized (Vector &v, size_t n, VectorCategories::DenseVectorTag<Trait> tag)
			{ v.resize (n); }

		template <class Vector, class Trait>
		inline void ensureDimSpecialized (Vector &v, size_t n, VectorCategories::SparseSequenceVectorTag<Trait> tag)
			{}

		template <class Vector, class Trait>
		inline void ensureDimSpecialized (Vector &v, size_t n, VectorCategories::SparseAssociativeVectorTag<Trait> tag)
			{}

		template <class Vector>
		inline void ensureDim (Vector &v, size_t n) 
			{ ensureDimSpecialized (v, n, VectorTraits<Vector>::VectorCategory()); }
	}

} // namespace LinBox

//@} Vector traits
#endif // __VECTOR_TRAITS_H
