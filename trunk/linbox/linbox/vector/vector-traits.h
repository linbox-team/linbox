/* linbox/vector/vector-traits.h
 * Copyright (C) 1999-2001 William J Turner,
 *               2001-2002 Bradford Hovinen
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
 * see COPYING for license details
 *
 * ------------------------------------ 
 */

#ifndef __LINBOX_vector_traits_H
#define __LINBOX_vector_traits_H

#include <vector>	// STL vectors
#include <list>		// STL lists
#include <deque>	// STL deques
#include <utility>      // STL pairs
#include <functional>   // STL functions
#include <map>          // STL maps
#include <algorithm>    // STL algorithms

#include "linbox/field/archetype.h"
#include "linbox/field/rebind.h"

namespace LinBox
{

/** @name Vector traits.
 * Vector traits are use to allow template specialization to choose different
 * code for dense and sparse vectors.
	\ingroup vector
 */
//@{

	/** \brief List of vector categories.

	 * This structure contains three structures: one relating to dense vectors,
	 * one relating to sparse vectors implemented as sequences of pairs, and 
	 * one relating to sparse vectors implemented as associative containers.
	 * These types allow us to use template specialization to use different 
	 * code for different types of vectors.
	 */
        struct VectorCategories
        {
            struct GenericVectorTag { 
                friend std::ostream& operator<< (std::ostream& o, 
                                                 const GenericVectorTag& ) {
                    return o << "GenericVectorTag"; 
                } 
            };
            
                // These are valid for GF2 only

            struct SparseZeroOneVectorTag : public GenericVectorTag { 
                friend std::ostream& operator<< (std::ostream& o, 
                                                 const SparseZeroOneVectorTag& ) { 
                    return o << "SparseZeroOneVectorTag"; 
                } 
            };
                // These are valid for all fields
            struct DenseVectorTag : public SparseZeroOneVectorTag { 
                    // Inherits from SparseZeroOneVectorTag:
                    // This simplifies gf2 management of vectors
                    // and makes Dense vector also a generic vector
                friend std::ostream& operator<< (std::ostream& o, 
                                                 const DenseVectorTag& ) { 
                    return o << "DenseVectorTag"; 
                } 
            };
            
            struct SparseVectorTag : public GenericVectorTag { 
                friend std::ostream& operator<< (std::ostream& o, 
                                                 const SparseVectorTag& ) { 
                    return o << "SparseVectorTag"; 
                } 
            };
            
            struct SparseSequenceVectorTag : public SparseVectorTag { 
                friend std::ostream& operator<< (std::ostream& o, 
                                                 const SparseSequenceVectorTag& ) { 
                    return o << "SparseSequenceVectorTag"; 
                } 
            };
            struct SparseAssociativeVectorTag : public SparseVectorTag { 
                friend std::ostream& operator<< (std::ostream& o, 
                                                 const SparseAssociativeVectorTag& ) { 
                    return o << "SparseAssociativeVectorTag"; 
                } 
            };
            struct SparseParallelVectorTag : public SparseVectorTag { 
                friend std::ostream& operator<< (std::ostream& o, 
                                                 const SparseParallelVectorTag& ) { 
                    return o << "SparseParallelVectorTag"; 
                } 
            };

            struct DenseZeroOneVectorTag : public DenseVectorTag { 
                friend std::ostream& operator<< (std::ostream& o, 
                                                 const DenseZeroOneVectorTag& ) { 
                    return o << "DenseZeroOneVectorTag"; 
                } 
            };

        };

	// Helper structure used for various STL's sorts (std::list::sort and std::stable_sort) 
	// for comparison of two pairs of elements (by their first elements)
	template<class Element>
	struct SparseSequenceVectorPairLessThan :
		public std::binary_function<const std::pair<size_t, Element>&, const std::pair<size_t, Element>&, bool >
	{
		bool operator() (const std::pair<size_t, Element>& p1, const std::pair<size_t, Element>& p2)
			{ return p1.first < p2.first; }
	};


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
		typedef typename VectorCategories::DenseVectorTag VectorCategory; 
	};

	// Specialization for STL vectors of pairs of size_t and elements
	template <class Element> 
	struct VectorTraits< std::vector< std::pair<size_t, Element> > >
	{ 
		typedef std::vector< std::pair<size_t, Element> > VectorType;
		typedef typename VectorCategories::SparseSequenceVectorTag VectorCategory; 

		static void sort (VectorType& v) { std::stable_sort(v.begin(), v.end(), SparseSequenceVectorPairLessThan<Element>()); }
	};

	// Specialization for STL lists of pairs of size_t and elements
	template <class Element> 
	struct VectorTraits< std::list< std::pair<size_t, Element> > >
	{ 
		typedef std::list< std::pair<size_t, Element> > VectorType;
		typedef typename VectorCategories::SparseSequenceVectorTag VectorCategory; 

		static void sort (VectorType& v) { v.sort(SparseSequenceVectorPairLessThan<Element>()); }
	};

	// Specialization for STL singly linked lists of pairs of size_t and elements
	template <class Element> 
	struct VectorTraits< std::deque< std::pair<size_t, Element> > >
	{ 
		typedef std::deque< std::pair<size_t, Element> > VectorType;
		typedef typename VectorCategories::SparseSequenceVectorTag VectorCategory; 

		static void sort (VectorType& v) { std::stable_sort(v.begin, v.end(), SparseSequenceVectorPairLessThan<Element>()); }
	};
  
	// Specialization for STL maps of size_t and elements
	template <class Element> 
	struct VectorTraits< std::map<size_t, Element> >
	{ 
		typedef std::map<size_t, Element> VectorType;
		typedef typename VectorCategories::SparseAssociativeVectorTag VectorCategory; 
	};

	// Specialization for an STL pair of an STL vector of size_t's and an STL vector of elements
	template <class Element> 
	struct VectorTraits< std::pair<std::vector<size_t>, std::vector<Element> > >
	{ 
		typedef std::pair<std::vector<size_t>, std::vector<Element> > VectorType;
		typedef typename VectorCategories::SparseParallelVectorTag VectorCategory; 
	};

	// Namespace containing some useful generic functions

	namespace VectorWrapper 
	{
		template <class T>
		class CompareSparseEntries
		{
		    public:
                        template<typename PairType>
			inline bool operator () (const PairType &i, const size_t j) const
				{ return i.first < j; }
		};

		template <class Field, class Vector>
		inline typename Field::Element &refSpecialized
			(Vector &v, size_t i, VectorCategories::DenseVectorTag)
			{ return v[i]; }

		template <class Field, class Vector>
		inline typename Field::Element &refSpecialized
			(Vector &v, size_t i, VectorCategories::SparseSequenceVectorTag)
		{
			static typename Field::Element zero;
			typename Vector::iterator j;

			if (v.size () == 0) {
				v.push_back (std::pair <size_t, typename Field::Element> (i, zero));
				return v[0].second;
			}

			j = std::lower_bound (v.begin (), v.end (), i, CompareSparseEntries<typename Field::Element> ());

			if (j == v.end () || j->first != i)
				j = v.insert (j, std::pair <size_t, typename Field::Element> (i, zero));

			return j->second;
		}

		template <class Field, class Vector>
		inline typename Field::Element &refSpecialized
			(Vector &v, size_t i, VectorCategories::SparseAssociativeVectorTag)
			{ return v[i]; }

		template <class Field, class Vector>
		inline typename Field::Element &refSpecialized
			(Vector &v, size_t i, VectorCategories::SparseParallelVectorTag)
		{
			static typename Field::Element zero;
			typename Vector::first_type::iterator j_idx;
			typename Vector::second_type::iterator j_elt;

			if (v.first.size () == 0) {
				v.first.push_back (i);
				v.second.push_back (zero);
				return v.second.front ();
			}

			j_idx = std::lower_bound (v.first.begin (), v.first.end (), i);
			j_elt = v.second.begin () + (j_idx - v.first.begin ());

			if (j_idx == v.first.end () || *j_idx != i) {
				v.first.insert (j_idx, i);
				j_elt = v.second.insert (j_elt, zero);
			}

			return *j_elt;
		}

		template <class Field, class Vector>
		inline typename Field::Element &ref (Vector &v, size_t i) 
			{ return refSpecialized<Field, Vector> (v, i, typename VectorTraits<Vector>::VectorCategory()); }

		template <class Field, class Vector>
		inline const typename Field::Element &constRefSpecialized
			(Vector &v, size_t i, VectorCategories::DenseVectorTag)
			{ return v[i]; }

		template <class Field, class Vector>
		inline const typename Field::Element &constRefSpecialized
			(Vector &v, size_t i, VectorCategories::SparseSequenceVectorTag)
		{
			static typename Field::Element zero;
			typename Vector::const_iterator j;

			if (v.size () == 0)
				return zero;

			j = std::lower_bound (v.begin (), v.end (), i, CompareSparseEntries<typename Field::Element> ());

			if (j == v.end () || j->first != i)
				return zero;
			else
				return j->second;
		}

		template <class Field, class Vector>
		inline const typename Field::Element &constRefSpecialized
			(Vector &v, size_t i, VectorCategories::SparseAssociativeVectorTag)
			{ return v[i]; }

		template <class Field, class Vector>
		inline typename Field::Element &constRefSpecialized
			(Vector &v, size_t i, VectorCategories::SparseParallelVectorTag)
		{
			static typename Field::Element zero;
			typename Vector::first_type::iterator j_idx;
			typename Vector::second_type::iterator j_elt;

			if (v.first.size () == 0)
				return zero;

			j_idx = std::lower_bound (v.first.begin (), v.first.end (), i);

			if (j_idx == v.first.end () || *j_idx != i)
				return zero;
			else {
				j_elt = v.second.begin () + (j_idx - v.first.begin ());
				return *j_elt;
			}
		}

		template <class Vector>
		inline void ensureDimSpecialized (Vector &v, size_t n, VectorCategories::DenseZeroOneVectorTag)
			{ v.resize (n); }

		template <class Vector>
		inline void ensureDimSpecialized (Vector &v, size_t n, VectorCategories::SparseZeroOneVectorTag)
			{}

		template <class Field, class Vector>
		inline const typename Field::Element &constRef (Vector &v, size_t i) 
			{ return constRefSpecialized<Field, Vector> (v, i, typename VectorTraits<Vector>::VectorCategory()); }

		template <class Vector>
		inline void ensureDimSpecialized (Vector &v, size_t n, VectorCategories::DenseVectorTag)
			{ v.resize (n); }

		template <class Vector>
		inline void ensureDimSpecialized (Vector &v, size_t n, VectorCategories::SparseSequenceVectorTag)
			{}

		template <class Vector>
		inline void ensureDimSpecialized (Vector &v, size_t n, VectorCategories::SparseAssociativeVectorTag)
			{}

		template <class Vector>
		inline void ensureDimSpecialized (Vector &v, size_t n, VectorCategories::SparseParallelVectorTag)
			{}

		template <class Vector>
		inline void ensureDim (Vector &v, size_t n) 
			{ ensureDimSpecialized (v, n, typename VectorTraits<Vector>::VectorCategory()); }
	} // Namespace VectorWrapper

	// Now we create some "canonical" vector types, so that users don't 
	// always have to typedef everything

	/** Canonical vector types
	 *
	 * This class includes some typedefs that avoid the necessity to typedef
	 * the vector type whenever it is used. In a typical case, one would say
	 * Vector<Field>::Dense for a dense vector and Vector<Field>::Sparse for
	 * a sparse vector.
	 */

	template <class Element>
	struct RawVector 
	{
            typedef std::vector<Element> Dense;
            typedef std::pair<std::vector<size_t>, std::vector<Element> > Sparse;
            typedef std::vector<std::pair<size_t, Element> > SparseSeq;
            typedef std::map<size_t, Element> SparseMap;
            typedef std::pair<std::vector<size_t>, std::vector<Element> > SparsePar;
            
            template<class VType> static size_t size(const VType& d) { return d.size(); }        
            static size_t size(const Sparse& d) { return d.first.size(); }        
//             static size_t size(const SparsePar& d) { return d.first.size(); }      

            template<class VType>
            static VType& convert(VType& u, const VType& v) { return u = v; }
            
            static Dense& convert(Dense& u, const SparseSeq& v) {
                for(typename SparseSeq::const_iterator sit = v.begin(); sit != v.end(); ++sit) {
                    if (sit->first > u.size()) u.resize(sit->first);
                    u[sit->first] = sit->second;
                }
                return u;
            }
            static Sparse& convert(Sparse& u, const SparseSeq& v) {
                u.first.resize(v.size());
                u.second.resize(v.size());
                typename SparseSeq::const_iterator ssit = v.begin();
                std::vector<size_t>::iterator vit = u.first.begin();
                typename std::vector<Element>::iterator eit = u.second.begin();
                for(; ssit != v.end(); ++ssit, ++vit, ++eit) {
                    *vit = ssit->first;
                    *eit = ssit->second;
                }
                return u;
            }
             static SparseSeq& convert(SparseSeq& u, const Sparse& v) {
                u.resize(v.first.size());
                typename SparseSeq::iterator ssit = u.begin();
                std::vector<size_t>::const_iterator vit = v.first.begin();
                typename std::vector<Element>::const_iterator eit = v.second.begin();
                for(; ssit != u.end(); ++ssit, ++vit, ++eit) {
                    ssit->first = *vit;
                    ssit->second = *eit;
                }
                return u;
            }
                  
                // TODO : other convertions
            
  	};

	template <class Field>
	struct Vector : public RawVector<typename Field::Element>
	{
            typedef typename Field::Element Element;
            
            template<class U>
            struct rebind {
                typedef Vector<U> other;
            };
	};

 
    	template<class T, class U>
    	struct Rebind< std::vector<T>, U > {
        	typedef typename Vector<U>::Dense other;
    	};

    
	template<class T, class U>
	struct Rebind< std::pair<std::vector<size_t>, std::vector<T> >,U > {
	    typedef typename Vector<U>::Sparse other;
	};

	template<class T, class U>
	struct Rebind< std::vector<std::pair<size_t, T> >, U > {
	    typedef typename Vector<U>::SparseSeq other;
	};

	template<class T, class U>
	struct Rebind< std::map<size_t, T>, U > {
	    typedef typename Vector<U>::SparseMap other;
	};

   
//@} Vector traits

} // namespace LinBox

#endif // __LINBOX_vector_traits_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
