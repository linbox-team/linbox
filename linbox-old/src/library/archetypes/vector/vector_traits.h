/* File: src/library/archetypes/vector/vector_traits.h
 * Author: William Turner for the LinBox group
 */

#ifndef _VECTOR_TRAITS_
#define _VECTOR_TRAITS_

#include <vector>	// STL vectors
#include <list>		// STL lists
#include <deque>	// STL deques
#include <utility>	// STL pairs
#include <map>		// STL maps

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
  struct vector_categories
  {
    struct dense_vector_tag {};
    struct sparse_sequence_vector_tag {};
    struct sparse_associative_vector_tag {};

  }; // struct vector_categories

  /** Vector traits template structure.
   * By default, it tries to take all information from the vector class,
   * but it cannot usually do this.  For example, the vector_category
   * type is not defined in STL types, so this must be done through 
   * template specialization.
   * @param Vector \Ref{LinBox} dense or sparse vector.
   */
  template <class Vector> struct vector_traits
  {
    typedef typename Vector::vector_category vector_category;

    // These are defined for all STL vectors and sequence containers.

  }; // template <class Vector> struct vector_traits

  // Specialization for STL vectors
  template <class Element> struct vector_traits< std::vector<Element> >
  { typedef typename vector_categories::dense_vector_tag vector_category; };

  // Specialization for STL vectors of pairs of size_t and elements
  template <class Element> 
    struct vector_traits< std::vector< std::pair<size_t, Element> > >
  { typedef typename vector_categories::sparse_sequence_vector_tag vector_category; };

  // Specialization for STL lists of pairs of size_t and elements
  template <class Element> 
    struct vector_traits< std::list< std::pair<size_t, Element> > >
  { typedef typename vector_categories::sparse_sequence_vector_tag vector_category; };

  // Specialization for STL singly linked lists of pairs of size_t and elements
  template <class Element> 
    struct vector_traits< std::deque< std::pair<size_t, Element> > >
  { typedef typename vector_categories::sparse_sequence_vector_tag vector_category; };
  
  // Specialization for STL maps of size_t and elements
  template <class Element> 
    struct vector_traits< std::map<size_t, Element> >
  { typedef typename vector_categories::sparse_associative_vector_tag vector_category; };

} // namespace LinBox

//@} Vector traits
#endif // _VECTOR_TRAITS_
