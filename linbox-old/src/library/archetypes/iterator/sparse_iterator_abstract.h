/* File: src/library/archetypes/vector/sparse_iterator_abstract.h
 * Author: William J Turner for the LinBox group
 */
 
#ifndef _SPARSE_ITERATOR_ABSTRACT_
#define _SPARSE_ITERATOR_ABSTRACT_

#include <utility>   // STL pair
#include "integer.h"

namespace LinBox // LinBox namespace in which all code resides
{

  /** Abstract sparse iterator base class
   * Used to implement sparse vector iterator archetype.
   * Actual sparse vector iterator classes are derived from this and must
   * implement all purely virtual methods.
   */
  template <class T>
  class Sparse_iterator_abstract
  {
  public:

    // Virtual object management
    virtual Sparse_iterator_abstract* clone(void) const = 0; // virtual copy constructor

    // Forward iterators
    virtual Sparse_iterator_abstract& operator=(const Sparse_iterator_abstract& iter) = 0;   // assignment
    virtual bool operator==(const Sparse_iterator_abstract& iter) const = 0; // equality
    virtual bool operator!=(const Sparse_iterator_abstract& iter) const = 0; // inequality
    virtual const pair<integer, T>& operator*(void) const = 0;        // dereference
    // operator->(void);                         // access operator
    virtual Sparse_iterator_abstract& operator++(void) = 0;             // prefix increment
    virtual const Sparse_iterator_abstract& operator++(int) = 0;  // postfix increment

    // Bidirectional iterators

    virtual Sparse_iterator_abstract& operator--(void) = 0;             // prefix decrement
    virtual const Sparse_iterator_abstract& operator--(int) = 0;        // postfix decrement

  protected:

    // Constructors are protected because they are needed by derived 
    // classes, but they should never be called without a derived class.
    Sparse_iterator_abstract(void) {}  // default constructor

  }; // class Sparse_iterator_abstract

} // namespace LinBox

#endif // #ifndef _SPARSE_ITERATOR_ABSTRACT_

