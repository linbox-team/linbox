/* File: src/library/archetypes/iterator/iterator_archetype.h
 * Author: William J Turner for the LinBox group
 */
 
#ifndef _ITERATOR_ARCHETYPE_
#define _ITERATOR_ARCHETYPE_

#include "integer.h"
#include "iterator_abstract.h"

/// Namespace in which all LinBox library code resides
namespace LinBox
{

  /// Forward declaration
  template <class T> class Vector_archetype;

  /** Iterator archetype.
   * Archetype for the iterator iterator common object interface to LinBox.
   *
   * Like the LinBox field common object interface, the common object 
   * interface for LinBox iterators requires the copy constructor. 
   * This means the base object cannot be the archetype.
   *
   * The class contains a private pointer to an abstract base class iterators.
   * This abstract base class, in turn, uses virtual member functions to define 
   * all operations through derived classes.
   * 
   * Because of their use of virtual member funtions, these archetypes can be 
   * inefficient.
   */
  template <class T>
  class Iterator_archetype
  {
  public:

    /** @name Common object interface.
     * These methods, members, and types are required of all LinBox iterators.
     */
    //@{
    
    /** Copy constructor.
     * Constructs Iterator_archetype by copying the iterator.
     * In this implementation, iterator to which iter_ptr points is copied.
     * @param iter constant reference to a Iterator_archetype object
     */
    Iterator_archetype(const Iterator_archetype& iter) 
    { iter_ptr = iter.iter_ptr->clone(); }

    /** Assignment operator.
     * Assigns iterator iter to iterator.
     * In this implementation, iterator to which iter_ptr points is copied.
     * @return Iterator_archetype reference to self
     * @param iter constant Iterator_archetype reference
     */
    Iterator_archetype& operator= (const Iterator_archetype& iter)
    {
      if (this != &iter) // quard against self-assignment
      {
	delete iter_ptr;
	iter_ptr = iter.iter_ptr->clone();
      }

      return *this;
    }

    /** Equality operator.
     * Tests the equality of two iterators.
     * In this implementation, tests the equality of the two iterators
     * to which iter_ptr points.
     * @return bool true if equal, false if not
     * @param iter constant Vector_abstract reference
     */
    bool operator== (const Iterator_archetype& iter) const
    { return *iter_ptr == *iter.iter_ptr; }

    /** Inequality operator.
     * Tests the inequality of two iterators.
     * In this implementation, tests the inequality of the two iterators
     * to which iter_ptr points.
     * @return bool false if equal, true if not
     * @param iter constant Vector_abstract reference
     */
    bool operator!= (const Iterator_archetype& iter) const
    { return *iter_ptr != *iter.iter_ptr; }

    //@} Public object management

    /** Forward iterator operators.
     * These operators, together with the object management operators, are
     * required of all STL forward iterators.
     */
    //@{
    
    /** Derefernce operator.
     * This operator returns a constant reference to the value stored in the
     * given position.  The constant reference is required because the vector
     * may not actually store the value in memory.  For example, a dense
     * representation of a vector may not store all of the zero elements.
     * @return constant reference to element
     */
    const T& operator*(void) const { return *(*iter_ptr); } 

    /** Access operator.
     */
    Iterator_archetype* operator->(void) { return this; }

    /** Pre-fix increment operator.
     * Increments iterator before return reference to self.
     * @return reference to iterator after increment
     */
    Iterator_archetype<T>& operator++(void)
    {
      ++(*iter_ptr);
      return *this;
    }
     
    /** Post-fix increment operator.
     * Increments iterator after return reference to self.
     * Required by abstract base class.
     * @return reference to abstract iterator before increment
     * @param int marks operator as post-fix
     */
    const Iterator_archetype<T>& operator++(int)
    {
      Iterator_archetype* temp_ptr = new Iterator_archetype(*this);
      ++(*iter_ptr);
      return *temp_ptr;
    }
 
    //@} Forward iterator operators

    /** Bidirectional iterator operators.
     * These operators, together with the object management and forward
     * iterator operators, are required of all STL bidirectional iterators.
     */
    //@{

    /** Pre-fix decrement operator.
     * Decrements iterator before return reference to self.
     * Required by abstract base class.
     * @return reference to abstract iterator after decrement
     */
    Iterator_archetype<T>& operator--(void)
    {
      --(*iter_ptr);
      return *this;
    }
 
    /** Post-fix decrement operator.
     * Decrements iterator after return reference to self.
     * Required by abstract base class.
     * @return reference to abstract iterator before decrement
     * @param int marks operator as post-fix
     */
    const Iterator_archetype<T>& operator--(int)
    {
      Iterator_archetype* temp_ptr = newIterator_archetype(*this);
      --(*iter_ptr);
      return *temp_ptr;
    }

    //@} Bidirectional iterators

    /** Random access iterator operators.
     * These operators, together with the object management, forward
     * iterator and bidirectional iterator operators, are required of 
     * all STL random acess iterators.
     *
     * These are not implemented yet because of problems related to the size
     * and difference types.
     *
     * These all have implementation problems because there is no guarantee
     * that integer can be converted to the correct type.
     */
    //@{

/*
    Iterator_archetype<T>& operator+=(integer n);
    const Iterator_archetype<T>& operator+(integer n) const;
    Iterator_archetype<T>& operator-=(integer n);
    const Iterator_archetype<T>& operator-(integer n) const;
    integer operator-(const Iterator_archetype<T>& iter) const;
    T operator[](integer n);
    bool operator<(const Iterator_archetype<T>& iter) const;
    bool operator>(const Iterator_archetype<T>& iter) const;
    bool operator>=(const Iterator_archetype<T>& iter) const;
    bool operator<=(const Iterator_archetype<T>& iter) const;
*/
    //@} Random access iterator operators

    //@} Common object interface

    /** @name Implementation specific methods.
     * These methods are not required of all LinBox vectors, but are only
     * included here for this particular implementation of the archetype.
     */
    //@{

    /** Constructor from pointer to Iterator_abstract object.
     * Copies Iterator_abstract by calling its clone function, thus
     * copying the derived class.
     * @param vect pointer to Iterator_abstract object
     */
    Iterator_archetype(Iterator_abstract<T>* init_ptr) 
    { iter_ptr = init_ptr->clone(); }

    //@} Implementation specific methods

  private:

    /// Pointer to Iterator_abstract through which iterator is implemented.
    Iterator_abstract<T>* iter_ptr;

    /// Friend declaration for Vector_archetype
    friend Vector_archetype<T>;
    
  }; // class Iterator_archetype

} // namespace LinBox

#endif // _ITERATOR_ARCHETYPE

