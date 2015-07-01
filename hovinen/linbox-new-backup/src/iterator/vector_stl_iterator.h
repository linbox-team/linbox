/* File: src/library/objects/vector/vector_stl_iterator.h
 * Author: William Turner for the LinBox group
 */

#ifndef _ITERATOR_STL_
#define _ITERATOR_STL_

#include <iostream>
#include <vector>
#include "integer.h"
#include "iterator_abstract.h"

/// Namespace in which all LinBox library code resides
namespace LinBox
{

  /// Forward declaration
  template <class T> class Vector_stl; // forward declaration

  /** LinBox iterator wrapper for STL vector iterators.
   * Wraps an STL vector iterator as a LinBox vector and a subclass of
   * Iterator_abstract so it can be used in LinBox code as a vector iterator.
   * This class implements all purely virtual functions from the class
   * Iterator_abstract.
   */
  template <class T>
  class Vector_stl_iterator : public Iterator_abstract<T>
  {
  public:

    /** @name Public object management.
     * Methods and operators used for object management.
     * Constructors are protected because they are required by
     * the derived class, but they should never be called without one.
     */
    //@{

    /** Default constructor.
     * Constructors are not derived from Iterator_abstract because they cannot
     * be virtual, but they are still required.  They use the protected
     * default constructor of the abstract base class.
     */
    Vector_stl_iterator(void) {}

    /** Copy constructor.
     * Constructors are not derived from Iterator_abstract because they cannot
     * be virtual, but they are still required.  They use the protected
     * default constructor of the abstract base class.
     * @param iter constant reference to a Vector_stl_iterator object
     */
    Vector_stl_iterator(const Vector_stl_iterator& iter) 
      : _iterator(iter._iterator) {}
    
    /** Constructor from stl vector iterator.
     * Constructors are not derived from Iterator_abstract because they cannot
     * be virtual, but they are still required.  They use the protected
     * default constructor of the abstract base class.
     * @param vect constant reference to an STL vector iterator
     */
    Vector_stl_iterator(const std::vector<T>::iterator& iter) 
      : _iterator(iter) {}

 
    /** Virtual copy constructor.
     * Required because copy constructors cannot be virtual.
     * Required by abstract base class.
     * @return Iterator_abstract pointer
     */
    Iterator_abstract<T>* clone(void) const 
    { return new Vector_stl_iterator(*this); }

    /** Assignment operator.
     * Required by abstract base class.
     * @return Iterator_abstract reference to self
     * @param vect constant Iterator_abstract reference
     */
    Iterator_abstract<T>& operator=(const Iterator_abstract<T>& iter)
    {
      if (this != &iter) // guard against self-assignment
	_iterator = 
	  static_cast<const Vector_stl_iterator&>(iter)._iterator;
      return *this;
    }
 
    /** Equality operator.
     * Required by abstract base class.
     * @return bool true if equal, false if not
     * @param vect constant Iterator_abstract reference
     */
    bool operator==(const Iterator_abstract<T>& iter) const
    { 
      return _iterator == 
	static_cast<const Vector_stl_iterator&>(iter)._iterator;
    }

    /** Inequality operator.
     * Required by abstract base class.
     * @return bool false if equal, true if not
     * @param vect constant Iterator_abstract reference
     */
    bool operator!=(const Iterator_abstract<T>& iter) const
    { 
      return _iterator != 
	static_cast<const Vector_stl_iterator&>(iter)._iterator;
    }

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
     * Required by abstract base class.
     * @return constant reference to element
     */
    const T& operator*(void) const { return *_iterator; } 

    /** Access operator.
     * Required by abstract base class.
     */
    Vector_stl_iterator* operator->(void) { return this; }

    /** Pre-fix increment operator.
     * Increments iterator before return reference to self.
     * Required by abstract base class.
     * @return reference to abstract iterator after increment
     */
    Iterator_abstract<T>& operator++(void)
    {
      ++_iterator;
      return *this;
    }
     
    /** Post-fix increment operator.
     * Increments iterator after return reference to self.
     * Required by abstract base class.
     * @return reference to abstract iterator before increment
     * @param int marks operator as post-fix
     */
    const Iterator_abstract<T>& operator++(int)
    {
      Vector_stl_iterator* temp_ptr = new Vector_stl_iterator(*this);
      ++_iterator;
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
    Iterator_abstract<T>& operator--(void)
    {
      --_iterator;
      return *this;
    }
 
    /** Post-fix decrement operator.
     * Decrements iterator after return reference to self.
     * Required by abstract base class.
     * @return reference to abstract iterator before decrement
     * @param int marks operator as post-fix
     */
    const Iterator_abstract<T>& operator--(int)
    {
      Vector_stl_iterator* temp_ptr  = new Vector_stl_iterator(*this);
      --_iterator;
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
     */
    //@{

/*
    Iterator_abstract<T>& operator+=(integer n);
    const Iterator_abstract<T>& operator+(integer n) const;
    Iterator_abstract<T>& operator-=(integer n);
    const Iterator_abstract<T>& operator-(integer n) const;
    integer operator-(const Iterator_abstract<T>& iter) const;
    T operator[](integer n);
    bool operator<(const Iterator_abstract<T>& iter) const;
    bool operator>(const Iterator_abstract<T>& iter) const;
    bool operator>=(const Iterator_abstract<T>& iter) const;
    bool operator<=(const Iterator_abstract<T>& iter) const;
*/
    //@} Random access iterator operators

  private:

    /// Wrapped stl vector iterator
    std::vector<T>::iterator _iterator;

    /// Declare Vector_stl<T> as a friend
    friend Vector_stl<T>;

  }; // class Vector_stl_iterator

} // namespace LinBox

#endif // _VECTOR_STL_ITERATOR_

