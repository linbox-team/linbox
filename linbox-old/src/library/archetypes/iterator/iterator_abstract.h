/* File: src/library/archetypes/vector/iterator_abstract.h
 * Author: William J Turner for the LinBox group
 */
 
#ifndef _ITERATOR_ABSTRACT_
#define _ITERATOR_ABSTRACT_

#include "LinBox/integer.h"

/// Namespace in which all LinBox library code resides
namespace LinBox
{

  /** Abstract iterator base class.
   * Used to implement vector iterator archetype.
   * Like the LinBox field common object interface, the common object 
   * interface for LinBox iterators requires the copy constructor. 
   * This means the base object cannot be the archetype.
   * Actual iterator classes are derived from this and must implement all 
   * purely virtual methods.
   */
  template <class T>
  class Iterator_abstract
  {
  public:

    /** @name Public object management.
     * Methods and operators used for object management.
     * Constructors are protected because they are required by
     * the derived class, but they should never be called without one.
     */
    //@{

    /** Virtual copy constructor.
     * Required because copy constructors cannot be virtual.
     * Purely virtual.
     * @return Iterator_abstract pointer
     */
    virtual Iterator_abstract* clone(void) const = 0;

     
    /** Destructor.
     * Required to be virtual because virtual member functions exist.
     * Not purely virtual.
     */
    virtual ~Iterator_abstract(void) {}

    /** Assignment operator.
     * Purely virtual.
     * @return Iterator_abstract reference to self
     * @param vect constant Iterator_abstract reference
     */
    virtual Iterator_abstract& operator=(const Iterator_abstract& iter) = 0;

    /** Equality operator.
     * Purely virtual.
     * @return bool true if equal, false if not
     * @param vect constant Iterator_abstract reference
     */
    virtual bool operator==(const Iterator_abstract& iter) const = 0;

    /** Inequality operator.
     * Purely virtual.
     * @return bool false if equal, true if not
     * @param vect constant Iterator_abstract reference
     */
    virtual bool operator!=(const Iterator_abstract& iter) const = 0;

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
     * Purely virtual.
     * @return constant reference to element
     */
    virtual const T& operator*(void) const = 0;

    /** Access operator.
     * Purely virtual.
     */
    // operator->(void);

    /** Pre-fix increment operator.
     * Increments iterator before return reference to self.
     * Purely virtual.
     * @return reference to abstract iterator after increment
     */
    virtual Iterator_abstract& operator++(void) = 0;

    /** Post-fix increment operator.
     * Increments iterator after return reference to self.
     * Purely virtual.
     * @return reference to abstract iterator before increment
     * @param int marks operator as post-fix
     */
    virtual const Iterator_abstract& operator++(int) = 0;

    //@} Forward iterator operators

    /** Bidirectional iterator operators.
     * These operators, together with the object management and forward
     * iterator operators, are required of all STL bidirectional iterators.
     */
    //@{

    /** Pre-fix decrement operator.
     * Decrements iterator before return reference to self.
     * Purely virtual.
     * @return reference to abstract iterator after decrement
     */
    virtual Iterator_abstract& operator--(void) = 0;

    /** Post-fix decrement operator.
     * Decrements iterator after return reference to self.
     * Purely virtual.
     * @return reference to abstract iterator before decrement
     * @param int marks operator as post-fix
     */
    virtual const Iterator_abstract& operator--(int) = 0;

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
    virtual Iterator_abstract& operator+=(integer n) = 0;             //  
    virtual const Iterator_abstract& operator+(integer n) const = 0;
    virtual Iterator_abstract& operator-=(integer n) = 0;
    virtual const Iterator_abstract& operator-(integer n) const = 0;
    virtual integer operator-(const Iterator_abstract& iter) const = 0;
    virtual T operator[](integer n) = 0;
    virtual bool operator<(const Iterator_abstract& iter) const = 0;
    virtual bool operator>(const Iterator_abstract& iter) const = 0;
    virtual bool operator>=(const Iterator_abstract& iter) const = 0;
    virtual bool operator<=(const Iterator_abstract& iter) const = 0;
*/   
    //@} Random access iterator operators

  protected:

    /** @name Protected object management
     * Constructors are protected because they are required by
     * the derived class, but they should never be called without one.
     */
    //@{

    /// Default constructor.
    Iterator_abstract(void) {}

    //@} Protected object management

  }; // class Iterator_abstract

} // namespace LinBox

#endif // #ifndef _ITERATOR_ABSTRACT_

