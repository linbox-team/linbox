/* File: src/library/archetypes/localring/localring_archetype.h
 * Author: BDS for the LinBox group (buildt from field_archetype.h
 */

#ifndef _LOCALRING_ARCHETYPE_
#define _LOCALRING_ARCHETYPE_

#include <iostream>
#include "LinBox/element_archetype.h"
#include "LinBox/localring_abstract.h"
#include "LinBox/element_abstract.h"
#include "LinBox/randiter_abstract.h"
#include "LinBox/integer.h"

// Namespace in which all LinBox library code resides
namespace LinBox
{
  // Forward declarations
  class RandIter_archetype;

  /** localring Archetype.
   * Archetype for the localring common object interface to \Ref{LinBox}.
   *
   * The \Ref{Localring_archetype} and its encapsulated
   * element class contain pointers to the \Ref{Localring_abstract}
   * and its encapsualted localring element, respectively.
   * \Ref{Localring_abstract} then uses virtual member functions to
   * define operations on its encapsulated localring element.  This localring 
   * element has no knowledge of the localring properties being used on it 
   * which means the localring object must supply these operations.
   *
   * It does not contain elements zero and one because they can be created 
   * whenever necessary, although it might be beneficial from an efficiency
   * stand point to include them.  However, because of archetype use three,
   * the elements themselves cannot be contained, but rather pointers to them.
   */
  class Localring_archetype
  {
  public:
    
    /** @name Common Object Interface for a LinBox localring.
     * These methods are required of all \Ref{LinBox} localrings.
     */
    //@{
    
    /// Element type.
    typedef Element_archetype element;

    /// Random iterator generator type.
    typedef RandIter_archetype randIter;
    
    /** @name Object Management
     * x <- convert(y)
     */
    //@{
    
    /** Copy constructor.
     * Constructs Localring_archetype object by copying the localring.
     * This is required to allow localring objects to be passed by value
     * into functions.
     * In this implementation, this means copying the localring to
     * which F._localring_ptr points, the element to which F._elem_ptr points, 
     * and the random element generator to which F._randIter_ptr points.
     * @param  F Localring_archetype object.
     */
    Localring_archetype(const Localring_archetype& F) 
    { 
      if (F._localring_ptr != 0) _localring_ptr = F._localring_ptr->clone(); 
      if (F._elem_ptr != 0) _elem_ptr = F._elem_ptr->clone();
      if (F._randIter_ptr != 0) _randIter_ptr = F._randIter_ptr->clone();
    }

    /** Destructor.
     * This destructs the localring object, but it does not destroy the localring 
     * element objects.  The destructor for each localring element must also 
     * be called.
     * In this implementation, this destroys localring by deleting localring 
     * object to which _localring_ptr points, the localring element to which 
     * _elem_ptr points, and the random element generator to which 
     * _randIter_ptr points.
     */
    ~Localring_archetype(void) 
    {
      if (_localring_ptr != 0) delete _localring_ptr;
      if (_elem_ptr != 0) delete _elem_ptr; 
      if (_randIter_ptr != 0) delete _randIter_ptr;
    }
    
    /** Assignment operator.
     * Assigns Localring_archetype object F to localring.
     * In this implementation, this means copying the localring to
     * which F._localring_ptr points, the element to which F._elem_ptr points, 
     * and the random element generator to which F._randIter_ptr points.
     * @param  F Localring_archetype object.
     */
    Localring_archetype& operator=(const Localring_archetype& F)
    {
      if (this != &F) // guard against self-assignment
      {
        if (_localring_ptr != 0) delete _localring_ptr;
        if (_elem_ptr != 0) delete _elem_ptr;
        if (_randIter_ptr != 0) delete _randIter_ptr;
        if (F._localring_ptr != 0) _localring_ptr = F._localring_ptr->clone(); 
        if (F._elem_ptr != 0) _elem_ptr = F._elem_ptr->clone();
        if (F._randIter_ptr != 0) _randIter_ptr = F._randIter_ptr->clone();
      }
      return *this;
    }
    
    /** Initialization of localring element from an integer.
     * Behaves like C++ allocator construct.
     * This function assumes the output localring element x has already been 
     * constructed, but that it is not necessarily already initialized.
     * In this implementation, this means the _elem_ptr of x exists, but
     * that it may be the null pointer.
     * @return reference to localring element.
     * @param x localring element to contain output (reference returned).
     * @param y constant reference to integer.
     */
    element& init(element& x, const integer& y) const
    {
      if (x._elem_ptr != 0) delete x._elem_ptr;
      x._elem_ptr = _elem_ptr->clone();
      _localring_ptr->init(*x._elem_ptr, y);
      return x;
    }
  
    /** Conversion of localring element to an integer.
     * This function assumes the output localring element x has already been 
     * constructed, but that it is not already initialized.
     * In this implementation, this means the _elem_ptr of y exists, and
     * that it is not the null pointer.
     * @return reference to integer.
     * @param x reference to integer to contain output (reference returned).
     * @param y constant reference to localring element.
     */
    integer& convert(integer& x, const element& y) const
    {
      _localring_ptr->convert(x, *y._elem_ptr);
      return x;
    }
    
    /** Assignment of one localring element to another.
     * This function assumes both localring elements have already been 
     * constructed and initialized.
     * In this implementation, this means for both x and y, 
     * _elem_ptr exists and does not point to null.
     * @return reference to x
     * @param  x localring element (reference returned).
     * @param  y localring element.
     */
    element& assign(element& x, const element& y) const
    {
      _localring_ptr->assign(*x._elem_ptr, *y._elem_ptr);
      return x;
    }
    
    /** Cardinality.
     * Return integer representing cardinality of the localring.
     * Returns a non-negative integer for all localrings with finite
     * cardinality, and returns -1 to signify a localring of infinite 
     * cardinality.
     * @return constant reference to integer representing cardinality 
     *	       of the localring
     */
    integer& cardinality(integer& c) const 
    { return _localring_ptr->cardinality(c); }
    
    /** Characteristic.
     * Return integer representing characteristic of the localring.
     * Returns a positive integer to all localrings with finite characteristic,
     * and returns 0 to signify a localring of infinite characteristic.
     * @return constant reference to integer representing characteristic 
     * 	       of the localring.
     */
    integer& characteristic(integer& c) const
    { return _localring_ptr->characteristic(c); }
    
    //@} Object Management
    
    /** @name Arithmetic Operations 
     * x <- y op z; x <- op y
     * These operations require all elements, including x, to be initialized
     * before the operation is called.  Uninitialized localring elements will
     * give undefined results.
     */
    //@{
    
    /** Equality of two elements.
     * This function assumes both localring elements have already been 
     * constructed and initialized.
     * In this implementation, this means for both x and y, 
     * _elem_ptr exists and does not point to null.
     * @return boolean true if equal, false if not.
     * @param  x localring element
     * @param  y localring element
     */
    bool areEqual(const element& x, const element& y) const
    { return _localring_ptr->areEqual(*x._elem_ptr, *y._elem_ptr); }
    
    /** Addition.
     * x = y + z
     * This function assumes all the localring elements have already been 
     * constructed and initialized.
     * In this implementation, this means for x, y, and z, 
     * _elem_ptr exists and does not point to null.
     * @return reference to x.
     * @param  x localring element (reference returned).
     * @param  y localring element.
     * @param  z localring element.
     */
    element& add(element& x, const element& y, const element& z) const
    {
      _localring_ptr->add(*x._elem_ptr, *y._elem_ptr, *z._elem_ptr);
      return x;
    }
    
    /** Subtraction.
     * x = y - z
     * This function assumes all the localring elements have already been 
     * constructed and initialized.
     * In this implementation, this means for x, y, and z, 
     * _elem_ptr exists and does not point to null.
     * @return reference to x.
     * @param  x localring element (reference returned).
     * @param  y localring element.
     * @param  z localring element.
     */
    element& sub(element& x, const element& y, const element& z) const
    {
      _localring_ptr->sub(*x._elem_ptr, *y._elem_ptr, *z._elem_ptr);
      return x;
    }
    
    /** Multiplication.
     * x = y * z
     * This function assumes all the localring elements have already been 
     * constructed and initialized.
     * In this implementation, this means for x, y, and z, 
     * _elem_ptr exists and does not point to null.
     * @return reference to x.
     * @param  x localring element (reference returned).
     * @param  y localring element.
     * @param  z localring element.
     */
    element& mul(element& x, const element& y, const element& z) const
    {
      _localring_ptr->mul(*x._elem_ptr, *y._elem_ptr, *z._elem_ptr);
      return x;
    }
    
    /** Divisibility
     * @return boolean true iff there exists ring element z such that xz = y.
     * @param  x localring element.
     * @param  y localring element.
     */
    bool divides(const element& x, const element& y) const 
    { return _localring_ptr->divides(*x._elem_ptr, *y._elem_ptr); }

    /** Division.
     * x = y / z
     * Division is only defined for y,z such that divides(z, y).
     * This function assumes all the localring elements have already been 
     * constructed and initialized.
     * In this implementation, this means for x, y, and z, 
     * _elem_ptr exists and does not point to null.
     * @return reference to x.
     * @param  x localring element (reference returned). 
     * Some x such that xz = y is returned.
     * All divisions of y by z will return the same value.
     * @param  y localring element.
     * @param  z localring element.  Must divide y.  
     */
    element& div(element& x, const element& y, const element& z) const
    {
      _localring_ptr->div(*x._elem_ptr, *y._elem_ptr, *z._elem_ptr);
      return x;
    }
    
    /** Additive Inverse (Negation).
     * x = - y
     * This function assumes both localring elements have already been 
     * constructed and initialized.
     * In this implementation, this means for both x and y 
     * _elem_ptr exists and does not point to null.
     * @return reference to x.
     * @param  x localring element (reference returned).
     * @param  y localring element.
     */
    element& neg(element& x, const element& y) const
    {
      _localring_ptr->neg(*x._elem_ptr, *y._elem_ptr);
      return x;
    }
    
    /** Multiplicative Inverse.
     * x = 1 / y
     * This function assumes both localring elements have already been 
     * constructed and initialized.
     * In this implementation, this means for both x and y 
     * _elem_ptr exists and does not point to null.
     * @return reference to x.
     * @param  x localring element (reference returned).
     * @param  y localring element which is a unit.
     */
    element& inv(element& x, const element& y) const
    {
      _localring_ptr->inv(*x._elem_ptr, *y._elem_ptr);
      return x;
    }
    
    //@} Arithmetic Operations
    
    /** @name Inplace Arithmetic Operations 
     * x <- x op y; x <- op x
     * These operations require all elements, including x, to be initialized
     * before the operation is called.  Uninitialized localring elements will
     * give undefined results.
     */
    //@{
    
    /** Zero equality.
     * Test if localring element is equal to zero.
     * This function assumes the localring element has already been 
     * constructed and initialized.
     * In this implementation, this means the _elem_ptr of x
     * exists and does not point to null.
     * @return boolean true if equals zero, false if not.
     * @param  x localring element.
     */
    bool isZero(const element& x) const 
    { return _localring_ptr->isZero(*x._elem_ptr); }
    
    /** One equality.
     * Test if localring element is equal to one.
     * This function assumes the localring element has already been 
     * constructed and initialized.
     * In this implementation, this means the _elem_ptr of x
     * exists and does not point to null.
     * @return boolean true if equals one, false if not.
     * @param  x localring element.
     */
    bool isOne(const element& x) const 
    { return _localring_ptr->isOne(*x._elem_ptr); }
    
    /** Invertibility
     * Test of localring element is a unit, i.e. has an inverse in the ring.
     * @return boolean true if x is a unit, false if not.
     * @param  x localring element.
     */
    bool isUnit(const element& x) const 
    { return _localring_ptr->isUnit(*x._elem_ptr); }

    /** Inplace Addition.
     * x += y
     * This function assumes both localring elements have already been 
     * constructed and initialized.
     * In this implementation, this means for both x and y 
     * _elem_ptr exists and does not point to null.
     * @return reference to x.
     * @param  x localring element (reference returned).
     * @param  y localring element.
     */
    element& addin(element& x, const element& y) const
    {
      _localring_ptr->addin(*x._elem_ptr, *y._elem_ptr);
      return x;
    }
    
    /** Inplace Subtraction.
     * x -= y
     * This function assumes both localring elements have already been 
     * constructed and initialized.
     * In this implementation, this means for both x and y 
     * _elem_ptr exists and does not point to null.
     * @return reference to x.
     * @param  x localring element (reference returned).
     * @param  y localring element.
     */
    element& subin(element& x, const element& y) const
    {
      _localring_ptr->subin(*x._elem_ptr, *y._elem_ptr);
      return x;
    }
 
    /** Inplace Multiplication.
     * x *= y
     * This function assumes both localring elements have already been 
     * constructed and initialized.
     * In this implementation, this means for both x and y 
     * _elem_ptr exists and does not point to null.
     * @return reference to x.
     * @param  x localring element (reference returned).
     * @param  y localring element.
     */
    element& mulin(element& x, const element& y) const
    {
      _localring_ptr->mulin(*x._elem_ptr, *y._elem_ptr);
      return x;
    }
    
    /** Inplace Division.
     * x /= y
     * This function assumes both localring elements have already been 
     * constructed and initialized.
     * In this implementation, this means for both x and y 
     * _elem_ptr exists and does not point to null.
     * @return reference to x.
     * @param  x localring element (reference returned).
     * @param  y localring element.
     */
    element& divin(element& x, const element& y) const
    {
      _localring_ptr->divin(*x._elem_ptr, *y._elem_ptr);
      return x;
    }
    
    /** Inplace Additive Inverse (Inplace Negation).
     * x = - x
     * This function assumes the localring element has already been 
     * constructed and initialized.
     * In this implementation, this means the _elem_ptr of x
     * exists and does not point to null.
     * @return reference to x.
     * @param  x localring element (reference returned).
     */
    element& negin(element& x) const
    {
      _localring_ptr->negin(*x._elem_ptr);
      return x;
    }
    
    /** Inplace Multiplicative Inverse.
     * x = 1 / x
     * This function assumes the localring elementhas already been 
     * constructed and initialized.
     * In this implementation, this means the _elem_ptr of x
     * exists and does not point to null.
     * @return reference to x.
     * @param  x localring element (reference returned).
     */
    element& invin(element& x) const
    {
      _localring_ptr->invin(*x._elem_ptr);
      return x;
    }
    
    //@} Inplace Arithmetic Operations
    
    /** @name Input/Output Operations */
    //@{
    
    /** Print localring.
     * @return output stream to which localring is written.
     * @param  os  output stream to which localring is written.
     */
    ostream& write(ostream& os) const { return _localring_ptr->write(os); }
    
    /** Read localring.
     * @return input stream from which localring is read.
     * @param  is  input stream from which localring is read.
     */
    istream& read(istream& is) { return _localring_ptr->read(is); }
    
    /** Print localring element.
     * This function assumes the localring element has already been 
     * constructed and initialized.
     * In this implementation, this means for the _elem_ptr for x 
     * exists and does not point to null.
     * @return output stream to which localring element is written.
     * @param  os  output stream to which localring element is written.
     * @param  x   localring element.
     */
    ostream& write(ostream& os, const element& x) const 
    { return _localring_ptr->write(os, *x._elem_ptr); }
    
    /** Read localring element.
     * This function assumes the localring element has already been 
     * constructed and initialized.
     * In this implementation, this means for the _elem_ptr for x 
     * exists and does not point to null.
     * @return input stream from which localring element is read.
     * @param  is  input stream from which localring element is read.
     * @param  x   localring element.
     */
    istream& read(istream& is, element& x) const
    { return _localring_ptr->read(is, *x._elem_ptr); }
    
    //@} Input/Output Operations
    
    //@} Common Object Interface
    
    /** @name Implementation-Specific Methods.
     * These methods are not required of all \Ref{LinBox localrings}
     * and are included only for this implementation of the archetype.
     */
    //@{
    
    /** Constructor.
     * Constructs localring from pointer to \Ref{Localring_abstract} and its
     * encapsulated element and random element generator.
     * Not part of the interface.
     * Creates new copies of localring, element, and random iterator generator
     * objects in dynamic memory.
     * @param  localring_ptr pointer to \Ref{Localring_abstract}.
     * @param  elem_ptr  pointer to \Ref{Element_abstract}, which is the
     *                   encapsulated element of \Ref{Localring_abstract}.
     * @param  randIter_ptr  pointer to \Ref{RandIter_abstract}, which is the
     *                       encapsulated random iterator generator
     *                       of \Ref{Localring_abstract}.
     */
    Localring_archetype(Localring_abstract* localring_ptr,
		    Element_abstract* elem_ptr = 0,
		    RandIter_abstract* randIter_ptr = 0)
      : _localring_ptr(localring_ptr->clone()), 
        _elem_ptr(elem_ptr->clone())//, 
//	_randIter_ptr(randIter_ptr->clone())
    {
      if (randIter_ptr != 0) _randIter_ptr = randIter_ptr->clone();
    }
    
    //@} Implementation-Specific Methods
    
  private:
    
    friend Localring_archetype::element;
    friend Localring_archetype::randIter;
    
    /** Pointer to Localring_abstract object.
     * Not part of the interface.
     * Included to allow for archetype use three.
     */
    mutable Localring_abstract* _localring_ptr;
    
    /** Pointer to Element_abstract object.
     * Not part of the interface.
     * Included to allow for archetype use three.
     */
    mutable Element_abstract* _elem_ptr;
    
    /** Pointer to RandIter_abstract object.
     * Not part of the interface.
     * Included to allow for archetype use three.
     */
    mutable RandIter_abstract* _randIter_ptr;
    
  }; // class Localring_archetype
  
} // namespace LinBox

#include "LinBox/randiter_archetype.h"

#endif // _LOCALRING_ARCHETYPE_
