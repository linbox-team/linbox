/* File: src/library/objects/field/abstract_double.h
 * Author: William Turner for the LinBox group
 */

#ifndef _ABSTRACT_DOUBLE_
#define _ABSTRACT_DOUBLE_

#include <iostream>
#include <cmath>
#include "LinBox/integer.h"
#include "LinBox/abstract_double_element.h"
#include "LinBox/field_abstract.h"
#include "LinBox/element_abstract.h"
#include "LinBox/randiter_abstract.h"

// Namespace in which all LinBox code resides
namespace LinBox 
{ 

  // Forward declarations
  class abstract_double_randIter;

  /** Abstract double LinBox field.
   * Derived class used to implement the field archetype to minimize
   * code bloat.  This class implements all purely virtual member functions
   * of the abstract base class.  This class implements a fields of doubles 
   * through a private double members of the element class.
   */
  class abstract_double : public Field_abstract
  {
  public:

    /** Element type.
     * It is derived from the class Element_abstract.
     */
    typedef abstract_double_element element;

    /** Random iterator generator type.
     * It is derived from the class RandIter_abstract.
     */
    typedef abstract_double_randIter randIter;

    /** @name Object Management
     */
    //@{
 
    /** Default constructor.
     */
    abstract_double(void) {}

    /** Copy constructor.
     * Constructs abstract_double object by copying the field.
     * This is required to allow field objects to be passed by value
     * into functions.
     * @param  F abstract_double object.
     */
    abstract_double(const abstract_double& F) {}
 
    /** Virtual copy constructor.
     * Required because constructors cannot be virtual.
     * Passes construction on to derived classes.
     * This function is not part of the common object interface.
     * @return pointer to new object in dynamic memory.
     */
    Field_abstract* clone(void) const { return new abstract_double(*this); }

    /** Assignment operator.
     * Required by abstract base class.
     * @return reference to Field_abstract object for self
     * @param F constant reference to Field_abstract object
     */
    Field_abstract& operator= (const Field_abstract& F)
    { return *this; }

    /** Initialization of field base element from an integer.
     * Behaves like C++ allocator construct.
     * This function assumes the output field base element x has already been
     * constructed, but that it is not already initialized.
     * This is not a specialization of the template function because
     * such a specialization is not allowed inside the class declaration.
     * @return reference to field base element.
     * @param x field base element to contain output (reference returned).
     * @param y integer.
     */
    Element_abstract& init(Element_abstract& x, const integer& y) const
    {
      static_cast<abstract_double_element&>(x)._elem = static_cast<double>(y);
      return x;
    }
 
   /** Conversion of field base element to a template class T.
     * This function assumes the output field base element x has already been
     * constructed, but that it is not already initialized.
     * @return reference to template class T.
     * @param x template class T to contain output (reference returned).
     * @param y constant field base element.
     */
    integer& convert(integer& x, const Element_abstract& y) const
    { 
      return 
	x = static_cast<integer>(floor(static_cast<const abstract_double_element&>(y)._elem)); 
    }
 
    /** Assignment of one field base element to another.
     * This function assumes both field base elements have already been
     * constructed and initialized.
     * @return reference to x
     * @param  x field base element (reference returned).
     * @param  y field base element.
     */
    Element_abstract& assign(Element_abstract& x, 
			     const Element_abstract& y) const
    {
      static_cast<abstract_double_element&>(x)._elem 
	= static_cast<const abstract_double_element&>(y)._elem;
      return x;
    }

    /** Cardinality.
     * Return integer representing cardinality of the domain.
     * Returns a non-negative integer for all domains with finite
     * cardinality, and returns -1 to signify a domain of infinite
     * cardinality.
     * @return integer representing cardinality of the domain
     */
    const integer& cardinality(void) const { return *(new integer(-1)); }
 
    /** Characteristic.
     * Return integer representing characteristic of the domain.
     * Returns a positive integer to all domains with finite characteristic,
     * and returns 0 to signify a domain of infinite characteristic.
     * @return integer representing characteristic of the domain.
     */
    const integer& characteristic(void) const { return *(new integer(0)); }

    //@} Object Management

    /** @name Arithmetic Operations
     * x <- y op z; x <- op y
     * These operations require all elements, including x, to be initialized
     * before the operation is called.  Uninitialized field base elements will
     * give undefined results.
     */
    //@{

    /** Equality of two elements.
     * This function assumes both field base elements have already been
     * constructed and initialized.
     * @return boolean true if equal, false if not.
     * @param  x field base element
     * @param  y field base element
     */
    bool areEqual(const Element_abstract& x, const Element_abstract& y) const
    {
      return static_cast<const abstract_double_element&>(x)._elem
		== static_cast<const abstract_double_element&>(y)._elem;
    }

    /** Addition.
     * x = y + z
     * This function assumes all the field base elements have already been
     * constructed and initialized.
     * @return reference to x.
     * @param  x field base element (reference returned).
     * @param  y field base element.
     * @param  z field base element.
     */
    Element_abstract& add(Element_abstract& x,
  			  const Element_abstract& y,
  			  const Element_abstract& z) const
    {
      static_cast<abstract_double_element&>(x)._elem
	= static_cast<const abstract_double_element&>(y)._elem
		+ static_cast<const abstract_double_element&>(z)._elem;
      return x;
    }
 
    /** Subtraction.
     * x = y - z
     * This function assumes all the field base elements have already been
     * constructed and initialized.
     * @return reference to x.
     * @param  x field base element (reference returned).
     * @param  y field base element.
     * @param  z field base element.
     */
    Element_abstract& sub(Element_abstract& x,
  			  const Element_abstract& y,
  			  const Element_abstract& z) const
    {
      static_cast<abstract_double_element&>(x)._elem
	= static_cast<const abstract_double_element&>(y)._elem
		- static_cast<const abstract_double_element&>(z)._elem;
      return x;
    }
 
    /** Multiplication.
     * x = y * z
     * This function assumes all the field base elements have already been
     * constructed and initialized.
     * @return reference to x.
     * @param  x field base element (reference returned).
     * @param  y field base element.
     * @param  z field base element.
     */
    Element_abstract& mul(Element_abstract& x,
  			  const Element_abstract& y,
  			  const Element_abstract& z) const
    {
      static_cast<abstract_double_element&>(x)._elem
	= static_cast<const abstract_double_element&>(y)._elem
		* static_cast<const abstract_double_element&>(z)._elem;
      return x;
    }
 
    /** Division.
     * x = y / z
     * This function assumes all the field base elements have already been
     * constructed and initialized.
     * @return reference to x.
     * @param  x field base element (reference returned).
     * @param  y field base element.
     * @param  z field base element.
     */
    Element_abstract& div(Element_abstract& x,
  			  const Element_abstract& y,
  			  const Element_abstract& z) const
    {
      static_cast<abstract_double_element&>(x)._elem
	= static_cast<const abstract_double_element&>(y)._elem
		/ static_cast<const abstract_double_element&>(z)._elem;
      return x;
    }
 
    /** Additive Inverse (Negation).
     * x = - y
     * This function assumes both field base elements have already been
     * constructed and initialized.
     * @return reference to x.
     * @param  x field base element (reference returned).
     * @param  y field base element.
     */
    Element_abstract& neg(Element_abstract& x, const Element_abstract& y) const
    {
      static_cast<abstract_double_element&>(x)._elem
	= - static_cast<const abstract_double_element&>(y)._elem;
      return x;
    }
 
    /** Multiplicative Inverse.
     * x = 1 / y
     * This function assumes both field base elements have already been
     * constructed and initialized.
     * @return reference to x.
     * @param  x field base element (reference returned).
     * @param  y field base element.
     */
    Element_abstract& inv(Element_abstract& x, const Element_abstract& y) const
    {
      static_cast<abstract_double_element&>(x)._elem
	= 1 / static_cast<const abstract_double_element&>(y)._elem;
      return x;
    }

    //@} Arithmetic Operations
 
    /** @name Inplace Arithmetic Operations
     * x <- x op y; x <- op x
     */
    //@{

    /** Zero equality.
     * Test if field base element is equal to zero.
     * This function assumes the field base element has already been
     * constructed and initialized.
     * @return boolean true if equals zero, false if not.
     * @param  x field base element.
     */
    bool isZero(const Element_abstract& x) const
    { return static_cast<const abstract_double_element&>(x)._elem == 0; }
 
    /** One equality.
     * Test if field base element is equal to one.
     * This function assumes the field base element has already been
     * constructed and initialized.
     * @return boolean true if equals one, false if not.
     * @param  x field base element.
     */
    bool isOne(const Element_abstract& x) const
    { return static_cast<const abstract_double_element&>(x)._elem == 1; }

    /** Inplace Addition.
     * x += y
     * This function assumes both field base elements have already been
     * constructed and initialized.
     * @return reference to x.
     * @param  x field base element (reference returned).
     * @param  y field base element.
     */
    Element_abstract& addin(Element_abstract& x, 
			    const Element_abstract& y) const
    {
      static_cast<abstract_double_element&>(x)._elem
	+= static_cast<const abstract_double_element&>(y)._elem;
      return x;
    }
 
    /** Inplace Subtraction.
     * x -= y
     * This function assumes both field base elements have already been
     * constructed and initialized.
     * @return reference to x.
     * @param  x field base element (reference returned).
     * @param  y field base element.
     */
    Element_abstract& subin(Element_abstract& x, 
			    const Element_abstract& y) const
    {
      static_cast<abstract_double_element&>(x)._elem
	-= static_cast<const abstract_double_element&>(y)._elem;
      return x;
    }
 
    /** Inplace Multiplication.
     * x *= y
     * This function assumes both field base elements have already been
     * constructed and initialized.
     * @return reference to x.
     * @param  x field base element (reference returned).
     * @param  y field base element.
     */
    Element_abstract& mulin(Element_abstract& x, 
			    const Element_abstract& y) const
    {
      static_cast<abstract_double_element&>(x)._elem
	*= static_cast<const abstract_double_element&>(y)._elem;
      return x;
    }
 
    /** Inplace Division.
     * x /= y
     * This function assumes both field base elements have already been
     * constructed and initialized.
     * @return reference to x.
     * @param  x field base element (reference returned).
     * @param  y field base element.
     */
    Element_abstract& divin(Element_abstract& x, 
			    const Element_abstract& y) const
    {
      static_cast<abstract_double_element&>(x)._elem
	/= static_cast<const abstract_double_element&>(y)._elem;
      return x;
    }
 
    /** Inplace Additive Inverse (Inplace Negation).
     * x = - x
     * This function assumes the field base element has already been
     * constructed and initialized.
     * @return reference to x.
     * @param  x field base element (reference returned).
     */
    Element_abstract& negin(Element_abstract& x) const
    {
      static_cast<abstract_double_element&>(x)._elem
	= - static_cast<abstract_double_element&>(x)._elem;
      return x;
    }
 
    /** Inplace Multiplicative Inverse.
     * x = 1 / x
     * This function assumes the field base elementhas already been
     * constructed and initialized.
     * @return reference to x.
     * @param  x field base element (reference returned).
     */
    Element_abstract& invin(Element_abstract& x) const
    {
      static_cast<abstract_double_element&>(x)._elem
	= 1 / static_cast<abstract_double_element&>(x)._elem;
      return x;
    }

    //@} Inplace Arithmetic Operations

    /** @name Input/Output Operations */
    //@{

    /** Print field.
     * @return output stream to which field is written.
     * @param  os  output stream to which field is written.
     */
    ostream& write(ostream& os) const { return os << "Abstract double field"; }
 
    /** Read field.
     * @return input stream from which field is read.
     * @param  is  input stream from which field is read.
     */
    istream& read(istream& is) { return is; }

    /** Print field base element.
     * This function assumes the field base element has already been
     * constructed and initialized.
     * @return output stream to which field base element is written.
     * @param  os  output stream to which field base element is written.
     * @param  x   field base element.
     */
    ostream& write(ostream& os, const Element_abstract& x) const
    { return os << static_cast<const abstract_double_element&>(x)._elem; }
 
    /** Read field base element.
     * This function assumes the field base element has already been
     * constructed and initialized.
     * @return input stream from which field base element is read.
     * @param  is  input stream from which field base element is read.
     * @param  x   field base element.
     */
    istream& read(istream& is, Element_abstract& x) const
    { return is >> static_cast<abstract_double_element&>(x)._elem; }

    //@}

  }; // class abstract_double

} // namespace LinBox

#include "LinBox/abstract_double_randiter.h"

#endif // _ABSTRACT_DOUBLE_
