/* File: src/library/objects/field/abstract_modular.h
 * Author: William Turner for the LinBox group
 */

#ifndef _ABSTRACT_MODULAR_
#define _ABSTRACT_MODULAR_

#include <iostream>
#include "LinBox/integer.h"
#include "LinBox/abstract_modular_element.h"
#include "LinBox/field_abstract.h"
#include "LinBox/element_abstract.h"
#include "LinBox/randiter_abstract.h"

// Namespace in which all LinBox code resides
namespace LinBox 
{ 

  // Forward declarations
  class abstract_modular_randIter;

  /** Abstract unparameterized field modulo prime number.
   * The primality of the modulus will not be checked, so 
   * it is the programmer's responsibility to supply a prime modulus.
   * Derived class used to implement the field archetype to minimize
   * code bloat.  This class implements all purely virtual member functions
   * of the abstract base class.  This class implements 
   * a field of unparamterized integers modulo a prime integer.
   * Field has static member to contain modulus of field.
   */
  class abstract_modular : public Field_abstract
  {
  public:

    /** Element type.
     * It is derived from the class Element_abstract, and it contains 
     * an integer _residue for the residue class.
     */
    typedef abstract_modular_element element;

    /** Random iterator generator type.
     * It is derived from the class RandIter_abstract.
     */
    typedef abstract_modular_randIter randIter;

    /** @name Object Management
     */
    //@{
 
    /** Default constructor.
     */
    abstract_modular(void) {}

    /** Constructor from an integer.
     * Sets the modulus of the field throug the static member of the 
     * element type.
     * @param value constant reference to integer prime modulus
     */
    abstract_modular(const integer& value) { _modulus = value; }

    /** Copy constructor.
     * Constructs abstract_modular object by copying the field.
     * This is required to allow field objects to be passed by value
     * into functions.
     * @param  F abstract_modular object.
     */
    abstract_modular(const abstract_modular& F) {}
 
    /** Virtual copy constructor.
     * Required because constructors cannot be virtual.
     * Passes construction on to derived classes.
     * This function is not part of the common object interface.
     * @return pointer to new object in dynamic memory.
     */
    Field_abstract* clone(void) const { return new abstract_modular(*this); }

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
      integer residue = y % _modulus;
      if (residue < 0) residue += _modulus;
      static_cast<abstract_modular_element&>(x)._residue = residue;
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
    { return x = static_cast<const abstract_modular_element>(y)._residue; }
 
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
      static_cast<abstract_modular_element&>(x)._residue
	= static_cast<const abstract_modular_element&>(y)._residue;
      return x;
    }

    /** Cardinality.
     * Return integer representing cardinality of the domain.
     * Returns a non-negative integer for all domains with finite
     * cardinality, and returns -1 to signify a domain of infinite
     * cardinality.
     * @return integer representing cardinality of the domain
     */
    const integer& cardinality(void) const 
    { return *(new integer(_modulus)); }
 
    /** Characteristic.
     * Return integer representing characteristic of the domain.
     * Returns a positive integer to all domains with finite characteristic,
     * and returns 0 to signify a domain of infinite characteristic.
     * @return integer representing characteristic of the domain.
     */
    const integer& characteristic(void) const
    { return *(new integer(_modulus)); }

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
      return static_cast<const abstract_modular_element&>(x)._residue
		== static_cast<const abstract_modular_element&>(y)._residue;
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
      static_cast<abstract_modular_element&>(x)._residue
	= ( static_cast<const abstract_modular_element&>(y)._residue
	    + static_cast<const abstract_modular_element&>(z)._residue )
	  % _modulus;
		
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
      static_cast<abstract_modular_element&>(x)._residue
	= ( static_cast<const abstract_modular_element&>(y)._residue
	    - static_cast<const abstract_modular_element&>(z)._residue )
	  % _modulus;

      if (static_cast<abstract_modular_element&>(x)._residue < 0)
	static_cast<abstract_modular_element&>(x)._residue += _modulus;

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
      static_cast<abstract_modular_element&>(x)._residue
	= ( static_cast<const abstract_modular_element&>(y)._residue
	    * static_cast<const abstract_modular_element&>(z)._residue )
	  % _modulus;

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
      element temp;
      inv(temp, static_cast<const abstract_modular_element&>(z));
    
      mul(x, y, temp);
      
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
      static_cast<abstract_modular_element&>(x)._residue
	= _modulus - static_cast<const abstract_modular_element&>(y)._residue;
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
      // The extended Euclidean algoritm
      integer x_int, y_int, q, tx, ty, temp;
      x_int = _modulus; 
      y_int = static_cast<const abstract_modular_element&>(y)._residue;
      tx = 0; 
      ty = 1;

      while(y_int != 0)
      {
	// always: gcd(modulus,residue) = gcd(x_int,y_int)
        //         sx*modulus + tx*residue = x_int
        //         sy*modulus + ty*residue = y_int
        q = x_int / y_int; // integer quotient
        temp = y_int;  y_int  = x_int  - q*y_int;  x_int  = temp;
        temp = ty; ty = tx - q*ty; tx = temp;
      }

      // now x_int = gcd(modulus,residue)
      integer residue = tx;
      if (residue < 0) residue += _modulus;
      
      static_cast<abstract_modular_element&>(x)._residue = residue;

      return x;
    } // Element_abstract& inv(Element_abstract&, const Element_abstract&) const

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
    { return static_cast<const abstract_modular_element&>(x)._residue == 0; }
 
    /** One equality.
     * Test if field base element is equal to one.
     * This function assumes the field base element has already been
     * constructed and initialized.
     * @return boolean true if equals one, false if not.
     * @param  x field base element.
     */
    bool isOne(const Element_abstract& x) const
    { return static_cast<const abstract_modular_element&>(x)._residue == 1; }

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
      add(static_cast<abstract_modular_element&>(x),
	  static_cast<const abstract_modular_element&>(x),
	  static_cast<const abstract_modular_element&>(y));
		
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
      sub(static_cast<abstract_modular_element&>(x),
	  static_cast<const abstract_modular_element&>(x),
	  static_cast<const abstract_modular_element&>(y));
		
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
      mul(static_cast<abstract_modular_element&>(x),
	  static_cast<const abstract_modular_element&>(x),
	  static_cast<const abstract_modular_element&>(y));
		
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
      div(static_cast<abstract_modular_element&>(x),
	  static_cast<const abstract_modular_element&>(x),
	  static_cast<const abstract_modular_element&>(y));
		
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
      neg(static_cast<abstract_modular_element&>(x),
	  static_cast<abstract_modular_element&>(x));

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
      inv(static_cast<abstract_modular_element&>(x),
	  static_cast<abstract_modular_element&>(x));

      return x;
    }

    //@} Inplace Arithmetic Operations

    /** @name Input/Output Operations */
    //@{

    /** Print field.
     * @return output stream to which field is written.
     * @param  os  output stream to which field is written.
     */
    ostream& write(ostream& os) const 
    { return os << " mod " << _modulus; }
 
    /** Read field.
     * @return input stream from which field is read.
     * @param  is  input stream from which field is read.
     */
    istream& read(istream& is) { return is >> _modulus; }

    /** Print field base element.
     * This function assumes the field base element has already been
     * constructed and initialized.
     * @return output stream to which field base element is written.
     * @param  os  output stream to which field base element is written.
     * @param  x   field base element.
     */
    ostream& write(ostream& os, const Element_abstract& x) const
    { 
      return os << static_cast<const abstract_modular_element&>(x)._residue
	        << " mod " << _modulus;
    }
 
    /** Read field base element.
     * This function assumes the field base element has already been
     * constructed and initialized.
     * @return input stream from which field base element is read.
     * @param  is  input stream from which field base element is read.
     * @param  x   field base element.
     */
    istream& read(istream& is, Element_abstract& x) const
    { 
      integer x_int;
      is >> x_int;

      x_int %= _modulus;
      if (x_int < 0) x_int += _modulus;

      static_cast<abstract_modular_element&>(x)._residue = x_int;
      return is; 
    }

    //@}

  private:

    /// Private static integer for modulus
    static integer _modulus;

  }; // class abstract_modular

  integer abstract_modular::_modulus;  // declare static member

} // namespace LinBox

#include "LinBox/abstract_modular_randiter.h"

#endif // _ABSTRACT_MODULAR_
