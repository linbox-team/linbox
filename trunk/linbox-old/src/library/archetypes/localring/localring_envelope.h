/* File: src/library/archetypes/localring/localring_envelope.h
 * Author: William Turner for the LinBox group
 */

#ifndef _LOCALRING_ENVELOPE_
#define _LOCALRING_ENVELOPE_

#include <iostream>
#include "LinBox/integer.h"
#include "LinBox/element_envelope.h"
#include "LinBox/localring_abstract.h"
#include "LinBox/element_abstract.h"
#include "LinBox/randiter_abstract.h"

// Namespace in which all LinBox code resides
namespace LinBox 
{ 
  // Forward declarations
  template <class Localring> class RandIter_envelope;

  /** Localring Envelope Template.
   * Derived class used to implement the localring archetype to minimize
   * code bloat.  This class implements all purely virtual member functions
   * of the abstract base class.  This class is used to wrap a
   * \Ref{LinBox}
   * localring so that it might be used with the Localring archetype.
   */
  template <class Localring> class Localring_envelope : public Localring_abstract
  {
  public:

    /** Element type.
     * It is derived from the class Element_abstract, and it must contain
     * a wrapped localring element.
     */
    typedef Element_envelope<Localring> element;

    /** Random iterator generator type.
     * It is derived from the class RandIter_abstract, and it must contain
     * a wrapped localring random iterator generator.
     */
    typedef RandIter_envelope<Localring> randIter;

    /** @name Object Management
     */
    //@{
 
    /** Default constructor.
     * In this implementation, this means copying the localring E._localring.
     */
    Localring_envelope(void) {}

    /** Constructor from localring to be wrapped.
     * @param F Localring object to be wrapped.
     */
    Localring_envelope(const Localring& F) : _localring(F) {}
 
    /** Copy constructor.
     * Constructs Localring_envelope object by copying the localring.
     * This is required to allow localring objects to be passed by value
     * into functions.
     * In this implementation, this means copying the localring E._localring.
     * @param  E Localring_envelope object.
     */
    Localring_envelope(const Localring_envelope& E) : _localring(E._localring) {}
 
    /** Virtual copy constructor.
     * Required because constructors cannot be virtual.
     * Passes construction on to derived classes.
     * This function is not part of the common object interface.
     * @return pointer to new object in dynamic memory.
     */
    Localring_abstract* clone() const
    { return new Localring_envelope(*this); }

    /** Assignment operator.
     * Required by abstract base class.
     * @return reference to Localring_abstract object for self
     * @param F constant reference to Localring_abstract object
     */
    Localring_abstract& operator= (const Localring_abstract& F)
    {
      if (this != &F) // guard against self-assignment
  	_localring = static_cast<const Localring_envelope&>(F)._localring;

      return *this;
    }

    /** Initialization of localring base element from an integer.
     * Behaves like C++ allocator construct.
     * This function assumes the output localring base element x has already been
     * constructed, but that it is not already initialized.
     * This is not a specialization of the template function because
     * such a specialization is not allowed inside the class declaration.
     * @return reference to localring base element.
     * @param x localring base element to contain output (reference returned).
     * @param y integer.
     */
    Element_abstract& init(Element_abstract& x, const integer& y) const
    {
      _localring.init(static_cast<Element_envelope<Localring>&>(x)._elem, y);
      return x;
    }
 
   /** Conversion of localring base element to a template class T.
     * This function assumes the output localring base element x has already been
     * constructed, but that it is not already initialized.
     * @return reference to template class T.
     * @param x template class T to contain output (reference returned).
     * @param y constant localring base element.
     */
    integer& convert(integer& x, const Element_abstract& y) const
    {
      _localring.convert(x, static_cast<const Element_envelope<Localring>&>(y)._elem);
      return x;
    }
 
    /** Assignment of one localring base element to another.
     * This function assumes both localring base elements have already been
     * constructed and initialized.
     * @return reference to x
     * @param  x localring base element (reference returned).
     * @param  y localring base element.
     */
    Element_abstract& assign(Element_abstract& x, const Element_abstract& y) const
    {
      _localring.assign(static_cast<Element_envelope<Localring>&>(x)._elem,
  		    static_cast<const Element_envelope<Localring>&>(y)._elem);
      return x;
    }

    /** Cardinality.
     * Return integer representing cardinality of the domain.
     * Returns a non-negative integer for all domains with finite
     * cardinality, and returns -1 to signify a domain of infinite
     * cardinality.
     * @return integer representing cardinality of the domain
     */
    integer& cardinality(integer& c) const { return _localring.cardinality(c); }
 
    /** Characteristic.
     * Return integer representing characteristic of the domain.
     * Returns a positive integer to all domains with finite characteristic,
     * and returns 0 to signify a domain of infinite characteristic.
     * @return integer representing characteristic of the domain.
     */
    integer& characteristic(integer& c) const { return _localring.characteristic(c); }

    //@} Object Management

    /** @name Arithmetic Operations
     * x <- y op z; x <- op y
     * These operations require all elements, including x, to be initialized
     * before the operation is called.  Uninitialized localring base elements will
     * give undefined results.
     */
    //@{

    /** Equality of two elements.
     * This function assumes both localring base elements have already been
     * constructed and initialized.
     * @return boolean true if equal, false if not.
     * @param  x localring base element
     * @param  y localring base element
     */
    bool areEqual(const Element_abstract& x, const Element_abstract& y) const
    {
      return
  	_localring.areEqual(static_cast<const Element_envelope<Localring>&>(x)._elem,
  			static_cast<const Element_envelope<Localring>&>(y)._elem);
    }

    /** Addition.
     * x = y + z
     * This function assumes all the localring base elements have already been
     * constructed and initialized.
     * @return reference to x.
     * @param  x localring base element (reference returned).
     * @param  y localring base element.
     * @param  z localring base element.
     */
    Element_abstract& add(Element_abstract& x,
  			  const Element_abstract& y,
  			  const Element_abstract& z) const
    {
      _localring.add(static_cast<Element_envelope<Localring>&>(x)._elem,
  		 static_cast<const Element_envelope<Localring>&>(y)._elem,
  		 static_cast<const Element_envelope<Localring>&>(z)._elem);
      return x;
    }
 
    /** Subtraction.
     * x = y - z
     * This function assumes all the localring base elements have already been
     * constructed and initialized.
     * @return reference to x.
     * @param  x localring base element (reference returned).
     * @param  y localring base element.
     * @param  z localring base element.
     */
    Element_abstract& sub(Element_abstract& x,
  			  const Element_abstract& y,
  			  const Element_abstract& z) const
    {
      _localring.sub(static_cast<Element_envelope<Localring>&>(x)._elem,
  		 static_cast<const Element_envelope<Localring>&>(y)._elem,
  		 static_cast<const Element_envelope<Localring>&>(z)._elem);
      return x;
    }
 
    /** Multiplication.
     * x = y * z
     * This function assumes all the localring base elements have already been
     * constructed and initialized.
     * @return reference to x.
     * @param  x localring base element (reference returned).
     * @param  y localring base element.
     * @param  z localring base element.
     */
    Element_abstract& mul(Element_abstract& x,
  			  const Element_abstract& y,
  			  const Element_abstract& z) const
    {
      _localring.mul(static_cast<Element_envelope<Localring>&>(x)._elem,
  		 static_cast<const Element_envelope<Localring>&>(y)._elem,
  		 static_cast<const Element_envelope<Localring>&>(z)._elem);
      return x;
    }
 
     /** Divisibility
       * @return boolean true iff there exists ring element z such that xz = y.
       * @param  x localring element.
       * @param  y localring element.
       * @param  z localring element.
       */
    bool divides(const Element_abstract& x, const Element_abstract& y) const
    { return _localring.divides(
        static_cast<const Element_envelope<Localring>&>(x)._elem, 
        static_cast<const Element_envelope<Localring>&>(y)._elem); 
    }

    /** Division.
     * x = y / z
     * This function assumes all the localring base elements have already been
     * constructed and initialized.
     * @return reference to x.
     * @param  x localring base element (reference returned).
     * @param  y localring base element.
     * @param  z localring base element. z divides y.
     */
    Element_abstract& div(Element_abstract& x,
                          const Element_abstract& y,
  			  const Element_abstract& z) const
    {
      _localring.div(static_cast<Element_envelope<Localring>&>(x)._elem,
  		 static_cast<const Element_envelope<Localring>&>(y)._elem,
  		 static_cast<const Element_envelope<Localring>&>(z)._elem);
      return x;
    }
 
    /** Additive Inverse (Negation).
     * x = - y
     * This function assumes both localring base elements have already been
     * constructed and initialized.
     * @return reference to x.
     * @param  x localring base element (reference returned).
     * @param  y localring base element.
     */
    Element_abstract& neg(Element_abstract& x, const Element_abstract& y) const
    {
      _localring.neg(static_cast<Element_envelope<Localring>&>(x)._elem,
  		 static_cast<const Element_envelope<Localring>&>(y)._elem);
      return x;
    }
 
    /** Multiplicative Inverse.
     * x = 1 / y
     * This function assumes both localring base elements have already been
     * constructed and initialized.
     * @return reference to x.
     * @param  x localring base element (reference returned).
     * @param  y localring base element.
     */
    Element_abstract& inv(Element_abstract& x, const Element_abstract& y) const
    {
      _localring.inv(static_cast<Element_envelope<Localring>&>(x)._elem,
  		 static_cast<const Element_envelope<Localring>&>(y)._elem);
      return x;
    }

    //@} Arithmetic Operations
 
    /** @name Inplace Arithmetic Operations
     * x <- x op y; x <- op x
     */
    //@{

    /** Zero equality.
     * Test if localring base element is equal to zero.
     * This function assumes the localring base element has already been
     * constructed and initialized.
     * @return boolean true if equals zero, false if not.
     * @param  x localring base element.
     */
    bool isZero(const Element_abstract& x) const
    {
      return _localring.isZero(static_cast<const Element_envelope<Localring>&>(x)._elem);
    }
 
    /** One equality.
     * Test if localring base element is equal to one.
     * This function assumes the localring base element has already been
     * constructed and initialized.
     * @return boolean true if equals one, false if not.
     * @param  x localring base element.
     */
    bool isOne(const Element_abstract& x) const
    {
      return _localring.isOne(static_cast<const Element_envelope<Localring>&>(x)._elem);
    }

    /** Invertibility
     * Test of localring element is a unit, i.e. has an inverse in the ring.
     * @return boolean true if x is a unit, false if not.
     * @param  x localring element.
     */
    bool isUnit(const element& x) const
    { return _localring.isUnit(static_cast<const Element_envelope<Localring>&>(x)._elem);
    }

    /** Inplace Addition.
     * x += y
     * This function assumes both localring base elements have already been
     * constructed and initialized.
     * @return reference to x.
     * @param  x localring base element (reference returned).
     * @param  y localring base element.
     */
    Element_abstract& addin(Element_abstract& x, const Element_abstract& y) const
    {
      _localring.addin(static_cast<Element_envelope<Localring>&>(x)._elem,
  		   static_cast<const Element_envelope<Localring>&>(y)._elem);
      return x;
    }
 
    /** Inplace Subtraction.
     * x -= y
     * This function assumes both localring base elements have already been
     * constructed and initialized.
     * @return reference to x.
     * @param  x localring base element (reference returned).
     * @param  y localring base element.
     */
    Element_abstract& subin(Element_abstract& x, const Element_abstract& y) const
    {
      _localring.subin(static_cast<Element_envelope<Localring>&>(x)._elem,
  		   static_cast<const Element_envelope<Localring>&>(y)._elem);
      return x;
    }
 
    /** Inplace Multiplication.
     * x *= y
     * This function assumes both localring base elements have already been
     * constructed and initialized.
     * @return reference to x.
     * @param  x localring base element (reference returned).
     * @param  y localring base element.
     */
    Element_abstract& mulin(Element_abstract& x, const Element_abstract& y) const
    {
      _localring.mulin(static_cast<Element_envelope<Localring>&>(x)._elem,
  		   static_cast<const Element_envelope<Localring>&>(y)._elem);
      return x;
    }
 
    /** Inplace Division.
     * x /= y
     * This function assumes both localring base elements have already been
     * constructed and initialized.
     * @return reference to x.
     * @param  x localring base element (reference returned).
     * @param  y localring base element.
     */
    Element_abstract& divin(Element_abstract& x, const Element_abstract& y) const
    {
      _localring.divin(static_cast<Element_envelope<Localring>&>(x)._elem,
  		   static_cast<const Element_envelope<Localring>&>(y)._elem);
      return x;
    }
 
    /** Inplace Additive Inverse (Inplace Negation).
     * x = - x
     * This function assumes the localring base element has already been
     * constructed and initialized.
     * @return reference to x.
     * @param  x localring base element (reference returned).
     */
    Element_abstract& negin(Element_abstract& x) const
    {
      _localring.negin(static_cast<Element_envelope<Localring>&>(x)._elem);
      return x;
    }
 
    /** Inplace Multiplicative Inverse.
     * x = 1 / x
     * This function assumes the localring base elementhas already been
     * constructed and initialized.
     * @return reference to x.
     * @param  x localring base element (reference returned).
     */
    Element_abstract& invin(Element_abstract& x) const
    {
      _localring.invin(static_cast<Element_envelope<Localring>&>(x)._elem);
      return x;
    }

    //@} Inplace Arithmetic Operations

    /** @name Input/Output Operations */
    //@{

    /** Print localring.
     * @return output stream to which localring is written.
     * @param  os  output stream to which localring is written.
     */
    ostream& write(ostream& os) const { return _localring.write(os); }
 
    /** Read localring.
     * @return input stream from which localring is read.
     * @param  is  input stream from which localring is read.
     */
    istream& read(istream& is) { return _localring.read(is); }

    /** Print localring base element.
     * This function assumes the localring base element has already been
     * constructed and initialized.
     * @return output stream to which localring base element is written.
     * @param  os  output stream to which localring base element is written.
     * @param  x   localring base element.
     */
    ostream& write(ostream& os, const Element_abstract& x) const
    { return _localring.write(os, static_cast<const Element_envelope<Localring>&>(x)._elem); }
 
    /** Read localring base element.
     * This function assumes the localring base element has already been
     * constructed and initialized.
     * @return input stream from which localring base element is read.
     * @param  is  input stream from which localring base element is read.
     * @param  x   localring base element.
     */
    istream& read(istream& is, Element_abstract& x) const
    { return _localring.read(is, static_cast<Element_envelope<Localring>&>(x)._elem); }

    //@}

  private:

    friend RandIter_envelope<Localring>;

    /// Wrapped localring.
    Localring _localring;

  }; // class Localring_envelope

} // namespace LinBox

#include "LinBox/randiter_envelope.h"

#endif // _LOCALRING_ENVELOPE_
