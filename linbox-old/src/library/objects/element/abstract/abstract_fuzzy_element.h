/* File: src/library/objects/element/abstract_fuzzy_element.h
 * Author: William Turner for the LinBox group
 */

#ifndef _ABSTRACT_FUZZY_ELEMENT_
#define _ABSTRACT_FUZZY_ELEMENT_

#include <iostream>
#include "LinBox/integer.h"
#include "LinBox/element_abstract.h"

// Namespace in which all LinBox code resides
namespace LinBox 
{ 
  // Forward declarations
  class abstract_fuzzy;
  
  /** Abstract "fuzzy" double LinBox field element.
   * Derived class used to implement the field element archetype to minimize
   * code bloat.  This class implements all purely virtual member functions
   * of the abstract base class.  This class implements the element class
   * of a field of unparamterized integers modulo a prime integer.
   *
   * @see abstract_fuzzy
   */
  class abstract_fuzzy_element : public Element_abstract
  {
  public:

    /** Default Constructor.
     */
    abstract_fuzzy_element(void) {}

    /** Copy constructor.
     * Constructs abstract_fuzzy_element object by copying the element
     * it wraps.
     * This is required to allow element objects to be passed by value
     * into functions.
     * In this implementation, this means copying the element E._residue.
     * @param  E Field_envelope object.
     */
    abstract_fuzzy_element(const Element_abstract& E)
      : _residue(static_cast<const abstract_fuzzy_element&>(E)._residue) 
    {}
  
    /** Constructor from an integer.
     * Sets modulus to value supplied.
     * @param value constant reference to integer
     */
    abstract_fuzzy_element(const integer& value) 
    {
      _residue = value % _modulus;
      if (_residue < 0)
        _residue += _modulus;
    }

    /** Virtual copy constructor.
     * Required because constructors cannot be virtual.
     * Passes construction on to derived classes.
     * @return pointer to new element object in dynamic memory.
     */
    Element_abstract* clone(void) const 
    { return new abstract_fuzzy_element(*this); }

    /** Assignment operator.
     * @return reference to self
     * @param  x parameterized field base element
     */
    Element_abstract& operator=(const Element_abstract& E)
    {
      if (this != &E) // guard against self-assignment
        _residue = static_cast<const abstract_fuzzy_element&>(E)._residue;

      return *this;
    }

    /** Destructor.
     */
    ~abstract_fuzzy_element(void) {}

  private:

    // Friend declarations
    friend abstract_fuzzy;

    /// Private integer for residue class
    integer _residue;

    /// Private static integer for modulus
    static integer _modulus;

  }; // class abstract_fuzzy
  
  integer abstract_fuzzy_element::_modulus;  // declare static member


} // namespace LinBox

#endif // _ABSTRACT_FUZZY_ELEMENT_
