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
  class abstract_fuzzy_randIter;
  
  /** Abstract fuzzy LinBox field element.
   * Derived class used to implement the field element archetype to minimize
   * code bloat.  This class implements all purely virtual member functions
   * of the abstract base class.  This class implements the elements of
   * a field of doubles with two doubles considered to be equal if
   * their difference is less than a static fuzz value.
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
     * In this implementation, this means copying the element E._value.
     * @param  E Field_envelope object.
     */
    abstract_fuzzy_element(const Element_abstract& E)
      : _value(static_cast<const abstract_fuzzy_element&>(E)._value) {}
  
    /** Constructor from an integer.
     * Sets residue to value supplied.
     * @param value constant reference to integer
     */
    abstract_fuzzy_element(const integer& value) 
    { _value = static_cast<const double&>(value); }

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
        _value = static_cast<const abstract_fuzzy_element&>(E)._value;

      return *this;
    }

    /** Destructor.
     */
    ~abstract_fuzzy_element(void) {}

  private:

    // Friend declarations
    friend abstract_fuzzy;
    friend abstract_fuzzy_randIter;

    /// Private integer for residue class
    double _value;

  }; // class abstract_fuzzy

} // namespace LinBox

#endif // _ABSTRACT_FUZZY_ELEMENT_
