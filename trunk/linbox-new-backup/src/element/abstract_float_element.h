/* File: src/library/objects/element/abstract_float_element.h
 * Author: William Turner for the LinBox group
 */

#ifndef _ABSTRACT_FLOAT_ELEMENT_
#define _ABSTRACT_FLOAT_ELEMENT_

#include <iostream>
#include "LinBox/integer.h"
#include "LinBox/element_abstract.h"

// Namespace in which all LinBox code resides
namespace LinBox 
{ 
  // Forward declarations
  class abstract_float;
  class abstract_float_randIter;
  
  /** Abstract float LinBox field element.
   * Derived class used to implement the field element archetype to minimize
   * code bloat.  This class implements all purely virtual member functions
   * of the abstract base class.  This class implements the element class 
   * of a field of floats 
   * through a private float member.
   *
   * @see abstract_float
   */
  class abstract_float_element : public Element_abstract
  {
  public:
    
    /** Default Constructor.
     */
    abstract_float_element(void) {}

    /** Copy constructor.
     * Constructs abstract_float_element object by copying the element
     * it wraps.
     * This is required to allow element objects to be passed by value
     * into functions.
     * In this implementation, this means copying the element E._elem.
     * @param  E Field_envelope object.
     */
    abstract_float_element(const Element_abstract& E)
      : _elem(static_cast<const abstract_float_element&>(E)._elem) 
    {}
  
    /** Constructor from a float.
     * Sets private float member to value supplied.
     * @param elem float
     */
    abstract_float_element(float elem) : _elem(elem) {}

    /** Virtual copy constructor.
     * Required because constructors cannot be virtual.
     * Passes construction on to derived classes.
     * @return pointer to new element object in dynamic memory.
     */
    Element_abstract* clone(void) const 
    { return new abstract_float_element(*this); }

    /** Assignment operator.
     * @return reference to self
     * @param  x parameterized field base element
     */
    Element_abstract& operator=(const Element_abstract& E)
    {
      if (this != &E) // guard against self-assignment
        _elem = static_cast<const abstract_float_element&>(E)._elem;

      return *this;
    }

    /** Destructor.
     */
    ~abstract_float_element(void) {}

  private:

    // Friend declarations
    friend abstract_float;
    friend abstract_float_randIter;

    float _elem;

  }; // class abstract_float

} // namespace LinBox

#endif // _ABSTRACT_FLOAT_ELEMENT_
