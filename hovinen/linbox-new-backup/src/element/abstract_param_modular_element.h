/* File: src/library/objects/element/abstract_param_modular_element.h
 * Author: William Turner for the LinBox group
 */

#ifndef _ABSTRACT_PARAM_MODULAR_ELEMENT_
#define _ABSTRACT_PARAM_MODULAR_ELEMENT_

#include <iostream>
#include "LinBox/integer.h"
#include "LinBox/element_abstract.h"

// Namespace in which all LinBox code resides
namespace LinBox 
{ 
  // Forward declarations
  class abstract_param_modular;
  class abstract_param_modular_randIter;
  
  /** Abstract param_modular LinBox field element.
   * Derived class used to implement the field element archetype to minimize
   * code bloat.  This class implements all purely virtual member functions
   * of the abstract base class.  This class implements the element class
   * of a field of unparamterized integers modulo a prime integer.
   *
   * @see abstract_param_modular
   */
  class abstract_param_modular_element : public Element_abstract
  {
  public:

    /** Default Constructor.
     */
    abstract_param_modular_element(void) {}

    /** Copy constructor.
     * Constructs abstract_param_modular_element object by copying the element
     * it wraps.
     * This is required to allow element objects to be passed by value
     * into functions.
     * In this implementation, this means copying the element E._residue.
     * @param  E Field_envelope object.
     */
    abstract_param_modular_element(const Element_abstract& E)
      : _residue(static_cast<const abstract_param_modular_element&>(E)._residue) {}
  
    /** Constructor from an integer.
     * Sets residue to value supplied.
     * @param value constant reference to integer
     */
    abstract_param_modular_element(const integer& value) { _residue = value; }

    /** Virtual copy constructor.
     * Required because constructors cannot be virtual.
     * Passes construction on to derived classes.
     * @return pointer to new element object in dynamic memory.
     */
    Element_abstract* clone(void) const 
    { return new abstract_param_modular_element(*this); }

    /** Assignment operator.
     * @return reference to self
     * @param  x parameterized field base element
     */
    Element_abstract& operator=(const Element_abstract& E)
    {
      if (this != &E) // guard against self-assignment
        _residue = static_cast<const abstract_param_modular_element&>(E)._residue;

      return *this;
    }

    /** Destructor.
     */
    ~abstract_param_modular_element(void) {}

  private:

    // Friend declarations
    friend abstract_param_modular;
    friend abstract_param_modular_randIter;

    /// Private integer for residue class
    integer _residue;

  }; // class abstract_param_modular

} // namespace LinBox

#endif // _ABSTRACT_PARAM_MODULAR_ELEMENT_
