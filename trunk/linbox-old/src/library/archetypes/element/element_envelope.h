/* File: src/library/archetypes/element/element_envelope.h
 * Author: William Turner for the LinBox group
 */

#ifndef _ELEMENT_ENVELOPE_
#define _ELEMENT_ENVELOPE_

#include <iostream>
#include "LinBox/element_abstract.h"

// Namespace in which all LinBox code resides
namespace LinBox 
{ 
  // Forward declarations
  template <class Field> class Field_envelope;
  template <class Field> class RandIter_envelope;

  /** Element envelope template.
   * Encapsulated class of Field_envelope class.
   * This element has no knowledge of the field to which it belongs, 
   * so all operations and functions requiring knolwedge of the field,
   * such as addition and other arithmetic operations, must be supplied
   * by the field and not the element.
   */
  template <class Field>
  class Element_envelope : public Element_abstract
  {
  public:

    /** Default Constructor.
     */
    Element_envelope() {}

    /** Constructor from the Field element to be wrapped.
     * @param elem Field element object to be wrapped.
     */
    Element_envelope(const typename Field::element& elem) : _elem(elem) {}

    /** Copy constructor.
     * Constructs Element_envelope object by copying the element
     * it wraps.
     * This is required to allow element objects to be passed by value
     * into functions.
     * In this implementation, this means copying the element E._elem.
     * @param  E Field_envelope object.
     */
    Element_envelope(const Element_abstract& E)
      : _elem(static_cast<const Element_envelope&>(E)._elem) {}
  
    /** Virtual copy constructor.
     * Required because constructors cannot be virtual.
     * Passes construction on to derived classes.
     * @return pointer to new element object in dynamic memory.
     */
    Element_abstract* clone(void) const { return new Element_envelope(*this); }

    /** Assignment operator.
     * @return reference to self
     * @param  x parameterized field base element
     */
    Element_abstract& operator=(const Element_abstract& E)
    {
      if (this != &E) // guard against self-assignment
	_elem = static_cast<const Element_envelope&>(E)._elem;
      return *this;
    }

    /** Destructor.
     */
    ~Element_envelope() {}

  private:

    // Friend declarations
    friend Field_envelope<Field>;
    friend RandIter_envelope<Field>;

    typename Field::element _elem;

  }; // class Element_envelope

} // namespace LinBox

#endif // _ELEMENT_ENVELOPE_
