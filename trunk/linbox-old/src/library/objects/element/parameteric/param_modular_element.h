/* File: src/library/objects/element/param_modular_element.h
 * Author: William J Turner for the LinBox group
 */

#ifndef _PARAM_MODULAR_ELEMENT_
#define _PARAM_MODULAR_ELEMENT_

// Namespace in which all LinBox library code resides
namespace LinBox
{

  // forward declaration
  param_modular;

  /** Element of parameterized modular field \Ref{param_modular}.
    * The field element contains the residue but little else.
    * It has no knowledge of the modulus being employed and so
    * the field object must provide all operations and funtions
    * in which it is used.
    */
  class param_modular_element
  {
  public:
    /** @name Common Object Interface for LinBox Field Elements.
      * These methods are required of all LinBox field elements.
      */
    //@{

    /// Default constructor
    param_modular_element(void) {}

    /** Constructor from field and integer.
      * Assigns the representative between 0 and (modulus - 1) 
      * of the residue class to which start belongs to the residue.
      * @param  F parameteric modular field object
      * @param  start integer (default = 0)
      */
    param_modular_element(const param_modular& F, const integer& start = 0) 
    {
      residue = start % F.modulus;
      if (residue < 0)
      residue = residue + F.modulus;
    }

    /** Assignment operator.
      * @param  x parameteric modular field element
      * @return reference to self
      */
    param_modular_element& operator=(const param_modular_element& x)
    { 
      residue = x.residue;
      return *this;
    }

    /** Destructor.
      */
    ~param_modular_element() {}

    //@} Common Object Interface

  private:

    friend param_modular;

    /** Residue class representative.
      * No necessarily between 0 and (modulus - 1) because field 
      * element
      * has no knoledge of modulus.
      */
    integer residue; 

  }; // class param_modular_element

} // namespace LinBox

#endif // _PARAM_MODULAR_ELEMENT_
