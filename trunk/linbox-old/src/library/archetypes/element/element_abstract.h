/* File: src/library/archetypes/element/element_abstract.h
 * Author: William Turner for the LinBox group
 */

#ifndef _ELEMENT_ABSTRACT_
#define _ELEMENT_ABSTRACT_

// Namespace in which all LinBox code resides
namespace LinBox 
{ 

  /** Abstract element base class.
   * Element class of \Ref{Field_abstract}.
   * This element has no knowledge of the field to which it belongs, 
   * so all operations and functions requiring knolwedge of the field,
   * such as addition and other arithmetic operations, must be supplied
   * by the field and not the element.
   */
  class Element_abstract 
  {
  public:
    
    /** Virtual copy constructor.
     * Required because constructors cannot be virtual.
     * Passes construction on to derived classes.
     * Purely virtual.
     * @return pointer to new Element_abstract object in dynamic memory.
     */
    virtual Element_abstract* clone(void) const = 0;

    /** Assignment operator.
     * Purely virtual.
     * @param  x constant reference to Element_abstract object
     * @return reference to self
     */
    virtual Element_abstract& operator=(const Element_abstract& x) // = 0;
    // This should be purely virtual, but for some reason the compiler
    // has linking errors for abstract_double and abstract_float.
    { cout << "called Element_abstract::operator=" << endl; return *this; }

    /** Destructor.
     */
    virtual ~Element_abstract(void) {}

  protected:

    /** Default Constructor.
     * Required by derived classes, but protected because this class should
     * never be constructed by itself.
     */
    Element_abstract(void) {}

  }; // class Element_abstract

} // namespace LinBox

#endif // _ELEMENT_ABSTRACT_

