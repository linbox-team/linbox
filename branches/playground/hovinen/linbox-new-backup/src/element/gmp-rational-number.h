/* File: src/library/objects/element/gmp-rational-number.h
 * Author: Bradford Hovinen for the LinBox group
 */

#ifndef _GMP_RATIONAL_NUMBER_
#define _GMP_RATIONAL_NUMBER_

#include "LinBox/integer.h"

#include <gmp.h>

// Namespace in which all LinBox library code resides
namespace LinBox
{

  // Forward declarations
  class GMP_Rational_Field;
  class GMP_Rational_Random;

  /** Element archetype.
   * Archetype for the element common object interface for \Ref{LinBox}.
   *
   * This class must contain a default constructor, a copy constructor, 
   * an assignment operator, and a destructor.  This is to allow field 
   * elements to be primitive C++ types such as double and int.  The copy
   * constructor is also used to allow elements to be passed by value to 
   * a function.
   */
  class GMP_Rational_Number
  {
  public:

    /** @name Common Object Interface for LinBox Field Elements.
     * These methods are required of all \Ref{LinBox} 
     * {@link Fields field} elements.
     */
   //@{

    /** Default constructor.
     * This constructor is required to allow 
     * {@link Fields field} elements to be primitive C++ types.
     * Because constructor does not know what {@link Fields field} 
     * the element belongs to, it cannot actually construct the element.
     * In this implementation, the constructor it sets _elem_ptr
     * to the null pointer.  Initialization of the element is done through
     * the field function init where the field is known.
     */
    GMP_Rational_Number(void) { mpq_init (rep); }

    /** Copy constructor.
     * This constructor is required to allow 
     * {@link Fields field} elements to be primitive C++ types, 
     * and to allow field elements to be passed by value into 
     * functions.
     * Constructs {@link Fields field} element by copying the 
     * {@link Fields field} element.
     * In this implementation, this means copying the element to
     * which a._elem_ptr points.
     * @param  a field element.
     */
    GMP_Rational_Number(const GMP_Rational_Number& a) 
    { mpq_init (rep); mpq_set (rep, a.rep); }

    /** Destructor.
     * In this implementation, this destroys element by deleting field 
     * element to which _elem_ptr points.
     */
    ~GMP_Rational_Number() { mpq_clear (rep); }

    /** Assignment operator.
     * Assigns element a to element.  
     * In this implementation, this is done 
     * by copying field element to which _elem_ptr points.
     * @param  a field element.
     */
    GMP_Rational_Number& operator=(const GMP_Rational_Number& a)
    {
	 if (this != &a) { // guard against self-assignment
	      mpq_set (rep, a.rep);
	 }
	 return *this;
    }

    //@} Common Object Interface

    /** @name Implementation-Specific Methods.
     * These methods are not required of all LinBox field elements
     * and are included only for this implementation of the archetype.
     */
    //@{

    /** Constructor.
     * Constructs field element from an mpq_t
     * Not part of the interface.
     * Creates new copy of element object in dynamic memory.
     * @param  elem_ptr  pointer to \Ref{Element_abstract}
     */
    GMP_Rational_Number (mpq_t _rep) {
	 mpq_init (rep);
	 mpq_set (rep, _rep);
    }

    /** Constructor
     * Initialize from numerator and denominator
     */
    GMP_Rational_Number (const integer &num, const integer &den) 
	 {
	      mpq_init (rep);
	      mpz_set_si (mpq_numref (rep), num);
	      mpz_set_si (mpq_denref (rep), den);
	 }

    //@}p
    
  private:

    friend class GMP_Rational_Field;
    friend class GMP_Rational_Random;

    /** @name Implementation-Specific Data.
     * This data is not required of all LinBox field elements
     * and is included only for this implementation of the archetype.
     */
    //@{
    
    /** Pointer to parameterized field element.
     * Not part of the common object interface for \Ref{LinBox} field elements.
     * Included to avoid code bloat.
     */
    mutable mpq_t rep;
    
    //@} Non-Interface

  }; // class element

} // namespace LinBox

#endif // _GMP_RATIONAL_NUMBER_
