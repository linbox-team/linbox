/* File: src/wrappers/by_scope/field/unparam_field.h
 * Author: William J Turner for the LinBox group
 */
 
#ifndef __UNPARAM_FIELD_WRAPPER_
#define __UNPARAM_FIELD_WRAPPER_

#include <string>
#include <algorithm>
#include "LinBox/integer.h"
#include "LinBox/unparam_randiter.h"

// Namespace in which all LinBox library code resides
namespace LinBox
{
  /** Unparameterized field template.
   * Implements LinBox field common object interface for unparameterized 
   * fields.
   * Used to generate efficient field classes for unparameterized fields.
   * Constructs LinBox unparameterized fields from field types K.
   * In particular, constructs LinBox fields for
   * unparameterized fields from field types that
   * adhere to the operations for double, for
   * example unparam_field< float >.
   * Can be used as a pattern to write a particular
   * field interface, such as, unparam_field< SaclibQ > as
   * a template specialization.
   * @param  K unparameterized field class
   */
  template <class K> class unparam_field
  {
  public:
    
    /** @name Common Object Interface for a LinBox Field.
     * These methods are required of all LinBox fields.
     */
    //@{
    
    /** Field element type.
     * The field element must contain a default constructor, 
     * a copy constructor, a destructor, and an assignment operator.
     */
    typedef K element;    

    /// Random field element generator type.
    typedef unparam_randIter<K> randIter;
    
    /** @name Object Management.
     * x <- convert(y)
       */
    //@{
    
    /** Copy constructor.
     * Constructs unparam_field object by copying the field.
     * This is required to allow field objects to be passed by value
     * into functions.
     * @param  F unparam_field object.
       */
    unparam_field(const unparam_field& F) {}
    
    /** Destructor.
     * This destructs the field object, but it does not destroy the field 
     * element objects.  The destructor for each field element must also 
     * be called.
     * _elem_ptr points.
     */
    ~unparam_field() {}
    
    /** Assignment operator.
      * Assigns unparam_field object F to field.
      * @param  F unparam_field object.
      */
    unparam_field& operator=(const unparam_field& F) { return *this; }
    
    /** Initialization of field element from an integer.
     * Behaves like C++ allocator construct.
     * This function assumes the output field element x has already been 
     * constructed, but that it is not already initialized.
     * For now, this is done by converting the integer type to a C++ 
     * long and then to the element type through the use of static casts.
     * This, of course, assumes such static casts are possible.
     * This function should be changed in the future to avoid using long.
     * @return reference to field element.
     * @param x field element to contain output (reference returned).
     * @param y integer.
     */
    element& init(element& x, const integer& y=0) const 
    { return x = static_cast<const element&>(static_cast<const long&>(y)); }
    
    /** Conversion of field element to an integer.
     * This function assumes the output field element x has already been 
     * constructed, but that it is not already initialized.
     * For now, this is done by converting the element type to a C++ 
     * long and then to the integer type through the use of static casts.
     * This, of course, assumes such static casts are possible.
     * This function should be changed in the future to avoid using long.
     * @return reference to integer.
     * @param x reference to integer to contain output (reference returned).
     * @param y constant reference to field element.
     */
    integer& convert(integer& x, const element& y) const 
    { 
      element temp(y);
      return x = static_cast<long>(temp); 
    }
    
    /** Assignment of one field element to another.
     * This function assumes both field elements have already been 
     * constructed and initialized.
     * @return reference to x
     * @param  x field element (reference returned).
     * @param  y field element.
     */
    element& assign(element& x, const element& y) const { return x = y; }
    
    /** Cardinality.
     * Return integer representing cardinality of the field.
     * Returns a non-negative integer for all fields with finite
     * cardinality, and returns -1 to signify a field of infinite 
     * cardinality.
     * The default behavior of this function is to return -1 to signify
     * a field of infinite cardinality.  This should be changed via
     * partial template specialization for fields of other cardinalities.
     * @return integer representing cardinality of the field
     */
    const integer& cardinality(void) const { return *(new integer(-1)); }
    
    /** Characteristic.
     * Return integer representing characteristic of the field.
     * Returns a positive integer to all fields with finite characteristic,
     * and returns 0 to signify a field of infinite characteristic.
     * The default behavior of this function is to return 0 to signify
     * a field of infinite characteristic.  This should be changed via
     * partial template specialization for fields of other characteristics.
     * @return integer representing characteristic of the field.
     */
    const integer& characteristic(void) const { return *(new integer(0)); }
    
    //@} Object Management
    
    /** @name Arithmetic Operations 
     * x <- y op z; x <- op y
     * These operations require all elements, including x, to be initialized
     * before the operation is called.  Uninitialized field elements will
     * give undefined results.
     */
    //@{
    
    /** Equality of two elements.
     * This function assumes both field elements have already been 
     * constructed and initialized.
     * @return boolean true if equal, false if not.
     * @param  x field element
     * @param  y field element
     */
    bool areEqual(const element& x, const element& y) const { return x == y; }
    
    /** Addition.
     * x = y + z
     * This function assumes all the field elements have already been 
     * constructed and initialized.
     * @return reference to x.
     * @param  x field element (reference returned).
     * @param  y field element.
     * @param  z field element.
     */
    element& add(element& x, const element& y, const element& z) const
      { return x = y + z; }
    
    /** Subtraction.
     * x = y - z
     * This function assumes all the field elements have already been 
     * constructed and initialized.
     * @return reference to x.
     * @param  x field element (reference returned).
     * @param  y field element.
     * @param  z field element.
     */
    element& sub(element& x, const element& y, const element& z) const
      { return x = y - z; }
    
    /** Multiplication.
     * x = y * z
     * This function assumes all the field elements have already been 
     * constructed and initialized.
     * @return reference to x.
     * @param  x field element (reference returned).
     * @param  y field element.
     * @param  z field element.
     */
    element& mul(element& x, const element& y, const element& z) const
      { return x = y * z; }
    
    /** Division.
     * x = y / z
     * This function assumes all the field elements have already been 
     * constructed and initialized.
     * @return reference to x.
     * @param  x field element (reference returned).
     * @param  y field element.
     * @param  z field element.
     */
    element& div(element& x, const element& y, const element& z) const
      { return x = y / z; }
    
    /** Additive Inverse (Negation).
     * x = - y
     * This function assumes both field elements have already been 
     * constructed and initialized.
     * @return reference to x.
     * @param  x field element (reference returned).
     * @param  y field element.
     */
    element& neg(element& x, const element& y) const { return x = - y; }
    
    /** Multiplicative Inverse.
     * x = 1 / y
     * This function assumes both field elements have already been 
     * constructed and initialized.
     * @return reference to x.
     * @param  x field element (reference returned).
     * @param  y field element.
     */
    element& inv(element& x, const element& y) const 
      { return x = element(1) / y; }
    
    /** Natural AXPY.
     * r  = a * x + y
     * This function assumes all field elements have already been 
     * constructed and initialized.
     * @return reference to r.
     * @param  r field element (reference returned).
     * @param  a field element.
     * @param  x field element.
     * @param  y field element.
     */
    element& axpy(element& r, 
		  const element& a, 
		  const element& x, 
		  const element& y) const
    { return r = a * x + y; }
 
    //@} Arithmetic Operations
    
    /** @name Inplace Arithmetic Operations 
     * x <- x op y; x <- op x
     * These operations require all elements, including x, to be initialized
     * before the operation is called.  Uninitialized field elements will
     * give undefined results.
     */
    //@{
    
    /** Zero equality.
     * Test if field element is equal to zero.
     * This function assumes the field element has already been 
     * constructed and initialized.
     * @return boolean true if equals zero, false if not.
     * @param  x field element.
     */
    bool isZero(const element& x) const { return x == element(0); }
    
    /** One equality.
     * Test if field element is equal to one.
     * This function assumes the field element has already been 
     * constructed and initialized.
     * @return boolean true if equals one, false if not.
     * @param  x field element.
     */
    bool isOne(const element& x) const { return x == element(1); }
    
    /** Inplace Addition.
     * x += y
     * This function assumes both field elements have already been 
     * constructed and initialized.
     * @return reference to x.
     * @param  x field element (reference returned).
     * @param  y field element.
     */
    element& addin(element& x, const element& y) const { return x += y; }
    
    /** Inplace Subtraction.
     * x -= y
     * This function assumes both field elements have already been 
     * constructed and initialized.
     * @return reference to x.
     * @param  x field element (reference returned).
     * @param  y field element.
     */
    element& subin(element& x, const element& y) const { return x -= y; }
    
    /** Inplace Multiplication.
     * x *= y
     * This function assumes both field elements have already been 
     * constructed and initialized.
     * @return reference to x.
     * @param  x field element (reference returned).
     * @param  y field element.
     */
    element& mulin(element& x, const element& y) const { return x *= y; }
    
    /** Inplace Division.
     * x /= y
     * This function assumes both field elements have already been 
     * constructed and initialized.
     * @return reference to x.
     * @param  x field element (reference returned).
     * @param  y field element.
     */
    element& divin(element& x, const element& y) const { return x /= y; }
    
    /** Inplace Additive Inverse (Inplace Negation).
     * x = - x
     * This function assumes the field element has already been 
     * constructed and initialized.
     * @return reference to x.
     * @param  x field element (reference returned).
     */
    element& negin(element& x) const { return x = - x; }
    
    /** Inplace Multiplicative Inverse.
     * x = 1 / x
     * This function assumes the field elementhas already been 
     * constructed and initialized.
     * @return reference to x.
     * @param  x field element (reference returned).
     */
    element& invin(element& x) const { return x = element(1) / x; }
    
    /** Inplace AXPY.
     * r  += a * x
     * This function assumes all field elements have already been 
     * constructed and initialized.
     * @return reference to r.
     * @param  r field element (reference returned).
     * @param  a field element.
     * @param  x field element.
     */
    element& axpyin(element& r, const element& a, const element& x) const
    { return r += a * x; }
 
    //@} Inplace Arithmetic Operations
    
    /** @name Input/Output Operations */
    //@{
    
    /** Print field.
     * @return output stream to which field is written.
     * @param  os  output stream to which field is written.
     */
    ostream& write(ostream& os) const { return os << "unparamterized field"; }
    
    /** Read field.
     * @return input stream from which field is read.
     * @param  is  input stream from which field is read.
     */
    istream& read(istream& is) const { return is; }
    
    /** Print field element.
     * @return output stream to which field element is written.
     * @param  os  output stream to which field element is written.
     * @param  x   field element.
     */
    ostream& write(ostream& os, const element& x) const { return os << x; }
    
    /** Read field element.
     * @return input stream from which field element is read.
     * @param  is  input stream from which field element is read.
     * @param  x   field element.
     */
    istream& read(istream& is, element& x) const { return is >> x; }
    
    //@}
    
    //@} Common Object Interface
    
    /** @name Implementation-Specific Methods.
     * These methods are not required of all LinBox fields
     * and are included only for the implementation of this field
     * template.
     */
    //@{
    
    /// Default constructor
    unparam_field(void) {}
    
    /** Constructor from field object.
     * @param  A unparameterized field object
     */
    unparam_field(const K& A) {} 
    
    /** Constant access operator.
     * @return constant reference to field object
     */
    const K& operator () (void) const { return element(); }
    
    /** Access operator.
     * @return reference to field object
     */
    K& operator () (void) { return element(); }
    
    //@} Implementation-Specific Methods
    
  }; // template <class K> class unparam_field
  
} // namespace LinBox

#include "LinBox/unparam_randiter.h"

#endif // _UNPARAM_FIELD_WRAPPER_
