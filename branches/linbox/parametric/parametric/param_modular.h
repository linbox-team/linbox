/* File: src/library/objects/field/param_modular.h
 * Author: William J Turner for the LinBox group
 */

#ifndef _PARAM_MODULAR_
#define _PARAM_MODULAR_

// Namespace in which all LinBox library code resides
namespace LinBox
{

/** Parameterized field modulo small (integer) prime number.
  * Implementation of parametric field modulo a prime integer.
  * Field object contains modulus and all arithmetic and other functions
  * in which the modulus is needed.
  * The encapsulated element class contains residue and constructors,
  * but has no knowledge of the modulus in use.
  * The absense of a static data member allows for many fields
  * to be implemented at the same time.
  *  
  * This parametric class is from the work of Mark Giesbrecht.
  * @see modular
  */
class param_modular
{
public:
  /** @name Common Object Interface for a LinBox Field.
    * These methods are required of all LinBox fields.
    */
  //@{
  
  /// Element type.
  typedef param_modular_element element;

  /// Random iterator generator type.
  typedef param_modular_randIter randIter;
  
  /** @name Object Management */
  //@{
  
  /** Initialization of parametric modular field element.
    * @return reference to parametric modular field element.
    * @param  x parametric modular field element to be initialized
    *           (Reference returned)
    * @param  y parametric modular field element from which to 
    *           initialize x
    */
  element& init(element& x, const element& y) const
  {
    x.residue = y.residue;
    return x;
  }
 
  /** Initialization of parametric modular field element.
    * @return reference to parametric modular field element.
    * @param  x parametric modular field element to be initialized
    *           (Reference returned)     
    * @param  y integer from which to initialize x
    */
  element& init(element& x, const integer& y) const
  {
    integer temp = y % modulus;
    if (temp < 0) temp += modulus;
    x.residue = temp;
    return x;
  }
 
  /** Access to representation of parametric modular field element.
    * @return integer representing field element.
    */
  const integer& access(const element& x) const { return x.residue; }
 
  /** Assignment of parametric field element.
    * @return reference to parametric modular field element.
    * @param  x parametric modular field element to be initialized
    *           (Reference returned)
    * @param  y parametric modular field element from which to 
    *           assign x
    */
  element& assign(element& x, const element& y) const
  {
    x.residue = y.residue;
    return x;
  }
 
  //@} Object Management
 
  /** @name Arithmetic Operators */
  //@{
 
  /** Equality.
    * @return boolean true if elements are equal, false if not.
    * @param  x parametric modular field element.
    * @param  y parametric modular field element.
    */
  bool isequal(const element& x, const element& y) const
    { return x.residue == y.residue; }

  /** Addition.
    * x = y + z
    * @return reference to x.
    * @param  x parametric modular field element 
    *           (Reference returned).
    * @param  y parametric modular field element.
    * @param  z parametric modular field element.
    */
  element& add(element& x, const element& y, const element& z) const
    { return assign(x, y.residue + z.residue); }
 
  /** Subtraction.
    * x = y - z
    * @return reference to x.
    * @param  x parametric modular field element 
    *           (Reference returned).
    * @param  y parametric modular field element.
    * @param  z parametric modular field element.
    */
  element& sub(element& x, const element& y, const element& z) const
    { return assign(x, y.residue - z.residue); }
 
  /** Multiplication.
    * x = y * z
    * @return reference to x.
    * @param  x parametric modular field element 
    *           (Reference returned).
    * @param  y parametric modular field element.
    * @param  z parametric modular field element.
    */
  element& mul(element& x, const element& y, const element& z) const
    { return assign(x, y.residue * z.residue); }
 
  /** Division.
    * x = y / z
    * @return reference to x.
    * @param  x parametric modular field element 
    *           (Reference returned).
    * @param  y parametric modular field element.
    * @param  z parametric modular field element.
    */
  element& div(element& x, const element& y, const element& z) const
  { 
    element temp(*this, 1);
    return mul(x, y, inv(temp, z));
  }
 
  /** Additive Inverse (Negation).
    * x = - y
    * @return reference to x.
    * @param  x parametric modular field element 
    *           (Reference returned).
    * @param  y parametric modular field element.
    */
  element& neg(element& x, const element& y) const
    { return assign(x, - y.residue); }
 
  /** Multiplicative Inverse.
    * x = 1 / y
    * @return reference to x.
    * @param  x parametric modular field element 
    *           (Reference returned).
    * @param  y parametric modular field element.
    */
  element& inv(element& x, const element& y) const;

  //@} Arithmetic Operations
 
  /** @name Inplace Arithmetic Operations
    * x <- x op y; x <- op x
    */
  //@{
 
  /** Zero equality.
    * Test if parametric modular field element is equal to zero.
    * @return boolean true if equals zero, false if not.
    * @param  x parametric modular field element.
    */
  bool iszero(const element& x) const { return (x.residue == 0); }
 
  /** One equality.
    * Test if parametric modular field element is equal to one.
    * @return boolean true if equals one, false if not.
    * @param  x parametric modular field element.
    */                                                       
  bool isone(const element& x) const { return (x.residue == 1); }
 
  /** Inplace Addition.
    * x = x + y
    * @return reference to x.
    * @param  x parametric modular field element 
    *           (Reference returned).
    * @param  y parametric modular field element.
    */
  element& addin(element& x, const element& y) const { return add(x, x, y); }

  /** Inplace Subtraction.
    * x = x - y
    * @return reference to x.
    * @param  x parametric modular field element 
    *           (Reference returned).
    * @param  y parametric modular field element.
    */
  element& subin(element& x, const element& y) const { return sub(x, x, y); }

  /** Inplace Multiplication.
    * x = x * y
    * @return reference to x.
    * @param  x parametric modular field element 
    *           (Reference returned).
    * @param  y parametric modular field element.
    */
  element& mulin(element& x, const element& y) const { return mul(x, x, y); }

  /** Inplace Division.
    * x = x / y
    * @return reference to x.
    * @param  x parametric modular field element 
    *           (Reference returned).
    * @param  y parametric modular field element.
    */
  element& divin(element& x, const element& y) const { return div(x, x, y); }

  /** Inplace Additive Inverse (Inplace Negation).
    * x = - x
    * @return reference to x.                                
    * @param  x parametric modular field element 
    *           (Reference returned).
    */
  element& negin(element& x) const { return neg(x, x); }
 
  /** Inplace Multiplicative Inverse.
    * x = 1 / x
    * @return reference to x.
    * @param  x parametric modular field element 
    *           (Reference returned).
    */
  element& invin(element& x) const { return inv(x, x); }
 
  //@} Inplace Arithmetic Operations

  /** @name Input/Output Operations
    */
  //@{
 
  /** Print field object.
    * @param  os  output stream to which field is printed.
    */
  void print(ostream& os) const
    { os << "Derived parametric field modulo " << modulus; }
 
  /** Read field object.
    * @param  is  input stream from which field is read.
    */
  void read(istream& is) { is >> modulus; }
 
  /** Print field element.
    * @param  os  output stream to which element is printed.
    * @param  x   parametric modular field element.
    */
  void print(ostream& os, const element& x) const
    { os << x.residue << " mod " << modulus; }
 
  /** Read field element.
    * @param  is  input stream from which element is read.
    * @param  x   parametric modular field element.
    */                                                     
  void read(istream& is, element& x) const
  {
    integer i;
    is >> i;
    init(x, i);
  }
 
  //@} Input/Output Operations
 
  //@} Common Object Interface

  /** @name Implementation-Specific Methods.
    * These methods are not required of all \Ref{LinBox Fields}
    * and are included only for this implementation of the field.
    */
  //@{
  
  /** @name Object Management */
  //@{

  /** Constructor.
    * Sets modulus to integer.
    * @param  m 15-bit prime number (default = 2)
    */
  param_modular (const integer& m = 2) : modulus(m) {}

  /** Assignment operator.
    * @param  x parametric modular field
    * @return reference to self                             
    */
  param_modular& operator=(const param_modular& x)
  {
    modulus = x.modulus;
    return *this;
  }
 
  /** Destructor.
    */
  ~param_modular() {}

  /** Retreive prime modulus.
    * @return prime modulus
    */
  const integer& get_modulus() { return modulus; }
 
  /** Change prime modulus.
    * Changes modulus of field to prime integer.
    * @param  start prime modulus
    */
  void put_modulus(const const integer& start) { modulus = start; }
  
  //@} Object management

  //@} Implementation-Specific Methods
 
private:

  friend param_modular::element;

  /** @name Implementation-Specific Data.
    * Member data not required of all \Ref{LinBox Fields} and are 
    * only included for this implementation of the field.
    */
  //@{

  /** Prime modulus.
    * Modulus is stored in the field object and not the field elements.
    */
  integer modulus;

  //@} Implementation-Specific Data

}; // class param_modular

param_modular::element& 
param_modular::inv(param_modular::element& a, 
                   const param_modular::element& z) const 
{
  integer x, y, q, tx, ty, temp;
  x = modulus; y = z.residue;
  tx = 0; ty = 1;
  
  while(y != 0) {
    // always: gcd(modulus, residue) = gcd(x,y)
    //         sx*modulus + tx*residue = x
    //         sy*modulus + ty*residue = y
    q = x / y; // integer quotient
    temp = y;  y  = x  - q*y;  x  = temp;
    temp = ty; ty = tx - q*ty; tx = temp;
  }; // while
  
  if (tx < 0) tx += modulus;
  
  a.residue = tx;

  return a;
} // inv()

#endif // _PARAM_MODULAR_FIELD_
