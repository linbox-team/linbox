/* File: src/library/objects/field/NTL_lzz_Field.h
 * Author: Li Chen for the LinBox group
 */

#ifndef _NTL_lzz_Field_
#define _NTL_lzz_Field_

#include <iostream>
#include <fstream>
#include <NTL/lzz_p.h>
#include "LinBox/integer.h"

// Namespace in which all LinBox code resides
namespace LinBox 
{ 

  /** NTL Field lzz_p wrapper class .
   */
  class  NTL_lzz_Field
  {
  public:

    /// Element type.
    typedef zz_p element;

    /// Random iterator generator type.
    /* typedef RandIter_abstract randIter; */
 
    /** @name Object Management
     * There are no public constructors for this class.
     * It should only be used in tandem with \Ref{Field_archetype}.
     */
    //@{

    /** Destructor.
     */
    ~NTL_lzz_Field(void) {}

    /**  copy constructor.
     * @return pointer to new object in dynamic memory.
     */
    NTL_lzz_Field* clone() const{
      NTL_lzz_Field* x;
      x = new(NTL_lzz_Field)();
      return x;
    }

    /** Assignment operator.
     */
    NTL_lzz_Field& operator= (const NTL_lzz_Field& F){return *this;}

    /** Initialization of field element from an integer.
     * Behaves like C++ allocator construct.
     * This function assumes the output field element x has already been
     * constructed, but that it is not already initialized.
     * @return reference to field element.
     * @param x field element to contain output (reference returned).
     * @param y integer.
     */
    element& init(element& x, const integer& y) const{return x = y;}
 
   /** Conversion of field element to an integer.
     * This function assumes the output field element x has already been
     * constructed, but that it is not already initialized.
     * @return reference to integer.
     * @param x reference to interger to contain output (reference returned).
     * @param y constant field element.
     */
    integer& convert(integer& x, const element& y) const{return x = rep(y);}
 
    /** Assignment of one field element to another.
     * This function assumes both field elements have already been
     * constructed and initialized.
     * @return reference to x
     * @param  x field element (reference returned).
     * @param  y field element.
     */
    element& assign(element& x, const element& y) const{return x = y;}

    /** Cardinality.
     * Return integer representing cardinality of the domain.
     * Returns a non-negative integer for all domains with finite
     * cardinality, and returns -1 to signify a domain of infinite
     * cardinality.
     * @return integer representing cardinality of the domain
     */
    const integer& cardinality(void) const{
      return zz_pInfo->p;
    }
 
    /** Characteristic.
     * Return integer representing characteristic of the domain.
     * Returns a positive integer to all domains with finite characteristic,
     * and returns 0 to signify a domain of infinite characteristic.
     * @return integer representing characteristic of the domain.
     */
    const integer& characteristic(void) const{
      return zz_pInfo->p;
    }
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
    bool areEqual(const element& x, const element& y) const{return (x == y);}

    /** Addition.
     * x = y + z
     * This function assumes all the field elements have already been
     * constructed and initialized.
     * @return reference to x.
     * @param  x field element (reference returned).
     * @param  y field element.
     * @param  z field element.
     */
    element& add(element& x, const element& y, const element& z) const{return x = y + z;}
 
    /** Subtraction.
     * x = y - z
     * This function assumes all the field elements have already been
     * constructed and initialized.
     * @return reference to x.
     * @param  x field element (reference returned).
     * @param  y field element.
     * @param  z field element.
     */
    element& sub(element& x, const element& y, const element& z) const{return x = y - z;}
 
    /** Multiplication.
     * x = y * z
     * This function assumes all the field elements have already been
     * constructed and initialized.
     * @return reference to x.
     * @param  x field element (reference returned).
     * @param  y field element.
     * @param  z field element.
     */

    element& mul(element& x, const element& y, const element& z) const{return x = y * z;}
 
    /** Division.
     * x = y / z
     * This function assumes all the field elements have already been
     * constructed and initialized.
     * @return reference to x.
     * @param  x field element (reference returned).
     * @param  y field element.
     * @param  z field element.
     */
    element& div(element& x, const element& y, const element& z) const{return x = y / z;}
 
    /** Additive Inverse (Negation).
     * x = - y
     * This function assumes both field elements have already been
     * constructed and initialized.
     * @return reference to x.
     * @param  x field element (reference returned).
     * @param  y field element.
     */
    element& neg(element& x, const element& y) const{return x = -y;}
 
    /** Multiplicative Inverse.
     * x = 1 / y
     * This function assumes both field elements have already been
     * constructed and initialized.
     * @return reference to x.
     * @param  x field element (reference returned).
     * @param  y field element.
     */
    element& inv(element& x, const element& y) const{return x = 1 / y;}

    //@} Arithmetic Operations
 
    /** @name Inplace Arithmetic Operations
     * x <- x op y; x <- op x
     */
    //@{

    /** Zero equality.
     * Test if field element is equal to zero.
     * This function assumes the field element has already been
     * constructed and initialized.
     * @return boolean true if equals zero, false if not.
     * @param  x field element.
     */
    bool isZero(const element& x) const{return x == 0;}
 
    /** One equality.
     * Test if field element is equal to one.
     * This function assumes the field element has already been
     * constructed and initialized.
     * @return boolean true if equals one, false if not.
     * @param  x field element.
     */
    bool isOne(const element& x) const{return x == 1;}

    /** Inplace Addition.
     * x += y
     * This function assumes both field elements have already been
     * constructed and initialized.
     * @return reference to x.
     * @param  x field element (reference returned).
     * @param  y field element.
     */
    element& addin(element& x, const element& y) const{return x += y;}
 
    /** Inplace Subtraction.
     * x -= y
     * This function assumes both field elements have already been
     * constructed and initialized.
     * @return reference to x.
     * @param  x field element (reference returned).
     * @param  y field element.
     */
    element& subin(element& x, const element& y) const{return x -= y;}
 
    /** Inplace Multiplication.
     * x *= y
     * This function assumes both field elements have already been
     * constructed and initialized.
     * @return reference to x.
     * @param  x field element (reference returned).
     * @param  y field element.
     */
    element& mulin(element& x, const element& y) const{return x *= y;}
 
    /** Inplace Division.
     * x /= y
     * This function assumes both field elements have already been
     * constructed and initialized.
     * @return reference to x.
     * @param  x field element (reference returned).
     * @param  y field element.
     */
    element& divin(element& x, const element& y) const{return x /= y;}
 
    /** Inplace Additive Inverse (Inplace Negation).
     * x = - x
     * This function assumes the field element has already been
     * constructed and initialized.
     * @return reference to x.
     * @param  x field element (reference returned).
     */
    element& negin(element& x) const{return x = -x;}
 
    /** Inplace Multiplicative Inverse.
     * x = 1 / x
     * This function assumes the field elementhas already been
     * constructed and initialized.
     * @return reference to x.
     * @param  x field element (reference returned).
     */
    element& invin(element& x) const{return x = 1 / x;} 

    //@} Inplace Arithmetic Operations

    /** @name Input/Output Operations */
    //@{

    /** Print field.
     * @return output stream to which field is written.
     * @param  os  output stream to which field is written.
     */
    ostream& write(ostream& os) const{return os;}
 
    /** Read field.
     * @return input stream from which field is read.
     * @param  is  input stream from which field is read.
     */
    istream& read(istream& is) {return is;}

    /** Print field element.
     * This function assumes the field element has already been
     * constructed and initialized.
     * @return output stream to which field element is written.
     * @param  os  output stream to which field element is written.
     * @param  x   field element.
     */
    ostream& write(ostream& os, const element& x) const{return os << rep(x);}
 
    /** Read field element.
     * This function assumes the field element has already been
     * constructed and initialized.
     * @return input stream from which field element is read.
     * @param  is  input stream from which field element is read.
     * @param  x   field element.
     */
    istream& read(istream& is, element& x) const{ 
      long m;
      is >> m;
      x = m;
      return is;
    }
    //@}

    // size function is to return the size of element, it is not required  in Field archetype.
    long size(){return 1;}

    /** Default Constructor.
     */
    NTL_lzz_Field () {}
 
  private:

    /// Field_archetype is friend.
    friend class Field_archetype;


  }; // class NTL_lzz_Field

} // namespace LinBox

#endif // _NTL_lzz_Field_
