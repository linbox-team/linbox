/* File: src/wrappers/by_library/ntl.h
 * Author: William J Turner for the LinBox group
 */

#ifndef _NTL_
#define _NTL_

#include <NTL/RR.h>
#include <NTL/ZZ_p.h>
#include <NTL/lzz_p.h>
#include "LinBox/unparam_field.h"
#include "LinBox/unparam_randiter.h"

// Namespace in which all LinBox library code resides
namespace LinBox
{
  
  /** @name NTL Fields.
   * Partial template instantiantions for wrapping objects from
   * \URL[Victor Shoup]{http://www.shoup.net/index.html}'s
   * package \URL[NTL]{http://www.shoup.net/ntl/index.html} 5.0b
   * objects to comply with the common object interface for
   * \Ref{LinBox} fields.
   *
   * See \URL{http://www.shoup.net/ntl/doc/tour.html} for a tour of the
   * \URL[NTL]{http://www.shoup.net/ntl/index.html} library.
   */
  //@{

  /** @name class RR.
   * Arbitrary precision floating point numbers.
   * These specializations allow the \Ref{unparam_field} template class to be
   * used to wrap NTL's RR class as a LinBox field.
   */
  //@{

  /** Initialization of field element from an integer.
   * Behaves like C++ allocator construct.
   * This function assumes the output field element x has already been
   * constructed, but that it is not already initialized.
   * For now, this is done by converting the integer type to a C++
   * long and then to the element type through the use of static cast and
   * NTL's to_RR function.
   * This, of course, assumes such static casts are possible.
   * This function should be changed in the future to avoid using long.
   * @return reference to field element.
   * @param x field element to contain output (reference returned).
   * @param y integer.
   */
  template <>
  NTL::RR& unparam_field<NTL::RR>::init(NTL::RR& x, const integer& y) const
  { return x = NTL::to_RR(static_cast<const long&>(y)); }

  /** Conversion of field element to an integer.
   * This function assumes the output field element x has already been
   * constructed, but that it is not already initialized.
   * For now, this is done by converting the element type to a C++
   * long and then to the integer type through the use of static cast and
   * NTL's to_long function.
   * This, of course, assumes such static casts are possible.
   * This function should be changed in the future to avoid using long.
   * @return reference to integer.
   * @param x reference to integer to contain output (reference returned).
   * @param y constant reference to field element.
   */
  template <>
  integer& unparam_field<NTL::RR>::convert(integer& x, const NTL::RR& y) const
  { return x = static_cast<integer>(to_long(y)); }

  /** Multiplicative Inverse.
   * x = 1 / y
   * This function assumes both field elements have already been
   * constructed and initialized.
   * @return reference to x.
   * @param  x field element (reference returned).
   * @param  y field element.
   */
  template <> 
  NTL::RR& unparam_field<NTL::RR>::inv(NTL::RR& x, const NTL::RR& y) const
  { return x = NTL::inv(y); }
 
  /** Zero equality.
   * Test if field element is equal to zero.
   * This function assumes the field element has already been
   * constructed and initialized.
   * In this specialization, NTL's IsZero function is called.
   * @return boolean true if equals zero, false if not.
   * @param  x field element.
   */
  template <> bool unparam_field<NTL::RR>::isZero(const NTL::RR& x) const
  { return static_cast<bool>(IsZero(x)); }

  /** One equality.
   * Test if field element is equal to one.
   * This function assumes the field element has already been
   * constructed and initialized.
   * In this specialization, NTL's IsOne function is called.
   * @return boolean true if equals one, false if not.
   * @param  x field element.
   */
  template <> bool unparam_field<NTL::RR>::isOne(const NTL::RR& x) const
  { return static_cast<bool>(IsOne(x)); }

  /** Inplace Multiplicative Inverse.
   * x = 1 / x
   * This function assumes both field elements have already been
   * constructed and initialized.
   * @return reference to x.
   * @param  x field element (reference returned).
   */
  template <> NTL::RR& unparam_field<NTL::RR>::invin(NTL::RR& x) const
  { return x = NTL::inv(x); }

	/** Print field.
	 * @return output stream to which field is written.
	 * @param  os  output stream to which field is written.
	 */
	 template <> ostream& unparam_field<NTL::RR>::write(ostream& os) const 
	{ return os << "unparamterized field NTL::RR"; }

  /** Random field element creator.
   * This returns a random field element from the information supplied
   * at the creation of the generator.
   * This generator uses the built-in C++ random number generator instead of
   * NTL's random function because the NTL function does not allow as much
   * control over the sampling size as the generic LinBox template.  This 
   * specialization is included only to allow conversion to an NTL 
   * object.
   * @return random field element
   */
  template <> NTL::RR& unparam_randIter<NTL::RR>::operator() (void)
  {
    NTL::RR* temp_ptr = new NTL::RR();
    double temp;
    
    // Create new random elements
    if (_size == 0)
      temp = rand();
    else
      temp = static_cast<long>((double(rand())/RAND_MAX)*double(_size));

    *temp_ptr = temp;

#ifdef TRACE
    cout << "random double = " << temp 
         << "    random NTL::RR = " << *temp_ptr << endl;
#endif // TRACE

    return *(temp_ptr);
    
  } // element& operator() (void)

  //@} // class RR

  /** @name class ZZ_p.
   * Arbitrary precision integers modulus a positive integer.
   * While NTL allows any integer to serve as the modulus, only prime
   * moduli yield fields.  Therefore, while arthmetic operations may be
   * valid for any modulus, only prime moduli are supported in this
   * implementation.  The primality of the modulus will not be checked, so
   * it is the programmer's responsibility to supply a prime modulus.
   * These specializations allow the \Ref{unparam_field} template class to be
   * used to wrap NTL's ZZ_p class as a LinBox field.
   */
  //@{

  /** Initialization of field element from an integer.
   * Behaves like C++ allocator construct.
   * This function assumes the output field element x has already been
   * constructed, but that it is not already initialized.
   * For now, this is done by converting the integer type to a C++
   * long and then to the element type through the use of static cast and
   * NTL's to_ZZ_p function.
   * This, of course, assumes such static casts are possible.
   * This function should be changed in the future to avoid using long.
   * @return reference to field element.
   * @param x field element to contain output (reference returned).
   * @param y integer.
   */
  template <>
  NTL::ZZ_p& unparam_field<NTL::ZZ_p>::init(NTL::ZZ_p& x, const integer& y) const
  { return x = NTL::to_ZZ_p(static_cast<const long&>(y)); }

  /** Conversion of field element to an integer.
   * This function assumes the output field element x has already been
   * constructed, but that it is not already initialized.
   * For now, this is done by converting the element type to a C++
   * long and then to the integer type through the use of static cast and
   * NTL's to_long function.
   * This, of course, assumes such static casts are possible.
   * This function should be changed in the future to avoid using long.
   * @return reference to integer.
   * @param x reference to integer to contain output (reference returned).
   * @param y constant reference to field element.
   */
  template <>
  integer& unparam_field<NTL::ZZ_p>::convert(integer& x, const NTL::ZZ_p& y) const
  { return x = static_cast<integer>(to_long(rep(y))); }

  /** Cardinality.
   * Return integer representing cardinality of the field.
   * Returns the modulus of the field, which should be prime.
   * @return integer representing cardinality of the field
   */
  template <> 
  integer& unparam_field<NTL::ZZ_p>::cardinality(integer& c) const
  { return c = static_cast<integer>(to_long(NTL::ZZ_p::modulus())); }

  /** Characteristic.
   * Return integer representing characteristic of the field.
   * Returns the modulus of the field, which should be prime.
   * @return integer representing characteristic of the field.
   */
  template <> 
  integer& unparam_field<NTL::ZZ_p>::characteristic(integer& c) const
  //FIXME we shouldn't go thru long here as p may be larger than that.
  // check if NTL has cast ZZp to gmp integers.
  { return c = static_cast<integer>(to_long(NTL::ZZ_p::modulus())); }

   /** Multiplicative Inverse.
   * x = 1 / y
   * This function assumes both field elements have already been
   * constructed and initialized.
   * @return reference to x.
   * @param  x field element (reference returned).
   * @param  y field element.
   */
  template <> NTL::ZZ_p& 
  unparam_field<NTL::ZZ_p>::inv(NTL::ZZ_p& x, const NTL::ZZ_p& y) const
  { return x = NTL::inv(y); }
 
  /** Zero equality.
   * Test if field element is equal to zero.
   * This function assumes the field element has already been
   * constructed and initialized.
   * In this specialization, NTL's IsZero function is called.
   * @return boolean true if equals zero, false if not.
   * @param  x field element.
   */
  template <> bool unparam_field<NTL::ZZ_p>::isZero(const NTL::ZZ_p& x) const
  { return static_cast<bool>(IsZero(x)); }

  /** One equality.
   * Test if field element is equal to one.
   * This function assumes the field element has already been
   * constructed and initialized.
   * In this specialization, NTL's IsOne function is called.
   * @return boolean true if equals one, false if not.
   * @param  x field element.
   */
  template <> bool unparam_field<NTL::ZZ_p>::isOne(const NTL::ZZ_p& x) const
  { return static_cast<bool>(IsOne(x)); }

  /** Inplace Multiplicative Inverse.
   * x = 1 / x
   * This function assumes both field elements have already been
   * constructed and initialized.
   * @return reference to x.
   * @param  x field element (reference returned).
   */
  template <> NTL::ZZ_p& unparam_field<NTL::ZZ_p>::invin(NTL::ZZ_p& x) const
  { return x = NTL::inv(x); }

	/** Print field.
	 * @return output stream to which field is written.
	 * @param  os  output stream to which field is written.
	 */
	 template <> ostream& unparam_field<NTL::ZZ_p>::write(ostream& os) const 
	{ 
		return os << "unparamterized field NTL::ZZ_p with p = " 
			<< NTL::ZZ_p::modulus(); 
	}

  /** Random field element creator.
   * This returns a random field element from the information supplied
   * at the creation of the generator.
   * This generator uses the built-in C++ random number generator instead of
   * NTL's random function because the NTL function does not allow as much
   * control over the sampling size as the generic LinBox template.  This 
   * specialization is included only to allow conversion to an NTL 
   * object.
   * @return random field element
   */
  template <> NTL::ZZ_p& unparam_randIter<NTL::ZZ_p>::operator() (void)
  {
    NTL::ZZ_p* temp_ptr = new NTL::ZZ_p();
    long temp;
    
    // Create new random elements
    temp = static_cast<long>((double(rand())/RAND_MAX)*double(_size));

    *temp_ptr = temp;

#ifdef TRACE
    cout << "random long = " << temp 
         << "    random NTL::ZZ_p = " << *temp_ptr << endl;
#endif // TRACE

    return *(temp_ptr);
    
  } // element& operator() (void)

  //@} // class ZZ_p

  /** @name class zz_p.
   * Arbitrary precision integers modulus a positive integer.
   * While NTL allows any integer to serve as the modulus, only prime
   * moduli yield fields.  Therefore, while arthmetic operations may be
   * valid for any modulus, only prime moduli are supported in this
   * implementation.  The primality of the modulus will not be checked, so
   * it is the programmer's responsibility to supply a prime modulus.
   * These specializations allow the \Ref{unparam_field} template class to be
   * used to wrap NTL's zz_p class as a LinBox field.
   */
  //@{

  /** Initialization of field element from an integer.
   * Behaves like C++ allocator construct.
   * This function assumes the output field element x has already been
   * constructed, but that it is not already initialized.
   * For now, this is done by converting the integer type to a C++
   * long and then to the element type through the use of static cast and
   * NTL's to_zz_p function.
   * This, of course, assumes such static casts are possible.
   * This function should be changed in the future to avoid using long.
   * @return reference to field element.
   * @param x field element to contain output (reference returned).
   * @param y integer.
   */
  template <>
  NTL::zz_p& unparam_field<NTL::zz_p>::init(NTL::zz_p& x, const integer& y) const
  { return x = NTL::to_zz_p(static_cast<const long&>(y)); }

  /** Conversion of field element to an integer.
   * This function assumes the output field element x has already been
   * constructed, but that it is not already initialized.
   * For now, this is done by converting the element type to a C++
   * long and then to the integer type through the use of static cast and
   * NTL's to_long function.
   * This, of course, assumes such static casts are possible.
   * This function should be changed in the future to avoid using long.
   * @return reference to integer.
   * @param x reference to integer to contain output (reference returned).
   * @param y constant reference to field element.
   */
  template <>
  integer& unparam_field<NTL::zz_p>::convert(integer& x, const NTL::zz_p& y) const
  { return x = static_cast<integer>(rep(y)); }

  /** Cardinality.
   * Return integer representing cardinality of the field.
   * Returns the modulus of the field, which should be prime.
   * @return integer representing cardinality of the field
   */
  template <> 
  integer& unparam_field<NTL::zz_p>::cardinality(integer& c) const
  { return c = static_cast<integer>(NTL::zz_p::modulus()); }

  /** Characteristic.
   * Return integer representing characteristic of the field.
   * Returns the modulus of the field, which should be prime.
   * @return integer representing characteristic of the field.
   */
  template <> 
  integer& unparam_field<NTL::zz_p>::characteristic(integer& c) const
  { return c = static_cast<integer>(NTL::zz_p::modulus()); }

   /** Multiplicative Inverse.
   * x = 1 / y
   * This function assumes both field elements have already been
   * constructed and initialized.
   * @return reference to x.
   * @param  x field element (reference returned).
   * @param  y field element.
   */
  template <> NTL::zz_p& 
  unparam_field<NTL::zz_p>::inv(NTL::zz_p& x, const NTL::zz_p& y) const
  { return x = NTL::inv(y); }
 
  /** Zero equality.
   * Test if field element is equal to zero.
   * This function assumes the field element has already been
   * constructed and initialized.
   * In this specialization, NTL's IsZero function is called.
   * @return boolean true if equals zero, false if not.
   * @param  x field element.
   */
  template <> bool unparam_field<NTL::zz_p>::isZero(const NTL::zz_p& x) const
  { return static_cast<bool>(IsZero(x)); }

  /** One equality.
   * Test if field element is equal to one.
   * This function assumes the field element has already been
   * constructed and initialized.
   * In this specialization, NTL's IsOne function is called.
   * @return boolean true if equals one, false if not.
   * @param  x field element.
   */
  template <> bool unparam_field<NTL::zz_p>::isOne(const NTL::zz_p& x) const
  { return static_cast<bool>(IsOne(x)); }

  /** Inplace Multiplicative Inverse.
   * x = 1 / x
   * This function assumes both field elements have already been
   * constructed and initialized.
   * @return reference to x.
   * @param  x field element (reference returned).
   */
  template <> NTL::zz_p& unparam_field<NTL::zz_p>::invin(NTL::zz_p& x) const
  { return x = NTL::inv(x); }

	/** Print field.
	 * @return output stream to which field is written.
	 * @param  os  output stream to which field is written.
	 */
	 template <> ostream& unparam_field<NTL::zz_p>::write(ostream& os) const 
	{ 
		return os << "unparamterized field NTL::zz_p with p = " 
			<< NTL::zz_p::modulus(); 
	}

  /** Random field element creator.
   * This returns a random field element from the information supplied
   * at the creation of the generator.
   * This generator uses the built-in C++ random number generator instead of
   * NTL's random function because the NTL function does not allow as much
   * control over the sampling size as the generic LinBox template.  This 
   * specialization is included only to allow conversion to an NTL 
   * object.
   * @return random field element
   */
  template <> NTL::zz_p& unparam_randIter<NTL::zz_p>::operator() (void)
  {
    NTL::zz_p* temp_ptr = new NTL::zz_p();
    long temp;
    
    // Create new random elements
    temp = static_cast<long>((double(rand())/RAND_MAX)*double(_size));

    *temp_ptr = temp;

#ifdef TRACE
    cout << "random long = " << temp 
         << "    random NTL::zz_p = " << *temp_ptr << endl;
#endif // TRACE

    return *(temp_ptr);
    
  } // element& operator() (void)

  

  //@} // class zz_p

  //@} NTL Fields

} // namespace LinBox

#endif // _NTL_
