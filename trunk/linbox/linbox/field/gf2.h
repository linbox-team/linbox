/* linbox/field/gf2.h
 * Copyright (C) 2003-2007 The LinBox group
 *
 * Authors : B. Hovinen, JG Dumas, C. Pernet
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __LINBOX_field_gf2_H
#define __LINBOX_field_gf2_H

#include <iostream>
#include <climits>
#include <cmath>

#include "linbox/integer.h"
#include "linbox/field/field-interface.h"
#include "linbox/util/debug.h"
#include "linbox/vector/bit-vector.h"
#include "linbox/linbox-config.h"
#include "linbox/field/field-traits.h"

// Namespace in which all LinBox code resides
namespace LinBox 
{

class GF2RandIter;

/** 
 * \brief Integers modulo 2
 *
 * This is a tuned implementation of the field of integers modulo
 * 2. In particular, when one constructs a VectorDomain object over
 * this field, highly optimized bit operations will be used to make
 * vector arithmetic very fast.
 \ingroup field
 */

template <class Ring>
struct ClassifyRing;

class GF2;

template<>
struct ClassifyRing<GF2> {
	typedef RingCategories::ModularTag categoryTag;
};

class GF2 : public FieldInterface
{
    public:
    const bool zero,one,mone;
    

	/** Element type
	 */
	typedef bool Element;

	/** Random iterator generator type.
	 * It must meet the common object interface of random element generators
	 * as given in the the archetype RandIterArchetype.
	 */
	typedef GF2RandIter RandIter;

	/** @name Object Management
	 */
	//@{
 
	/** Default constructor.
	 */
	GF2 () : zero(false),one(true),mone(true) {}
	GF2 (int p, int exp = 1) : zero(false),one(true),mone(true) {
		if(p != 2) throw PreconditionFailed(__FUNCTION__,__LINE__,"modulus must be 2");
		if(exp != 1) throw PreconditionFailed(__FUNCTION__,__LINE__,"exponent must be 1");
	}

	/** Copy constructor.
	 * Constructs Modular object by copying the field.
	 * This is required to allow field objects to be passed by value
	 * into functions.
	 * @param  F Modular object.
	 */
	GF2 (const GF2 &) : zero(false),one(true),mone(true) {}

	/** Assignment operator
	 * Required by the archetype
	 *
	 * @param F constant reference to Modular object
	 * @return reference to Modular object for self
	 */
	const GF2 &operator = (const GF2 &) 
		{ return *this; }

	/** Initialization of field base element from an integer.
	 * Behaves like C++ allocator construct.
	 * This function assumes the output field base element x has already been
	 * constructed, but that it is not already initialized.
	 * This is not a specialization of the template function because
	 * such a specialization is not allowed inside the class declaration.
	 * @return reference to field base element.
	 * @param x field base element to contain output (reference returned).
	 * @param y integer.
	 */
	Element &init (Element &x, const int &y = 0) const
		{ return x = y & 1; }

	Element &init (Element &x, const unsigned int &y = 0) const
		{ return x = y & 1; }

	Element &init (Element &x, const long &y = 0) const
		{ return x = y & 1; }

	Element &init (Element &x, const unsigned long &y = 0) const
		{ return x = y & 1; }

	Element &init (Element &x, const float &y) const
		{ return x = static_cast<unsigned char>(y) & 1; }

	Element &init (Element &x, const double &y) const
		{ return x = static_cast<unsigned char>(y) & 1; }

	Element &init (Element &x, const integer &y) const
		{ return x = static_cast<long>(y) & 1; }

	BitVector::reference init (BitVector::reference x, const integer &y = 0) const
		{ return x = long (y) & 1; }

	std::_Bit_reference init (std::_Bit_reference x, const integer &y = 0) const
		{ return x = long (y) & 1; }

	/** Conversion of field base element to a template class T.
	 * This function assumes the output field base element x has already been
	 * constructed, but that it is not already initialized.
	 * @return reference to template class T.
	 * @param x template class T to contain output (reference returned).
	 * @param y constant field base element.
	 */
	integer &convert (integer &x, Element y) const
		{ return x = y; }
 
        std::_Bit_reference convert (std::_Bit_reference x, Element y) const {
            return x = y;
        }

    	template<class XXX> 
        XXX& convert (XXX& x, Element y) const {
            return x = static_cast<XXX>(y);
        }
            
// 	unsigned int &convert (unsigned int &x, Element y) const
// 		{ return x = static_cast<unsigned int>(y); }
 
// 	int &convert (int &x, Element y) const
// 		{ return x = static_cast<int>(y); }
 
// 	unsigned long &convert (unsigned long &x, Element y) const
// 		{ return x = static_cast<unsigned long>(y); }
 
// 	long &convert (long &x, Element y) const
// 		{ return x = static_cast<int>(y); }
 
// 	float &convert (float &x, Element y) const
// 		{ return x = static_cast<float>(y); }
 
// 	double &convert (double &x, Element y) const
// 		{ return x = static_cast<double>(y); }
 
	/** Assignment of one field base element to another.
	 * This function assumes both field base elements have already been
	 * constructed and initialized.
	 * @return reference to x
	 * @param  x field base element (reference returned).
	 * @param  y field base element.
	 */
	Element &assign (Element &x, Element y) const
		{ return x = y; }

        BitVector::reference assign (BitVector::reference x, Element y) const
		{ return x = y; }

        std::_Bit_reference assign (std::_Bit_reference x, Element y) const
		{ return x = y; }
    
	/** Cardinality.
	 * Return integer representing cardinality of the domain.
	 * Returns a non-negative integer for all domains with finite
	 * cardinality, and returns -1 to signify a domain of infinite
	 * cardinality.
	 * @return integer representing cardinality of the domain
	 */
	integer &cardinality (integer &c) const
		{ return c = 2; }

	/** Characteristic.
	 * Return integer representing characteristic of the domain.
	 * Returns a positive integer to all domains with finite characteristic,
	 * and returns 0 to signify a domain of infinite characteristic.
	 * @return integer representing characteristic of the domain.
	 */
	integer &characteristic (integer &c) const
		{ return c = 2; }

	//@} Object Management

	/** @name Arithmetic Operations
	 * x <- y op z; x <- op y
	 * These operations require all elements, including x, to be initialized
	 * before the operation is called.  Uninitialized field base elements will
	 * give undefined results.
	 */
	//@{

	/** Equality of two elements.
	 * This function assumes both field base elements have already been
	 * constructed and initialized.
	 * @return boolean true if equal, false if not.
	 * @param  x field base element
	 * @param  y field base element
	 */
	bool areEqual (Element x, Element y) const
		{ return x == y; }

	/** Zero equality.
	 * Test if field base element is equal to zero.
	 * This function assumes the field base element has already been
	 * constructed and initialized.
	 * @return boolean true if equals zero, false if not.
	 * @param  x field base element.
	 */
	bool isZero (Element x) const
		{ return !x; }
 
	/** One equality.
	 * Test if field base element is equal to one.
	 * This function assumes the field base element has already been
	 * constructed and initialized.
	 * @return boolean true if equals one, false if not.
	 * @param  x field base element.
	 */
	bool isOne (Element x) const
		{ return x; }

	//@} Arithmetic Operations

	/** @name Input/Output Operations */
	//@{

	/** Print field.
	 * @return output stream to which field is written.
	 * @param  os  output stream to which field is written.
	 */
	std::ostream &write (std::ostream &os) const 
		{ return os << "integers mod 2"; }

	/** Read field.
	 * @return input stream from which field is read.
	 * @param  is  input stream from which field is read.
	 */
	std::istream &read (std::istream &is)
		{ return is; }

	/** Print field base element.
	 * This function assumes the field base element has already been
	 * constructed and initialized.
	 * @return output stream to which field base element is written.
	 * @param  os  output stream to which field base element is written.
	 * @param  x   field base element.
	 */
	std::ostream &write (std::ostream &os, Element x) const
		{ return os << x; }
 
	/** Read field base element.
	 * This function assumes the field base element has already been
	 * constructed and initialized.
	 * @return input stream from which field base element is read.
	 * @param  is  input stream from which field base element is read.
	 * @param  x   field base element.
	 */
	std::istream &read (std::istream &is, Element &x) const
		{ is >> x; return is; }

	std::istream &read (std::istream &is, BitVector::reference x) const
		{ is >> x; return is; }

	std::istream &read (std::istream &is, std::_Bit_reference x) const
		{ bool a; is >> a; x=a; return is; }

	//@}

	/** @name Arithmetic Operations
	 * x <- y op z; x <- op y
	 * These operations require all elements, including x, to be initialized
	 * before the operation is called.  Uninitialized field base elements will
	 * give undefined results.
	 */
	//@{

	/** Addition.
	 * x = y + z
	 * This function assumes all the field base elements have already been
	 * constructed and initialized.
	 * @return reference to x.
	 * @param  x field base element (reference returned).
	 * @param  y field base element.
	 * @param  z field base element.
	 */
	Element &add (Element &x, Element y, Element z) const
		{ return x = y ^ z; }

	BitVector::reference add (BitVector::reference x, Element y, Element z) const
		{ return x = y ^ z; }
 
	std::_Bit_reference add (std::_Bit_reference x, Element y, Element z) const
		{ return x = y ^ z; }
 
	/** Subtraction.
	 * x = y - z
	 * This function assumes all the field base elements have already been
	 * constructed and initialized.
	 * @return reference to x.
	 * @param  x field base element (reference returned).
	 * @param  y field base element.
	 * @param  z field base element.
	 */
	Element &sub (Element &x, Element y, Element z) const
		{ return x = y ^ z; }

	BitVector::reference sub (BitVector::reference x, Element y, Element z) const
		{ return x = y ^ z; }
 
	std::_Bit_reference sub (std::_Bit_reference x, Element y, Element z) const
		{ return x = y ^ z; }
 
	/** Multiplication.
	 * x = y * z
	 * This function assumes all the field base elements have already been
	 * constructed and initialized.
	 * @return reference to x.
	 * @param  x field base element (reference returned).
	 * @param  y field base element.
	 * @param  z field base element.
	 */
	Element &mul (Element &x, Element y, Element z) const
		{ return x = y & z; }

	BitVector::reference mul (BitVector::reference x, Element y, Element z) const
		{ return x = y & z; }
 
	std::_Bit_reference mul (std::_Bit_reference x, Element y, Element z) const
		{ return x = y & z; }
 
	/** Division.
	 * x = y / z
	 * This function assumes all the field base elements have already been
	 * constructed and initialized.
	 * @return reference to x.
	 * @param  x field base element (reference returned).
	 * @param  y field base element.
	 * @param  z field base element.
	 */
	Element &div (Element &x, Element y, Element ) const // z is unused!
		{ return x = y; }

	BitVector::reference div (BitVector::reference x, Element y, Element ) const
		{ return x = y; }
 
	std::_Bit_reference div (std::_Bit_reference x, Element y, Element ) const
		{ return x = y; }
 
	/** Additive Inverse (Negation).
	 * x = - y
	 * This function assumes both field base elements have already been
	 * constructed and initialized.
	 * @return reference to x.
	 * @param  x field base element (reference returned).
	 * @param  y field base element.
	 */
	Element &neg (Element &x, Element y) const
		{ return x = y; }

	BitVector::reference neg (BitVector::reference x, Element y) const
		{ return x = y; }
 
	std::_Bit_reference neg (std::_Bit_reference x, Element y) const
		{ return x = y; }
 
	/** Multiplicative Inverse.
	 * x = 1 / y
	 * This function assumes both field base elements have already been
	 * constructed and initialized.
	 * @return reference to x.
	 * @param  x field base element (reference returned).
	 * @param  y field base element.
	 */
	Element &inv (Element &x, Element y) const
		{ return x = y; }

	BitVector::reference inv (BitVector::reference x, Element y) const
		{ return x = y; }

	std::_Bit_reference inv (std::_Bit_reference x, Element y) const
		{ return x = y; }

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
	BitVector::reference axpy (BitVector::reference r, 
				   Element a, 
				   Element x, 
				   Element y) const
		{ return r = (a & x) ^ y; }

	std::_Bit_reference axpy (std::_Bit_reference r, 
				   Element a, 
				   Element x, 
				   Element y) const
		{ return r = (a & x) ^ y; }

	Element &axpy (Element &r, Element a, Element x, Element y) const
		{ return r = (a & x) ^ y; }

	//@} Arithmetic Operations
 
	/** @name Inplace Arithmetic Operations
	 * x <- x op y; x <- op x
	 */
	//@{

	/** Inplace Addition.
	 * x += y
	 * This function assumes both field base elements have already been
	 * constructed and initialized.
	 * @return reference to x.
	 * @param  x field base element (reference returned).
	 * @param  y field base element.
	 */
	Element &addin (Element &x, Element y) const
		{ return x ^= y; }

	BitVector::reference addin (BitVector::reference x, Element y) const
		{ return x ^= y; }
 
	std::_Bit_reference addin (std::_Bit_reference x, Element y) const
        	{ return x = x ^ y; }
 
	/** Inplace Subtraction.
	 * x -= y
	 * This function assumes both field base elements have already been
	 * constructed and initialized.
	 * @return reference to x.
	 * @param  x field base element (reference returned).
	 * @param  y field base element.
	 */
	Element &subin (Element &x, Element y) const
		{ return x ^= y; }

	BitVector::reference subin (BitVector::reference x, Element y) const
		{ return x ^= y; }
 
	std::_Bit_reference subin (std::_Bit_reference x, Element y) const
		{ return x = x ^ y; }
 
	/** Inplace Multiplication.
	 * x *= y
	 * This function assumes both field base elements have already been
	 * constructed and initialized.
	 * @return reference to x.
	 * @param  x field base element (reference returned).
	 * @param  y field base element.
	 */
	Element &mulin (Element &x, Element y) const
		{ return x &= y; }

	BitVector::reference mulin (BitVector::reference x, Element y) const
		{ return x &= y; }
    
	Element& mulin (std::_Bit_reference& x, Element y) const
		{ return mulin((bool&)x,y); }
 
	/** Inplace Division.
	 * x /= y
	 * This function assumes both field base elements have already been
	 * constructed and initialized.
	 * @return reference to x.
	 * @param  x field base element (reference returned).
	 * @param  y field base element.
	 */
	Element &divin (Element &x, Element ) const //y is unsued !
		{ return x; }

	BitVector::reference divin (BitVector::reference x, Element ) const
		{ return x; }
 
	std::_Bit_reference divin (std::_Bit_reference x, Element ) const
		{ return x; }
 
	/** Inplace Additive Inverse (Inplace Negation).
	 * x = - x
	 * This function assumes the field base element has already been
	 * constructed and initialized.
	 * @return reference to x.
	 * @param  x field base element (reference returned).
	 */
	Element &negin (Element &x) const
		{ return x; }

	BitVector::reference negin (BitVector::reference x) const
		{ return x; }
 
	std::_Bit_reference negin (std::_Bit_reference x) const
		{ return x; }
 
	/** Inplace Multiplicative Inverse.
	 * x = 1 / x
	 * This function assumes the field base elementhas already been
	 * constructed and initialized.
	 * @return reference to x.
	 * @param  x field base element (reference returned).
	 */
	Element &invin (Element &x) const
		{ return x; }

	BitVector::reference invin (BitVector::reference x) const
		{ return x; }

	std::_Bit_reference invin (std::_Bit_reference x) const
		{ return x; }

	/** Inplace AXPY.
	 * r  += a * x
	 * This function assumes all field elements have already been 
	 * constructed and initialized.
	 * Purely virtual
	 * @return reference to r.
	 * @param  r field element (reference returned).
	 * @param  a field element.
	 * @param  x field element.
	 */
	Element &axpyin (Element &r, Element a, Element x) const
		{ return r ^= a & x; }

	BitVector::reference axpyin (BitVector::reference r, Element a, Element x) const
		{ return r ^= a & x; }

	std::_Bit_reference axpyin (std::_Bit_reference r, Element a, Element x) const
		{ return r = r ^ (a & x); }

	Element &axpyin (Element &r, const std::_Bit_reference a, Element x) const
		{ return r ^= a & x; }

	std::_Bit_reference axpyin (std::_Bit_reference r, const std::_Bit_reference a, Element x) const
		{ return r = r ^ (a & x); }

	Element &axpyin (Element &r, Element a, const std::_Bit_reference x) const
		{ return r ^= a & static_cast<bool>(x); }

	std::_Bit_reference axpyin (std::_Bit_reference r, Element a, const std::_Bit_reference x) const
		{ return r = r ^ (a & static_cast<bool>(x)); }

	Element &axpyin (Element &r, const std::_Bit_reference a, const std::_Bit_reference x) const
		{ return r ^= a & static_cast<bool>(x); }

	std::_Bit_reference axpyin (std::_Bit_reference r, const std::_Bit_reference a, const std::_Bit_reference x) const
		{ return r = r ^ (a & static_cast<bool>(x)); }

	//@} Inplace Arithmetic Operations

	static inline int getMaxModulus() { return 2; }

}; // class GF2

} // namespace LinBox


// Specialization of GivaroField for GF2
#include "linbox/field/givaro-field.h"
namespace LinBox 
{

  /** 
  \brief give LinBox fields an allure of Givaro Fields
  \ingroup field

   *  This class adds the necessary requirements allowing 
   *  the construction of an extension of a LinBox field.
   */ 
    template<>
    struct GivaroField<LinBox::GF2> : public LinBox::GF2
    {
        typedef LinBox::GF2 BaseField;
        typedef BaseField::Element TT;
        typedef Signed_Trait<TT>::unsigned_type UTT;
        typedef TT Rep;
        typedef GivaroField<BaseField> Self_t;
        typedef Rep Element;
        typedef UTT Residu_t;

        Element zero, one;
        GivaroField(const BaseField& bf) : BaseField(bf) {
            this->init(zero,0UL);
            this->init(one, 1UL);
        }


            // -- amxy: r <- c - a * b mod p
        Rep& amxy (Rep& r, const Rep a, const Rep b, const Rep c) const {
            Rep tmp;
            this->mul(tmp, a, b);
            this->assign(r,c);
            return this->subin(r,tmp);
        }
        std::_Bit_reference amxy (std::_Bit_reference r, const Rep a, const Rep b, const Rep c) const {
            Rep tmp;
            this->mul(tmp, a, b);
            this->assign(r,c);
            return this->subin(r,tmp);
        }


            // -- maxpy: r <- y - a * x
        Rep& maxpy (Rep& r, const Rep a, const Rep x, const Rep y) const {
            Rep tmp; this->mul(tmp, a, x);
            return this->sub(r,y,tmp);
        }
        std::_Bit_reference maxpy (std::_Bit_reference r, const Rep a, const Rep x, const Rep y) const {
            Rep tmp; this->mul(tmp, a, x);
            return this->sub(r,y,tmp);
        }
            // -- axmyin: r <- r - a * x 
        Rep& axmyin (Rep& r, const Rep a, const Rep x) const {
            Rep tmp; this->mul(tmp, a, x);
            return this->subin(r,tmp);
        }
        std::_Bit_reference axmyin (std::_Bit_reference r, const Rep a, const Rep x) const {
            Rep tmp; this->mul(tmp, a, x);
            return this->subin(r,tmp);
        }
            // -- maxpyin: r <- r - a * x
        Rep& maxpyin (Rep& r, const Rep a, const Rep x) const {
	    return axmyin(r,a,x);
        }
        std::_Bit_reference maxpyin (std::_Bit_reference r, const Rep a, const Rep x) const {
	    return axmyin(r,a,x);
        }

 

        bool areNEqual ( const Rep a, const Rep b) const {
            return ! this->areEqual(a,b);
        }

            // Access to the modulus, characteristic, size, exponent
        UTT residu() const { integer c; BaseField::characteristic(c); return UTT(c); }
        UTT characteristic() const  { integer c; BaseField::characteristic(c); return UTT(c); }
        UTT cardinality() const  { integer c; BaseField::cardinality(c); return UTT(c); }
        UTT exponent() const { return 1; }
        UTT size() const  { integer c; BaseField::cardinality(c); return UTT(c); }


            // ----- random generators
        template<class RandIter> Rep& random(RandIter& g, Rep& r) const { return r = g() ; }
        template<class RandIter> Rep& random(RandIter& g, Rep& r, long s) const { return r = g() ; }
        template<class RandIter> Rep& random(RandIter& g, Rep& r, const Rep& b) const { return r = g() ; }
        template<class RandIter> Rep& nonzerorandom(RandIter& g, Rep& r) const { return r = g() ; }
        template<class RandIter> Rep& nonzerorandom(RandIter& g, Rep& r, long s) const { return r = g() ; }
        template<class RandIter> Rep& nonzerorandom(RandIter& g, Rep& r, const Rep& b) const { return r = g() ; }

        template<class RandIter> std::_Bit_reference random(RandIter& g, std::_Bit_reference r) const { return r = g() ; }
        template<class RandIter> std::_Bit_reference random(RandIter& g, std::_Bit_reference r, long s) const { return r = g() ; }
        template<class RandIter> std::_Bit_reference random(RandIter& g, std::_Bit_reference r, const std::_Bit_reference b) const { return r = g() ; }
        template<class RandIter> std::_Bit_reference nonzerorandom(RandIter& g, std::_Bit_reference r) const { return r = g() ; }
        template<class RandIter> std::_Bit_reference nonzerorandom(RandIter& g, std::_Bit_reference r, long s) const { return r = g() ; }
        template<class RandIter> std::_Bit_reference nonzerorandom(RandIter& g, std::_Bit_reference r, const Rep& b) const { return r = g() ; }
        template<class RandIter> std::_Bit_reference nonzerorandom(RandIter& g, std::_Bit_reference r, const std::_Bit_reference b) const { return r = g() ; }

    };
    
}



// Specialization of homomorphism for basefield
#include "linbox/field/hom.h"
#include "linbox/field/givaro-extension.h"
namespace LinBox 
{

    template <>
    class Hom<GF2,GF2> {
        
    public:
        typedef GF2 Target;
        typedef GF2 Source;
        typedef Source::Element SrcElt;
        typedef Target::Element Elt;
	
        Hom(const Source& S, const Target& ) : _source (S){}
        Elt& image(Elt& t, const SrcElt& s) {
            return _source.assign (t, s);
        }
        SrcElt& preimage(SrcElt& s, const Elt& t) {
            return _source.assign (s, t);
        }
        const Source& source() { return _source;}
        const Target& target() { return _source;}
        
    protected:
        Source _source;
    };

    template<class Target > 
    class Hom<GF2, Target > 
    {   public:
        typedef GF2 Source;
        typedef typename GF2::Element SrcElt;
        typedef typename Target::Element Elt;

        Hom(const Source& S, const Target& T) : _source(S), _target(T){ }
        Elt& image(Elt& t, const SrcElt& s) {
            return _source.convert(t,s);
        }
        SrcElt& preimage(SrcElt& s, const Elt& t) {
            return _target.convert(s,t);
        }
        std::_Bit_reference preimage(std::_Bit_reference s, const Elt& t) const {
            int ts;
            return s = _target.convert(ts, t);
        }

        const Source& source() { return _source;}
        const Target& target() { return _target;}

    private:
        Source _source;
        Target _target;
    }; // end Hom 



    template<>
    class Hom < GF2, GivaroExtension<GF2> >
    {
        typedef GF2 Source;
        typedef GivaroExtension<GF2> Target;
    public:
        typedef Source::Element SrcElt;
        typedef Target::Element Elt;

            //Hom(){}
            /**
             * Construct a homomorphism from a specific source ring S and target 
             * field T with Hom(S, T).  The default behaviour is error.  
             * Specializations define all actual homomorphisms.
             */
        Hom(const Source& S, const Target& T) : _source(S), _target(T){}

            /** 
             * image(t, s) implements the homomorphism, assigning the 
             * t the value of the image of s under the mapping.
             *
             * The default behaviour is a no-op.
             */
        Elt& image(Elt& t, const SrcElt& s) const {return _target.assign(t,s);}
           
            /** If possible, preimage(s,t) assigns a value to s such that 
             * the image of s is t.  Otherwise behaviour is unspecified.
             * An error may be thrown, a conventional value may be set, or
             * an arb value set.
             *
             * The default behaviour is a no-op.
             */
        SrcElt& preimage(SrcElt& s, const Elt& t) const {
            return _target.convert(s, t);
        }
        std::_Bit_reference preimage(std::_Bit_reference s, const Elt& t) const {
            bool ts;
            return s = _target.convert(ts, t);
        }

        const Source& source() const { return _source;}
        const Target& target() const { return _target;}

    private:
        Source _source;
        Target _target;
    }; // end Hom 
}

// #include <bits/stl_bvector.h>
namespace std 
{
//! @todo JGD 05.11.2009 : it should be in bits/stl_bvector.h  ...
    inline void swap(_Bit_reference __x, _Bit_reference __y)
    {
      bool __tmp = __x;
      __x = __y;
      __y = __tmp;
    }
}


#include "linbox/randiter/gf2.h"
#include "linbox/field/gf2.inl"

#endif // __LINBOX_field_gf2_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
