/* linbox/field/gmp-rational.h
 * Copyright (C) 2001, 2002 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * ----------------------------------------
 *
 * See COPYING for license information
 */

#ifndef __LINBOX_field_gmp_rational_H
#define __LINBOX_field_gmp_rational_H

#include <iostream>
#include <cctype>

#include <gmp.h>

#include "linbox/integer.h"
#include "linbox/field/field-interface.h"
#include "linbox/element/gmp-rational.h"
#include "linbox/linbox-config.h"
#include "linbox/util/debug.h"
#include <linbox/field/field-traits.h>

// Namespace in which all LinBox library code resides
namespace LinBox
{

// Forward declarations
class GMPRationalRandIter;;

/**
 * \brief Field of rational numbers using GMP
 \ingroup field
 *
 * This is a wrapper for the GMP rational number facility, built to the
 * interface of the field archetype. 
 */

template <class Ring>
struct ClassifyRing;

class GMPRationalField;

template<>
struct ClassifyRing<GMPRationalField> {
	typedef RingCategories::RationalTag categoryTag;
};

class GMPRationalField : public FieldInterface
{
    private:

	const integer _cardinality;
	const integer _characteristic;

	const integer _zero;
	const integer _one;
	const integer _neg_one;

    public:

	/** @name Common Object Interface for a LinBox Field.
	 * These methods are required of all \ref{LinBox} fields.
	 */
	//@{
    
	/// element type.
	typedef GMPRationalElement Element;

	/// Random iterator generator type.
	typedef GMPRationalRandIter RandIter;
    
	const Element zero;
	const Element one;
	const Element neg_one;

	/** @name Object Management
	 * x <- convert (y)
	 */
	//@{
    
	/** Copy constructor.
	 *
	 * Vacuous, since this field is unparametric so there is no need to
	 * construct multiple field objects
	 */

	GMPRationalField (const GMPRationalField &) 
		: _cardinality (0), _characteristic (0), _zero (0), _one (1), _neg_one (-1),
		  zero (_zero, _one), one (_one, _one), neg_one (_neg_one, _one)
	{}

	/** Destructor.
	 * 
	 * Also vacuous, since there is no de-initialization system
	 */
	~GMPRationalField (void) 
	{}
    
	/** Assignment operator.
	 * 
	 * Also vacuous
	 */
	GMPRationalField &operator= (const GMPRationalField &)
	{ return *this; }
    
	/** Initialization of field element from an integer.
	 * Behaves like C++ allocator construct.
	 * This function assumes the output field element x has already been 
	 * constructed, but that it is not necessarily already initialized.
	 * In this implementation, this means the _elem_ptr of x exists, but
	 * that it may be the null pointer.
	 * @return reference to field element.
	 * @param x field element to contain output (reference returned).
	 * @param y constant reference to integer.
	 */
	Element &init (Element &x, const integer &y = 0) const
	{
		mpq_set_z (x. rep, SpyInteger::get_mpz(const_cast<integer&>(y)));
		//mpq_set_si (x.rep, (signed long) y, 1L);
		//mpq_canonicalize (x.rep);
		return x;
	}

	/*
	 * aniau@astronet.pl: 06/2009 Initialization of field element from numerator and denominator
	 * */
	Element &init (Element &x, const integer &num, const integer &den) const
	{
		init(x,num);
		Element y;
	        init(y,den);
	        divin(x,y);
	        return x;
	}
  
	/** Conversion of field element to an integer.
	 * This function assumes the output field element x has already been 
	 * constructed, but that it is not already initialized.
	 * In this implementation, this means the _elem_ptr of y exists, and
	 * that it is not the null pointer.
	 *
	 * Returns floor (numerator (y) / denominator (y))
	 *
	 * @return reference to integer.
	 * @param x reference to integer to contain output (reference returned).
	 * @param y constant reference to field element.
	 */
	integer &convert (integer &x, const Element &y = 0) const
	{
		mpz_t n, d;

		mpz_init (n);
		mpz_init (d);
		mpq_get_num (n, y.rep);
		mpq_get_den (d, y.rep);

		mpz_divexact (x.get_mpz(), n, d);

		/* Shouldn't there be something like this? We'll assume integer is gmp integers.
		x.set_mpz(n);
		*/
// 		x = integer::zero;

// 		// Really bad, but I know of no other general way to do this
// 		while (mpz_sgn (n) != 0) {
// 			// We need to be ready for multiple word sizes and so on here...
// 			x = (x << (sizeof (unsigned long) << 3)) + mpz_get_ui (n);
// 			mpz_tdiv_q_2exp (n, n, sizeof (unsigned long) << 3);
// 		}

		return x;
	}
    
	/** Assignment of one field element to another.
	 * This function assumes both field elements have already been 
	 * constructed and initialized.
	 * In this implementation, this means for both x and y, 
	 * _elem_ptr exists and does not point to null.
	 * @return reference to x
	 * @param  x field element (reference returned).
	 * @param  y field element.
	 *
	 * FIXME: Is this x := y? I am assuming so.
	 */
	Element &assign (Element &x, const Element &y) const
	{
		mpq_set (x.rep, y.rep);
		return x;
	}
    
	/// is infinite (signified by -1 here)
	integer &cardinality (integer &c) const 
	{ return c = _cardinality; }

	/// of the rationals is 0.
	integer &characteristic (integer &c) const
	{ return c = _characteristic; }

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
	 * In this implementation, this means for both x and y, 
	 * _elem_ptr exists and does not point to null.
	 * @return boolean true if equal, false if not.
	 * @param  x field element
	 * @param  y field element
	 */
	bool areEqual (const Element &x, const Element &y) const
	{ return mpq_equal (x.rep, y.rep); }

	/** Addition.
	 * x = y + z
	 * This function assumes all the field elements have already been 
	 * constructed and initialized.
	 * In this implementation, this means for x, y, and z, 
	 * _elem_ptr exists and does not point to null.
	 * @return reference to x.
	 * @param  x field element (reference returned).
	 * @param  y field element.
	 * @param  z field element.
	 */
	Element &add (Element &x, const Element &y, const Element &z) const
	{
		mpq_add (x.rep, y.rep, z.rep);
		return x;
	}
    
	/** Subtraction.
	 * x = y - z
	 * This function assumes all the field elements have already been 
	 * constructed and initialized.
	 * In this implementation, this means for x, y, and z, 
	 * _elem_ptr exists and does not point to null.
	 * @return reference to x.
	 * @param  x field element (reference returned).
	 * @param  y field element.
	 * @param  z field element.
	 */
	Element &sub (Element &x, const Element &y, const Element &z) const
	{
		mpq_sub (x.rep, y.rep, z.rep);
		return x;
	}
    
	/** Multiplication.
	 * x = y * z
	 * This function assumes all the field elements have already been 
	 * constructed and initialized.
	 * In this implementation, this means for x, y, and z, 
	 * _elem_ptr exists and does not point to null.
	 * @return reference to x.
	 * @param  x field element (reference returned).
	 * @param  y field element.
	 * @param  z field element.
	 */
	Element &mul (Element &x, const Element &y, const Element &z) const
	{
		mpq_mul (x.rep, y.rep, z.rep);
		return x;
	}
    
	/** Division.
	 * x = y / z
	 * This function assumes all the field elements have already been 
	 * constructed and initialized.
	 * In this implementation, this means for x, y, and z, 
	 * _elem_ptr exists and does not point to null.
	 * @return reference to x.
	 * @param  x field element (reference returned).
	 * @param  y field element.
	 * @param  z field element.
	 */
	Element &div (Element &x, const Element &y, const Element &z) const
	{
		mpq_div (x.rep, y.rep, z.rep);
		return x;
	}

	/** Additive Inverse (Negation).
	 * x = - y
	 * This function assumes both field elements have already been 
	 * constructed and initialized.
	 * In this implementation, this means for both x and y 
	 * _elem_ptr exists and does not point to null.
	 * @return reference to x.
	 * @param  x field element (reference returned).
	 * @param  y field element.
	 */
	Element &neg (Element &x, const Element &y) const
	{
		mpq_neg (x.rep, y.rep);
		return x;
	}

	/** Multiplicative Inverse.
	 * x = 1 / y
	 * This function assumes both field elements have already been 
	 * constructed and initialized.
	 * In this implementation, this means for both x and y 
	 * _elem_ptr exists and does not point to null.
	 * @return reference to x.
	 * @param  x field element (reference returned).
	 * @param  y field element.
	 */
	Element &inv (Element &x, const Element &y) const
	{
		mpq_inv (x.rep, y.rep);
		return x;
	}

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
	 * In this implementation, this means the _elem_ptr of x
	 * exists and does not point to null.
	 * @return boolean true if equals zero, false if not.
	 * @param  x field element.
	 */
	bool isZero (const Element &x) const 
	{ return mpq_sgn (x.rep) == 0; }
    
	/** One equality.
	 * Test if field element is equal to one.
	 * This function assumes the field element has already been 
	 * constructed and initialized.
	 * In this implementation, this means the _elem_ptr of x
	 * exists and does not point to null.
	 * @return boolean true if equals one, false if not.
	 * @param  x field element.
	 */
	bool isOne (const Element &x) const 
	{ return mpq_cmp_ui (x.rep, 1L, 1L) == 0; }
    
	/** Inplace Addition.
	 * x += y
	 * This function assumes both field elements have already been 
	 * constructed and initialized.
	 * In this implementation, this means for both x and y 
	 * _elem_ptr exists and does not point to null.
	 * @return reference to x.
	 * @param  x field element (reference returned).
	 * @param  y field element.
	 */
	Element &addin (Element &x, const Element &y) const
	{
		mpq_add (x.rep, x.rep, y.rep);
		return x;
	}

	/** Inplace Subtraction.
	 * x -= y
	 * This function assumes both field elements have already been 
	 * constructed and initialized.
	 * In this implementation, this means for both x and y 
	 * _elem_ptr exists and does not point to null.
	 * @return reference to x.
	 * @param  x field element (reference returned).
	 * @param  y field element.
	 */
	Element &subin (Element &x, const Element &y) const
	{
		mpq_sub (x.rep, x.rep, y.rep);
		return x;
	}
 
	/** Inplace Multiplication.
	 * x *= y
	 * This function assumes both field elements have already been 
	 * constructed and initialized.
	 * In this implementation, this means for both x and y 
	 * _elem_ptr exists and does not point to null.
	 * @return reference to x.
	 * @param  x field element (reference returned).
	 * @param  y field element.
	 */
	Element &mulin (Element &x, const Element &y) const
	{
		mpq_mul (x.rep, x.rep, y.rep);
		return x;
	}

	Element &axpy (Element &r, const Element &a, const Element &x, const Element &y) const
	{
		mpq_mul (r.rep, a.rep, x.rep);
		mpq_add (r.rep, r.rep, y.rep);
		return r;
	}

	Element &axpyin (Element &r, const Element &a, const Element &x) const
	{
		Element tmp;
		mpq_mul (tmp.rep, a.rep, x.rep);
		mpq_add (r.rep, r.rep, tmp.rep);
		return r;
	}

	/** Inplace Division.
	 * x /= y
	 * This function assumes both field elements have already been 
	 * constructed and initialized.
	 * In this implementation, this means for both x and y 
	 * _elem_ptr exists and does not point to null.
	 * @return reference to x.
	 * @param  x field element (reference returned).
	 * @param  y field element.
	 */
	Element &divin (Element &x, const Element &y) const
	{
		mpq_div (x.rep, x.rep, y.rep);
		return x;
	}
    
	/** Inplace Additive Inverse (Inplace Negation).
	 * x = - x
	 * This function assumes the field element has already been 
	 * constructed and initialized.
	 * In this implementation, this means the _elem_ptr of x
	 * exists and does not point to null.
	 * @return reference to x.
	 * @param  x field element (reference returned).
	 */
	Element &negin (Element &x) const
	{
		mpq_neg (x.rep, x.rep);
		return x;
	}

	/** Inplace Multiplicative Inverse.
	 * x = 1 / x
	 * This function assumes the field elementhas already been 
	 * constructed and initialized.
	 * In this implementation, this means the _elem_ptr of x
	 * exists and does not point to null.
	 * @return reference to x.
	 * @param  x field element (reference returned).
	 */
	Element &invin (Element &x) const
	{
		mpq_inv (x.rep, x.rep);
		return x;
	}
    
	//@} Inplace Arithmetic Operations

	/** @name Input/Output Operations */
	//@{
    
	/** Print field.
	 * @return output stream to which field is written.
	 * @param  os  output stream to which field is written.
	 *
	 * This does not do much...
	 */
	std::ostream &write (std::ostream &os) const 
	{ 
		os << "GMP rational numbers"; 
		return os;
	}
    
	/** Read field.
	 * @return input stream from which field is read.
	 * @param  is  input stream from which field is read.
	 *
	 * This does not do much either...
	 *
	 * FIXME: Read the same thing written above, and throw an exception if the
	 * strings do not match.
	 */
	std::istream &read (std::istream &is) { return is; }
    
	/** Print field element.
	 * This function assumes the field element has already been 
	 * constructed and initialized.
	 * In this implementation, this means for the _elem_ptr for x 
	 * exists and does not point to null.
	 * @return output stream to which field element is written.
	 * @param  os  output stream to which field element is written.
	 * @param  x   field element.
	 */
	std::ostream &write (std::ostream &os, const Element &x) const 
	{
		char *str;

		str = new char[mpz_sizeinbase (mpq_numref (x.rep), 10) + 2];
		mpz_get_str (str, 10, mpq_numref (x.rep));
		os << str;
		delete[] str;

		if (mpz_cmp_ui (mpq_denref (x.rep), 1L) != 0) {
			str = new char[mpz_sizeinbase (mpq_denref (x.rep), 10) + 2];
			mpz_get_str (str, 10, mpq_denref (x.rep));
			os << '/' << str;
			delete[] str;
		}

		return os;
	}

	/** Read field element.
	 * This function assumes the field element has already been 
	 * constructed and initialized.
	 * In this implementation, this means for the _elem_ptr for x 
	 * exists and does not point to null.
	 * @return input stream from which field element is read.
	 * @param  is  input stream from which field element is read.
	 * @param  x   field element.
	 *
	 * FIXME: Avoid the magical limit on size here
	 * FIXME: Right now it skips over everything until it finds something that
	 * looks like a number. Is this really the correct policy?
	 *
	 * aniau@astronet.pl: 06/2009: supports scientific E/e notation for decimal fractions
	 */
	std::istream &read (std::istream &is, Element &x) const
	{
		char buffer[65535], endc;
		bool found_space = false;
		int i = 0;			

		do {
			is.get (endc);
		} while (is && !isdigit (endc) && endc != '-' && endc != '.' &&  endc !='e' && endc != 'E');		

		buffer[i]=endc;
		
		while ((buffer[i] == '-' || isdigit (buffer[i])) && i < 65535)  {
			i++;
			is.get (buffer[i]);
		}

		endc = buffer[i];
	       	buffer[i] = '\0';

		if (i > 0)
			mpz_set_str (mpq_numref (x.rep), buffer, 10);
		else
			mpq_set_si (x.rep, 0L, 1L);

		if (endc == ' ') {
			found_space = true;
			while (endc == ' ') is >> endc;
		}

		if (endc == '/') {
			i = 0;

			is.get (endc);
			while (isspace (endc)) is.get (endc);
			is.putback (endc);

			do {
				is.get (buffer[i++]);
			} while (isdigit (buffer[i - 1]) && i < 65536);

			is.putback (buffer[i - 1]);
			buffer[i - 1] = '\0';

			mpz_set_str (mpq_denref (x.rep), buffer, 10);
			mpq_canonicalize (x.rep);
			return is;
		} else {
			 mpz_set_si (mpq_denref (x.rep), 1L);
		}

		if (endc == '.' && !found_space) {
			Element decimal_part;

			mpz_set_si (mpq_denref (x.rep), 1L);
			mpq_set_si (decimal_part.rep, 1L, 1L);
			mpz_set_si (mpq_denref (decimal_part.rep), 1L);

			i = 0;

			do {
				is.get (buffer[i++]);
				if (isdigit (buffer[i - 1]))
					mpz_mul_ui (mpq_denref (decimal_part.rep),
						    mpq_denref (decimal_part.rep), 10L);
			} while (isdigit (buffer[i - 1]) && i < 65536);

			is.putback (buffer[i - 1]);
			buffer[i - 1] = '\0';

			mpz_set_str (mpq_numref (decimal_part.rep), buffer, 10);
			mpq_canonicalize (decimal_part.rep);

			mpq_add (x.rep, x.rep, decimal_part.rep);

			do {
				is.get (endc);
			} while (is && endc == ' ') ;
		}
		
		if ((endc == 'e') || (endc == 'E')) {
	                is.get(endc);
	                bool minus = false;
                        if (endc == '-') {
                                minus = true;
                        } else if (endc == '+') {
                                minus = false;
                        } else {
	                        is.putback(endc);
                        }

	                i=0;
	                do {
		                is.get (buffer[i++]);
		        } while (isdigit (buffer[i-1]) && i < 65536);
		        is.putback(buffer[i-1]);
		        buffer[i-1] = '\0';

			integer pow(buffer), powten=1;

			for (integer it=0; it< pow; ++it) powten *=10;
                        if (minus) {
	                        div(x,x,powten);
	                } else {
	                        mul(x,x,powten);
	                }
                }
		else {
			is.putback (endc);
		//	mpz_set_si (mpq_denref (x.rep), 1L);
		}

		mpq_canonicalize (x.rep);

		return is;
	}

	//@} Input/Output Operations
    

	//@} Common Object Interface

	GMPRationalField (int p = 0, int exp = 1)
		: _cardinality (0), _characteristic (0), _zero (0), _one (1), _neg_one (-1),
		  zero (_zero, _one), one (_one, _one), neg_one (_neg_one, _one)
	{
		if(p != 0) throw PreconditionFailed(__FUNCTION__,__LINE__,"modulus must be 0 (no modulus)");
		if(exp != 1) throw PreconditionFailed(__FUNCTION__,__LINE__,"exponent must be 1");
	}
    
	static inline int getMaxModulus() { return 0; } // no modulus
	
	// x = numerator of y
	integer& get_num (integer& x, const Element& y)  const{
		mpq_get_num (SpyInteger::get_mpz(x), y. rep);
		return x;

	}

	// x = denominator of y
	integer& get_den (integer& x, const Element& y) const {
		mpq_get_den (SpyInteger::get_mpz(x), y. rep);
		return x;
	}

	int sign (const Element& x) const {
		return mpq_sgn (x. rep);
	}
	
	//bitsize as sum of bitsizes of numerator and denominator
        integer& bitsize(integer& bs, const Element q) const {
                integer y; get_den(y,q);
                integer x; get_num(x,q);
                bs = x.bitsize() + y.bitsize();
                return bs;
        }


}; // class GMPRationalField

/*  use GMPRationalField::read() and GMPRationalField::write(), not these operators.
std::ostream &operator << (std::ostream &os, GMPRationalElement &elt)
{
	GMPRationalField field;

	field.write (os, elt);
	return os;
}

std::istream &operator >> (std::istream &is, GMPRationalElement &elt)
{
	GMPRationalField field;

	field.read (is, elt);
	return is;
}
*/

} // namespace LinBox

#include "linbox/randiter/gmp-rational.h"

#endif // __LINBOX_field_gmp_rational_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
