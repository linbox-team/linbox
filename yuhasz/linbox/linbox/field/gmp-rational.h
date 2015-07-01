/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/field/gmp-rational.h
 * Copyright (C) 2001, 2002 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * ----------------------------------------
 *
 * See COPYING for license information
 */

#ifndef __FIELD_GMP_RATIONAL_H
#define __FIELD_GMP_RATIONAL_H

#include <iostream>
#include <cctype>

#include <gmp.h>

#include "linbox/integer.h"
#include "linbox/field/field-interface.h"
#include "linbox/element/gmp-rational.h"
#include "linbox-config.h"

#ifdef __LINBOX_XMLENABLED

#include "linbox/util/xml/linbox-reader.h"
#include "linbox/util/xml/linbox-writer.h"

using LinBox::Reader;
using LinBox::Writer;

#include <string>

using std::istream;
using std::ostream;
using std::string;

#endif

// Namespace in which all LinBox library code resides
namespace LinBox
{

// Forward declarations
class GMPRationalRandIter;;

/** @name GMP Rational field
 * @memo Field of rational numbers using GMP
 *
 * @doc
 * This is a wrapper for the GMP rational number facility, built to the
 * interface of the field archetype. 
 */
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
	 * These methods are required of all \Ref{LinBox} fields.
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

	GMPRationalField (const GMPRationalField &F) 
		: _cardinality (0), _characteristic (0), _zero (0), _one (1), _neg_one (-1),
		  zero (_zero, _one), one (_one, _one), neg_one (_neg_one, _one)
	{}

#ifdef __LINBOX_XMLENABLED
	// XML Reader constructor
	GMPRationalField(Reader &R) : _cardinality(0), _characteristic(0), _zero(0), _one(1), _neg_one(-1)
	{
		if(!R.expectTagName("field") || !R.expectChildTag()) return;
		R.traverseChild();
		if(!R.expectTagName("rational")) return;

		R.upToParent();
	}
#endif

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
	GMPRationalField &operator= (const GMPRationalField &F)
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
		mpq_set_si (x.rep, (signed long) y, 1L);
		mpq_canonicalize (x.rep);
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

		mpz_div (n, n, d);

		x = integer::zero;

		// Really bad, but I know of no other general way to do this
		while (mpz_sgn (n) != 0) {
			// We need to be ready for multiple word sizes and so on here...
			x = (x << (sizeof (unsigned long) << 3)) + mpz_get_ui (n);
			mpz_tdiv_q_2exp (n, n, sizeof (unsigned long) << 3);
		}

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
    
	/** Cardinality.
	 * Return integer representing cardinality of the field.
	 * Returns a non-negative integer for all fields with finite
	 * cardinality, and returns -1 to signify a field of infinite 
	 * cardinality.
	 * @return constant reference to integer representing cardinality 
	 *	       of the field
	 */
	integer &cardinality (integer &c) const 
	{ return c = _cardinality; }

	/** Characteristic.
	 * Return integer representing characteristic of the field.
	 * Returns a positive integer to all fields with finite characteristic,
	 * and returns 0 to signify a field of infinite characteristic.
	 * @return constant reference to integer representing characteristic 
	 * 	       of the field.
	 */
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

#ifndef __LINBOX_XMLENABLED
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
		delete str;

		if (mpz_cmp_ui (mpq_denref (x.rep), 1L) != 0) {
			str = new char[mpz_sizeinbase (mpq_denref (x.rep), 10) + 2];
			mpz_get_str (str, 10, mpq_denref (x.rep));
			os << '/' << str;
			delete str;
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
	 */
	std::istream &read (std::istream &is, Element &x) const
	{
		char buffer[65535], endc;
		bool found_space = false;
		int i = 0;

		do {
			is.get (endc);
		} while (is && !isdigit (endc) && endc != '-' && endc != '.');


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
		}
		else if (endc == '.' && !found_space) {
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
		}
		else {
			is.putback (endc);
			mpz_set_si (mpq_denref (x.rep), 1L);
		}

		mpq_canonicalize (x.rep);

		return is;
	}

	//@} Input/Output Operations
    
#else

	// XML field writer
	ostream &GMPRationalField::write(ostream &os) const
	{
		Writer W;
		if( toTag(W) )
			W.write(os);

		return os;
	}

	// XML toTag function
	bool GMPRationalField::toTag(Writer &W) const
	{
		W.setTagName("field");
		W.setAttribute("implDetail", "gmp-rational");
		W.setAttribute("cardinality", "0");

		W.addTagChild();
		W.setTagName("rational");
		W.upToParent();

		return true;
	}

	// XML element writer method
	ostream &GMPRationalField::write(ostream &os, const Element &e) const 
	{
		Writer W;
		if( toTag(W, e))
			W.write(os);

		return os;
	}

	// This code may look familiar, since I lifted the heart of it
	// from Bradford's write method above (sorry Bradford :-))
	//
	bool GMPRationalField::toTag(Writer &W, const Element &e) const 
	{
		string s;
		char *str;

		W.setTagName("cn");

		str = new char[mpz_sizeinbase (mpq_numref (e.rep), 10) + 2];
		mpz_get_str (str, 10, mpq_numref (e.rep));
		s = str;
		delete [] str;

		if (mpz_cmp_ui (mpq_denref (e.rep), 1L) != 0) {
			str = new char[mpz_sizeinbase (mpq_denref (e.rep), 10) + 2];
			mpz_get_str (str, 10, mpq_denref (e.rep));
			s += '/';
			s += str;
			delete [] str;
			W.addDataChild(s);
		}
		else {
			W.addDataChild(s);
		}

		return true;
	}

	// XML element Read function
	istream &GMPRationalField::read(istream &is, Element &e) const 
	{
		Reader R(is);
		if(!fromTag(R, e)) {
			is.setstate(istream::failbit);
			if(!R.initalized())
				is.setstate(istream::badbit);
		}
		
		return is;
	}

	// (real) XML element read method, initalizes the element using the
	// Reader.  As I've stolen a bunch of code from Bradford, so as above,
	// it is assumed that this field element has already been constructed
	// when it is put to this function
	//
	bool GMPRationalField::fromTag(Reader &R, Element &e) const
	{
		string s, temp;
		size_t i;
		if(!R.expectTagName("cn") || !R.expectChildTextString(s)) return false;
		// first get the location of the /, if there is one
		i = s.find_first_of("/");

		if(i == string::npos) { // there isn't a denominator
			if(!R.isNum(s)) {
				R.setErrorString("Tried to convert the string \"" + s + "\" to a GMP Rational type, but couldn't.");
				R.setErrorCode(Reader::OTHER);
				return false;
			}
			mpz_set_str (mpq_numref (e.rep), s.c_str(), 10);
			mpz_set_str (mpq_denref (e.rep), "1", 10);
		}
		else {
			if(!R.isNum(s.substr(0, i))) {
				R.setErrorString("Tried to convert the string \"" + s + "\" to a GMP Rational type, but couldn't.");
				R.setErrorCode(Reader::OTHER);
				return false;
			}
			mpz_set_str(mpq_numref (e.rep), s.substr(0, i).c_str(), 10);

			if(!R.isNum(s.substr(i + 1, s.length() - i))) {
				R.setErrorString("Tried to convert the string \"" + s + "\" to a GMP Rational type, but couldn't.");
				R.setErrorCode(Reader::OTHER);
				return false;
			}
			mpz_set_str(mpq_denref (e.rep), s.substr(i + 1, s.length() - i).c_str(), s.length() - (2*i +1));
		}

		return true;
	}

#endif


	//@} Common Object Interface

	GMPRationalField ()
		: _cardinality (0), _characteristic (0), _zero (0), _one (1), _neg_one (-1),
		  zero (_zero, _one), one (_one, _one), neg_one (_neg_one, _one)
	{}
    
}; // class GMPRationalField

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

} // namespace LinBox

#include "linbox/randiter/gmp-rational.h"

#endif // __FIELD_GMP_RATIONAL_H
