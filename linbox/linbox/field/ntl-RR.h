/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/field/ntl-RR.h
 * Copyright (C) 1999-2002 William J Turner,
 *               2001 Bradford Hovinen
 *
 * Written by William J Turner <wjturner@math.ncsu.edu>,
 *            Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __FIELD_NTL_RR_H
#define __FIELD_NTL_RR_H
#include <NTL/tools.h>

using namespace NTL;
#include <NTL/RR.h>

#include "linbox/field/unparametric.h"
#include "linbox-config.h"

#ifdef __LINBOX_XMLENABLED

#include "linbox/util/xml/linbox-reader.h"
#include "linbox/util/xml/linbox-writer.h"

using LinBox::Reader;
using LinBox::Writer;

#include <iostream>
#include <string>

using std::istream;
using std::ostream;
using std::string;

#endif

// Namespace in which all LinBox library code resides
namespace LinBox
{
  
	///
	//class NTL_RR: public UnparametricField<NTL::RR>, public FieldInterface{};
	//typedef UnparametricField<NTL::RR> NTL_RR;

	/** @name class RR.
	 * Arbitrary precision floating point numbers.
	 * These specializations allow the \Ref{UnparametricField} template class to be
	 * used to wrap NTL's RR class as a LinBox field.
	 */
	//@{

#ifdef __LINBOX_XMLENABLED
	
	// This function is the XML Reader constructor.  Take in a Reader
	// and attempt to "initalize" the field.  This really does nothing
	// except check that the field matches the XML given
	//
	template <> UnparametricField<NTL::RR>::UnparametricField(Reader &R)
	{
		if(!R.expectTagName("field") || !R.expectChildTag()) return;
		R.traverseChild();
		if(!R.expectTagName("real")) return;
		R.upToParent();

		// set the proper characteristic & cardinality
		_p = Integer(0);
		_card = Integer(-1);

		return;
	}
#endif



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
		NTL::RR& UnparametricField<NTL::RR>::init(NTL::RR& x, const integer& y) const
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
		integer& UnparametricField<NTL::RR>::convert(integer& x, const NTL::RR& y) const
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
		NTL::RR& UnparametricField<NTL::RR>::inv(NTL::RR& x, const NTL::RR& y) const
		{ return x = NTL::inv(y); }
 
	/** Zero equality.
	 * Test if field element is equal to zero.
	 * This function assumes the field element has already been
	 * constructed and initialized.
	 * In this specialization, NTL's IsZero function is called.
	 * @return boolean true if equals zero, false if not.
	 * @param  x field element.
	 */
	template <> bool UnparametricField<NTL::RR>::isZero(const NTL::RR& x) const
		{ return static_cast<bool>(IsZero(x)); }

	/** One equality.
	 * Test if field element is equal to one.
	 * This function assumes the field element has already been
	 * constructed and initialized.
	 * In this specialization, NTL's IsOne function is called.
	 * @return boolean true if equals one, false if not.
	 * @param  x field element.
	 */
	template <> bool UnparametricField<NTL::RR>::isOne(const NTL::RR& x) const
		{ return static_cast<bool>(IsOne(x)); }

	/** Inplace Multiplicative Inverse.
	 * x = 1 / x
	 * This function assumes both field elements have already been
	 * constructed and initialized.
	 * @return reference to x.
	 * @param  x field element (reference returned).
	 */
	template <> NTL::RR& UnparametricField<NTL::RR>::invin(NTL::RR& x) const
		{ return x = NTL::inv(x); }

#ifndef __LINBOX_XMLENABLED // <- old writer

	/** Print field.
	 * @return output stream to which field is written.
	 * @param  os  output stream to which field is written.
	 */
	template <> std::ostream& UnparametricField<NTL::RR>::write(std::ostream& os) const 
		{ return os << "unparamterized field NTL::RR"; }


#else // <- new writer / reader methods

	// output field to Writer
	template <> bool UnparametricField<NTL::RR>::toTag(Writer &W) const
	{

		W.setTagName("field");
		W.setAttribute("implDetail", "ntl-RR");
		W.setAttribute("cardinality", "-1");

		W.addTagChild();
		W.setTagName("real");
		W.upToParent();

		return true;
	}

	 // print field
        template <> ostream &UnparametricField<NTL::RR>::write(ostream &os) const
        {
                Writer W;
                if( toTag(W))
                        W.write(os);

                return os;
        }


	// write field element to writer
	template <> bool UnparametricField<NTL::RR>::toTag(Writer &W, const Element &e) const
	{
		string s;
		W.setTagName("cn");

		// note, this call is supported because
		// the RR class in the NTL has an operator<< 
		// which the template version of numToString uses
		//
		W.addDataChild(Writer::numToString(s, e));

		return true;
	}

        // output field element
        template <> ostream &UnparametricField<NTL::RR>::write(ostream &os, const Element &e) const
        {

                Writer W;
                if( toTag(W, e))
                        W.write(os);

                return os;
        }


	// read field element using Reader
	template <> bool UnparametricField<NTL::RR>::fromTag(Reader &R, Element &e) const
	{
		// This method uses an overloaded istream operator<<, which
		// NTL properly provides.  So this call is correct
		//
		return R.expectTagNum(e);
	}

	// read field element from istream
        template <> istream &UnparametricField<NTL::RR>::read(istream &is, Element &e) const
        {
                Reader R(is);
                if( !fromTag(R, e)) {
                        is.setstate(istream::failbit);
                        if(!R.initalized())
                                is.setstate(istream::badbit);
                }

                return is;
        }

#endif 

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
	template <> NTL::RR& UnparametricRandIter<NTL::RR>::random(NTL::RR &elt) const
		{
			// Create new random elements
			if (_size == 0)
				elt = rand();
			else
				elt = static_cast<long>((double(rand())/RAND_MAX)*double(_size));

#ifdef TRACE
			double temp = elt;
			cout << "random double = " << temp 
			     << "    random NTL::RR = " << elt << endl;
#endif // TRACE

			return elt;
    
		} // element& operator() (void)



	//@} 
} // namespace LinBox

#endif // __FIELD_NTL_RR_H
