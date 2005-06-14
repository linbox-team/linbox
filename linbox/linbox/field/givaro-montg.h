/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/field/givaro-gfq.h
 * Copyright (C) 2004 Jean-Guillaume Dumas
 *
 * Written by Jean-Guillaume Dumas <Jean-Guillaume.Dumas@imag.fr>
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */

/* WARNING this wrapper works only with an improved version of Givaro.
 * This version of givaro won't be available for public yet.
 * But it is available on my web page.
 * You can send me a mail to get it or for others details.
 */

#ifndef __FIELD_GIVARO_MONTGOMERY
#define __FIELD_GIVARO_MONTGOMERY


#include <linbox/integer.h>
#include <linbox/field/field-interface.h>
#include <linbox/util/debug.h>
#include "linbox-config.h"
#include <linbox/field/field-traits.h>

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

//------------------------------------
// Files of Givaro library


#include <givaro/givmontg32.h>
#include <givaro/giv_randiter.h>
//------------------------------------

// Namespace in which all LinBox code resides
namespace LinBox 
{ 


	template <class Ring>
	struct ClassifyRing;

	class GivaroMontg;

	template<>
	struct ClassifyRing<GivaroMontg> {
		typedef RingCategories::ModularTag categoryTag;
	};
	/** This template class is define just to be in phase with the LinBox
	 *  archetype.
	 *  Most of all methods are inherited from Montgomery<Std32>  class
	 *  of Givaro.
	 *  this class is a modular representation with a Montgomery reduction
	 */   
	class GivaroMontg : public Montgomery<Std32>, public FieldInterface
	{
 
	public:

		/** Element type.
		 *  This type is inherited from the Givaro class Montgomery<Std32>
		 */
		typedef  Montgomery<Std32>::Rep Element;
    
		/** RandIter type
		 *  This type is inherited from the Givaro class Montgomery<Std32>
		 */	
		typedef GIV_randIter< Montgomery<Std32>, LinBox::integer >  RandIter;

		/** Constructor from an integer
		 *  this constructor use the ZpzDom<TAG> constructor
		 */
		GivaroMontg(const integer& p) :
			Montgomery<Std32>(static_cast<uint32>(long(p))) { }

		/** Constructor from an integer (takes degree of extension as 2nd parameter, must be 1)
		 *  this constructor use the ZpzDom<TAG> constructor
		 */
	  	GivaroMontg(const integer& p, const integer& k) :
			Montgomery<Std32>(static_cast<uint32>(long(p))) { 
			
			if (k!=1)
				throw PreconditionFailed(__FUNCTION__,__LINE__,"exponent must be 1");
		}
    
#ifdef __LINBOX_XMLENABLED
		// XML Reader constructor
		GivaroMontg(Reader &R)
		{
			integer p, n;

			if(!R.expectTagName("field") || !R.expectChildTag()) return;
			R.traverseChild();

			if(!R.expectTagName("finite") || !R.expectChildTag()) return;
			R.traverseChild();

			if(!R.expectTagName("characteristic") || !R.expectChildTag()) return;
			R.traverseChild();
			if(!R.expectTagNum(p));
			R.upToParent();

			if(R.getNextChild()) {
				if(!R.expectChildTag()) return;
				R.traverseChild();

				if(!R.expectTagName("extension") || !R.expectChildTag()) return;
				R.traverseChild();
				if(!R.expectTagNum(n)) return;
				R.upToParent();

				R.upToParent();
				R.getPrevChild();
			}
			else {
				n = Integer(1);
			}
			R.upToParent();
			R.upToParent();

			// now try building using the above constructor.  Note, NO
			// ATTEMPT IS MADE TO CATCH THE ERROR THIS METHOD CAN THROW, 
			// IT IS ALLOWED TO PASS THROUGH
			//
			GivaroMontg oth(p, n);
			*this = oth;

			return;
		}
#endif


		/** Characteristic.
		 * Return integer representing characteristic of the domain.
		 * Returns a positive integer to all domains with finite characteristic,
		 * and returns 0 to signify a domain of infinite characteristic.
		 * @return integer representing characteristic of the domain.
		 */
		integer& characteristic(integer& c) const
		{return c=integer(static_cast<long>(Montgomery<Std32>::characteristic()));}

		long characteristic() const
		{return static_cast<long>(Montgomery<Std32>::characteristic());}
    
      
		/** Cardinality. 
		 * Return integer representing cardinality of the domain.
		 * Returns a non-negative integer for all domains with finite
		 * cardinality, and returns -1 to signify a domain of infinite
		 * cardinality.
		 * @return integer representing cardinality of the domain
		 */
		integer& cardinality(integer& c) const
		{ return c=integer(static_cast<long>(Montgomery<Std32>::size()));}
 

		/** Initialization of field base Element from an integer.
		 * Behaves like C++ allocator construct.
		 * This function assumes the output field base Element x has already been
		 * constructed, but that it is not already initialized.
		 * We assume that the type of Element is short int.
		 * this methos is just a simple cast.
		 * @return reference to field base Element.
		 * @param x field base Element to contain output (reference returned).
		 * @param y integer.
		 */  
		Element& init(Element& x , const integer& y=0) const
		{ 			
			return Montgomery<Std32>::init( x,long(y % (integer)_p));
		}
      
		Element& init(Element& x , const double y=0.0) const
		{ return Montgomery<Std32>::init( x, y);}

		/** Conversion of field base element to an integer.
		 * This function assumes the output field base element x has already been
		 * constructed, but that it is not already initialized.
		 * @return reference to an integer.
		 * @param x integer to contain output (reference returned).
		 * @param y constant field base element.
		 */
		integer& convert(integer& x, const Element& y) const
		{
			long tmp;	
			//	return x = *(new integer(Montgomery<Std32>::convert(tmp,y)));
			return x = integer(Montgomery<Std32>::convert(tmp,y));
		}
		double& convert(double& x, const Element& y) const
		{
			return Montgomery<Std32>::convert( x, y);
		}

		//bool isZero(const Element& x) const { return Montgomery<Std32>::isZero(x); }


#ifdef __LINBOX_XMLENABLED

		ostream &write(ostream &os) const
		{
			Writer W;
			if( toTag(W) )
				W.write(os);

			return os;
		}

		bool toTag(Writer &W) const
		{
			string s;
			long card = Montgomery<Std32>::size();
			size_t i = 0;

			W.setTagName("field");
			W.setAttribute("implDetail", "givaro-gfq");
			W.setAttribute("cardinality", Writer::numToString(s, card));

			W.addTagChild();
			W.setTagName("finite");

			W.addTagChild();
			W.setTagName("characteristic");
			W.addNum(Montgomery<Std32>::characteristic());
			W.upToParent();
			W.addTagChild();
			W.setTagName("extension");

			while(card > 1) {
				card /= Montgomery<Std32>::characteristic();
				++i;
			}
			W.addNum(i);
			W.upToParent();

			W.upToParent();

			return true;
		}


		// Special Note:  In LinBox, all elements of a field will be written
		// in the following manner:  for e in ZZp[x] with 
		// e = a0 + a1x + a2x^2 + ..., e is represented as:
		// "<cn>n</cn>" where n = a0 + a1 * p + a2 * p^2 + ...
		// 
		ostream &write(ostream &os, const Element &e) const
		{
			Writer W;
			if( toTag(W, e))
				W.write(os);

			return os;
		}

		bool toTag(Writer &W, const Element &e) const
		{
			string s;
			long rep = _log2pol[ (unsigned long) e];

			W.setTagName("cn");
			W.addDataChild(Writer::numToString(s, rep));
		  
			return true;
		}

		istream &read(istream &is, Element &e) const
		{
			Reader R(is);
			if( !fromTag(R, e)) {
				is.setstate(istream::failbit);
				if(!R.initalized())
					is.setstate(istream::badbit);
			}

			return is;
		}

		bool fromTag(Reader &R, Element &e) const
		{
			unsigned long i;

			if(!R.expectTagName("cn") || !R.expectChildTextNum(i))
				return false;

			e = _pol2log[i];
			return true;
		}
			  

#endif

		static inline int getMaxModulus() { return 40504; }

	}; // class GivaroMontg
 


} // namespace LinBox

#endif // __FIELD_GIVARO_MONTGOMERY
