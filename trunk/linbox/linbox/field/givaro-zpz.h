/* -*- mode: C++; style: linux -*- */

/* linbox/field/givaro-zpz.h
 * Copyright (C) 2002 Pascal Giorgi
 *
 * Written by Pascal Giorgi <pascal.giorgi@ens-lyon.fr>
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

#ifndef __FIELD_GIVARO_ZPZ
#define __FIELD_GIVARO_ZPZ


#include "linbox/integer.h"

//-------------------------------------
// Files of Givaro library
#include <givzpz16std.h>
#include <givzpz32std.h>
#include <givzpz16table1.h>
#include <giv_randiter.h>

//--------------------------------------



// Namespace in which all LinBox code resides
namespace LinBox 
{ 

	/*  This wrappers allows to use three sorts of givaro fiels :
	 *  Elements represent by a 32 bits integer
	 *  Elements represent by a 16 bits integer
	 *  Elements represent in Zech log representation by a 16 bits integer
	 *
	 *  To use this fields with the wrapper below just replace the template
	 *  parameter by the Tag appropriated.
	 *  "Std16"  for 16 bits integer
	 *  "Std32"  for 32 bits integer
	 *  "Log16"  for Zech log representation in 16 bits
	 */




	/** This template class is define just to be in phase with the LinBox
	 *  archetype. Read the archetype to know all functions are available.
	 *  Most of all methods are inherited from ZpzDom<Std16>, ZpzDom<Std32>
	 *  and ZpzDom<log16> class of Givaro.
	 *  these class allow to construct only finite field with a prime modulus.
	 */   

	template <class TAG> class GivaroZpz : public ZpzDom<TAG>
	{
	    public:

		//typedef integer Integer;

		/** Element type.
		 *  This type is inherited from the Givaro class ZpzDom<TAG>
		 */
		typedef typename ZpzDom<TAG>::Rep Element;

		/** RandIter type
		 *  This type is inherited from the Givaro class ZpzDom<TAG>
		 */	
		typedef GIV_randIter< ZpzDom<TAG>, Element > RandIter;

		/** Constructor from an integer
		 *  this constructor use the ZpzDom<TAG> constructor
		 */
		GivaroZpz (const integer &p) : ZpzDom<TAG> (static_cast<int> (p)) {}

		/** Characteristic.
		 * Return integer representing characteristic of the domain.     
		 * @return integer representing characteristic of the domain.
		 */
		integer &characteristic (integer &c) const
			{ return c = integer (static_cast<int> (ZpzDom<TAG>::size ())); }

		/** Cardinality. 
		 * Return integer representing cardinality of the domain.     
		 * @return integer representing cardinality of the domain
		 */
		integer &cardinality (integer &c) const
			{ return c = integer (static_cast<int> (ZpzDom<TAG>::size ())); }

		/** Conversion of field base element to an integer.
		 * This function assumes the output field base element x has already been
		 * constructed, but that it is not already initialized.
		 * @return reference to an integer.
		 * @param x integer to contain output (reference returned).
		 * @param y constant field base element.
		 */
		integer &convert (integer &x, const Element &y) const
			{ return x = integer (static_cast<int> (y)); }

		/** Initialization of field base element from an integer.
		 * Behaves like C++ allocator construct.
		 * This function assumes the output field base element x has already been
		 * constructed, but that it is not already initialized.     
		 * @return reference to field base element.
		 * @param x field base element to contain output (reference returned).
		 * @param y integer.
		 */  
		Element &init (Element &x , const integer &y = 0) const
			{ return ZpzDom<TAG>::init (x, long (y)); }

	}; // class GivaroZpz<TAG>
 
	/** Specialisation of the convert function for the zech log representation
	 *	of givaro-zpz (GivaroZpz<Log16>.
	 *  this function translates the internal representation to the real 
	 *	value of the element.
	 *	This can have no sense but can be usefull
	 *  NB : the init function for this specialisation does the same thing.
	 *  the function transaltes the values to her internal representation.
	 */ 
	integer& GivaroZpz<Log16>::convert(integer& x, const Element& y) const
	{     
		int tmp = _tab_rep2value[y];
		return x = integer (tmp);
	}

} // namespace LinBox

#endif // __FIELD_GIVARO_ZPZ
