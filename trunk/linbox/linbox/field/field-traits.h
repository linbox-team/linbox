/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/field/field-traits.h
 * Copyright (C) June 2004 Dan Roche
 * see COPYING for permissions etc.
 */

#ifndef __FIELD_TRAITS_H
#define __FIELD_TRAITS_H

#include <linbox/integer.h>

// Namespace in which all LinBox library code resides
namespace LinBox {

/* FieldTraits will hold some basic information about each field.  It will try
 * to take the information from the field when possible, and use defaults
 * otherwise.
 * maxModulus returns the greatest modulus that is usable for a given field, -1
 * for infinite. goodModulus returns takes an integer and returns true if that
 * integer is a valid modulus for the given field.
 * maxExponent and goodExponent do the same for the prime power.
 */

class RingCategories {

public:

	//generica ring.
	struct GenericTag{};
	//If it is isomorphic to Z/mZ, for some m or its extensions.
	struct ModularTag : public virtual GenericTag{};
	//If it is isomorphic to Z
	struct IntegerTag : public virtual GenericTag{};
	//If it is isomorphic to Q
	struct RationalTag : public virtual GenericTag{};
};

template <class Field>
struct ClassifyRing {
typedef	RingCategories::GenericTag categoryTag;
};

template <class Field>
struct FieldTraits
{
	typedef ClassifyRing<Field>::categoryTag categoryTag;
	
	static integer& maxModulus( integer& i ) {
		return i = static_cast<integer>(Field::getMaxModulus());
	}
	static bool goodModulus( const integer& i ) {
		integer max;
		maxModulus( max );
		if( max == -1 ) return ( i >= 2 );
		else if( max == 0 ) return ( i == 0 );
		else return ( i >= 2 && i <= max );
	}

	static integer& maxExponent( integer& i ) { return i = 1; }
	static bool goodExponent( const integer& i ) {
		integer max;
                maxExponent( max );
                if( max == -1 ) return ( i >= 1 );
                else return ( i >= 1 && i <= max );
	}
};

} // Namespace LinBox

#endif // __FIELD_TRAITS_H
  
