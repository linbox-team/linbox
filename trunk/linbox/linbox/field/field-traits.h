/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/field/field-traits.h
 * Copyright (C) June 2004 Dan Roche
 * see COPYING for permissions etc.
 */

#ifndef __FIELD_TRAITS_H
#define __FIELD_TRAITS_H
#include <linbox/ntl.h>

// Namespace in which all LinBox library code resides
namespace LinBox {

/* FieldTraits will hold various methods for each field that do not need to be put
 * in the field definition itself.
 * maxModulus returns the greatest modulus that is usable for a given field, -1
 * for infinite.
 * goodModulus returns takes an integer and returns true if that integer is a valid
 * modulus for the given field.
 */

template <class Field>
struct FieldTraits
{
	static integer& maxModulus( integer& );
	static bool goodModulus( const integer& i ) {
		integer max;
		return ( i >= 2 && i <= maxModulus(max) );
	}
};

// Specializations for FieldTraits

template< class E > class Modular;
template< class E > class PIRModular;
class NTL_zz_p;

integer& FieldTraits< Modular< int32 > >::maxModulus( integer& i )
	{ return i = 1073741789; }

integer& FieldTraits< Modular< int16 > >::maxModulus( integer& i )
	{ return i = 32749; }

integer& FieldTraits< Modular< int8 > >::maxModulus( integer& i )
	{ return i = 127; }

integer& FieldTraits< Modular< double > >::maxModulus( integer& i )
	{ return i = 94906249; }

integer& FieldTraits< PIRModular<int32> >::maxModulus( integer& i )
	{ return FieldTraits< Modular<int32> >::maxModulus( i ); }

integer& FieldTraits< NTL_zz_p >::maxModulus( integer& i )
	{ return i = NTL_SP_BOUND; }


} // namespace Linbox

#endif // __FIELD_TRAITS_H
