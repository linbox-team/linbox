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
 * maxExponent and goodExponent do the same for the prime power.
 */

template <class Field>
struct FieldTraits
{
	static integer& maxModulus( integer& );
	static bool goodModulus( const integer& i ) {
		integer max;
		maxModulus( max );
		if( max == -1 ) return ( i >= 2 );
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

// Class prototypes

template< class E > class Modular;
template< class E > class PIRModular;
class PIR_ntl_ZZ_p;
struct NTL_zz_p;
class NTL_ZZ_p;
class GF2;
class GMPRationalField;
struct Local2_32;
class NTL_ZZ;
class LidiaGfq;
struct NTL_PID_zz_p;

// Specializations for various fields

integer& FieldTraits< Modular< int32 > >::maxModulus( integer& i )
	{ return i = 1073741823; } // 2^30 - 1

integer& FieldTraits< Modular< int16 > >::maxModulus( integer& i )
	{ return i = 32767; } // 2^15 - 1

integer& FieldTraits< Modular< int8 > >::maxModulus( integer& i )
	{ return i = 127; } // 2^7-1

integer& FieldTraits< Modular< double > >::maxModulus( integer& i )
	{ return i = 94906265; } // floor(2^26.5)

integer& FieldTraits< PIRModular<int32> >::maxModulus( integer& i )
	{ return FieldTraits< Modular<int32> >::maxModulus( i ); }

integer& FieldTraits< NTL_zz_p >::maxModulus( integer& i )
	{ return i = NTL_SP_BOUND; }

integer& FieldTraits< NTL_zz_p >::maxExponent( integer& i )
	{ return i = integer( "4294967295" ); } // 2^32 - 1

integer& FieldTraits< NTL_PID_zz_p >::maxModulus( integer& i )
	{ return i = FieldTraits<NTL_zz_p>::maxModulus(i); }

integer& FieldTraits< NTL_ZZ_p >::maxModulus( integer& i )
	{ return i = -1; }

integer& FieldTraits< NTL_ZZ_p >::maxExponent( integer& i )
	{ return i = integer( "4294967295" ); } // 2^32 - 1

integer& FieldTraits< PIR_ntl_ZZ_p >::maxModulus( integer& i )
	{ return FieldTraits< NTL_ZZ_p >::maxModulus( i ); }

integer& FieldTraits< GF2 >::maxModulus( integer& i )
	{ return i = 2; }

integer& FieldTraits< GMPRationalField >::maxModulus( integer& i )
	{ return i = 0; }

integer& FieldTraits< Local2_32 >::maxModulus( integer& i )
	{ return i = 2; }

integer& FieldTraits< Local2_32 >::maxExponent( integer& i )
	{ return i = 32; }

bool FieldTraits< Local2_32 >::goodExponent( const integer& i )
	{ return (i == 32); }

integer& FieldTraits< NTL_ZZ >::maxModulus( integer& i )
	{ return i = 0; }

integer& FieldTraits< LidiaGfq >::maxModulus( integer& i )
	{ return i = integer( "9007199254740881" ); } // prevprime(2^53)

integer& FieldTraits< LidiaGfq >::maxExponent( integer& i )
	{ return i = 2147483647; } // 2^31-1

} // namespace Linbox

#endif // __FIELD_TRAITS_H
