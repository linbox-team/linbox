#ifndef _GIV_COMPGMP_H_
#define _GIV_COMPGMP_H_
// ==========================================================================
// $Source$
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Author: T. Gautier
// $Id$
// ==========================================================================

#ifdef __MWERKS__
#ifdef __cplusplus
// -- CW release 10 do not define __STDC__ when compiling using C++
#define __STDC__
#endif
#endif

// Useful macro
#define SGN(l)    ((l) <0 ? -1 : (l >0 ? 1 : 0))
#define ABS(l)     ((l) <0 ? -l : l)
#define MIN(a,b)   ((a) < (b) ? (a) : (b))
// #define SZ_REP(r)  (ABS(r.size))
#define SZ_REP(r)  ( mpz_size( (mpz_ptr)&(r) ) )
inline long MAX(long a, long b){ return (a<b ? b: a ); }   

// --- Compatibility Gmp1.3.2 -> Gmp1.906
// #define mpz_ptr MP_INT*
// #define __STDC__
extern "C" {
#include "gmp.h"
}

#endif
