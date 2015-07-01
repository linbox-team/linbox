#ifndef _GMPplusplus_H_
#define _GMPplusplus_H_
// ========================================================================
// LinBox version of gmp++.h
// Time-stamp: <19 Dec 06 10:49:28 Jean-Guillaume.Dumas@imag.fr> 
// ========================================================================

#ifndef __DONOTUSE_64__
#define __USE_64_bits__
#endif

#if !defined(GMP_VERSION_3) && !defined(GMP_NO_CXX)
#include <gmpxx.h>
#endif

// If GMP is at least version 4, do not need extern
#ifdef GMP_VERSION_3
extern "C" {
#endif

#include "gmp.h"

#ifdef GMP_VERSION_3
}
#endif

#include <gmp++/gmp++_int.h>


#ifdef LinBoxSrcOnly
#include <linbox/util/gmp++/gmp++_int.C>
#endif

#endif
