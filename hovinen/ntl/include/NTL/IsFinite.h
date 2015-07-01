
#ifndef NTL_IsFinite__H
#define NTL_IsFinite__H

#include <NTL/config.h>

#if (defined(__cplusplus) && !defined(NTL_CPLUSPLUS_ONLY))
extern "C" {
#endif

long IsFinite(double *p);
/* This forces a double into memory, and tests if it is "normal";
   that means, not NaN, not +/- infinity, not denormalized, etc.
   Forcing into memory is sometimes necessary on machines 
   with "extended" double precision registers (e.g., Intel x86s)
   to force the standard IEEE format. */

void ForceToMem(double *p);
/* This is do-nothing routine that has the effect of forcing
   a double into memory (see comment above). */

#if (defined(__cplusplus) && !defined(NTL_CPLUSPLUS_ONLY))
}
#endif

#endif

