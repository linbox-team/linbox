#include <NTL/config.h>

#include <time.h>

#ifndef CLOCKS_PER_SEC
#define CLOCKS_PER_SEC (1)
#endif


#if (defined(__cplusplus) && !defined(NTL_CPLUSPLUS_ONLY))
extern "C" double GetTime();
#endif


double GetTime(void)
{
   return clock()/((double)CLOCKS_PER_SEC);
}

