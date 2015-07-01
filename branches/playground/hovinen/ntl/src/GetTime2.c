#include <NTL/config.h>

#include <sys/time.h>
#include <sys/resource.h>


#if (defined(__cplusplus) && !defined(NTL_CPLUSPLUS_ONLY))
extern "C" double GetTime();
#endif


double GetTime(void)
{
   struct rusage used;

   getrusage(RUSAGE_SELF, &used);
   return (used.ru_utime.tv_sec + used.ru_stime.tv_sec +
      (used.ru_utime.tv_usec + used.ru_stime.tv_usec) / 1e6);
}

