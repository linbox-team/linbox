
#include <NTL/IsFinite.h>

double IsFinite__local;
double *IsFinite__ptr1 = &IsFinite__local;
double *IsFinite__ptr2 = &IsFinite__local;

long IsFinite(double *p)
{
   *IsFinite__ptr1 = *p;
   *IsFinite__ptr2 = (*IsFinite__ptr2 - *p);
   if (*IsFinite__ptr1 != 0.0) return 0;
   return 1;
}

void ForceToMem(double *p)
{ }
