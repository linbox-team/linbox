#ifndef NTL_PrimeSource__H
#define NTL_PrimeSource__H

#include <NTL/vec_long.h>


// The following are used to manage lists of random single-precision primes.
// Usage:  PrimeSource src; ... src >> p; ...
// This assigns successive random primes to p.
// Primes within a source are independent, but sources
// are "recycled".


struct PrimeSourceRep {
vec_long list;
PrimeSourceRep *link;
};

class PrimeSource {
private:
PrimeSourceRep *rep;
long n;

void operator=(const PrimeSource&);  // disabled
PrimeSource(const PrimeSource&);     // disabled

static PrimeSourceRep *FreeList;

public:
PrimeSource();
~PrimeSource();

friend PrimeSource& operator>>(PrimeSource& src, long& p);
};


#endif

