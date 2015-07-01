
#ifndef NTL_FacVec__H
#define NTL_FacVec__H

#include <NTL/vector.h>

struct IntFactor {
   IntFactor() { }
   ~IntFactor() { }

   long q;
   long a;
   long val;
   long link;
};


NTL_vector_decl(IntFactor,vec_IntFactor)
typedef vec_IntFactor FacVec;

void FactorInt(FacVec& fvec, long n);

#endif
