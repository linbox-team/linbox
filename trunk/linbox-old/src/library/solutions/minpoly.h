#ifndef __LINBOX_MINPOLY_H__
#define __LINBOX_MINPOLY_H__

#include "LinBox/lin_rand.h"                  // Random Iterator
#include "LinBox/lin_bbit.h"       // BB iterator
#include "LinBox/lin_massey.C"                // massey reccuring sequence solver

#include "LinBox/lin_methods.h"

namespace LinBox
{
  typedef double Pbound;  // an upper bound for the probability of erroneous result.

  /** p <-- minpoly(A), a factor of the minpoly of A.
   *  $\P ( p \neq \minpoly(A) ) \le p$.
   *  Separate calls are independent trials.
   */
template <class Polynomial, class BB, class MT>
Polynomial& minpoly(Polynomial& P, Pbound& p, const BB& A, MT mt = 0)
{ // FIXME mt ignored for now
  Random generator;
  unsigned long deg;

/// A is supposed to be symmetric
  BB_Container< BB > TF( &A, generator);
  MasseyDom< BB >  WD(&TF);
  WD.pseudo_minpoly(P, deg);
  p = 1; // until we compute a val.
  return P;
};

} // namespace LinBox
#endif


