//--------------------------------------------------------------------------
//          Test for the LUdivine LUP factorisation
//--------------------------------------------------------------------------
// Clement Pernet, Morgan Brassel
// 17/06/2002
// Last modified 02/04/2003
// usage: Test_LUdivine p A n, for n LUP fatorisation of A over Z/pZ
//-------------------------------------------------------------------------

#define COMPUTE_L 0
//-------------------------------------------------------------------------
#define DEBUG 0
// Debug option  0: no debug
//               1: check A = LUP 
//               2: trace the recursive algorithm
//-------------------------------------------------------------------------

//#define GLOBALTIMINGS

#include <stdio.h>
#include <stdlib.h>
#include "givinit.h"
#include "givgfq.h"
#include "givtimer.h"
#include "FFLAP.h"
#include "Matio.h"
#include "givzpz32std.h"



using namespace std;

//----------------------------------------------------------------------------
// Choice of the finite field representation
//typedef GFqDom<long> GFqDomain;
typedef ZpzDom<Std32> GFqDomain;
//----------------------------------------------------------------------------


int main(int argc, char** argv){
  long i,j,id,R;
  int r=atoi(argv[1]);
  int k=2*r+1;
  int p=2*k+1;
  Timer tim, tim1;
  GFqDomain F(2*k+1);
  GFqDomain::element sq,tmp, *Ai;
  tim.clear();
  tim.start();

  GFqDomain::element * A = new GFqDomain::element[k*k];
  for ( i = 0; i < k*k; i++)
    F.neg(*(A+i), F.one);
  for ( i = 0; i< k; i++){
    F.init(tmp, i+1);
    F.mul(sq, tmp, tmp);
    F.convert(id, sq);
    if (id>k)
    {
      id = 2*k-id;
      Ai = A+ k-1+k*(k-id);
    }
    else
      Ai = A + id-1;
    for (j=id; j>0; j--, Ai+=k-1){
      *Ai = F.one;
    }
  }
  tim.stop();
#if DEBUG
  write_field(F,cerr,A,k,k,k);
#endif

   tim1.clear();
   tim1.start();
   R = FFLAP::Rank( F, k, k, A, k);
   tim1.stop();
  
  if (R<k){
    cerr<<"No  [r="<<r<<", k="<<k<<", p="<<p<<" Rank="<<R<<"] Const. of A: "
	<<tim.usertime()<<"s Rank: "<<tim1.usertime()<<"s"<<endl;
    return 0;}
  else{
    cerr<<"Yes [r="<<r<<", k="<<k<<", p="<<p<<" Rank="<<R<<"] Construction of A: "
	<<tim.usertime()<<"s Rank computation: "<<tim1.usertime()<<"s"<<endl;
    return 1;}
}







