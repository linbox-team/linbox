//--------------------------------------------------------------------------
//          Test for the CharPolylaswp computation 
//--------------------------------------------------------------------------
// Clement Pernet
// 15/04/2003
//-------------------------------------------------------------------------

//-----------------------------------------------------
#define DEBUG 2
// Option de debuggage 0: rien
//-----------------------------------------------------
using namespace std;

#include <iostream>
#include <fstream>
#include <vector>
#include <cstdio>
#include "givgfq.h"
#include "Field_laswp.h"
#include "Matio.h"

//----------------------------------------------------------------------------
// Choice of the finite field representation
typedef GFqDom<long> GFqDomain;
//typedef ZpzDom<Std32> GFqDomain;
//----------------------------------------------------------------------------

int main(int argc, char** argv){
  int m,n,i;

  const GFqDomain F(atoi(argv[1]));
  cerr<<"Test of Flaswp:"<<endl
      <<"Reading Matrix...";
  long * A = read_field(F, argv[2],&m,&n);
  cerr<<"Ok"<<endl;

  int r = atoi(argv[3]);
  int k = atoi(argv[4]);
  int P[n];
  for(i=0;i<n;i++)
    P[i]=n-i-1;
  
  GFqDomain::element * U = new  GFqDomain::element[n*n];
  for ( i=0;i<n*n;i++)
    F.init(U[i],F.zero);
  cerr<<"A="<<endl;
  write_field(F,cerr,A,n,n,n);
  
  cerr<<"U="<<endl;
  write_field(F,cerr,U,n,n,n);

  cerr<<"Call to flaswp_trsp_half_copy"<<endl;
  flaswp_trsp_copy( F, n , A, n, 0, r, P, 1, U, n, k);

  cerr<<"U="<<endl;
  write_field(F,cerr,U,n,n,n);
  
  return 0;
}














