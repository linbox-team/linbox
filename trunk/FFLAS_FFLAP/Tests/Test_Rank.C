//--------------------------------------------------------------------------
//          Test for the Rank computation
//--------------------------------------------------------------------------
// Clement Pernet, Morgan Brassel
// 07/05/2002
// usage: Test_Rank p A n, for n Rank computations of A over Z/pZ
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
#define DEBUG 0
// Debug option  0: no debug
//               2: trace the recursive algorithm
//-------------------------------------------------------------------------

#include <iostream.h>
#include "givinit.h"
#include "givgfq.h"
#include "givtimer.h"
#include "Matio.h"
#include "givzpz32std.h"
#include "FFLAP.h"


using namespace std;

//----------------------------------------------------------------------------
// Choice of the finite field representation
typedef GFqDom<long> GFqDomain;
//typedef ZpzDom<Std32> GFqDomain;
//----------------------------------------------------------------------------

int main(int argc, char** argv){
  int i,j,m,n,nbf,R;
  GFqDomain F(atoi(argv[1]));
  GFqDomain::element * A;
  
  A = read_field(F,argv[2],&m,&n);
  
  nbf = atoi(argv[3]);
  
  Timer tim,tim1;
  tim.clear();

#if DEBUG==2
  cerr<<"A="<<endl;
  write_field(F,cerr,A,m,n,n);
#endif
  for ( i=0;i<nbf;i++){
    if (i) A = read_field(F,argv[2],&m,&n);
    tim1.clear();      
    tim1.start();	
    R = FFLAP::Rank( F, m, n, A, n);
    tim1.stop();
    tim+=tim1;
#if DEBUG == 0
    delete[] A;
#endif
  }

  double t = tim.usertime();
  double numops = m*m/1000.0*(m+3*n)/6.0;

   cerr
   	    <<m
   	    <<"x" 
  	    << n 
   	    << " : rank = " << R << "  ["
   	    << ((double)nbf/1000.0*(double)numops / t) 
   	    << " MFops "
   	    << " in "
   	    << t<<"s]" 
   	    << endl;
  cout<<m<<" "<<((double)nbf/1000.0*(double)numops / t)<<endl;

  return 0;
}














