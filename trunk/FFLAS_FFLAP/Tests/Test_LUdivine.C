//--------------------------------------------------------------------------
//          Test for the LUdivine LUP factorisation
//--------------------------------------------------------------------------
// Clement Pernet, Morgan Brassel
// 17/06/2002
// Last modified 02/04/2003
// usage: Test_LUdivine p A n, for n LUP fatorisation of A over Z/pZ
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
#define DEBUG 1
// Debug option  0: no debug
//               1: check A = LUP 
//               2: trace the recursive algorithm
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
// Diplay positions of errors in LUP!=A
#define DISP_ERRORS 1
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
//typedef GFqDom<long> GFqDomain;
typedef ZpzDom<Std32> GFqDomain;
//----------------------------------------------------------------------------

int main(int argc, char** argv){
  int i,j,m,n,nbf,R;
  GFqDomain F(atoi(argv[1]));
  GFqDomain::element * A;
  
  A = read_field(F,argv[2],&m,&n);
  
  size_t *Permut = new size_t[n];
  for ( i=0;i<m;i++)
    Permut[i]=0;
#if DEBUG
  GFqDomain::element  * A_debug = new GFqDomain::element[m*n];
  for ( i=0;i<m;i++)
    for ( j=0;j<n;j++){
      A_debug[j+n*i] = A[j+n*i];
    }
#endif
  nbf = atoi(argv[3]);
  
  Timer tim,tim1;
  tim.clear();

#if DEBUG==1
  cerr<<"A="<<endl;
  write_field(F,cerr,A,m,n,n);
#endif
  for ( i=0;i<nbf;i++){
    if (i) A = read_field(F,argv[2],&m,&n);
    for (j=0;j<n;j++)
      Permut[j]=0;
    tim1.clear();      
    tim1.start();	
    R = FFLAP::LUdivine( F, FFLAS::FflasNonUnit, m, n, A, n, Permut, FFLAP:: FflapLSP);
    tim1.stop();
    tim+=tim1;
#if DEBUG == 0
    delete[] A;
#endif
  }

#if DEBUG  
  size_t MN=(m<n)?m:n;
  GFqDomain::element  U[m*n];
  GFqDomain::element  L[m*m];
  
  size_t pivot = 0;
  for ( i=0;i<n;i++){
    if (!F.iszero(A[i*n+pivot])){
       for (j=0;j<pivot;++j)
	 U[i*n+j]=F.zero;
       for (j=pivot;j<n;++j)
	 U[i*n+j]=A[i*n+j];
       pivot++;
    }
    else{
      for (j=0;j<=pivot;++j)
	U[i*n+j]=F.zero;
      for (j=pivot+1;j<n;++j)
	U[i*n+j]=A[i*n+j];
    }
  }
  
  size_t jl=0;
  pivot = 0;
  for ( j=0;j<R;++j){
    for ( i=0; i<pivot;++i){
      L[i*m+jl] = F.zero;
    }
    L[pivot*m+jl] = F.one;
    for ( i=pivot+1;i<m;++i){
      L[i*m+jl] = A[i*n+j];
    }
    // search for the next non zero pivot
    pivot++;
    jl++;
    while ( F.iszero ( A[j+1+n*pivot]) && jl<m){
      // inserting 0 column in L
      for ( i = 0; i<m; ++i )
	L[i*m+jl] = F.zero;
      pivot++;
      jl++;
    }
    
  }
  
#endif

  #if DEBUG==1
  cerr<<"L="<<endl;
  write_field(F,cerr,L,m,m,m);
  cerr<<"U="<<endl;
  write_field(F,cerr,U,m,n,n);
  #endif
#if DEBUG
  FFLAS::fgemm(F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, m,n,MN,F.one,L,n,U,n,
		 F.zero,U,n);
#endif
#if DEBUG==2
  cerr<<"L*U="<<endl;
  write_field(F,cerr,U,m,n,n);
#endif
  // Apply P
#if DEBUG
  int affiche=0;
  size_t nl = n;
  FFLAP::flaswp(F, m, U, nl, 0, R, Permut, -1);
  for ( i=0;i<m;++i){
    for ( j=0;j<n;++j){
      if (!F.areEqual(A_debug[n*i+j],U[n*i+j])){
	//cerr<< " error in i="<<i<<" j="<<j<<endl;
#if DISP_ERRORS
	cerr << "X ";
#endif
	affiche=1;
      }
      else{ 
#if DISP_ERRORS
	cerr << ". ";
#endif
      }
    }
#if DISP_ERRORS
    cerr<<endl;
#endif
  }
  if (affiche) cerr<<"Error occured during in computation: A!=LUP"<<endl;
#endif
#if DEBUG==1
  cerr<< "L*U*P ="<<endl;
  write_field(F, cerr, U,m,n,n);
  write_field(F,cerr, A,m,n,n);
#endif

  double t = tim.realtime();
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














