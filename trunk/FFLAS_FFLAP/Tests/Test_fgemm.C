//--------------------------------------------------------------------------
//                        Test for the  FFLAS_fgemm
//                  
//--------------------------------------------------------------------------
// Clement Pernet
// 24/04/2003
//-------------------------------------------------------------------------

#include <iomanip>
#include <iostream>
#include "givinit.h"
#include "givgfq.h"
#include "givtimer.h"
#include "Matio.h"
#include "FFLAS.h"

#define DEBUG 0

typedef GFqDom<long> GFqDomain;
typedef TTDom<double> DoubleDomain;


int main(int argc, char** argv){
  int m,n,k;
  int nbs=atoi(argv[4]);
  int nbf=atoi(argv[5]);
  long alpha = 1;
  if (argc>6){
    alpha=atoi(argv[6]);
    cerr<<"alpha="<<alpha<<endl;
  }
  long beta=0;
  if (argc>7){
    beta=atoi(argv[7]);
    cerr<<"beta="<<beta<<endl;
  }
  int domain=atoi(argv[1]);
  Timer tim;
  // 0 means double, >0 means finite field
  
  if (!domain)
    {
      DoubleDomain D;
      DoubleDomain::element* Ad;
      DoubleDomain::element* Bd;
      Ad = read_dbl(argv[2],&m,&k);
      Bd = read_dbl(argv[3],&k,&n);
      DoubleDomain::element* Cd = new DoubleDomain::element[m*n];
      DoubleDomain::element* Cd0 = new DoubleDomain::element[m*n];
      if (nbs==100)     // Tests for all Winograd's level
	{
	  for (int z=0;z<7;z++){
	    tim.clear(); tim.start();
	    for(int i = 0;i<nbf;++i)
		FFLAS::fgemm(D, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
			     m,n,k, 1.0, Ad, k, Bd, n, 0.0, Cd,n,z);
	    tim.stop();
	    cerr << z << " Winograd's level over double : t= "
		 << tim.usertime()/nbf 
		 <<  " s, Mflops: " << (2.0*m*k/1000000.0)*nbf*n/tim.usertime()
		 << endl;
	  }
	  tim.clear(); tim.start();
	  for(int i = 0;i<nbf;++i)
	    //	    FFLAS::fgemm(D, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
	    // m,n,k,1.0, Ad,k,Bd,n,0.0,Cd,n);
	  tim.stop();
	  cerr << " Automatic Winograd's level choice  over double : t= "
	       << tim.usertime()/nbf 
	       <<  " s, Mflops: " << (2.0*m*k/1000000.0)*nbf*n/tim.usertime()
	       << endl; 
	}
      else
	{
	  tim.clear(); tim.start();
	  for(int i = 0;i<nbf;++i)
	    //FFLAS::fgemm(D,FFLAS::FflasNoTrans,FFLAS::FflasNoTrans,
	    //	   m,n,k, 1.0, Ad, k,Bd, n, 0.0,Cd,n,nbs);
	  tim.stop();
	  cerr << nbs << " Winograd's level over double : t= "
	       << tim.usertime()/nbf 
	       << " s, Mflops = "<<(2.0*m*k/1000000.0)*nbf*n/tim.usertime() 
	       << endl;
	  
 
  	  //FFLAS::fgemm(D, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
	  //       m,n,k, 1.0, Ad,k, Bd,n, 0.0,Cd0,n,0);
	  for (int i=0;i<m;++i)
	    for (int j=0;j<n;++j)
	      if (Cd[n*i+j] != Cd0[n*i+j]) 
		cerr<< " error in i="<<i<<" j= "<<j<<endl;
	  cerr << "m, n, k = "<<m<<" "<< k<<" "<<n<<endl;
	} 
    }
  else
    {
      GFqDomain F(atoi(argv[1]));
      GFqDomain::element * A;
      GFqDomain::element * B;
      A = read_field(F,argv[2],&m,&k);
      B = read_field(F,argv[3],&k,&n);
      GFqDomain::element * C = new GFqDomain::element[m*n];
#if DEBUG
      GFqDomain::element * C0 = new GFqDomain::element[m*n];
#endif
      if (nbs==100)     // Tests for all Winograd's level

	{ 
	  for (int z=0;z<7;z++){
	    Timer tim; tim.clear(); tim.start();
	    for(int i = 0;i<nbf;++i)
		FFLAS::fgemm(F,FFLAS::FflasNoTrans,FFLAS::FflasNoTrans,
			     m,n,k,F.one, A,k, B,n,F.zero,C,n,z);
	    tim.stop();
	    cerr << z << " Winograd's level over finite Field : t= "
		 << tim.usertime()/nbf 
		 <<  " s, Mflops: " << (2.0*m*k/1000000.0)*nbf*n/tim.usertime()
		 << endl; 
	  }
	  tim.clear(); tim.start();
	  for(int i = 0;i<nbf;++i)
		FFLAS::fgemm(F,FFLAS::FflasNoTrans,FFLAS::FflasNoTrans,
			     m,n,k,F.one, A,k, B,n, F.zero,C,n);
	  tim.stop();
	  cerr << " Automatic Winograd's level choice  over finite Field : t= "
	       << tim.usertime()/nbf 
	       <<  " s, Mflops: " << (2.0*m*k/1000000.0)*nbf*n/tim.usertime()
	       << endl; 
}
      else
	{
	  GFqDomain::element Alp, Bet;
	  F.init(Alp,alpha);
	  F.init(Bet,beta);
	  if (beta)
	    C = read_field(F,argv[8],&m,&n);
	  Timer tim; tim.clear(); tim.start();
	  for(int i = 0;i<nbf;++i){
	    //cerr<<"m,n,k,Alp,Bet,nbs"<<m<<", "<<n<<", "<<k<<", "<<Alp
	    //	<<", "<<Bet<<", "<<nbs<<endl;
	    FFLAS::fgemm(F,FFLAS::FflasNoTrans,FFLAS::FflasNoTrans,
			 m,n,k,Alp, A,k, B,n, Bet,C,n,nbs);
	  }
	  tim.stop();
	  write_field(F,cerr,C,m,n,n);
	  cerr << nbs << " Winograd's level over finite Field : t= "
	       << tim.usertime()/nbf 
	       << " s, Mflops = "<<(2.0*m*k/1000000.0)*nbf*n/tim.usertime() 
	       << endl;
	  
	  F.init(Alp,alpha);
	  F.init(Bet,beta);
	  cerr<<"m,n,k,Alp,Bet,nbs"<<m<<", "<<n<<", "<<k<<", "<<Alp
	      <<", "<<Bet<<", "<<nbs<<endl;
#if DEBUG
	  FFLAS::fgemm(F,FFLAS::FflasNoTrans,FFLAS::FflasNoTrans,m,n,k,Alp,
		       A,k, B,n,Bet,C0,n,0);
	  write_field(F,cerr,C0,m,n,n);
	  cerr<<std::setiosflags(std::ios::fixed)<<std::setprecision(16);
	  for (int i=0;i<m;++i)
	    for (int j=0;j<n;++j)
	      if (C[n*i+j] != C0[n*i+j]) {
		cerr<< " error in i="<<i<<" j= "<<j<<endl;
		cerr<<"C[n*i+j]="<<C[n*i+j]<<" C0[n*i+j]="<<C0[n*i+j]<<endl;
	      }
	  cerr << "m, n, k = "<<m<<" "<< n<<" "<<k<<endl;
#endif
	} 
    }  
}














