/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
//--------------------------------------------------------------------------
//          Test for the CharPoly computation ( Krylov-LU )
//--------------------------------------------------------------------------
// Clement Pernet
// 17/04/2003
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
#define DEBUG 0
// Debug option:
//               0 = time +result
//               1 = factors + time + result
//               2 = matrix + details + time + result
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
#define CONSTRUCT 
// Choice for the minpoly algorithm:
//               ifdef  = computes the n Krylov vectors, then factorize U
//               ifndef = computes only the 2^k first Krylov vectors, where
//                   2^(k-1) <= deg(minpoly) < 2^k
//               second solution is better for matrices having a low degree
//               minpoly, but it may require more memory storage
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
//#define KGLU
// Select the Keller-Gehrig branching algorihtm
//-------------------------------------------------------------------------


#include <iostream>
#include <vector>
#include "givzpz32std.h"
#include "givgfq.h"
#include "givinit.h"
#include "givtimer.h"
#include "FFLAP.h"
#include "Matio.h"


//----------------------------------------------------------------------------
// Choice of the finite field representation
//typedef GFqDom<long> GFqDomain;
typedef ZpzDom<Std32> GFqDomain;
//----------------------------------------------------------------------------

template<class Field>
void 
print_poly(const Field & F, vector<typename Field::element> & P){
	int i = 0;
	long tmp;
	vector<typename Field::element>::iterator it = P.begin();
	F.convert(tmp,*it++);
	while (!tmp){
		F.convert(tmp,*it++); 
		i++;
	}
	cout<<tmp;
	if (i) cout<<"X^"<<i;
	i++;
	for (;it!=P.end(); it++,i++){
		F.convert(tmp,*it);
		if ( tmp ) cout<<" + "<<tmp<<"X^"<<i; 
	}
	cout<<endl;
}

template<class Field>
vector<typename Field::element> &
mulpoly(const Field & F, 
	vector<typename Field::element> &res, 
	const vector<typename Field::element> & P1,
	const vector<typename Field::element> & P2){
	int i,j;
	for (i=0;i<res.size();i++)
		F.init(res[i],0.0);
	for ( i=0;i<P1.size();i++)
		for ( j=0;j<P2.size();j++)
			F.axpyin(res[i+j],P1[i],P2[j]);
		
}

int main(int argc, char** argv){
  int m,n;
  Timer tim;
  list<vector<GFqDomain::element> > Charp;
  list<vector<GFqDomain::element> >::iterator it;
  GFqDomain F(atoi(argv[1]));
#if DEBUG
  cerr<<"Characteristic Polynomial Computation:"<<endl
      <<"Reading Matrix...";
#endif
  GFqDomain::element * A = read_field(F, argv[2],&m,&n);
  GFqDomain::element * U = new GFqDomain::element[n*(n+1)];
  vector<GFqDomain::element> prod(n+1);
  vector<GFqDomain::element> tmp(n+1);
#if DEBUG
  cerr<<"Ok"<<endl;
#endif
  for (int i=0;i<(n+1)*n;i++)
	  F.init(U[i],F.zero);

#if DEBUG==2
  cerr<<"A="<<endl;
  write_field(F,cerr,A,n,n,n);
#endif
#if DEBUG
  cerr<<"Starting computation of Charpoly:"<<endl;
#endif
  Charp.clear();
  tim.start();
  FFLAP::CharPoly( F, Charp, n, A, n, U, n );
  tim.stop();
  delete[] U;
  delete[] A;
#if DEBUG 
  cerr<<"Charpoly computed"<<endl;
#endif
  it=Charp.begin();
  F.init(prod[0],1.0);
  while(it!=Charp.end()){
#if DEBUG
	  print_poly(F, *it);
#endif
	  tmp = prod;
	  mulpoly(F, prod, tmp, *it);
	  it++;
  }

  cout<<"Charpoly(A) = ";
  print_poly(F, prod );

  //  cout<<"t="<<tim.usertime()<<endl;
  return 0;
}














