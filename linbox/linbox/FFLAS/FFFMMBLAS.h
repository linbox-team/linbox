/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/FFLAS/FFFMMBLAS.h
 *----------------------------------------------------------------------
 *                   Finite Field Fast Matrix Multiplication 
 *                using Basic Linear Algebra Subroutines
 *----------------------------------------------------------------------
 * Written by Clement PERNET ( clement.pernet@imag.fr )
 * 24/03/2003
 *----------------------------------------------------------------------
 *
 * Added and modified in LinBox by Pascal Giorgi <pascal.giorgi@ens-lyon.fr>
 *

#ifndef __FFFMMBLAS__
#define __FFFMMBLAS__

#include "linbox/FFLAS/lin_wrap_c++.h"
#include "linbox/integer.h"


typedef TTDom<double> DoubleDomain;

enum FFFMMBLAS_COEF {FffmmblasZero=0, FffmmblasOne=1, FffmmblasNone=-1};

class FFFMMBLAS {

  //num FFFMMBLAS_COEF ALPHA;
  //enum FFFMMBLAS_COEF BETA;
  int ALPHA;
  int BETA;


  //-----------------------------------------------------------------
  //            Matrix multiplications                               
  //-----------------------------------------------------------------

  // Classic multiplication
  template <class Field>
    void ClassicMatmul(const Field& F, 
		       typename Field::Element * A,
		       typename Field::Element * B,
		       typename Field::Element * C,
		       int ra, int rb, int rc,
		       int m, int k, int n);
  
  void ClassicMatmul(const DoubleDomain& F, 
		     double * Ad,
		     double * Bd,
		     double * Cd,
		     int ra, int rb, int rc,
		     int m, int k, int n);
  
  // Winograd Multiplication  A(n*k) * B(k*m) in C(n*m)

   // WinoCalc perform the 22 Winograd operations
   template<class Field>
     void WinoCalc( const Field& F, 
		    typename Field::Element * A,
		    typename Field::Element * B,
		    typename Field::Element * C,
		    typename Field::Element ** t_X1,
		    typename Field::Element ** t_X2,
		    typename Field::Element ** t_X3,
		    int ra, int rb, int rc,
		    int nr, int kr, int mr,
		    long long sf, int s);

   template<class Field>
     void WinoMain( const Field& F, 
		    typename Field::Element * A,
		    typename Field::Element * B,
		    typename Field::Element * C,
		    typename Field::Element ** t_X1,
		    typename Field::Element ** t_X2,
		    typename Field::Element ** t_X3,
		    int ra, int rb, int rc,
		    int m, int k, int n,
		    long long sf, int s);

    
   public:

   //-------------------------------------------------------------------
   // visible function: 
   //-------------------------------------------------------------------   
   
   // Compute C = alpha.A*B + beta.C over double
   double* 
     operator()( const DoubleDomain& F,
		 const int m,
		 const int n,
		 const int k,
		 const enum FFFMMBLAS_COEF alpha,
		 double* const A, const int lda,
		 double* const B, const int ldb, 
		 const enum FFFMMBLAS_COEF beta,
		 double* C, const int ldc,
		 const int nbe);
   
   // Compute C = alpha.A*B + beta.C over a Field
   template<class Field>
     typename Field::Element * 
     operator()( const Field& F,
		 const int m,
		 const int n,
		 const int k,
		 //const enum FFFMMBLAS_COEF alpha,
		 int alpha,
		 typename Field::Element* const A, const int lda,
		 typename Field::Element* const B, const int ldb, 
		 //const enum FFFMMBLAS_COEF beta,
		 int beta,
		 typename Field::Element* C, const int ldc,
		 const int nbe);
};


// Some conversion functions

// Finite Field matrix => double matrix
template<class Field>
void MatGFq2MatDouble(const Field& F,
		      int m, int n,
		      typename DoubleDomain::Element* S,
		      int lds,
		      typename Field::Element* E,
		      int lde){
  LinBox::integer tmp;
  for (int i= 0;i<m;++i)
    for (int j=0;j<n;++j){
      //F.convert(tmp,*(E+j+i*lde));
      //*(S+j+i*lds)=(double) tmp;
      *(S+j+i*lds)=(double)*(E+j+i*lde);
    }
}

// Special design for upper-triangular matrices
template<class Field>
void MatGFq2MatDouble_Triangular(const Field& F,
				 int m, int n,
				 typename DoubleDomain::Element* S,
				 int lds,
				 typename Field::Element* E,
				 int lde){
  LinBox::integer tmp;
  for (int i = 0; i<m;++i)    
    for (int j=i; j<n;++j){
      //F.convert(tmp,*(E+j+lde*i));
      //*(S+j+i*lds)=(double) tmp;
      *(S+j+i*lds)=(double)*(E+j+lde*i);
    }
  
}

// double matrix => Finite Field matrix
template<class Field>
void MatDouble2MatGFq(  const Field& F,
			int m, int n,
			typename Field::Element* S,
			int lds,
			typename DoubleDomain::Element* E,
			int lde){

  LinBox::integer charact;
  typename Field::Element tmp;
  F.characteristic(charact);
  double dchar = charact;
  for (int i = 0; i<m;++i) 
    for (int j=0; j<n;++j) {
      F.init(*(S+j+lds*i),*(E+j+lde*i));
    //  F.init(tmp,LinBox::integer(*(E+j+lde*i)));
    }   
      
}

#include "FFFMMBLAS.inl"

#endif
