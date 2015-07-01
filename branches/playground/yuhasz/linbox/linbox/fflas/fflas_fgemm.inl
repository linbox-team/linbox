/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/fflas/fflas_fgemm.inl
 * Copyright (C) 2003 Clement Pernet
 *
 * Written by Clement Pernet <Clement.Pernet@imag.fr>
 *
 * Warning, k.(p-1)^2 > 2^53 is still not implemented
 * ALPHA-BETA and TRANSPOSE options will force classic matrix multiplication
 *
 * See COPYING for license information.
 */

// Classic multiplication over a finite Field
// (Case where all Strassen operations have been performed over GFq)
#include "linbox/field/modular.h"
template <class Field>
void FFLAS::ClassicMatmul(const Field& F,  
			  const enum FFLAS_TRANSPOSE ta,
			  const enum FFLAS_TRANSPOSE tb,
			  const size_t m, const size_t n,const size_t k,
			  typename Field::Element ALPHA,
			  const typename Field::Element * A, 
			  const size_t lda,
			  const typename Field::Element * B, 
			  const size_t ldb,
			  typename Field::Element BETA,
			  typename Field::Element* C, const size_t ldc){
	static  typename Field::Element Mone;
	static  typename Field::Element one;
	F.init(Mone, -1);
	F.init(one, 1);
	size_t dlda,dldb;
	//if (k<kmax){
	// Double  matrices initialisation
	DoubleDomain::Element ALPHAd, BETAd;
	//cout<<"ALPHAd=";
	//	std::cout<<std::setiosflags(std::ios::fixed)
	//	 <<std::setprecision(16)
	//	 <<setw(30)<<ALPHAd<<"\n";
	//cout<<"BETAd=";
	//std::cout<<std::setiosflags(std::ios::fixed)
	//	 <<std::setprecision(16)
	//	 <<setw(30)<<BETAd<<"\n";


	DoubleDomain::Element * Add = new DoubleDomain::Element[m*k];
	DoubleDomain::Element * Bdd = new DoubleDomain::Element[k*n];
	DoubleDomain::Element * Cd = new DoubleDomain::Element[m*n];
	// Conversion finite Field => double
	if (ta == FflasTrans){
		dlda = m;
		MatF2MatD( F, Add, A, lda, k, m );
	}
	else{
		dlda = k;
		MatF2MatD( F, Add, A, lda, m, k );
	}
	
	if (tb == FflasTrans){
		dldb = k;
		MatF2MatD( F, Bdd, B, ldb, n, k );
 	}
	else{
		dldb = n;
		MatF2MatD( F, Bdd, B, ldb, k, n );
 	}
	if (!F.isZero(BETA)){
		MatF2MatD( F, Cd, C, ldc, m, n );
		if (F.areEqual(BETA, Mone))
			BETAd = -1.0;
		else
			F.convert( BETAd, BETA );
		
	}
	else
		F.convert( BETAd, BETA );
	if (F.areEqual(ALPHA, Mone))
		ALPHAd = -1.0;
	else
		F.convert( ALPHAd, ALPHA );

	// Call to the blas Multiplication 
	cblas_dgemm(CblasRowMajor, (enum CBLAS_TRANSPOSE) ta, 
		    (enum CBLAS_TRANSPOSE) tb,
		    m, n, k, (DoubleDomain::Element) ALPHAd,
		    Add, dlda, Bdd, dldb, (DoubleDomain::Element) BETAd,Cd, n);
	// Conversion double => Finite Field
	delete[] Add;
	delete[] Bdd;

	if ( !F.areEqual(ALPHA, Mone) && !F.isOne(ALPHA)){
		// Cd may contain approximations of C entries:
		// Perturbation of Cd to get C from truncation.
		double* Cdi=Cd;
		for (; Cdi<Cd+m*ldc; Cdi+=ldc)
			for(size_t j=0; j<n; ++j)
				if (*(Cdi+j)>=0){
					*(Cdi+j) += 0.1;
				}
				else
					*(Cdi+j) -= 0.1;
	}
	       
	MatD2MatF( F, C, ldc, Cd, m, n );
	delete[] Cd;
    /*}
 else{
    // NB reste a refaire le systeme de ALPHA BETA qui ne permet pas de faire C=A*B
    // en passant ici
    cerr<<"block computation over the field"<<endl;
    size_t nr = n/2;
    size_t kr = k/2;
    size_t mr = m/2;
    typename Field::Element * A12 = A+kr;
    typename Field::Element * A21 = A+mr*lda;
    typename Field::Element * A22 = A+kr+mr*lda;
    typename Field::Element * B12 = B+nr;
    typename Field::Element * B21 = B+kr*ldb;
    typename Field::Element * B22 = B+nr+kr*ldb;
    typename Field::Element * C12 = B+nr;
    typename Field::Element * C21 = B+mr*ldc;
    typename Field::Element * C22 = B+nr+mr*ldc;

    // C11 = A11.B11+A12.B21
    ClassicMatMul(F, A, B, C, lda, ldb, ldc, mr, kr, nr);
    ClassicMatMul(F, A12, B21, C, lda, ldb, ldc, mr, kr, nr);
    // C12 = A11.B12+A12.B22
    ClassicMatMul(F, A, B12, C12, lda, ldb, ldc, mr, kr, nr);
    ClassicMatMul(F, A12, B22, C12, lda, ldb, ldc, mr, kr, nr);
    // C21 = A21.B11+A22.B21
    ClassicMatMul(F, A21, B, C21, lda, ldb, ldc, mr, kr, nr);
    ClassicMatMul(F, A22, B21, C21, lda, ldb, ldc, mr, kr, nr);
    // C22 = A21.B12+A22.B22
    ClassicMatMul(F, A21, B12, C, lda, ldb, ldc, mr, kr, nr);
    ClassicMatMul(F, A22, B22, C, lda, ldb, ldc, mr, kr, nr);
    
    // Reste a traiter les cas oddsized
  }
*/
  
}
template <>
void FFLAS::ClassicMatmul(const Modular<double>& F,  
			  const enum FFLAS_TRANSPOSE ta,
			  const enum FFLAS_TRANSPOSE tb,
			  const size_t m, const size_t n,const size_t k,
			  double ALPHA,
			  const double * A, 
			  const size_t lda,
			  const double * B, 
			  const size_t ldb,
			  double BETA,
			  double* C, const size_t ldc){
	static  double Mone;
	F.init(Mone, -1);
	//size_t dlda,dldb;
	double* Ci=C;
	//cerr<<"PASSE EN MODULA DOUBLE"<<endl;
	// Call to the blas Multiplication 
	cblas_dgemm(CblasRowMajor, (enum CBLAS_TRANSPOSE) ta, 
		    (enum CBLAS_TRANSPOSE) tb, m, n, k, ALPHA,
		    A, lda, B, ldb,  BETA, C, ldc);
	
// 	if ( !F.areEqual(ALPHA, Mone) && !F.isOne(ALPHA)){
// 		// Cd may contain approximations of C entries:
// 		// Perturbation of Cd to get C from truncation.
// 		for (; Ci<C+m*ldc; Ci+=ldc)
// 			for(size_t j=0; j<lda; ++j){
// 				if (*(Ci+j)>=0){
// 					*(Ci+j) += 0.1;
// 				}
// 				else
// 					*(Ci+j) -= 0.1;
// 				F.init( *(Ci+j), *(Ci+j));
// 			}
// 	}
	size_t i=0, j;
	for (Ci=C ; i<m;++i, Ci+=ldc)
		for ( j=0; j<n;++j){
			F.init(*(Ci+j),*(Ci+j));
			//			*(Ci+j) = (*(Ci+j));
		}

    /*}
 else{
    // NB reste a refaire le systeme de ALPHA BETA qui ne permet pas de faire C=A*B
    // en passant ici
    cerr<<"block computation over the field"<<endl;
    size_t nr = n/2;
    size_t kr = k/2;
    size_t mr = m/2;
    double * A12 = A+kr;
    double * A21 = A+mr*lda;
    double * A22 = A+kr+mr*lda;
    double * B12 = B+nr;
    double * B21 = B+kr*ldb;
    double * B22 = B+nr+kr*ldb;
    double * C12 = B+nr;
    double * C21 = B+mr*ldc;
    double * C22 = B+nr+mr*ldc;

    // C11 = A11.B11+A12.B21
    ClassicMatMul(F, A, B, C, lda, ldb, ldc, mr, kr, nr);
    ClassicMatMul(F, A12, B21, C, lda, ldb, ldc, mr, kr, nr);
    // C12 = A11.B12+A12.B22
    ClassicMatMul(F, A, B12, C12, lda, ldb, ldc, mr, kr, nr);
    ClassicMatMul(F, A12, B22, C12, lda, ldb, ldc, mr, kr, nr);
    // C21 = A21.B11+A22.B21
    ClassicMatMul(F, A21, B, C21, lda, ldb, ldc, mr, kr, nr);
    ClassicMatMul(F, A22, B21, C21, lda, ldb, ldc, mr, kr, nr);
    // C22 = A21.B12+A22.B22
    ClassicMatMul(F, A21, B12, C, lda, ldb, ldc, mr, kr, nr);
    ClassicMatMul(F, A22, B22, C, lda, ldb, ldc, mr, kr, nr);
    
    // Reste a traiter les cas oddsized
  }
*/
  
}


// Classic Multiplication over double
template <>
inline  void FFLAS::ClassicMatmul(const DoubleDomain& F, 
				 const enum FFLAS_TRANSPOSE ta,
				 const enum FFLAS_TRANSPOSE tb,
				 const size_t m, const size_t n,const size_t k,
				 DoubleDomain::Element ALPHA,
				 const DoubleDomain::Element * Ad,
				 const size_t lda,
				 const DoubleDomain::Element * Bd,
				 const size_t ldb,
				 DoubleDomain::Element BETA,
				 DoubleDomain::Element * Cd, const size_t ldc){
	
	// Call to the blas multiplication
	cblas_dgemm(CblasRowMajor, (enum CBLAS_TRANSPOSE) ta, 
		    (enum CBLAS_TRANSPOSE) tb,
		    m, n, k, (DoubleDomain::Element) ALPHA,
		    Ad, lda, Bd, ldb, (DoubleDomain::Element) BETA,Cd, ldc);
}



// Winograd Multiplication  A(n*k) * B(k*m) in C(n*m)

// Computation of the 22 Winograd's operations
template<class Field>
inline  void FFLAS::WinoCalc( const Field& F, 
			      const size_t mr, const size_t nr, const size_t kr,
			      const typename Field::Element* A,const size_t lda,
			      const typename Field::Element* B,const size_t ldb,
			      typename Field::Element * C, const size_t ldc,
			      typename Field::Element ** t_X1,
			      typename Field::Element ** t_X2,
			      typename Field::Element ** t_X3,
			      long long kmax, size_t winostep){
	size_t i,j;
	typename Field::Element tmp;
	const typename Field::Element* d11,*d12,*d21,*d22;
	typename Field::Element* d11c,*d12c,*d21c,*d22c,*dx1,*dx2,*dx3;
	// Three temporary submatrices are required
	typename Field::Element* X1 = t_X1[winostep];
	typename Field::Element* X2 = t_X2[winostep];
	typename Field::Element* X3 = t_X3[winostep];
	
	// P2 = A12 * B21  in C11
	WinoMain( F, mr, nr, kr, A+kr, lda, (B+kr*ldb), ldb, C, ldc, t_X1, t_X2, t_X3,
		  kmax,winostep-1);
  
	// P1 = A11 * B11 in C12
	WinoMain( F, mr, nr, kr, A, lda, B, ldb, C+nr, ldc, t_X1,t_X2,t_X3,
		  kmax,winostep-1);
	
	// T3 = B22 - B12 in X3
	d12 = B + nr; d22 = d12 + kr*ldb; dx3 = X3;
	for ( i=0; i<kr; ++i, d12+=ldb, d22+=ldb, dx3+=nr){
		for (j=0;j<nr;++j)
			F.sub(*(dx3+j),*(d22 + j),*(d12 + j));
		
	}
	// S3 = A11 - A21 in X2 and  U1 = P2 + P1 in C11
	d11 = A; d21 = d11 + mr*lda; dx2 = X2; d11c = C;d12c = d11c + nr; 
	for (size_t i=0;i<mr;++i,d11+=lda,d11c+=ldc,d12c+=ldc, d21+=lda,dx2+=kr){
		for (size_t j=0;j<kr;++j)
			F.sub(*(dx2+j),*(d11 + j),*(d21 + j));
		for (size_t j=0;j<nr;++j)
			F.addin(*(d11c + j),*(d12c + j));
	}
	
	// P7 = S3*T3 in C21
	WinoMain( F, mr, nr, kr, X2, kr, X3, nr, C+mr*ldc, ldc, t_X1,t_X2,t_X3,
		  kmax,winostep-1);
	
	// S1 = A21 + A22 in X2
	d21 = A + mr*lda; d22 = d21 + kr; dx2 = X2;
	for (size_t i=0;i<mr;++i, d21+=lda, d22+=lda, dx2+=kr){
		for (size_t j=0;j<kr;++j)
			F.add(*(dx2+j),*( d21 + j),*(d22 + j));
	}
	
	// T1 = B12 - B11 in X3
	d11 = B; d12 = d11 + nr; dx3 = X3;
	for (size_t i=0;i<kr;++i, d11+=ldb, d12+=ldb,dx3+=nr){
		for (size_t j=0;j<nr;++j)
			F.sub(*(dx3+j),*(d12 + j),*(d11 + j));
	}
	
	// P5 = S1*T1 in C22
	WinoMain( F, mr, nr, kr, X2, kr, X3, nr, C+mr*ldc + nr, ldc, t_X1, t_X2, t_X3,
		  kmax,winostep-1);
	
	// T2 = B22 - T1 in X3
	d22 = B+kr*ldb+nr; dx3 = X3;
	for (size_t i=0;i<kr;++i,d22+=ldb,dx3+=nr){
		for (size_t j=0;j<nr;++j)
			F.sub(*(dx3+j),*(d22 + j),*(dx3+j));
	}
	
	// S2 = S1 - A11 in X2
	d11 = A; dx2 = X2;
	for (size_t i=0;i<mr;++i,d11+=lda,dx2+=kr){
		for (size_t j=0;j<kr;++j)
			F.subin(*(dx2+j),*(d11 + j));
	}
	
	// P6 = S2*T2 in X1
	WinoMain( F, mr, nr, kr, X2, kr, X3, nr, X1, nr, t_X1, t_X2, t_X3,
		  kmax,winostep-1);
	
	// T4 = T2 - B21 in X3
	d21 = B+kr*ldb;dx3=X3;
	for (size_t i=0;i<kr;++i,d21+=ldb,dx3+=nr){
		for (size_t j=0;j<nr;++j)
			F.subin(*(dx3+j),*( d21 + j));
	}
	
	// S4 = A12 -S2 in X2 and U2 = P1 + P6 in tmp and U3 = P7 + U2 in C21
	// and U4 = U2 + P5 in C12
	d12 = A+kr; dx2 = X2;d12c = C+nr;dx1=X1;d21c = C+mr*ldc; d22c = d21c+nr;
	for (size_t i=0;i<mr;++i, d12+=lda, dx2+=kr,d12c+=ldc,dx1+=nr,d21c+=ldc,d22c+=ldc){
		for (size_t j=0;j<kr;++j)
			F.sub(*(dx2+j),*(d12 + j),*(dx2+j));
		for (size_t j=0;j<nr;++j){
			F.add(tmp,*(d12c + j),*(dx1+j));// temporary U2
			F.addin(*(d21c + j),tmp);       // U3
			F.add(*(d12c + j),tmp,*(d22c +j));  // U4
		}
	}
	
	// P3 = S4*B22 in X1
	WinoMain( F, mr, nr, kr, X2, kr, B + kr*ldb + nr, ldb, X1, nr,
		  t_X1,t_X2,t_X3,
		  kmax,winostep-1);
	
	// U5 = U4 + P3 in C12
	d12c = C+nr;dx1=X1;
	for (size_t i=0;i<mr;++i, d12c+=ldc, dx1+=nr){
		for (size_t j=0;j<nr;++j)
			F.addin(*(d12c + j),*(dx1+j));
	}
	
	// P4 = A22*T4 in X1
	WinoMain( F, mr, nr, kr, A + mr*lda + kr, lda, X3, nr, X1, nr,
		  t_X1, t_X2, t_X3, kmax,winostep-1);
	
	// U7 = P5 + U3 in C22
	d21c = C+mr*ldc; d22c = d21c +nr;
	for (size_t i=0;i<mr;++i,d21c+=ldc,d22c+=ldc){
		for (size_t j=0;j<nr;++j)
			F.addin(*(d22c + j),*(d21c +j));
	}
	
	// U6 = U3 - P4 in C21
	d21c = C+mr*ldc; dx1=X1;
	for (size_t i=0;i<mr;++i, d21c+=ldc, dx1+=nr){
		for (size_t j=0;j<nr;++j)
			F.subin(*(d21c + j),*(dx1+j));
	}
	
}

// Control of the 2 cut-off criterias determining when to switch from Winograd's 
// to classic multiplication, and from finite field to double.
// Fix-up for odd-sized matrices using dynamic pealing ( coming soon...)
// for matrices over a finite Field
template <class Field>
inline  void FFLAS::WinoMain( const Field& F, 
			     const size_t m, const size_t n, const size_t k,
			     const typename Field::Element* A,const size_t lda,
			     const typename Field::Element* B,const size_t ldb,
			     typename Field::Element * C, const size_t ldc,
			     typename Field::Element ** t_X1,
			     typename Field::Element ** t_X2,
			     typename Field::Element ** t_X3,
			     long long kmax, size_t winostep){
	static typename Field::Element one,zero;
	F.init(one, 1);
	F.init(zero, 0);
	if (winostep<=0) // Winograd -> Classic
		ClassicMatmul( F, FflasNoTrans, FflasNoTrans, m, n, k,
			       one, A, lda, B, ldb, zero, C, ldc );
	else{
		size_t nr = n/2;
		size_t kr = k/2;
		size_t mr = m/2;
		if (kr < kmax){ // switch on double
			// Temporary double matrices
			DoubleDomain::Element* td_X1[winostep+1];
			DoubleDomain::Element* td_X2[winostep+1];
			DoubleDomain::Element* td_X3[winostep+1];
			for (size_t i=winostep;(i>=1);--i){
				td_X1[i] = new DoubleDomain::Element[mr*nr];
				td_X2[i] = new DoubleDomain::Element[mr*kr];
				td_X3[i] = new DoubleDomain::Element[kr*nr];
				nr/=2;
				mr/=2;
				kr/=2;
			}
			DoubleDomain::Element * Ad = new DoubleDomain::Element[m*k];
			DoubleDomain::Element * Bd = new DoubleDomain::Element[k*n];
			DoubleDomain::Element * Cd = new DoubleDomain::Element[m*n];
			// Conversion GFq => double
			MatF2MatD( F, Ad, A, lda, m, k);
			MatF2MatD( F, Bd, B, ldb, k, n); 
			//if (F.isnzero(BETA)){
			//	MatF2MatD( F, Cd, C, ldc, m, n ); 
			//}
			
			// recursive call
			WinoMain( DoubleDomain(), m, n, k, 
				  Ad, k, Bd, n, Cd, n, td_X1, td_X2, td_X3, 
				  kmax, winostep);
			// Conversion double => GFq
			MatD2MatF( F, C, ldc, Cd, m, n );
			// Temporary double matrices destruction
			delete[] Ad;
			delete[] Bd;
			delete[] Cd;
			for (size_t i=winostep;i>=1;--i){
				delete[] td_X1[i];
				delete[] td_X2[i];
				delete[] td_X3[i];
			}
			
		}
		else{
			WinoCalc( F, mr, nr, kr, A, lda, B, ldb, C, ldc,
				  t_X1, t_X2, t_X3, kmax,winostep);

			// Fix-up for odd-sized matrices  should be here
		}
	}
}


// Control of the cut-off criteria determing when to switch 
// to classic multiplication
// Fix-up for odd-sized matrices using dynamic pealing
// for matrices over double
template<>
inline  void FFLAS::WinoMain( const DoubleDomain& D, 
			      const size_t m, const size_t n, const size_t k,
			      const DoubleDomain::Element * A, const size_t lda,
			      const DoubleDomain::Element * B, const size_t ldb,
			      DoubleDomain::Element * C, const size_t ldc,
			      DoubleDomain::Element **t_X1,
			      DoubleDomain::Element **t_X2,
			      DoubleDomain::Element **t_X3,
			      long long kmax, size_t winostep){
	
	if (winostep<=0) // switch Winograd => Classic
		ClassicMatmul( D, FflasNoTrans, FflasNoTrans, m, n, k,
			       1.0, A, lda, B, ldb, 0.0, C, ldc);
	else{
		size_t nr = n/2;
		size_t kr = k/2;
		size_t mr = m/2;
		size_t mkn = (n & 0x1)+( (k & 0x1)<<1 )+ ( (m & 0x1)<<2 ); 
		// Waiting for a winograd implementation with alpha and beta
		const double ALPHAd = 1.0;
		const double BETAd = 0.0;
		WinoCalc( D, mr, nr, kr, A, lda, B, ldb, C, ldc,
			  t_X1, t_X2, t_X3,
			  kmax,winostep); 

		switch(mkn) { 
		case 1: // n oddsized
			cblas_dgemv( CblasRowMajor, CblasNoTrans, m, k,
				     (DoubleDomain::Element) ALPHAd, A, lda, 
				     B+n-1, ldb, 
				     (DoubleDomain::Element) BETAd, C+n-1,ldc);
			break;
      
		case 2: // k oddsized
			cblas_dger( CblasRowMajor, m, n,
				    (DoubleDomain::Element) ALPHAd,
				    A+k-1, lda, B+(k-1)*ldb, 1, C, ldc);
			break;
			
		case 3: // n, k oddsized
			cblas_dger( CblasRowMajor, m, n-1,
				    (DoubleDomain::Element) ALPHAd,
				    A+k-1, lda, B+(k-1)*ldb, 1, C, ldc);
			if (BETAd){
				cblas_daxpy( m, ALPHAd*BETAd*(*(B+(k-1)*ldb+n-1)), 
					     A+k-1, lda, C+n-1, ldc);
				cblas_dgemv( CblasRowMajor, CblasNoTrans, m, k-1,
					     (DoubleDomain::Element) ALPHAd, A, lda,
					     B+n-1, ldb, 
					     (DoubleDomain::Element) BETAd, C+n-1,ldc);
			}
			else{
				cblas_dcopy( m, A + k-1, lda, C+n-1, ldc);
				cblas_dscal( m, ALPHAd*(*(B+(k-1)*ldb+n-1)),
					     C + n-1, ldc);
				cblas_dgemv( CblasRowMajor, CblasNoTrans, m, k-1,
					     (DoubleDomain::Element) ALPHAd, 
					     A, lda, B+n-1, ldb, 1.0, C+n-1,ldc);
			}
			break;
			
		case 4: // m oddsized
			cblas_dgemv( CblasRowMajor, CblasTrans, k, n,
				     (DoubleDomain::Element) ALPHAd,
				     B, ldb, A+(m-1)*lda, 1,
				     (DoubleDomain::Element) BETAd, C+(m-1)*ldc, 1);
			break;
			
		case 5: // m, n oddsized
			cblas_dgemv( CblasRowMajor, CblasNoTrans, m-1, k,
				     (DoubleDomain::Element) ALPHAd,
				     A, lda, B+n-1, ldb,
				     (DoubleDomain::Element) BETAd, C+n-1, ldc);
			cblas_dgemv( CblasRowMajor, CblasTrans, k, n-1,
				     (DoubleDomain::Element) ALPHAd, 
				     B, ldb, A+(m-1)*lda, 1,
				     (DoubleDomain::Element) BETAd, C+(m-1)*ldc, 1);
			*(C+(m-1)*ldc+n-1) = ALPHAd*cblas_ddot(k, A+(m-1)*lda, 1,
							      B+n-1, ldb)
				+ BETAd*(*(C+(m-1)*ldc+n-1));
			break;
      
		case 6: // m, k oddsized
			cblas_dger( CblasRowMajor, m-1, n, 
				    (DoubleDomain::Element) ALPHAd, A+k-1, lda,
				    B+(k-1)*ldb, 1, C, ldc);
			if (BETAd){
				cblas_daxpy( n, ALPHAd*BETAd*(*(A+(m-1)*lda+k-1)),
					     B+(k-1)*ldb, 1, C+(m-1)*ldc, 1);
				cblas_dgemv( CblasRowMajor, CblasTrans, k-1, n,
					     (DoubleDomain::Element) ALPHAd,
					     B, ldb, A+(m-1)*lda, 1,
					     (DoubleDomain::Element) BETAd, 
					     C+(m-1)*ldc, 1);
			}
			else{
				cblas_dcopy( n, B + (k-1)*ldb, 1, C+(m-1)*ldc, 1);
				cblas_dscal( n, ALPHAd*(*(A+(m-1)*lda+k-1)),
					     C + (m-1)*ldc, 1);
				cblas_dgemv( CblasRowMajor, CblasTrans, k-1, n,
					    (DoubleDomain::Element) ALPHAd,
					     B, ldb, A+(m-1)*lda, 1,
					     1.0, C+(m-1)*ldc, 1);
			}
			break;
      
		case 7: // m, k, n oddsized
			// Block NW
			cblas_dger( CblasRowMajor, m-1, n-1,
				   (DoubleDomain::Element) ALPHAd, 
				    A+k-1, lda, B+(k-1)*ldb, 1, C, ldc);
			if (BETAd){
				// Block SW
				cblas_daxpy( n-1, ALPHAd*BETAd*(*(A+(m-1)*lda+k-1)), 
					     B+(k-1)*ldb, 1, C+(m-1)*ldc, 1);
				cblas_dgemv( CblasRowMajor, CblasTrans, k-1, n-1,
					     (DoubleDomain::Element) ALPHAd,
					     B, ldb, A+(m-1)*lda, 1,
					     (DoubleDomain::Element) BETAd,
					     C+(m-1)*ldc, 1);
				// Block NE
				cblas_daxpy( m-1, ALPHAd*BETAd*(*(B+(k-1)*ldb+n-1)), 
					     A+k-1, lda, C+n-1, ldc);
				cblas_dgemv( CblasRowMajor, CblasTrans, m-1, k-1,
					     (DoubleDomain::Element) ALPHAd,
					     A, lda, B+n-1, ldb,
					     (DoubleDomain::Element) BETAd,
					     C+n-1, ldc);
			}
			else{
				// Block SW
				cblas_dcopy( n-1, B + (k-1)*ldb, 1, C+(m-1)*ldc, 1);
				cblas_dscal( n-1, ALPHAd*(*(A+(m-1)*lda+k-1)), 
					     C + (m-1)*ldc, 1);
				cblas_dgemv( CblasRowMajor, CblasTrans, k-1, n-1,
					     (DoubleDomain::Element) ALPHAd,
					     B, ldb, A+(m-1)*lda, 1,
					     1.0, C+(m-1)*ldc, 1);
      
				// Block NE
				
				cblas_dcopy( m-1, A + k-1, lda, C+n-1, ldc);
				cblas_dscal( m-1, ALPHAd*(*(B+(k-1)*ldb+n-1)), 
					     C + n-1, ldc);
				cblas_dgemv( CblasRowMajor, CblasNoTrans, m-1, k-1,
					     (DoubleDomain::Element) ALPHAd,
					     A, lda, B+n-1, ldb,
					     1.0, C+n-1, ldc);
			}

			// Block SE
			*(C+(m-1)*ldc+n-1) = ALPHAd*cblas_ddot( k, A+(m-1)*lda, 1, 
							       B+n-1, ldb)
				+ BETAd*(*(C+(m-1)*ldc+n-1));
			break;
		}
		
	}
}

//-------------------------------------------------------------------
// visible function: operates C = A*B over double or a finite Field
//-------------------------------------------------------------------
template<>
inline  DoubleDomain::Element* 
FFLAS::fgemm( const DoubleDomain& D,
	      const enum FFLAS_TRANSPOSE ta,
	      const enum FFLAS_TRANSPOSE tb,
	      const size_t m, const size_t n, const size_t k,
	      const DoubleDomain::Element alpha,
	      const DoubleDomain::Element* A, const size_t lda,
	      const DoubleDomain::Element* B, const size_t ldb,
	      const DoubleDomain::Element beta,
	      DoubleDomain::Element* C, const size_t ldc,
	      const size_t winostep){
	size_t mr=m/2;
	size_t kr=k/2;
	size_t nr=n/2;
	DoubleDomain::Element* t_X1[winostep+1];
	DoubleDomain::Element* t_X2[winostep+1];
	DoubleDomain::Element* t_X3[winostep+1];
	// Temporary space allocation
	for (size_t i=winostep;(i>=1);--i){
		t_X1[i] = new DoubleDomain::Element[mr*nr];
		t_X2[i] = new DoubleDomain::Element[mr*kr];
		t_X3[i] = new DoubleDomain::Element[kr*nr];
		nr/=2;
		mr/=2;
		kr/=2;
	}
	if (winostep && (beta == 0.0) && (alpha == 1.0) && 
	    (ta == FflasNoTrans) && (tb == FflasNoTrans) ){
		WinoMain( D, m, n, k, A, lda, B, ldb, C, ldc, t_X1, t_X2, t_X3,
			  0, winostep);
	}
	else{
		ClassicMatmul( D, ta, tb, m, n, k, alpha, A, lda, B, ldb, beta, C,ldc);
	}
	kr=k/2;
	for (size_t i=winostep;(i>=1);--i){
		delete[] t_X1[i];
		delete[] t_X2[i];
		delete[] t_X3[i];
		kr/=2;
	}
	return C;
}

template<class Field>
inline  typename Field::Element* 
FFLAS::fgemm( const Field& F,
	      const enum FFLAS_TRANSPOSE ta,
	      const enum FFLAS_TRANSPOSE tb,
	      const size_t m, const size_t n, const size_t k,
	      const typename Field::Element alpha,
	      const typename Field::Element* A, const size_t lda,
	      const typename Field::Element* B, const size_t ldb, 
	      const typename Field::Element beta,
	      typename Field::Element* C, const size_t ldc,
	      const size_t winostep){

	if (!(m &&n && k)) return C;
	if ( !F.isZero(beta) || !F.isOne(alpha) ||
	     !winostep || ta==FflasTrans || tb==FflasTrans){
		ClassicMatmul( F, ta, tb,  m, n, k, alpha, A, lda, B, ldb,
			       beta, C,ldc);
	}
	else{
		typedef typename Field::Element elt;
		size_t mr=m/2;
		size_t kr=k/2;
		size_t nr=n/2;
		integer charac; //=F.characteristic();
		F.characteristic(charac);
		long long c = charac;
		// Threshold between GFq and double
		long long kmax = ((long long)1<<53)/(c*c);
		elt* t_X1[winostep+1];
		elt* t_X2[winostep+1];
		elt* t_X3[winostep+1];
		// Temporary space allocation
		for (size_t i=winostep;(i>=1)&&(kr>kmax);--i){
			t_X1[i] = new elt[mr*nr];
			t_X2[i] = new elt[mr*kr];
			t_X3[i] = new elt[kr*nr];
			nr/=2;
			mr/=2;
			kr/=2;
		}
		WinoMain(F, m, k, n, A, lda, B, ldb, C, ldc, t_X1, t_X2, t_X3, 
			 kmax, winostep);
		kr=k/2;
		for (size_t i=winostep;(i>=1)&&(kr>kmax);--i){
			delete[] t_X1[i];
			delete[] t_X2[i];
			delete[] t_X3[i];
			kr/=2;
		}
	}
	return C;
}

template<class Field>
inline  typename Field::Element*
FFLAS::fsquare( const Field& F,
		const enum FFLAS_TRANSPOSE ta,
		const size_t n, const typename Field::Element alpha,
		const typename Field::Element* A, const size_t lda,
		const typename Field::Element beta,
		typename Field::Element* C, const size_t ldc){
	double ALPHAd, BETAd;
	F.convert( ALPHAd, alpha );
	F.convert( BETAd, beta );
	// Double  matrices initialisation
	DoubleDomain::Element * Ad = new DoubleDomain::Element[n*n];
	DoubleDomain::Element * Cd = new DoubleDomain::Element[n*n];
	// Conversion finite Field => double
	MatF2MatD( F, Ad, A, lda, n, n);
	if (!F.isZero(beta))
		MatF2MatD( F, Cd, C, ldc, n, n); 
    
	// Call to the blas Multiplication 
	cblas_dgemm( CblasRowMajor, (enum CBLAS_TRANSPOSE)ta,
		     (enum CBLAS_TRANSPOSE)ta, n, n, n, 
		     (DoubleDomain::Element) ALPHAd, Ad, n, Ad, n,
		     (DoubleDomain::Element) BETAd, Cd, n);
	// Conversnion double => Finite Field
	delete[] Ad;
	MatD2MatF( F, C, ldc, Cd, n, n);
	delete[] Cd;
}

