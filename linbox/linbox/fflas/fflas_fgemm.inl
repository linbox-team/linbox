/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/fflas/fflas_fgemm.inl
 * Copyright (C) 2003 Clement Pernet
 *
 * Written by Clement Pernet <Clement.Pernet@imag.fr>
 *
 * Warning,
 * TRANSPOSE option will force classic matrix multiplication
 *
 * See COPYING for license information.
 */
#ifndef MAX
#define MAX(a,b) (a<b)?b:a
#endif

#include <linbox/field/modular-double.h>

// Classic multiplication over a finite Field
template <class Field>
void LinBox::FFLAS::ClassicMatmul( const Field& F,  
				   const enum FFLAS_TRANSPOSE ta,
				   const enum FFLAS_TRANSPOSE tb,
				   const size_t m, const size_t n,const size_t k,
				   const typename Field::Element alpha,
				   const typename Field::Element * A, const size_t lda,
				   const typename Field::Element * B, const size_t ldb,
				   const typename Field::Element beta,
				   typename Field::Element* C, const size_t ldc,
				   const size_t kmax ){
	static  typename Field::Element Mone;
	static  typename Field::Element one;
	F.init(Mone, -1);
	F.init(one, 1);
	
	if ( k < kmax ){
		size_t dlda,dldb;
		
		DoubleDomain::Element alphad, betad;
		
		DoubleDomain::Element * Add = new DoubleDomain::Element[m*k];
		DoubleDomain::Element * Bdd = new DoubleDomain::Element[k*n];
		DoubleDomain::Element * Cd = new DoubleDomain::Element[m*n];
		
#ifdef CHECK_MEMORY
		if (CUR_MEMORY+ (m*k+k*n+m*n)*sizeof(DoubleDomain::Element) > MAX_MEMORY)
			MAX_MEMORY=CUR_MEMORY+ (m*k+k*n+m*n)*sizeof(DoubleDomain::Element);
#endif
		// Conversion finite Field => double
		if (ta == FflasTrans){
			dlda = m;
			MatF2MatD( F, Add, A, lda, k, m );
		} else{
			dlda = k;
			MatF2MatD( F, Add, A, lda, m, k );
	}
		
		if (tb == FflasTrans){
			dldb = k;
			MatF2MatD( F, Bdd, B, ldb, n, k );
		} else{
			dldb = n;
			MatF2MatD( F, Bdd, B, ldb, k, n );
		}
		
		if (!F.isZero(beta)){
			MatF2MatD( F, Cd, C, ldc, m, n );
			if (F.areEqual(beta, Mone))
				betad = -1.0;
			else
				F.convert( betad, beta );
		}
		if (F.areEqual(alpha, Mone))
			alphad = -1.0;
		else
			F.convert( alphad, alpha );
		
		ClassicMatmul( DoubleDomain(), ta, tb, m, n, k, alphad, 
			       Add, dlda, Bdd, dldb, betad, Cd, n, kmax);
		// Conversion double => Finite Field
		delete[] Add;
		delete[] Bdd;

		if ( !F.areEqual(alpha, Mone) && !F.isOne(alpha)){
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
	}
	else{
		// NB reste a refaire le systeme de ALPHA BETA qui ne permet pas de faire C=A*B
		// en passant ici
		//std::cerr<<"block computation over the field"<<std::endl;
		size_t nr = n/2;
		size_t kr = k/2;
		size_t mr = m/2;
		const typename Field::Element * A12 = A+kr;
		const typename Field::Element * A21 = A+mr*lda;
		const typename Field::Element * A22 = A+kr+mr*lda;
		const typename Field::Element * B12 = B+nr;
		const typename Field::Element * B21 = B+kr*ldb;
		const typename Field::Element * B22 = B+nr+kr*ldb;
		typename Field::Element * C12 = C+nr;
		typename Field::Element * C21 = C+mr*ldc;
		typename Field::Element * C22 = C+nr+mr*ldc;
		
		// C11 = A11.B11+A12.B21
		ClassicMatmul(F, ta, tb, mr, nr, kr, alpha, A, lda,  B, ldb, beta, C, ldc, kmax );
		ClassicMatmul(F, ta, tb, mr, nr, kr, alpha, A12, lda, B21, ldb, one, C, ldc, kmax );
		// C12 = A11.B12+A12.B22
		ClassicMatmul(F, ta, tb, mr, nr, kr, alpha, A, lda,  B12, ldb, beta, C12, ldc, kmax );
		ClassicMatmul(F, ta, tb, mr, nr, kr, alpha, A12, lda, B22, ldb, one, C12, ldc, kmax );
		// C21 = A21.B11+A22.B21
		ClassicMatmul(F, ta, tb, mr, nr, kr, alpha, A21, lda,  B, ldb, beta, C21, ldc, kmax );
		ClassicMatmul(F, ta, tb, mr, nr, kr, alpha, A22, lda, B21, ldb, one, C21, ldc, kmax );
		// C22 = A21.B12+A22.B22
		ClassicMatmul(F, ta, tb, mr, nr, kr, alpha, A21, lda,  B12, ldb, beta, C22, ldc, kmax );
		ClassicMatmul(F, ta, tb, mr, nr, kr, alpha, A22, lda, B22, ldb, one, C22, ldc, kmax );
		
		// Reste a traiter les cas oddsized	
		DynamicPealing( F, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc, kmax );
	}
}

template <>
void LinBox::FFLAS::ClassicMatmul( const Modular<double>& F,  
				   const enum FFLAS_TRANSPOSE ta,
				   const enum FFLAS_TRANSPOSE tb,
				   const size_t m, const size_t n,const size_t k,
				   const double alpha, 
				   const double * A, const size_t lda,
				   const double * B, const size_t ldb,
				   const double beta,
				   double* C, const size_t ldc, 
				   const size_t kmax ){
	static  double Mone, one;
	F.init(Mone, -1);
	F.init(one, 1);

	if ( k < kmax ){
		double _alpha, _beta;
		double* Ci=C;

		if (F.areEqual( Mone, alpha )){
			_alpha = -1.0;
			F.convert( _beta, beta );
		}
		else{
			if (! F.areEqual( one, alpha) ){
				// Compute C = A*B + beta/alpha.C
				// and after C *= alpha
				F.div( _beta, beta, alpha );
			}
			else
				_beta = beta;
			_alpha = 1.0;
		}
		
		cblas_dgemm(CblasRowMajor, (enum CBLAS_TRANSPOSE) ta, 
			    (enum CBLAS_TRANSPOSE) tb, m, n, k, _alpha,
			    A, lda, B, ldb,  _beta, C, ldc);
		
		size_t i=0, j;
		for (Ci=C ; i<m;++i, Ci+=ldc)
			for ( j=0; j<n;++j)
				F.init(*(Ci+j),*(Ci+j));

		if ( !F.areEqual( one, alpha ) && !F.areEqual( Mone, alpha) ){
			for (double* Ci=C; Ci<C+m*ldc; Ci+=ldc)
				for ( size_t j=0; j<n; ++j ) 
					F.mulin( *( Ci + j ), alpha );
		}
		
	}
	else{
		//std::cerr<<"block computation over the field"<<std::endl;
		size_t nr = n/2;
		size_t kr = k/2;
		size_t mr = m/2;
		const double * A12 = A+kr;
		const double * A21 = A+mr*lda;
		const double * A22 = A+kr+mr*lda;
		const double * B12 = B+nr;
		const double * B21 = B+kr*ldb;
		const double * B22 = B+nr+kr*ldb;
		double * C12 = C+nr;
		double * C21 = C+mr*ldc;
		double * C22 = C+nr+mr*ldc;
		
		// C11 = A11.B11+A12.B21
		ClassicMatmul(F, ta, tb, mr, nr, kr, alpha, A, lda, B, ldb, beta, C, ldc, kmax );
		ClassicMatmul(F, ta, tb, mr, nr, kr, alpha, A12, lda, B21, ldb, one, C, ldc, kmax );
		// C12 = A11.B12+A12.B22
		ClassicMatmul(F, ta, tb, mr, nr, kr, alpha, A, lda, B12, ldb, beta, C12, ldc, kmax );
		ClassicMatmul(F, ta, tb, mr, nr, kr, alpha, A12, lda, B22, ldb, one, C12, ldc, kmax );
		// C21 = A21.B11+A22.B21
		ClassicMatmul(F, ta, tb, mr, nr, kr, alpha, A21, lda, B, ldb, beta, C21, ldc, kmax );
		ClassicMatmul(F, ta, tb, mr, nr, kr, alpha, A22, lda, B21, ldb, one, C21, ldc, kmax );
		// C22 = A21.B12+A22.B22
		ClassicMatmul(F, ta, tb, mr, nr, kr, alpha, A21, lda, B12, ldb, beta, C22, ldc, kmax );
		ClassicMatmul(F, ta, tb, mr, nr, kr, alpha, A22, lda, B22, ldb, one, C22, ldc, kmax );
		
		// Updates for oddsized dimensions
		DynamicPealing( F, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc, kmax );
			
	}
	
  
}


// Classic Multiplication over double
template <>
inline  void LinBox::FFLAS::ClassicMatmul(const DoubleDomain& F, 
					  const enum FFLAS_TRANSPOSE ta,
					  const enum FFLAS_TRANSPOSE tb,
					  const size_t m, const size_t n,const size_t k,
					  const DoubleDomain::Element alpha,
					  const DoubleDomain::Element * Ad,
					  const size_t lda,
					  const DoubleDomain::Element * Bd,
					  const size_t ldb,
					  const DoubleDomain::Element beta,
					  DoubleDomain::Element * Cd, const size_t ldc, 
					  const size_t kmax){
	
	// Call to the blas multiplication
	cblas_dgemm(CblasRowMajor, (enum CBLAS_TRANSPOSE) ta, 
		    (enum CBLAS_TRANSPOSE) tb,
		    m, n, k, (DoubleDomain::Element) alpha,
		    Ad, lda, Bd, ldb, (DoubleDomain::Element) beta,Cd, ldc);
}



// Winograd Multiplication  A(n*k) * B(k*m) in C(n*m)

// Computation of the 22 Winograd's operations
template<class Field>
inline  void LinBox::FFLAS::WinoCalc( const Field& F, 
				      const size_t mr, const size_t nr, const size_t kr,
				      const typename Field::Element alpha,
				      const typename Field::Element* A,const size_t lda,
				      const typename Field::Element* B,const size_t ldb,
				      const typename Field::Element beta,
				      typename Field::Element * C, const size_t ldc,
				      long long kmax, size_t winostep){
	static typename Field::Element zero, one;
	F.init ( zero, 0 );
	F.init ( one, 1 );
	typename Field::Element mbeta;
	F.neg(mbeta,beta);
	size_t i,j;
	typename Field::Element tmpU2, tmpU3;
	const typename Field::Element* d11,*d12,*d21,*d22;
	typename Field::Element* d11c,*d12c,*d21c,*d22c,*dx1,*dx2,*dx3;
	// Three temporary submatrices are required
	typename Field::Element* X1 = new typename Field::Element[mr*nr];//winostep];
	typename Field::Element* X2 = new typename Field::Element[mr*kr];//winostep];
	typename Field::Element* X3 = new typename Field::Element[(MAX(mr,kr))*nr];//winostep];
	
	// P2 = alpha . A12 * B21 + beta . C11  in C11
	WinoMain( F, mr, nr, kr, alpha, A+kr, lda, (B+kr*ldb), ldb, beta, C, ldc, kmax,winostep-1);
  
	// T3 = B22 - B12 in X3
	d12 = B + nr; d22 = d12 + kr*ldb; dx3 = X3;
	for ( i=0; i<kr; ++i, d12+=ldb, d22+=ldb, dx3+=nr){
		for (j=0;j<nr;++j)
			F.sub(*(dx3+j),*(d22 + j),*(d12 + j));
		
	}

	// S3 = A11 - A21 in X2 
	d11 = A; d21 = d11 + mr*lda; dx2 = X2; 
	for ( size_t i = 0; i < mr; ++i, d11 += lda, d21 += lda, dx2 += kr )
		for ( size_t j = 0;j < kr; ++j )
			F.sub( *(dx2+j), *(d11 + j), *(d21 + j) );

	// C22 = C22 - C12 if beta != 0
	if ( !F.isZero(beta) ){
		d12c = C + nr ;
		d22c = d12c + mr*ldc;
		for ( size_t i = 0; i< mr; ++i, d12c += ldc, d22c += ldc )
			for ( size_t j = 0; j < nr; ++j )
				F.subin( *(d22c + j), *(d12c + j) );

		// C21 = C21 - C22
		d21c = C + mr*ldc;
		d22c = d21c + nr ;
		for ( size_t i = 0; i< mr; ++i, d22c += ldc, d21c += ldc )
			for ( size_t j = 0; j < nr; ++j )
				F.subin( *(d21c + j), *(d22c + j) );
	}

	// P7 = alpha . S3 * T3 + beta . C22 in C22
	WinoMain( F, mr, nr, kr, alpha, X2, kr, X3, nr, beta, C+nr+mr*ldc, ldc, kmax,winostep-1);
	
	// T1 = B12 - B11 in X3
	d11 = B; d12 = d11 + nr; dx3 = X3;
	for ( size_t i = 0; i < kr; ++i, d11 += ldb, d12 += ldb, dx3 += nr ){
		for ( size_t j = 0; j < nr; ++j )
			F.sub( *(dx3 + j), *(d12 + j), *(d11 + j) );
	}

	// S1 = A21 + A22 in X2
	d21 = A + mr*lda; d22 = d21 + kr; dx2 = X2;
	for (size_t i=0;i<mr;++i, d21+=lda, d22+=lda, dx2+=kr){
		for (size_t j=0;j<kr;++j)
			F.add(*(dx2+j),*( d21 + j),*(d22 + j));
	}
	
	// P5 = alpha . S1*T1 + beta . C12 in C12
	WinoMain( F, mr, nr, kr, alpha, X2, kr, X3, nr, beta, C + nr, ldc, kmax,winostep-1);

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
	
	// P6 = alpha . S2 * T2 in X1
	WinoMain( F, mr, nr, kr, alpha, X2, kr, X3, nr, zero, X1, nr, kmax,winostep-1);
	
	// T4 = T2 - B21 in X3
	d21 = B+kr*ldb;dx3=X3;
	for (size_t i=0;i<kr;++i,d21+=ldb,dx3+=nr){
		for (size_t j=0;j<nr;++j)
			F.subin(*(dx3+j),*( d21 + j));
	}
	
	// S4 = A12 -S2 in X2 
	d12 = A+kr; dx2 = X2;
	for ( size_t i = 0;i < mr; ++i, d12 += lda, dx2 += kr ){
		for ( size_t j=0;j<kr;++j )
			F.sub( *(dx2+j), *(d12 + j), *(dx2+j) );
	}
	
	// P4 = alpha . A22 * T4 - beta . C21 in C21
	WinoMain( F, mr, nr, kr, alpha, A + mr*lda + kr, lda, X3, nr, mbeta, C + mr*ldc, ldc,
		  kmax,winostep-1);

	// P1 = alpha . A11 * B11 in X3
	WinoMain( F, mr, nr, kr, alpha, A, lda, B, ldb, zero, X3, nr, kmax,winostep-1);
	
	//  U1 = P2 + P1 in C11	
	d11c = C;dx3 = X3; 
	for ( size_t i = 0; i < mr; ++i, d11c += ldc, dx3 += nr )
		for ( size_t j = 0; j < nr; ++j )
			F.addin( *(d11c + j), *(dx3 + j) );
	
	// U2 = P1 + P6 in tmpU2  and
	// U3 = P7 + U2 in tmpU3  and 
	// U7 = P5 + U3 in C22    and
	// U4 = P5 + U2 in C12    and
	// U6 = U3 - P4 in C21    and
	d12c = C+nr; dx1=X1; dx3=X3; d21c = C+mr*ldc; d22c = d21c + nr; 
	for ( size_t i = 0; i < mr; 
	      ++i, d12c += ldc, dx1 += nr, dx3 += nr, d22c+=ldc, d21c += ldc ){
		for (size_t j=0;j<nr;++j){
			F.add( tmpU2, *(dx3 + j), *(dx1 + j));    // temporary U2 = P1 + P6
			F.add( tmpU3, tmpU2, *(d22c + j) );      // temporary U3 = U2 + P7
			F.add( *(d22c + j), *(d12c + j), tmpU3);  // U7 = P5 + U3 in C22
			F.addin( *(d12c + j), tmpU2);             // U4 = P5 + U2 in C12
			F.sub( *(d21c + j), tmpU3, *(d21c + j) ); // U6 = U3 - P4 in C21
		}
	}
	
	// U5 = P3 + U4 = alpha . S4*B22 + U4 in C12
	WinoMain( F, mr, nr, kr, alpha, X2, kr, B + kr*ldb + nr, ldb, one, C + nr, ldc,
		  kmax,winostep-1);

	delete[] X1;
	delete[] X2;
	delete[] X3;
}

// Control of the 2 cut-off criterias determining when to switch from Winograd's 
// to classic multiplication, and from finite field to double.
// Fix-up for odd-sized matrices using dynamic pealing ( coming soon...)
// for matrices over a finite Field
template <class Field>
inline  void LinBox::FFLAS::WinoMain( const Field& F, 
				      const size_t m, const size_t n, const size_t k,
				      const typename Field::Element alpha,
				      const typename Field::Element* A,const size_t lda,
				      const typename Field::Element* B,const size_t ldb,
				      const typename Field::Element beta,
				      typename Field::Element * C, const size_t ldc,
				      long long kmax, size_t winostep){
	static typename Field::Element one,zero,mone;
	F.init(one, 1);
	F.init(zero, 0);
	F.init(mone, -1);

	if (winostep<=0) // Winograd -> Classic
		ClassicMatmul( F, FflasNoTrans, FflasNoTrans, m, n, k,
			       alpha, A, lda, B, ldb, beta, C, ldc, kmax );
	else{
		size_t nr = n/2;
		size_t kr = k/2;
		size_t mr = m/2;
		if (kr < kmax){ // switch on double
			// Temporary double matrices
			DoubleDomain::Element alphad, betad;
			typename Field::Element _betabis;
			
 			if (F.areEqual( mone, alpha )){
				alphad = -1.0;
				F.convert( betad, beta );
			}
 			else{
				if (! F.areEqual( one, alpha) ){
					// Compute C = A*B + beta/alpha.C
					// and after C *= alpha
					F.div( _betabis, beta, alpha );
					F.convert( betad, _betabis );
				}
				else
					F.convert( betad, beta );
				alphad = 1.0;
			}
				
			DoubleDomain::Element * Ad = new DoubleDomain::Element[m*k];
			DoubleDomain::Element * Bd = new DoubleDomain::Element[k*n];
			DoubleDomain::Element * Cd = new DoubleDomain::Element[m*n];
			// Conversion GFq => double
			MatF2MatD( F, Ad, A, lda, m, k);
			MatF2MatD( F, Bd, B, ldb, k, n); 
			if ( !F.isZero(beta) )
				MatF2MatD( F, Cd, C, ldc, m, n ); 
			// recursive call
			WinoMain( DoubleDomain(), m, n, k, 
				  alphad, Ad, k, Bd, n, betad, Cd, n, kmax, winostep);
			// Conversion double => GFq
			MatD2MatF( F, C, ldc, Cd, m, n );
			
			if ( !F.areEqual( one, alpha ) && !F.areEqual( mone, alpha) ){
				// Fix-up: compute C *= alpha
				for (typename Field::Element* Ci=C; Ci<C+m*ldc; Ci+=ldc)
					for ( size_t j=0; j<n; ++j ) 
						F.mulin( *( Ci + j ), alpha );
			}
			// Temporary double matrices destruction
			delete[] Ad;
			delete[] Bd;
			delete[] Cd;
		}
		else{
			WinoCalc( F, mr, nr, kr, alpha, A, lda, B, ldb, beta, C, ldc,
				  kmax,winostep);

			// Updates for oddsized dimensions
			DynamicPealing( F, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc, kmax );
		}
	}
}


// Control of the cut-off criterion determing when to switch 
// to classic multiplication
// Fix-up for odd-sized matrices using dynamic pealing
// for matrices over double
template<>
inline  void LinBox::FFLAS::WinoMain( const DoubleDomain& D, 
				      const size_t m, const size_t n, const size_t k,
				      const DoubleDomain::Element alpha,
				      const DoubleDomain::Element * A, const size_t lda,
				      const DoubleDomain::Element * B, const size_t ldb,
				      const DoubleDomain::Element beta,
				      DoubleDomain::Element * C, const size_t ldc,
				      long long kmax, size_t winostep){
	
	if (winostep<=0) // switch Winograd => Classic
		ClassicMatmul( D, FflasNoTrans, FflasNoTrans, m, n, k,
			       alpha, A, lda, B, ldb, beta, C, ldc, kmax );
	else{
		WinoCalc( D, m/2, n/2, k/2, alpha, A, lda, B, ldb, beta, C, ldc,
			  kmax,winostep); 

		DynamicPealing( D, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc, kmax );
	}
}

// 		switch(mkn) { 
// 		case 1: // n oddsized
// 			cblas_dgemv( CblasRowMajor, CblasNoTrans, m, k,
// 				     (DoubleDomain::Element) alpha, A, lda, 
// 				     B+n-1, ldb, 
// 				     (DoubleDomain::Element) beta, C+n-1,ldc);
// 			break;
      
// 		case 2: // k oddsized
// 			cblas_dger( CblasRowMajor, m, n,
// 				    (DoubleDomain::Element) alpha,
// 				    A+k-1, lda, B+(k-1)*ldb, 1, C, ldc);
// 			break;
			
// 		case 3: // n, k oddsized
// 			cblas_dger( CblasRowMajor, m, n-1,
// 				    (DoubleDomain::Element) alpha,
// 				    A+k-1, lda, B+(k-1)*ldb, 1, C, ldc);
// 			cblas_dgemv( CblasRowMajor, CblasNoTrans, m, k-1,
// 				     (DoubleDomain::Element) alpha, A, lda,
// 				     B+n-1, ldb, 
// 				     (DoubleDomain::Element) beta, C+n-1,ldc);
// 			cblas_daxpy( m, alpha*(*(B+(k-1)*ldb+n-1)), 
// 				     A+k-1, lda, C+n-1, ldc);
// 			break;
			
// 		case 4: // m oddsized
// 			cblas_dgemv( CblasRowMajor, CblasTrans, k, n,
// 				     (DoubleDomain::Element) alpha,
// 				     B, ldb, A+(m-1)*lda, 1,
// 				     (DoubleDomain::Element) beta, C+(m-1)*ldc, 1);
// 			break;
			
// 		case 5: // m, n oddsized
// 			cblas_dgemv( CblasRowMajor, CblasNoTrans, m-1, k,
// 				     (DoubleDomain::Element) alpha,
// 				     A, lda, B+n-1, ldb,
// 				     (DoubleDomain::Element) beta, C+n-1, ldc);
// 			cblas_dgemv( CblasRowMajor, CblasTrans, k, n-1,
// 				     (DoubleDomain::Element) alpha, 
// 				     B, ldb, A+(m-1)*lda, 1,
// 				     (DoubleDomain::Element) beta, C+(m-1)*ldc, 1);
// 			*(C+(m-1)*ldc+n-1) = alpha*cblas_ddot(k, A+(m-1)*lda, 1, B+n-1, ldb)
// 				+ beta*(*(C+(m-1)*ldc+n-1));
// 			break;
      
// 		case 6: // m, k oddsized
// 			cblas_dger( CblasRowMajor, m-1, n, 
// 				    (DoubleDomain::Element) alpha, A+k-1, lda,
// 				    B+(k-1)*ldb, 1, C, ldc);
// 			cblas_dgemv( CblasRowMajor, CblasTrans, k-1, n,
// 				     (DoubleDomain::Element) alpha,
// 				     B, ldb, A+(m-1)*lda, 1,
// 				     (DoubleDomain::Element) beta, 
// 				     C+(m-1)*ldc, 1);
// 			cblas_daxpy( n, alpha*(*(A+(m-1)*lda+k-1)),
// 				     B+(k-1)*ldb, 1, C+(m-1)*ldc, 1);
// 			break;
      
// 		case 7: // m, k, n oddsized
// 			// Block NW
// 			cblas_dger( CblasRowMajor, m-1, n-1,
// 				    (DoubleDomain::Element) alpha, 
// 				    A+k-1, lda, B+(k-1)*ldb, 1, C, ldc);
// 			// Block SW
// 			cblas_dgemv( CblasRowMajor, CblasTrans, k-1, n-1,
// 				     (DoubleDomain::Element) alpha,
// 				     B, ldb, A+(m-1)*lda, 1,
// 				     (DoubleDomain::Element) beta,
// 				     C+(m-1)*ldc, 1);
			
// 			cblas_daxpy( n-1, alpha*(*(A+(m-1)*lda+k-1)), 
// 				     B+(k-1)*ldb, 1, C+(m-1)*ldc, 1);
// 			// Block NE
// 			cblas_dgemv( CblasRowMajor, CblasNoTrans, m-1, k-1,
// 				     (DoubleDomain::Element) alpha,
// 				     A, lda, B+n-1, ldb,
// 				     (DoubleDomain::Element) beta,
// 				     C+n-1, ldc);
// 			cblas_daxpy( m-1, alpha*(*(B+(k-1)*ldb+n-1)), 
// 				     A+k-1, lda, C+n-1, ldc);
				
// 			// Block SE
// 			*(C+(m-1)*ldc+n-1) = alpha*cblas_ddot(k, A+(m-1)*lda, 1, B+n-1, ldb)
// 				+ beta*(*(C+(m-1)*ldc+n-1));
// 			break;
// 		}
		
// 	}
// }


template <class Field>
inline void
LinBox::FFLAS::DynamicPealing( const Field& F, 
			       const size_t m, const size_t n, const size_t k,
			       const typename Field::Element alpha, 
			       const typename Field::Element* A, const size_t lda,
			       const typename Field::Element* B, const size_t ldb, 
			       const typename Field::Element beta,
			       typename Field::Element* C, const size_t ldc, 
			       const size_t kmax ){

	typename Field::Element tmp;
	size_t mkn = (n & 0x1)+( (k & 0x1)<<1 )+ ( (m & 0x1)<<2 ); 
	switch(mkn) { 
	case 1: // n oddsized
		fgemv( F, FflasNoTrans, m, k,
		       alpha, A, lda, B+n-1, ldb, beta, C+n-1,ldc);
		break;
      
	case 2: // k oddsized
		fger( F, m, n, alpha, A+k-1, lda, B+(k-1)*ldb, 1, C, ldc);
		break;
			
	case 3: // n, k oddsized
		fger( F, m, n-1, alpha, A+k-1, lda, B+(k-1)*ldb, 1, C, ldc);
		fgemv( F, FflasNoTrans, m, k-1,
		       alpha, A, lda, B+n-1, ldb, beta, C+n-1,ldc);
		F.mul( tmp, alpha, (*(B+(k-1)*ldb+n-1)) );
		faxpy( F, m, tmp, A+k-1, lda, C+n-1, ldc);
		break;
			
	case 4: // m oddsized
		fgemv( F, FflasTrans, k, n,
		       alpha, B, ldb, A+(m-1)*lda, 1, beta, C+(m-1)*ldc, 1);
		break;
			
	case 5: // m, n oddsized
		fgemv( F, FflasNoTrans, m-1, k,
		       alpha, A, lda, B+n-1, ldb, beta, C+n-1, ldc);
		fgemv( F, FflasTrans, k, n-1,
		       alpha, B, ldb, A+(m-1)*lda, 1, beta, C+(m-1)*ldc, 1);
		F.mulin( *(C+(m-1)*ldc+n-1), beta );
		F.axpyin( *(C+(m-1)*ldc+n-1), alpha, fdot( F, k, A+(m-1)*lda, 1, B+n-1, ldb) );
		break;
      
	case 6: // m, k oddsized
		fger( F, m-1, n, alpha, A+k-1, lda, B+(k-1)*ldb, 1, C, ldc);
		fgemv( F, FflasTrans, k-1, n,
		       alpha, B, ldb, A+(m-1)*lda, 1, beta, C+(m-1)*ldc, 1);
		F.mul( tmp, alpha, (*(A+(m-1)*lda+k-1)) ); 
		faxpy( F, n, tmp, B+(k-1)*ldb, 1, C+(m-1)*ldc, 1);
		break;
      
	case 7: // m, k, n oddsized
		// Block NW
		fger( F, m-1, n-1, alpha, A+k-1, lda, B+(k-1)*ldb, 1, C, ldc);
		// Block SW
		fgemv( F, FflasTrans, k-1, n-1,
		       alpha, B, ldb, A+(m-1)*lda, 1, beta, C+(m-1)*ldc, 1 );
		F.mul( tmp, alpha, *(A+(m-1)*lda+k-1) );
		faxpy( F, n-1, tmp, B+(k-1)*ldb, 1, C+(m-1)*ldc, 1);
		// Block NE
		fgemv( F, FflasNoTrans, m-1, k-1,
		       alpha, A, lda, B+n-1, ldb, beta, C+n-1, ldc);
		F.mul( tmp, alpha, *(B+(k-1)*ldb+n-1) );
		faxpy( F, m-1, tmp, A+k-1, lda, C+n-1, ldc);
				
		// Block SE
		F.mulin( *(C+(m-1)*ldc+n-1), beta);
		F.axpyin( *(C+(m-1)*ldc+n-1), alpha, fdot( F, k, A+(m-1)*lda, 1, B+n-1, ldb) );
		break;
	}
}

//-------------------------------------------------------------------
// visible function: operates C = A*B over double or a finite Field
//-------------------------------------------------------------------
template<class Field>
inline  typename Field::Element* 
LinBox::FFLAS::fgemm( const Field& F,
		      const enum FFLAS_TRANSPOSE ta,
		      const enum FFLAS_TRANSPOSE tb,
		      const size_t m, const size_t n, const size_t k,
		      const typename Field::Element alpha,
		      const typename Field::Element* A, const size_t lda,
		      const typename Field::Element* B, const size_t ldb, 
		      const typename Field::Element beta,
		      typename Field::Element* C, const size_t ldc,
		      const size_t winostep){
	if (!(m && n && k)) return C;
	
	integer charac;
	F.characteristic(charac);		
	long long c = charac-1;
	// Threshold between GFq and double
	long long kmax = ((long long)1<<53)/(c*c);
	//	long long kmax =6;
	if ( !winostep || ta==FflasTrans || tb==FflasTrans )
		ClassicMatmul( F, ta, tb,  m, n, k, alpha, A, lda, B, ldb,
			       beta, C,ldc, kmax );
	else
		WinoMain(F, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc, kmax, winostep);
	return C;
}

template<>
inline double*
LinBox::FFLAS::fsquare( const Modular<double>& F,
			const enum FFLAS_TRANSPOSE ta,
			const size_t n, const double alpha,
			const double* A, const size_t lda,
			const double beta,
			double* C, const size_t ldc){
	
	
	if ( C==A ){
		double * Ad = new double[n*n];
		for ( size_t i=0; i<n; ++i)
			fcopy( F, n,Ad+i*n, 1, A+i*lda, 1);
		cblas_dgemm( CblasRowMajor, (enum CBLAS_TRANSPOSE)ta,
			     (enum CBLAS_TRANSPOSE)ta, n, n, n, 
			     alpha, Ad, n, Ad, n, beta, C, ldc);
	}
	else
		cblas_dgemm( CblasRowMajor, (enum CBLAS_TRANSPOSE)ta,
			     (enum CBLAS_TRANSPOSE)ta, n, n, n, 
			     alpha, A, lda, A, lda, beta, C, ldc);
	// Conversion double => Finite Field
	size_t i;
	double *Ci;
	for ( i=0, Ci=C ; i<n;++i, Ci+=ldc)
		for ( size_t j=0; j<n;++j)
			F.init(*(Ci+j),*(Ci+j));
	return C;
}

template<class Field>
inline  typename Field::Element*
LinBox::FFLAS::fsquare( const Field& F,
			const enum FFLAS_TRANSPOSE ta,
			const size_t n, const typename Field::Element alpha,
			const typename Field::Element* A, const size_t lda,
			const typename Field::Element beta,
			typename Field::Element* C, const size_t ldc){
	double alphad, betad;
	F.convert( alphad, alpha );
	if ( F.areEqual( beta, mone ) )
		betad = -1.0;
	else
		F.convert( betad, beta );

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
		     (DoubleDomain::Element) alphad, Ad, n, Ad, n,
		     (DoubleDomain::Element) betad, Cd, n);
	// Conversnion double => Finite Field
	delete[] Ad;
	MatD2MatF( F, C, ldc, Cd, n, n);
	delete[] Cd;
	return C;
}
