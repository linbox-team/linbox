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
#ifndef MIN
#define MIN(a,b) (a>b)?b:a
#endif

#include <linbox/field/modular-double.h>
// Note:
// The domain is supposed to be a field since some divisions are required for efficiency purposes
// An alternative has to written for finite rings if necessary

// Classic Multiplication over double
template <>
inline  void LinBox::FFLAS::ClassicMatmul(const DoubleDomain& F, 
					  const enum FFLAS_TRANSPOSE ta,
					  const enum FFLAS_TRANSPOSE tb,
					  const size_t m, const size_t n,const size_t k,
					  const DoubleDomain::Element alpha,
					  const DoubleDomain::Element * Ad, const size_t lda,
					  const DoubleDomain::Element * Bd, const size_t ldb,
					  const DoubleDomain::Element beta,
					  DoubleDomain::Element * Cd, const size_t ldc, 
					  const size_t kmax){
	
	// Call to the blas multiplication
	cblas_dgemm(CblasRowMajor, (enum CBLAS_TRANSPOSE) ta, 
		    (enum CBLAS_TRANSPOSE) tb,
		    m, n, k, (DoubleDomain::Element) alpha,
		    Ad, lda, Bd, ldb, (DoubleDomain::Element) beta,Cd, ldc);
}


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
	typename Field::Element tmp;
	DoubleDomain::Element alphad, betad;
	size_t k2 = MIN(k,kmax-1);
	DoubleDomain::Element * Add = new DoubleDomain::Element[m*k2];
	DoubleDomain::Element * Bdd = new DoubleDomain::Element[k2*n];
	DoubleDomain::Element * Cd = new DoubleDomain::Element[m*n];
	
	// To ensure the initial computation with beta
	size_t nblock = k / (kmax-1);
	size_t remblock = k % (kmax-1);
	if (!remblock){
		remblock = kmax - 1 ;
		--nblock;
	}
	if ( F.areEqual( Mone, beta ) )
		betad = -1.0;
	else
		F.convert( betad, beta );
	
	if ( F.areEqual( Mone, alpha ) )
		alphad = -1.0;
	else{
		alphad = 1.0;
		if (! F.areEqual( one, alpha) ){
			// Compute y = A*x + beta/alpha.y
			// and after y *= alpha
			F.div( tmp, beta, alpha );
			F.convert( betad, tmp );
		}
	}
	
	size_t am, an, bm, bn, dlda, dldb;
	if (!F.isZero(beta)){
		MatF2MatD( F, Cd, n, C, ldc, m, n );
	}
	if (ta == FflasTrans){
		dlda = m; am = k2; an = m;
		MatF2MatD( F, Add, m, A+k2*nblock*lda, lda, remblock, m );
	}
	else{
		dlda = k2; am = m; an = k2;
		MatF2MatD( F, Add, k2, A+k2*nblock, lda, m, remblock);
	}
	if (tb == FflasTrans){
		dldb = k; bm = n; bn = k2; 
		MatF2MatD( F, Bdd, k2, B+k2*nblock, ldb, n, remblock );
	}
	else{
		dldb = n; bm = k2; bn = n;
		MatF2MatD( F, Bdd, n, B+k2*nblock*ldb, ldb, remblock, n );
	}
	
// 	std::cout<<"m,n,remblock,alphad,betad,k2,nblock"<<m<<" "<<n<<" "
//	 <<remblock<<" "<<alphad<<" "<<betad<<" "<<k2<<" "<<nblock<<std::endl;
// 	std::cout<<"A="<<std::endl;
// 	write_field(F, std::cout, A, m,k,lda);
// 	std::cout<<"Add="<<std::endl;
// 	write_field(DoubleDomain(), std::cout, Add, m,remblock,k2);
// 	std::cout<<"B="<<std::endl;
// 	write_field(F, std::cout, B, k,n,ldb);
// 	std::cout<<"Bdd="<<std::endl;
// 	write_field(DoubleDomain(), std::cout, Bdd, remblock,n,n);
// 	std::cout<<"C="<<std::endl;
// 	write_field(F, std::cout,C, m,n,ldc);
// 	std::cout<<"Cd="<<std::endl;
// 	write_field(DoubleDomain(), std::cout,Cd, m,n,n);

	ClassicMatmul( DoubleDomain(), ta, tb, m, n, remblock, alphad, Add, dlda,
		       Bdd, dldb, betad, Cd, n, kmax );
// 	std::cout<<"Avant conversion:"<<std::endl;
// 	std::cout<<"C="<<std::endl;
// 	write_field(F, std::cout,C, m,n,ldc);
// 	std::cout<<"Cd="<<std::endl;
// 	write_field(DoubleDomain(), std::cout,Cd, m,n,n);
	MatD2MatF( F, C, ldc, Cd, n, m, n );
	MatF2MatD( F, Cd, n, C, ldc, m, n );
// 	std::cout<<"Apres conversion:"<<std::endl;
// 	std::cout<<"C="<<std::endl;
// 	write_field(F, std::cout,C, m,n,ldc);
// 	std::cout<<"Cd="<<std::endl;
// 	write_field(DoubleDomain(), std::cout,Cd, m,n,n);

	for ( size_t i = 0; i < nblock; ++i ){

		if (ta == FflasTrans){
			dlda = m; am = k2; an = m;
			MatF2MatD( F, Add, m, A+k2*i*lda, lda, k2, m );
		}
		else{
			dlda = k2; am = m; an = k2;
			MatF2MatD( F, Add, k2,  A+k2*i, lda, m, k2);
		}
		if (tb == FflasTrans){
			dldb = k; bm = n; bn = k2; 
			MatF2MatD( F, Bdd, k2, B+k2*i, ldb, n, k2 );
		}
		else{
			dldb = n; bm = k2; bn = n;
			MatF2MatD( F, Bdd, n, B+k2*i*ldb, ldb, k2, n );
		}
		
// 		std::cout<<"Add="<<std::endl;
// 		write_field(DoubleDomain(), std::cout, Add, m,k2, k2);
// 		std::cout<<"Bdd="<<std::endl;
// 		write_field(DoubleDomain(), std::cout, Bdd, k2,n,n);
// 		std::cout<<"C="<<std::endl;
// 		write_field(F, std::cout,C, m,n,ldc);
// 		std::cout<<"Cd="<<std::endl;
// 		write_field(DoubleDomain(), std::cout,Cd, m,n,n);
		ClassicMatmul( DoubleDomain(), ta, tb, m, n, k2, alphad, Add, dlda,
			       Bdd, dldb, 1.0, Cd, n, kmax );

// 	std::cout<<"Avant conversion:"<<std::endl;
// 	std::cout<<"C="<<std::endl;
// 	write_field(F, std::cout,C, m,n,ldc);
// 	std::cout<<"Cd="<<std::endl;
// 	write_field(DoubleDomain(), std::cout,Cd, m,n,n);
	
	MatD2MatF( F, C, ldc, Cd, n, m, n );
	MatF2MatD( F, Cd, n, C, ldc, m, n );
	
// 	std::cout<<"Apres conversion:"<<std::endl;
// 	std::cout<<"C="<<std::endl;
// 	write_field(F, std::cout,C, m,n,ldc);
// 	std::cout<<"Cd="<<std::endl;
// 	write_field(DoubleDomain(), std::cout,Cd, m,n,n);
	}
	
	if ( (!F.areEqual( one, alpha )) && (!F.areEqual( Mone, alpha)) ){
		for ( typename Field::Element * Ci = C; Ci < C+m*ldc; Ci += ldc)
			for ( size_t j = 0; j < n; ++j ) 
				F.mulin( *( Ci + j ), alpha );
	}

	// Conversion double => Finite Field
	delete[] Add;
	delete[] Bdd;
// if ( !F.areEqual(alpha, Mone) && !F.isOne(alpha)){
// 		// Cd may contain approximations of C entries:
// 		// Perturbation of Cd to get C from truncation.
		
// 		for ( double* Cdi = Cd; ; Cdi<Cd+m*ldc; Cdi+=ldc)
// 			for(size_t j=0; j<n; ++j)
// 				if (*(Cdi+j)>=0){
// 					*(Cdi+j) += 0.1;
// 				}
// 				else
// 					*(Cdi+j) -= 0.1;
// 	}
	delete[] Cd;
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
	//std::cerr<<"Classic<modulardouble> k, kmax, m, n, k, alpha, beta="<<k<<" "<<kmax
	//	 <<" "<<m<<" "<<n<<" "<<k<<" "<<alpha<<" "<<beta<<std::endl;
	double _alpha, _beta;
	// To ensure the initial computation with beta
	size_t k2 = MIN(k,kmax-1);
	size_t nblock = k / (kmax-1);
	size_t remblock = k % (kmax-1);
	if (!remblock){
		remblock = kmax-1;
		--nblock;
	}
	if ( F.areEqual( Mone, beta ) )
		_beta = -1.0;
	else
		_beta = beta;

	if ( F.areEqual( Mone, alpha ) )
		_alpha = -1.0;
	else{
		_alpha = 1.0;
		if (! F.areEqual( one, alpha) ){
			// Compute y = A*x + beta/alpha.y
			// and after y *= alpha
			F.divin( _beta, alpha );
		}
	}
	
	ClassicMatmul( DoubleDomain(), ta, tb, m, n, remblock, _alpha, A+k2*nblock, lda,
		       B+k2*nblock*ldb, ldb, _beta, C, ldc, kmax );
	for ( double * Ci = C; Ci != C+m*ldc; Ci += ldc)
		for ( size_t j=0; j<n;++j)
			F.init(*(Ci+j),*(Ci+j));
	for ( size_t i = 0; i < nblock; ++i ){
		ClassicMatmul( DoubleDomain(), ta, tb, m, n, k2, _alpha, A+i*k2, lda,
			       B+i*k2*ldb, ldb, one, C, ldc, kmax );
		for ( double * Ci = C; Ci != C+m*ldc; Ci += ldc)
			for ( size_t j=0; j<n;++j)
				F.init(*(Ci+j),*(Ci+j));
	}
	
	if ( (!F.areEqual( one, alpha )) && (!F.areEqual( Mone, alpha)) ){
		for ( double * Ci = C; Ci < C+m*ldc; Ci += ldc)
			for ( size_t j = 0; j < n; ++j ) 
				F.mulin( *( Ci + j ), alpha );
	}
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

	if ( !F.isZero(beta) ){
		// C22 = C22 - C12 if beta != 0
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
	
	// P3 = alpha . S4*B22 in X1
	WinoMain( F, mr, nr, kr, alpha, X2, kr, B + kr*ldb + nr, ldb, zero, X1, nr,
		  kmax,winostep-1);
	// U5 = P3 + U4 in C12
	d12c = C+nr;dx1 = X1; 
	for ( size_t i = 0; i < mr; ++i, d12c += ldc, dx1 += nr )
		for ( size_t j = 0; j < nr; ++j )
			F.addin( *(d12c + j), *(dx1 + j) );

	delete[] X1;
	delete[] X2;
	delete[] X3;
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

		//std::cerr<<"avant dyn pealing sur double: alpha= "<<alpha<<std::endl;
		DynamicPealing( D, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc, kmax );
	}
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
			MatF2MatD( F, Ad, k, A, lda, m, k);
			MatF2MatD( F, Bd, n, B, ldb, k, n); 
			if ( !F.isZero(beta) )
				MatF2MatD( F, Cd, n, C, ldc, m, n ); 
			// recursive call
			WinoMain( DoubleDomain(), m, n, k, 
				  alphad, Ad, k, Bd, n, betad, Cd, n, kmax, winostep);
			// Conversion double => GFq
			MatD2MatF( F, C, ldc, Cd, n, m, n );
			
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

template <>
inline  void LinBox::FFLAS::WinoMain( const Modular<double>& F, 
				      const size_t m, const size_t n, const size_t k,
				      const double alpha,
				      const double* A,const size_t lda,
				      const double* B,const size_t ldb,
				      const double beta,
				      double * C, const size_t ldc,
				      long long kmax, size_t winostep){
	if (winostep<=0) // Winograd -> Classic
		ClassicMatmul( F, FflasNoTrans, FflasNoTrans, m, n, k,
			       alpha, A, lda, B, ldb, beta, C, ldc, kmax );
	else{
		size_t nr = n/2;
		size_t kr = k/2;
		size_t mr = m/2;
		if (kr < kmax){ // switch on double
			// Temporary double matrices
			DoubleDomain::Element _alpha, _beta;
			
			_beta = beta;
 			if (F.areEqual( -1.0, alpha )){
				_alpha = -1.0;
			}
 			else{
				if (! F.areEqual( 1.0, alpha) ){
					// Compute C = A*B + beta/alpha.C
					// and after C *= alpha
					F.divin( _beta, alpha );
				}
				_alpha = 1.0;
			}
				
			// recursive call
			WinoMain( DoubleDomain(), m, n, k, 
				  _alpha, A, lda, B, ldb, _beta, C, ldc, kmax, winostep);
			// Conversion double => GFq
			for ( double * Ci = C; Ci != C+m*ldc; Ci+=ldc )
				for (size_t j = 0; j < n; ++j )
					F.init( *(Ci + j), *(Ci + j) );
			
			if ( !F.areEqual( 1.0, alpha ) && !F.areEqual( -1.0, alpha) ){
				// Fix-up: compute C *= alpha
				for ( double* Ci=C; Ci<C+m*ldc; Ci+=ldc)
					for ( size_t j=0; j<n; ++j ) 
						F.mulin( *( Ci + j ), alpha );
			}
			// Temporary double matrices destruction
		}
		else{
			WinoCalc( F, mr, nr, kr, alpha, A, lda, B, ldb, beta, C, ldc,
				  kmax,winostep);

			// Updates for oddsized dimensions
			DynamicPealing( F, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc, kmax );
		}
	}
}


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
	size_t ex=1;
	for (size_t i=0;i<winostep; ++i)
		ex *= 3;
	long long c = (charac-1)*(1+ex)/2;
	// Threshold between GFq and double
	//	long long kmax = ((long long)1<<53)/(c*c);
	integer _beta;
	F.convert(_beta, beta);
	size_t kmax = ( (( (long long) 1<<53) )/(c*c) + 1)*(1<<winostep);
	//	size_t kmax = 2;
	if ( kmax == (size_t) (1<<winostep) )
		kmax=2;
	std::cout<<"kmax = "<<kmax<<"...";
	if ( !winostep || ta==FflasTrans || tb==FflasTrans ){
		ClassicMatmul( F, ta, tb,  m, n, k, alpha, A, lda, B, ldb,
			       beta, C,ldc, kmax );
	}
	else
		WinoMain( F, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc, kmax, winostep);
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
		fgemm( F, ta, ta, n, n, n, alpha, Ad, n, Ad, n,
		       beta, C, ldc );
		
// 		cblas_dgemm( CblasRowMajor, (enum CBLAS_TRANSPOSE)ta,
// 			     (enum CBLAS_TRANSPOSE)ta, n, n, n, 
//			     alpha, Ad, n, Ad, n, beta, C, ldc);
	}
	else
		fgemm( F, ta, ta, n, n, n, alpha, A, lda, A, lda,
		       beta, C, ldc );		
// 	cblas_dgemm( CblasRowMajor, (enum CBLAS_TRANSPOSE)ta,
// 			     (enum CBLAS_TRANSPOSE)ta, n, n, n, 
// 			     alpha, A, lda, A, lda, beta, C, ldc);
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
	static typename Field::Element mone;
	F.init( mone, -1 );
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
	MatF2MatD( F, Ad, n, A, lda, n, n);
	if (!F.isZero(beta))
		MatF2MatD( F, Cd, n, C, ldc, n, n); 
	
	// Call to the blas Multiplication 
	cblas_dgemm( CblasRowMajor, (enum CBLAS_TRANSPOSE)ta,
		     (enum CBLAS_TRANSPOSE)ta, n, n, n, 
		     (DoubleDomain::Element) alphad, Ad, n, Ad, n,
		     (DoubleDomain::Element) betad, Cd, n);
	// Conversnion double => Finite Field
	delete[] Ad;
	MatD2MatF( F, C, ldc, Cd, n, n, n);
	delete[] Cd;
	return C;
}
