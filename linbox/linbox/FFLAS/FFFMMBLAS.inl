/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
//----------------------------------------------------------------------
//                   Finite Field Fast Matrix Multiplication 
//                using Basic Linear Algebra Subroutines
//----------------------------------------------------------------------
// by Clement PERNET ( clement.pernet@imag.fr )
// 24/03/2003
//----------------------------------------------------------------------
// Note that the fix-up for odd-sized matrices is only implemented for 
// matrices over double.
//----------------------------------------------------------------------

extern "C" {
#include "cblas.h"
}
#include "linbox/integer.h"

// Classic Multiplication over double

inline void FFFMMBLAS::ClassicMatmul(const DoubleDomain& D, 
				     double * Ad,
				     double * Bd,
				     double * Cd,
				     int ra, int rb, int rc,
				     int m, int k, int n){
	// Call to the blas multiplication
	cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
		    m, n, k, (double) ALPHA,
		    Ad, ra,
		    Bd, rb,
		    (double) BETA,Cd, rc);
}

// Classic multiplication over a field
// (Case where all Winograd's operations have been performed over the field)
template <class Field>
inline void FFFMMBLAS::ClassicMatmul(const Field& F, 
				     typename Field::Element * A,
				     typename Field::Element * B,
				     typename Field::Element * C,
				     int ra, int rb, int rc,
				     int m, int k, int n){
	// Double  matrices initialisation
	double * Ad = new double[m*k];
	double * Bd = new double[k*n];
	double * Cd = new double[m*n];
	// Conversion finite Field => double

	MatGFq2MatDouble(F,m,k,Ad,k,A,ra);
	MatGFq2MatDouble(F,k,n,Bd,n,B,rb); 


	if (BETA)
		MatGFq2MatDouble(F,m,n,Cd,n,C,rc); 
		
    
	// Call to the blas Multiplication 
	cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
		    m, n, k, (double) ALPHA,
		    Ad, k, Bd, n,
		    (double) BETA,Cd, n);

	// Conversion double => Finite Field
	MatDouble2MatGFq(F,m,n,C,rc,Cd,n);
	delete[] Ad;
	delete[] Bd;
	delete[] Cd;
}


// Winograd Multiplication  A(n*k) * B(k*m) in C(n*m)

// Computation of the 22 Winograd's operations
template<class Field>
inline void FFFMMBLAS::WinoCalc( const Field& F, 
				 typename Field::Element * A,
				 typename Field::Element * B,
				 typename Field::Element * C,
				 typename Field::Element ** t_X1,
				 typename Field::Element ** t_X2,
				 typename Field::Element ** t_X3,
				 int ra, int rb, int rc,
				 int nr, int kr, int mr,
				 long long sf, int s){
  
	typename Field::Element tmp;
	typename Field::Element* d11,*d12,*d21,*d22,*d11c,*d12c,*dx1,*dx2,*dx3;
  
	// Three temporary submatrices are required
	typename Field::Element* X1 = t_X1[s];
	typename Field::Element* X2 = t_X2[s];
	typename Field::Element* X3 = t_X3[s];

	// P2 = alpha*A12 * B21 + betaC11 in C11
	WinoMain( F, 
		  A+kr,
		  (B+kr*rb), 
		  C,
		  t_X1,t_X2,t_X3,
		  ra,rb,rc,
		  nr,kr,mr,
		  sf,s-1);
  
	// P1 = A11 * B11 in C12
	WinoMain(F,A,B,C+mr,
		 t_X1,t_X2,t_X3,
		 ra,rb,rc,nr,kr,mr,sf,s-1);
  
	// T3 = B22 - B12 in X3
	d12 = B + mr; d22 = d12 + kr*rb; dx3 = X3;
	for (int i=0;i<kr;++i, d12+=rb, d22+=rb,dx3+=mr){
		for (int j=0;j<mr;++j)
			F.sub(*(dx3+j),*(d22 + j),*(d12 + j));
    
	}
	// S3 = A11 - A21 in X2 and  U1 = P2 + P1 in C11
	d11 = A; d21 = d11 + nr*ra; dx2 = X2; d11c = C;d12 = d11c + mr; 
	for (int i=0;i<nr;++i,d11+=ra,d11c+=rc,d12+=rc, d21+=ra,dx2+=kr){
		for (int j=0;j<kr;++j)
			F.sub(*(dx2+j),*(d11 + j),*(d21 + j));
		for (int j=0;j<mr;++j)
			F.addin(*(d11c + j),*(d12 + j));
	}
  
	// P7 = S3*T3 in C21
	WinoMain(F,X2,X3,C+nr*rc,
		 t_X1,t_X2,t_X3,
		 kr,mr,rc,nr,kr,mr,sf,s-1);
  
	// S1 = A21 + A22 in X2
	d21 = A + nr*ra; d22 = d21 + kr; dx2 = X2;
	for (int i=0;i<nr;++i, d21+=ra, d22+=ra, dx2+=kr){
		for (int j=0;j<kr;++j)
			F.add(*(dx2+j),*( d21 + j),*(d22 + j));
	}
  
	// T1 = B12 - B11 in X3
	d11 = B; d12 = d11 + mr; dx3 = X3;
	for (int i=0;i<kr;++i, d11+=rb, d12+=rb,dx3+=mr){
		for (int j=0;j<mr;++j)
			F.sub(*(dx3+j),*(d12 + j),*(d11 + j));
	}
  
	// P5 = S1*T1 in C22
	WinoMain(F,X2,X3,C+nr*rc + mr,
		 t_X1,t_X2,t_X3,
		 kr,mr,rc,nr,kr,mr,sf,s-1);

	// T2 = B22 - T1 in X3
	d22 = B+kr*rb+mr; dx3 = X3;
	for (int i=0;i<kr;++i,d22+=rb,dx3+=mr){
		for (int j=0;j<mr;++j)
			F.sub(*(dx3+j),*(d22 + j),*(dx3+j));
	}
  
	// S2 = S1 - A11 in X2
	d11 = A; dx2 = X2;
	for (int i=0;i<nr;++i,d11+=ra,dx2+=kr){
		for (int j=0;j<kr;++j)
			F.subin(*(dx2+j),*(d11 + j));
	}
  
	// P6 = S2*T2 in X1
	WinoMain(F,X2,X3,X1,
		 t_X1,t_X2,t_X3,
		 kr,mr,mr,nr,kr,mr,sf,s-1);
  
	// T4 = T2 - B21 in X3
	d21 = B+kr*rb;dx3=X3;
	for (int i=0;i<kr;++i,d21+=rb,dx3+=mr){
		for (int j=0;j<mr;++j)
			F.subin(*(dx3+j),*( d21 + j));
	}
  
	// S4 = A12 -S2 in X2 and U2 = P1 + P6 in tmp and U3 = P7 + U2 in C21
	// and U4 = U2 + P5 in C12
	d12 = A+kr; dx2 = X2;d12c = C+mr;dx1=X1;d21 = C+nr*rc; d22 = d21+mr;
	for (int i=0;i<nr;++i, d12+=ra, dx2+=kr,d12c+=rc,dx1+=mr,d21+=rc,d22+=rc){
		for (int j=0;j<kr;++j)
			F.sub(*(dx2+j),*(d12 + j),*(dx2+j));
		for (int j=0;j<mr;++j){
			F.add(tmp,*(d12c + j),*(dx1+j));// temporary U2
			F.addin(*(d21 + j),tmp);       // U3
			F.add(*(d12c + j),tmp,*(d22 +j));  // U4
		}
	}
  
	// P3 = S4*B22 in X1
	WinoMain(F,X2,B + kr*rb + mr,X1,
		 t_X1,t_X2,t_X3,
		 kr,rb,mr,nr,kr,mr,sf,s-1);
  
	// U5 = U4 + P3 in C12
	d12c = C+mr;dx1=X1;
	for (int i=0;i<nr;++i, d12c+=rc, dx1+=mr){
		for (int j=0;j<mr;++j)
			F.addin(*(d12c + j),*(dx1+j));
	}
  
	// P4 = A22*T4 in X1
	WinoMain(F,A + nr*ra + kr,X3,X1,
		 t_X1,t_X2,t_X3,
		 ra,mr,mr,nr,kr,mr,sf,s-1);
  
	// U7 = P5 + U3 in C22
	d21 = C+nr*rc; d22 = d21 +mr;
	for (int i=0;i<nr;++i,d21+=rc,d22+=rc){
		for (int j=0;j<mr;++j)
			F.addin(*(d22 + j),*(d21 +j));
	}
  
	// U6 = U3 - P4 in C21
	d21 = C+nr*rc; dx1=X1;
	for (int i=0;i<nr;++i, d21+=rc, dx1+=mr){
		for (int j=0;j<mr;++j)
			F.subin(*(d21 + j),*(dx1+j));
	}
  
}

// Control of the 2 cut-off criterias determining when to switch from Winograd's 
// to classic multiplication, and from finite field to double.
// Fix-up for odd-sized matrices using dynamic pealing ( coming soon...)
// for matrices over a finite Field
template <class Field>
inline void FFFMMBLAS::WinoMain( const Field& F, 
				 typename Field::Element * const A,
				 typename Field::Element * const B,
				 typename Field::Element * C,
				 typename Field::Element ** t_X1,
				 typename Field::Element ** t_X2,
				 typename Field::Element ** t_X3,
				 int ra, int rb, int rc,
				 int m, int k, int n,
				 long long sf, int s){
	if (s<=0) // Winograd -> Classic
		ClassicMatmul(F,A,B,C,ra,rb,rc,m,k,n);
	else{
		int nr = n/2;
		int kr = k/2;
		int mr = m/2;
		if (kr < sf){ // switch on double
			// Temporary double matrices
			double* td_X1[s+1];
			double* td_X2[s+1];
			double* td_X3[s+1];
			for (int i=s;(i>=1);--i){
				td_X1[i] = new double[mr*nr];
				td_X2[i] = new double[mr*kr];
				td_X3[i] = new double[kr*nr];
				nr/=2;
				mr/=2;
				kr/=2;
			}
			double * Ad = new double[m*k];
			double * Bd = new double[k*n];
			double * Cd = new double[m*n];
			// Conversion GFq => double
			MatGFq2MatDouble(F,m,k,Ad,k,A,ra);
			MatGFq2MatDouble(F,k,n,Bd,n,B,rb); 
      
			if (BETA)
				MatGFq2MatDouble(F,m,n,Cd,n,C,rc); 
      
			// recursive call
			WinoMain( DoubleDomain(),
				  Ad, Bd, Cd,
				  td_X1,td_X2,td_X3,
				  k, n, n, m, k, n, sf, s);
			// Conversion double => GFq
			MatDouble2MatGFq(F,m,n,C,rc,Cd,n);
			// Temporary double matrices destruction
			delete[] Ad;
			delete[] Bd;
			delete[] Cd;
			for (int i=s;i>=1;--i){
				delete[] td_X1[i];
				delete[] td_X2[i];
				delete[] td_X3[i];
			}
      
		}
		else
			WinoCalc(F,A,B,C,t_X1,t_X2,t_X3,ra,rb,rc,mr,kr,nr,sf,s);
	}
}


// Control of the cut-off criteria determing when to switch 
// to classic multiplication
// Fix-up for odd-sized matrices using dynamic pealing
// for matrices over double
template<>
inline void FFFMMBLAS::WinoMain( const DoubleDomain& D, 
				 double * A,
				 double * B,
				 double * C,
				 double **t_X1,
				 double **t_X2,
				 double **t_X3,
				 int ra, int rb, int rc,
				 int n, int k, int m,
				 long long sf, int s){
	if (s<=0) // switch Winograd => Classic
		ClassicMatmul(D,A,B,C,ra,rb,rc,n,k,m);
	else{
		int nr = n/2;
		int kr = k/2;
		int mr = m/2;
		int mkn = (m & 0x1)+( (k & 0x1)<<1 )+ ( (n & 0x1)<<2 ); 
		WinoCalc(D,A,B,C,t_X1,t_X2,t_X3,ra,rb,rc,nr,kr,mr,sf,s); 
    
		switch(mkn) { 
		case 1: // m oddsized
			cblas_dgemv(CblasRowMajor,CblasNoTrans,n,k,(double) ALPHA,
				    A,ra,B+m-1,rb,
				    (double) BETA,C+m-1,rc);
			break;
      
		case 2: // k oddsized
			cblas_dger(CblasRowMajor,n,m,(double) ALPHA,
				   A+k-1,ra,B+(k-1)*rb,1,C,rc);
			break;
      
		case 3: // m, k oddsized
			cblas_dger(CblasRowMajor,n,m-1,(double) ALPHA,
				   A+k-1,ra,B+(k-1)*rb,1,C,rc);
			if (BETA){
				cblas_daxpy(n,ALPHA*BETA*(*(B+(k-1)*rb+m-1)),A+k-1,ra,C+m-1,rc);
				cblas_dgemv(CblasRowMajor,CblasNoTrans,n,k-1,(double) ALPHA,
					    A,ra,B+m-1,rb,
					    (double) BETA, C+m-1,rc);
			}
			else{
				cblas_dcopy(n,A + k-1, ra, C+m-1, rc);
				cblas_dscal(n,ALPHA*(*(B+(k-1)*rb+m-1)),C + m-1, rc);
				cblas_dgemv(CblasRowMajor,CblasNoTrans,n,k-1,(double) ALPHA,
					    A,ra,B+m-1,rb,
					    1.0, C+m-1,rc);
			}
			break;
      
		case 4: // n oddsized
			cblas_dgemv(CblasRowMajor,CblasTrans,k,m,(double) ALPHA,
				    B,rb,A+(n-1)*ra,1,
				    (double) BETA,C+(n-1)*rc,1);
			break;
      
		case 5: // n, m oddsized
			cblas_dgemv(CblasRowMajor,CblasNoTrans,n-1,k,(double) ALPHA,
				    A,ra,B+m-1,rb,
				    (double) BETA, C+m-1,rc);
			cblas_dgemv(CblasRowMajor,CblasTrans,k,m-1,(double) ALPHA,
				    B,rb,A+(n-1)*ra,1,
				    (double) BETA, C+(n-1)*rc,1);
			*(C+(n-1)*rc+m-1) = ALPHA*cblas_ddot(k,
							     A+(n-1)*ra,1,
							     B+m-1,rb) + BETA*(*(C+(n-1)*rc+m-1));
			break;
      
		case 6: // n, k oddsized
			cblas_dger(CblasRowMajor,n-1,m,(double) ALPHA,
				   A+k-1,ra,B+(k-1)*rb,1,C,rc);
			if (BETA){
				cblas_daxpy(m,ALPHA*BETA*(*(A+(n-1)*ra+k-1)),B+(k-1)*rb,1,C+(n-1)*rc,1);
				cblas_dgemv(CblasRowMajor,CblasTrans,k-1,m,(double) ALPHA,
					    B,rb,A+(n-1)*ra,1,
					    (double) BETA,C+(n-1)*rc,1);
			}
			else{
				cblas_dcopy(m,B + (k-1)*rb,1,C+(n-1)*rc,1);
				cblas_dscal(m,ALPHA*(*(A+(n-1)*ra+k-1)),C + (n-1)*rc,1);
				cblas_dgemv(CblasRowMajor,CblasTrans,k-1,m,(double) ALPHA,
					    B,rb,A+(n-1)*ra,1,
					    1.0,C+(n-1)*rc,1);
			}
			break;
      
		case 7: // n, k, m oddsized
			// Block NW
			cblas_dger(CblasRowMajor,n-1,m-1,(double) ALPHA,
				   A+k-1,ra,B+(k-1)*rb,1,C,rc);
			if (BETA){
				// Block SW
				cblas_daxpy(m-1,ALPHA*BETA*(*(A+(n-1)*ra+k-1)),B+(k-1)*rb,1,C+(n-1)*rc,1);
				cblas_dgemv(CblasRowMajor,CblasTrans,k-1,m-1,(double) ALPHA,
					    B,rb,A+(n-1)*ra,1,
					    (double) BETA,C+(n-1)*rc,1);
				// Block NE
				cblas_daxpy(n-1,ALPHA*BETA*(*(B+(k-1)*rb+m-1)),A+k-1,ra,C+m-1,rc);
				cblas_dgemv(CblasRowMajor,CblasTrans,n-1,k-1,(double) ALPHA,
					    A,ra,B+m-1,rb,
					    (double) BETA, C+m-1,rc);
			}
			else{
				// Block SW
				cblas_dcopy(m-1,B + (k-1)*rb,1,C+(n-1)*rc,1);
				cblas_dscal(m-1, ALPHA*(*(A+(n-1)*ra+k-1)), C + (n-1)*rc, 1);
				cblas_dgemv(CblasRowMajor,CblasTrans,k-1,m-1,(double) ALPHA,
					    B,rb,A+(n-1)*ra,1,
					    1.0,C+(n-1)*rc,1);
      
				// Block NE

				cblas_dcopy(n-1,A + k-1, ra, C+m-1, rc);
				cblas_dscal(n-1,ALPHA*(*(B+(k-1)*rb+m-1)),C + m-1, rc);
				cblas_dgemv(CblasRowMajor,CblasNoTrans,n-1,k-1,(double) ALPHA,
					    A,ra,B+m-1,rb,
					    1.0, C+m-1,rc);
			}

			// Block SE
			*(C+(n-1)*rc+m-1) = ALPHA*cblas_ddot(k,A+(n-1)*ra,1,B+m-1,rb) + BETA*(*(C+(n-1)*rc+m-1));
			break;
		}

	}
}	

//-------------------------------------------------------------------
// visible function: operates C = A*B over double or a finite Field
//-------------------------------------------------------------------

inline double* FFFMMBLAS::operator ()( const DoubleDomain& D, 
				       const int m,
				       const int n,
				       const int k,
				       const enum FFFMMBLAS_COEF alpha,
				       double* const A, const int lda,
				       double* const B, const int ldb, 
				       const enum FFFMMBLAS_COEF beta,
				       double* C, const int ldc,
				       const int nbe){
  

	ALPHA = alpha;
	BETA = beta;
	int mr=m/2;
	int kr=k/2;
	int nr=n/2;
	double* t_X1[nbe+1];
	double* t_X2[nbe+1];
	double* t_X3[nbe+1];
	// Temporary space allocation
	for (int i=nbe;(i>=1);--i){
		t_X1[i] = new double[mr*nr];
		t_X2[i] = new double[mr*kr];
		t_X3[i] = new double[kr*nr];
		nr/=2;
		mr/=2;
		kr/=2;
	}
	WinoMain(D, A, B, C, t_X1, t_X2, t_X3, lda, ldb, ldc, m, k, n,0 , nbe);
	kr=k/2;
	for (int i=nbe;(i>=1);--i){
		delete[] t_X1[i];
		delete[] t_X2[i];
		delete[] t_X3[i];
		kr/=2;
	}
	return C;
}

template<class Field>
inline typename Field::Element* FFFMMBLAS::operator ()( const Field& F, 
							const int m,
							const int n,
							const int k,
							//const enum FFFMMBLAS_COEF alpha,
							int alpha,
							typename Field::Element* const A,
							const int lda,
							typename Field::Element* const B,
							const int ldb, 
							//const enum FFFMMBLAS_COEF beta,
							int beta,
							typename Field::Element* C,
							const int ldc,
							const int nbe){
	

	ALPHA = alpha;
	BETA = beta;
	typedef typename Field::Element elt;
	int mr=m/2;
	int kr=k/2;
	int nr=n/2;
	// Threshold between GFq and double
	LinBox::integer charact;
	F.characteristic(charact);
	long long chara= Integer2long(charact);
	long long sf = ((long long)1<<53)/(chara*chara);
	elt* t_X1[nbe+1];
	elt* t_X2[nbe+1];
	elt* t_X3[nbe+1];
	// Temporary space allocation
	for (int i=nbe;(i>=1)&&(kr>sf);--i){
		t_X1[i] = new elt[mr*nr];
		t_X2[i] = new elt[mr*kr];
		t_X3[i] = new elt[kr*nr];
		nr/=2;
		mr/=2;
		kr/=2;
	}
	WinoMain(F, A,
		 B,
		 C,
		 t_X1,t_X2,t_X3,
		 lda,ldb,ldc,m,k,n,sf,nbe);
	kr=k/2;
	for (int i=nbe;(i>=1)&&(kr>sf);--i){
		delete[] t_X1[i];
		delete[] t_X2[i];
		delete[] t_X3[i];
		kr/=2;
	}
	return C;
}
