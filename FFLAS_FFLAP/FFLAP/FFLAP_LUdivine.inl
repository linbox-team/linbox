/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
//----------------------------------------------------------------------
//                FFLAP: Finite Field Fast Linear Algebra Package
//                FFLAP_LUdivine: LUP factorisation
//----------------------------------------------------------------------
// by Clement PERNET ( clement.pernet@imag.fr )
// 30/04/2003
//----------------------------------------------------------------------

#ifndef MIN
#define MIN(a,b) (a<b)?a:b
#endif

//---------------------------------------------------------------------
// LUdivine: LUP factorisation of A 
// P is the permutation matrix stored in an array of indexes
//---------------------------------------------------------------------

template <class Field>
inline size_t 
FFLAP::LUdivine( const Field& F, const enum FFLAS_DIAG Diag,
		 const size_t M, const size_t N,		
		 typename Field::element * A, const size_t lda, size_t*P, 
		 const enum FFLAP_LUDIVINE_TAG LuTag){
	
	typedef typename Field::element elt;
	static elt Mone;
	F.neg(Mone, F.one);
	size_t MN = MIN(M,N);

	if (MN == 1){
		size_t ip=0;
		while (ip<N && iszero(*(A+ip))){ip++;}
		if (ip==N){ // current row is zero
			*P=0;
			return 0;
		}
		*P=ip;
		if (ip!=0){
			// swap the pivot
			typename Field::element tmp=*A;
			*A = *(A+ip);
			*(A+ip) = tmp;
		}
		if ( Diag == FflasUnit ){
			// Normalisation of the row
			for (size_t k=1; k<N; k++)
				F.divin(*(A+k), *A);
		}
		return 1;
	} 
	else{ // MN>1
		size_t Nup = MN>>1;
		size_t Ndown =  M - Nup;

		// Recursive call on NW
		size_t R = LUdivine(F, Diag, Nup, N, A, lda, P, LuTag);

		typename Field::element *Ar = A + Nup*lda; // SW
		typename Field::element *Ac = A + R;     // NE
		typename Field::element *An = Ar + R;    // SE

		if (  R==Nup || (LuTag != FflapLUP) ){ 

			if (R){
				// Apply the permutation on SW
				flaswp(F,Ndown,Ar,lda,0,R,P,1);
#if DEBUG==2
				cerr<<"Apres le premier LUdivine rec et le laswp"<<endl;
				write_field(F,cerr,A,M,N,lda);
#endif			
				// Triangular block inversion of NW and apply to SW
				// Ar <- Ar.U1^-1
				if ( LuTag == FflapLSP ){
					size_t ldt=N-R;
					elt * T = new elt[R*ldt]; // contiguous matrix
					TriangleCopy( F, Diag, R, T, ldt, A, lda); 
					cerr<<"Apres TriangleCopy T="<<endl;
					write_field(F,cerr,T,R,R,ldt);
					ftrsm( F, FflasRight, FflasUpper, 
					       FflasNoTrans, Diag, Ndown, R, 
					       F.one, T, ldt, Ar, lda);
#if DEBUG==2
					cerr<<"Apres le Ftrsm"<<endl;
					write_field(F,cerr,A,M,N,lda);
#endif
					// Update of SE
					// An <- An - Ar*Ac
					RectangleCopy( F, R, ldt, T, ldt, Ac,lda );
					fgemm( F, FflasNoTrans, FflasNoTrans,
					       Ndown, ldt, R, Mone, 
					       Ar, lda, T, ldt, F.one, An, lda);
#if DEBUG==2
					cerr<<"Apres le FFFMMBLAS"<<endl;
					write_field(F,cerr,A,M,N,lda);
#endif
					delete[] T;
				}
				else{
					ftrsm( F, FflasRight, FflasUpper, 
					       FflasNoTrans, Diag, Ndown, R, 
					       F.one, A, lda, Ar, lda);
#if DEBUG==2
					cerr<<"Apres le Ftrsm"<<endl;
					write_field(F,cerr,A,M,N,lda);
#endif
					// Update of SE
					// An <- An - Ar*Ac
					fgemm( F, FflasNoTrans, FflasNoTrans, Ndown, N-R, R,
					       Mone, Ar, lda, Ac, lda, F.one, An, lda);
#if DEBUG==2
					cerr<<"Apres le FFFMMBLAS"<<endl;
					write_field(F,cerr,A,M,N,lda);
#endif
				}
			}
			// Recursive call on SE
			
			size_t R2=LUdivine(F, Diag, Ndown, N-R, An, lda,P+R,LuTag);
			if (R2){
				for (size_t i=R;i!=R+R2;i++)
					P[i] += R;
				// Apply P on An
				flaswp(F, Nup, A, lda, R, R+R2, P, 1);
			}
#if DEBUG==2 
			cerr<<"Apres le deuxieme LUdivine rec et flaswp"<<endl;
			write_field(F,cerr,A,M,N,lda);
#endif

			// Non zero rows permutations
			if ( LuTag == FflapRank && (R<Nup)){
				for ( size_t i=R, i0=Nup; i0<(Nup+R2);i++,i0++){
					for (size_t j=0; j<N; j++){
						A[i*lda+j] = A[i0*lda+j];
					}
				}//FOR
#if DEBUG==2 
				cerr<<"Apres row permutation"<<endl;
				write_field(F,cerr,A,M,N,lda);
#endif
			}//if (R1<Nup)
			return R+=R2;
		}
		else{ return R;} 
		// Rank deficient matrices can only be factorized 
		// under the condition: the first R rows are linearly independent
		// If not, the lower block is never factorized as soon as the
		// upper block is rank defficient
	}
}


//---------------------------------------------------------------------
// LUdivine_construct: (Specialisation of LUdivine)
// LUP factorisation of X, the Krylov base matrix of A^t and v, in A.
// X contains the nRowX first vectors v, vA, .., vA^{nRowX-1}
// A contains the LUP factorisation of the nUsedRowX first row of X.
// When all rows of X have been factorized in A, and rank is full,
// then X is updated by the following scheme: X <= ( X; X.B ), where
// B = A^2^i.
// This enables to make use of Matrix multiplication, and stop computing
// Krylov vector, when the rank is not longer full.
// P is the permutation matrix stored in an array of indexes
//---------------------------------------------------------------------

template <class Field>
size_t
FFLAP::LUdivine_construct( const Field& F, const enum FFLAS_DIAG Diag,
			   const size_t M, const size_t N,
			   typename Field::element * B, const size_t ldb,
			   typename Field::element * X, const size_t ldx,
			   typename Field::element * A, const size_t lda, size_t* P,
			   size_t* nRowX, const size_t nRowXMax, size_t* nUsedRowX){

	static typename Field::element Mone;
	F.neg(Mone, F.one);
	size_t MN = MIN(M,N);

	if (MN == 1){ 
		size_t ip=0;
		while (ip<N && iszero(*(A+ip))){ip++;}
		if (ip==N){ // current row is zero
			*P=0;
			return 0;
		}
		(*nUsedRowX)++;
		*P=ip;
		if (ip!=0){
			// swap the pivot
			typename Field::element tmp=*A;
			*A = *(A+ip);
			*(A+ip) = tmp;
		}
		if ( Diag == FflasUnit ){
			// Normalisation of the row
			for (size_t k=1; k<N; k++)
				F.divin(*(A+k), *A);
		}
		if (N==1 && M>1){
			F.mul(*(A+lda),*A, *B);
		}
		return 1;
	}
	else{ // MN>1
		size_t Nup = MN>>1;
		size_t Ndown =  M - Nup;
		
		// Recursive call on NW
		size_t R = LUdivine_construct(F, Diag, Nup, N, B, ldb, X, ldx, A, lda,
					      P, nRowX, nRowXMax, nUsedRowX);
		size_t nNewRowX;
		if (R==Nup){
			typename Field::element * Xr = X+(*nRowX)*ldx;
			typename Field::element * Ar = A + Nup*lda; //  SW
			typename Field::element * Ac = A + Nup;     //  NE
			typename Field::element * An = Ar + Nup;    //  SE
			while (*nRowX < *nUsedRowX + Ndown){ 
				// All available rows have been used, compute new ones 
				// number of new lines to be computed : nrowX except if
			        nNewRowX = MIN( *nRowX,  nRowXMax - *nRowX);
				// B <= B*B
#if DEBUG==2
				cerr<<"nNewRowX ="<<nNewRowX<<" nRowX = "
				    <<(*nRowX)<<"  nRowXMax="<< nRowXMax<<endl;
#endif
				if (*nRowX > 1){
					fsquare(F, FflasNoTrans,
						       ldb, F.one, B, ldb, 
						       F.zero, B, ldb);
				}
				
				fgemm(F,FflasNoTrans,FflasTrans, nNewRowX , ldb,
				      ldb, F.one, X, ldx, B, ldb, F.zero, Xr, ldx);
#if DEBUG==2
				cerr<<"apres (X,BX), X="<<endl;
				write_field(F,cerr,X,nRowXMax,ldx,ldx);
#endif				
				// Copy of the Xr size_to Ar
				for ( typename Field::element* Ai = A+(*nRowX)*lda; 
				      Ai < A+(*nRowX+nNewRowX)*lda; 
				      Ai+=lda, Xr+=ldx){
					for (size_t j=0;j<ldb;j++){
						Ai[j] = Xr[j];
					}
				}
				*nRowX += nNewRowX;
			}
			
			// Apply the permutation on SW
			flaswp(F,Ndown,Ar,lda,0,R,P,1);
#if DEBUG==2
			cerr<<"Apres le premier LUdivinerec et le laswp"<<endl;
			write_field(F,cerr,A,M,N,lda);
#endif
			// Triangular block inversion of NW and apply to SW
			// Ar <- Ar.U1^-1
			ftrsm( F, FflasRight, FflasUpper, FflasNoTrans, Diag,
			       Ndown, R, F.one, A, lda, Ar, lda);
#if DEBUG==2
			cerr<<"Apres le Ftrsm"<<endl;
			write_field(F,cerr,A,M,N,lda);
#endif
			
			// Update of SE
			// An <- An - Ar*Ac
			fgemm( F, FflasNoTrans, FflasNoTrans, Ndown, N-Nup, Nup,
			       Mone, Ar, lda, Ac, lda, F.one, An, lda);
#if DEBUG==2
			cerr<<"Apres le FFFMMBLAS"<<endl;
			write_field(F,cerr,A,M,N,lda);
#endif
			// Recursive call on SE
			
			size_t R2 = LUdivine_construct(F, Diag, Ndown, N-Nup, B, ldb,
						       X, ldx, An, lda, P + Nup, 
						       nRowX, nRowXMax, nUsedRowX);
			for ( size_t i=Nup;i!=MN;i++) P[i] += Nup;
			
#if DEBUG==2
			cerr<<"avant d'appliquer le pivot: P=";
			cerr<<"Nup,R,R2="<<Nup<<", "<<R<<", "<<R2<<endl;
			for ( size_t i=0; i<Nup+R2;i++)
				cerr<<P[i]<<" ";
			cerr<<endl;
#endif
			flaswp(F, Nup, A, lda, Nup, Nup+R2, P, 1);
			
#if DEBUG==2
			cerr<<"Apres le deuxieme LSP rec et flaswp"<<endl;
			write_field(F,cerr,A,M,N,lda);
#endif
			
			return R+=R2;
		}
		else{ return R;}
		// Rank deficient matrices can only be factorized 
		// under the condition: the first R rows are linearly independent
		// If not, the lower block is never factorized as soon as the
		// upper block is rank defficient
	}
}
