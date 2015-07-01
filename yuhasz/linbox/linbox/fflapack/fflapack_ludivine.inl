/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/fflapack/fflapack_ludivine.inl
 * Copyright (C) 2003 Clement Pernet
 *
 * Written by Clement Pernet <Clement.Pernet@imag.fr>
 *
 * See COPYING for license information.
 */

#ifndef MIN
#define MIN(a,b) (a<b)?a:b
#endif
#ifndef MAX
#define MAX(a,b) (a<b)?b:a
#endif

//---------------------------------------------------------------------
// TURBO: rank computation algorithm 
//---------------------------------------------------------------------

template <class Field>
inline size_t 
FFLAPACK::TURBO( const Field& F, const size_t M, const size_t N,		
		 typename Field::Element * A, const size_t lda ){
	
	if ( !(M && N) ) return 0;
	typedef typename Field::Element elt;
	static elt Mone;
	F.init(Mone, -1);
	// Column permutation
	size_t * P = new size_t[N];
	// Row Permutation
	size_t * rowP = new size_t[M];
	
	size_t q1,q2,q3,q3b,q4;
	q1=q2=q3=q3b=q4=0;
	
	size_t mo2 = (M>>1);
	size_t no2 = (N>>1);
	
	elt *NW = A; 
	elt *NE = A + no2 ;
	elt *SW = A + mo2*lda ;
	elt *SE = SW + no2 ;
		
	// Step 1: NW = L1.Q1.U1.P1
	size_t mloc = mo2;
	size_t nloc = no2;
	q1 = LUdivine( F, FflasNonUnit, mloc, no2, NW, lda, P, FflapackLQUP, rowP );
		
	// B1 = L^-1.NE
	solveLB( F, mo2, N-no2, q1, NW, lda, rowP, NE, lda);
	
	// NE = Q^-1.NE
	applyrowp( F, mo2, N-no2, NE, lda, rowP );		
	
	// SW = SW.P1
	flaswp(F,M-mo2,SW,lda,0,q1,P,1);
	
	// N1 = SW_{1,q1} . U1^-1
	ftrsm( F, FflasRight, FflasUpper, FflasNoTrans, FflasNonUnit, M-mo2, q1, F.one, NW, lda , SW, lda );
	
	// I1 = SW_{q1+1,n} - N1.G1  
	fgemm(F, FflasNoTrans, FflasNoTrans, M-mo2,  no2-q1, q1, Mone, SW, lda, NW+q1, lda, F.one, SW+q1, lda);
			
	// E1 = SE - N1.B1_{1,q1}
	fgemm( F, FflasNoTrans, FflasNoTrans, M-mo2, N-no2, q1, Mone, SW, lda, NE, lda, F.one, SE, lda);

	//Step 2: E1 = L2.Q2.U2.P2
	mloc = M-mo2;
	nloc = N-no2;
	q2 = LUdivine( F, FflasNonUnit, mloc, nloc, SE, lda, P+no2, FflapackLQUP, rowP+mo2 );
		
	// Updating P
	for (size_t i=no2;i<N;++i)
		P[i] += no2;

	// [I2;F2] = L2^-1.I1
	solveLB( F, mloc, no2-q1, q2, SE, lda, rowP+mo2, SW+q1, lda);

	// I1 = Q2^-1.I1
	applyrowp( F, mloc, no2-q1, SW+q1, lda, rowP+mo2 );

	// B1 = B1.P2
	flaswp(F,mo2,A,lda,no2,no2+q2,P,1); 

	//alternative: de 0 a q2 avant
	// N2 = B1_{q1+1,mo2} . V2^-1
	ftrsm(F, FflasRight, FflasUpper,FflasNoTrans,FflasNonUnit, mo2-q1, q2, F.one, SE, lda, NE+q1*lda,lda);
		
	// H2 = B1_{q1+1,mo2;q2,N-no2} - N2.E2  
	fgemm(F, FflasNoTrans, FflasNoTrans, mo2-q1, N-no2-q2, q2, Mone, NE+q1*lda, lda, SE+q2, lda, F.one, NE+q1*lda+q2, lda);

	// O2 = NW_{q1+1,mo2;q1+1,N-no2} = - N2.I2  
	fgemm(F, FflasNoTrans, FflasNoTrans, mo2-q1, no2-q1, q2, Mone, NE+q1*lda, lda, SW+q1, lda, F.zero,
	      NW+q1*(lda+1), lda);

	//Step 3: F2 = L3.Q3.U3.P3
	mloc = M-mo2-q2;
	nloc = no2-q1;
	q3 = LUdivine( F, FflasNonUnit, mloc, nloc, SW+q2*lda+q1, lda, P+q1, FflapackLQUP, rowP+mo2+q2 );
	
	// Updating P
	for (size_t i=q1;i<no2;++i)
		P[i] += q1;
		
	//Step 3bis: H2 = L3b.Q3b.U3b.P3b
	mloc = mo2-q1;
	nloc = N-no2-q2;
	size_t * rP3b = new size_t[mo2-q1];
	for (int j=0;j<mo2-q1;++j)
		rP3b[j]=0;
	q3b = LUdivine( F, FflasNonUnit, mloc, nloc, NE+q1*lda+q2, lda, 
			P+no2+q2, FflapackLQUP, rP3b );
	
	// Updating P
	for (size_t i=no2+q2;i<N;++i)
		P[i] += no2+q2;
		
	if (( q3 < no2-q1) && (q3b<mo2-q1)){
			
		// [O3;_] = L3b^-1.O2
		if (q3b>0){
			if ( mo2-q1 < N-no2-q2+q1) 
				// L is expanded to a Lower triangular matrix
				solveLB( F, mloc, no2-q1, q3b, NE+q1*lda+q2 , lda, rP3b, NW+q1*(lda+1), lda);
			else{
				cerr<<"USING SOLVELB2"<<endl;
				//no modification of L
				solveLB2( F, mloc, no2-q1, q3b, NE+q1*lda+q2 , lda, rP3b, NW+q1*(lda+1), lda);
			}

			// O2 = Q3b^-1.O2
			applyrowp( F, mloc, no2-q1, NW+q1*(lda+1), lda, rP3b );

			//updating rowP
			size_t tmp;
			for (int j=0;j<mo2-q1;++j)
				if (rP3b[j]!=j){
					//	cerr<<"(rP3b["<<j<<"]="<<rP3b[j]<<endl;
					tmp = rowP[j+q1];
					rowP[j+q1] = rowP[rP3b[j]+q1];
					rowP[rP3b[j]+q1] = tmp;
				}
				
			// X2 = X2.P3
			// Si plusieurs niveaux rec, remplacer X2 par [NW;I2]
			flaswp(F,mo2-q1-q3b,NW+(q1+q3b)*lda,lda,q1,q1+q3,P,1); 
				
			
			// A faire si plusieurs niveaux recursifs
			// B2 = B2.P3b
			//flaswp(F,q1,NE,lda,no2+q2,no2+q2+q3b,P,1); 
			// E2 = E2.P3b
			//flaswp(F,q2,SE+q2,lda,no2+q2,no2+q2+q3b,P,1); 
		}
					
		// N3 = X2 . D3^-1
		ftrsm( F, FflasRight, FflasUpper, FflasNoTrans, FflasNonUnit, mo2-q1-q3b, q3, F.one, SW+q2*lda+q1, lda ,NW+(q1+q3b)*lda+q1,lda);

		// T2 = T2 - N3.F3
		fgemm( F, FflasNoTrans, FflasNoTrans, mo2-q1-q3b, no2-q1-q3,q3, Mone, NW+(q1+q3b)*lda+q1, lda, SW+q2*lda+q3+q1, lda, F.one, NW+(q1+q3b)*lda+q1+q3, lda );
				
		//Step 4: T2 = L4.Q4.U4.P4
		mloc = mo2-q1-q3b;
		nloc = no2-q1-q3;
		q4 = LUdivine( F, FflasNonUnit, mloc, nloc, NW+(q1+q3b)*lda+q1+q3, lda, P+q1+q3, FflapackLQUP, rowP+q1+q3b );

		// Updating P
		for (size_t i=q1+q3;i<no2;++i)
			P[i] += q1+q3;

		// A faire si plusieurs niveaux recursifs
		// [G1;O3] = [G1;O3].P4
		//flaswp(F,q1+q3b,NE,lda,no2+q2,no2+q2+q3b,P,1); 
		// [I2;F3] = [I2;F3].P4
		//flaswp(F,q2,SE+q2,lda,no2+q2,no2+q2+q3b,P,1); 
	}
		
	// Necessaire:
	// 1 traiter les flaswp manquants
	// Facultatif:
	// 2 permutations de lignes doivent etre coherentes
	// 3 effectuer les dernieres permutations lignes et colonnes
	//cerr<<q1<<" "<<q2<<" "<<q3<<" "<<q3b<<" "<<q4<<endl;
	return q1+q2+q3+q3b+q4;
}

//---------------------------------------------------------------------
// LUdivine: LUP factorisation of A 
// P is the permutation matrix stored in an array of indexes
//---------------------------------------------------------------------

template <class Field>
inline size_t 
FFLAPACK::LUdivine( const Field& F, const enum FFLAS_DIAG Diag,
		    const size_t M, const size_t N,		
		    typename Field::Element * A, const size_t lda, size_t*P, 
		    const enum FFLAPACK_LUDIVINE_TAG LuTag, size_t *rowP){
	
	if ( !(M && N) ) return 0;
	typedef typename Field::Element elt;
	static elt Mone;
	F.init(Mone, -1);

#if DEBUG==1
	cerr<<"Entering LUdivine with LUtag, M, N ="<<LuTag<<" "
	    <<M<<" "<<N<<endl;
#endif
	size_t MN = MIN(M,N);

	if (MN == 1){
		size_t ip=0;
		while (ip<N && iszero(*(A+ip))){ip++;}
		if (ip==N){ // current row is zero
			*P=0;
			if (N==1){
				while (ip<M && iszero(*(A+ip*lda))){
					rowP[ip]=ip;
					ip++;
						
				}
				if (ip==M){
					return 0;
				}
				else{
					rowP[ip]=0;
					rowP[0]=ip;
					return 1;
				}
			}
			else{
				if (LuTag == FflapackLQUP)
					*rowP=0;
				return 0;
			}
		}
		*P=ip;
		if (ip!=0){
			// swap the pivot
			typename Field::Element tmp=*A;
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
		size_t R = LUdivine(F, Diag, Nup, N, A, lda, P, LuTag, rowP);

		typename Field::Element *Ar = A + Nup*lda; // SW
		typename Field::Element *Ac = A + R;     // NE
		typename Field::Element *An = Ar + R;    // SE
		if ( !R ){
			if (LuTag == FflapackSingular ) 
				return 0;
		}
		else{ // Updating Ar and An 

			// Apply the permutation on SW
			// Ar <- Ar.P
			flaswp(F,Ndown,Ar,lda,0,R,P,1);
#if DEBUG==3
			cerr<<"Apres le premier LUdivine rec et le laswp"<<endl;
			write_field(F,cerr,A,M,N,lda);
#endif			
		      
			if ( LuTag == FflapackLSP ){
				// The triangle is not contiguous
				// => Copy into a temporary
				size_t ldt=MAX(N-R,R);
				elt * T = new elt[R*ldt]; // contiguous matrix
				TriangleCopy( F, FflasUpper, Diag, R, 
					      T, ldt, A, lda); 
#if DEBUG==3
				cerr<<"Apres TriangleCopy T="<<endl;
				write_field(F,cerr,T,R,R,ldt);
#endif
				// Triangular block inversion of NW and apply to SW
				// Ar <- Ar.U1^-1
				ftrsm( F, FflasRight, FflasUpper, 
				       FflasNoTrans, Diag, Ndown, R, 
				       F.one, T, ldt, Ar, lda);
#if DEBUG==3
				cerr<<"Apres le Ftrsm"<<endl;
				write_field(F,cerr,A,M,N,lda);
				
				// Update of SE
				// An <- An - Ar*Ac
				cerr<<"Avant Rectangle copy Ndown,R,ldt="
				    <<Ndown<<" "<<R<<" "<<ldt<<endl;
				cerr<<"Ac="<<Ac<<" *Ac="<<*Ac<<endl;
#endif
				RectangleCopy( F, R, ldt, T, ldt, Ac,lda );
#if DEBUG==3
				cerr<<"Ar="<<endl;
				write_field(F,cerr,Ar,Ndown,R,lda);
				cerr<<"T="<<endl;
				write_field(F,cerr,T,R,ldt,ldt);
				cerr<<"An="<<endl;
				write_field(F,cerr,An,Ndown,ldt,lda);
#endif					
				fgemm( F, FflasNoTrans, FflasNoTrans,
				       Ndown, ldt, R, Mone, 
				       Ar, lda, T, ldt, F.one, An, lda);
#if DEBUG==3
				cerr<<"Apres le FFFMMBLAS"<<endl;
				write_field(F,cerr,A,M,N,lda);
#endif
				delete[] T;
			}
			else{
				// The triangle is contiguous
				
				// Triangular block inversion of NW and apply to SW
				// Ar <- Ar.U1^-1
				ftrsm( F, FflasRight, FflasUpper, 
				       FflasNoTrans, Diag, Ndown, R, 
				       F.one, A, lda, Ar, lda);
#if DEBUG==3
				cerr<<"Apres le Ftrsm"<<endl;
				write_field(F,cerr,A,M,N,lda);
#endif
				// Updating SE
				// An <- An - Ar*Ac
				fgemm( F, FflasNoTrans, FflasNoTrans, Ndown, N-R, R,
				       Mone, Ar, lda, Ac, lda, F.one, An, lda);
#if DEBUG==3
				cerr<<"Apres le FFFMMBLAS"<<endl;
				write_field(F,cerr,A,M,N,lda);
#endif
			}
		}
		// Recursive call on SE
		
		size_t R2=LUdivine(F, Diag, Ndown, N-R, An, lda,P+R,LuTag, 
				   (LuTag == FflapackLQUP)?rowP+Nup:rowP);
		for (size_t i=R;i<N;++i)
			P[i] += R;

		if (R2){
			//if (LuTag == FflapackLQUP)
			// Apply P2 on An
			// An <- An.P2
			flaswp(F, Nup, A, lda, R, R+R2, P, 1);
		}
		else if( LuTag == FflapackSingular )
			return 0;
#if DEBUG==3
		cerr<<"Apres le deuxieme LUdivine rec et flaswp"<<endl;
		write_field(F,cerr,A,M,N,lda);
#endif
		
		// Non zero rows permutations
		if ( LuTag == FflapackRank && (R<Nup)){
			for ( size_t i=R, i0=Nup; i0<(Nup+R2);i++,i0++){
				for (size_t j=0; j<N; j++){
					A[i*lda+j] = A[i0*lda+j];
				}
			}//FOR
		}
		else if ( LuTag == FflapackLQUP ){
			for (size_t i=Nup;i<M;i++)
				rowP[i] += Nup;
			if (R<Nup){
				// Permutation of the 0 rows
				
				for ( size_t i=Nup, j=R ; i<Nup+R2;++i, ++j){
					fcopy( F, N-j,A+j*(lda+1),1, A+i*lda+j,1);
					for (typename Field::Element *Ai = A+i*lda+j;
					     Ai!=A+i*lda+N; ++Ai)
						F.assign(*Ai, F.zero);
					size_t t = rowP[j];
					rowP[j]=rowP[i];
					rowP[i] = t;
				}
				
			}
#if DEBUG==3
			cerr<<"Apres row permutation"<<endl;
			write_field(F,cerr,A,M,N,lda);
#endif
		}
		return R+=R2;
			
	
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
FFLAPACK::LUdivine_construct( const Field& F, const enum FFLAS_DIAG Diag,
			      const size_t M, const size_t N,
			      typename Field::Element * B, const size_t ldb,
			      typename Field::Element * X, const size_t ldx,
			      typename Field::Element * A, const size_t lda, size_t* P,
			      size_t* nRowX, const size_t nRowXMax, size_t* nUsedRowX){

	static typename Field::Element Mone;
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
			typename Field::Element tmp=*A;
			*A = *(A+ip);
			*(A+ip) = tmp;
		}
		if ( Diag == FflasUnit ){
			// Normalisation of the row
			for (size_t k=1; k<N; k++)
				F.divin(*(A+k), *A);
		}
		if (N==1 && M>1 && *nRowX<nRowXMax){
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
			typename Field::Element * Xr = X+(*nRowX)*ldx;
			typename Field::Element * Ar = A + Nup*lda; //  SW
			typename Field::Element * Ac = A + Nup;     //  NE
			typename Field::Element * An = Ar + Nup;    //  SE
			while (*nRowX < *nUsedRowX + Ndown){ 
				// All available rows have been used, compute new ones 
				// number of new lines to be computed : nrowX except if
			        nNewRowX = MIN( *nRowX,  nRowXMax - *nRowX);
				// B <= B*B
#if DEBUG==3
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
#if DEBUG==3
				cerr<<"apres (X,BX), X="<<endl;
				write_field(F,cerr,X,nRowXMax,ldx,ldx);
#endif				
				// Copy of the Xr size_to Ar
				for ( typename Field::Element* Ai = A+(*nRowX)*lda; 
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
#if DEBUG==3
			cerr<<"Apres le premier LUdivinerec et le laswp"<<endl;
			write_field(F,cerr,A,M,N,lda);
#endif
			// Triangular block inversion of NW and apply to SW
			// Ar <- Ar.U1^-1
			ftrsm( F, FflasRight, FflasUpper, FflasNoTrans, Diag,
			       Ndown, R, F.one, A, lda, Ar, lda);
#if DEBUG==3
			cerr<<"Apres le Ftrsm"<<endl;
			write_field(F,cerr,A,M,N,lda);
#endif
			
			// Update of SE
			// An <- An - Ar*Ac
			fgemm( F, FflasNoTrans, FflasNoTrans, Ndown, N-Nup, Nup,
			       Mone, Ar, lda, Ac, lda, F.one, An, lda);
#if DEBUG==3
			cerr<<"Apres le FFFMMBLAS"<<endl;
			write_field(F,cerr,A,M,N,lda);
#endif
			// Recursive call on SE
			
			size_t R2 = LUdivine_construct(F, Diag, Ndown, N-Nup, B, ldb,
						       X, ldx, An, lda, P + Nup, 
						       nRowX, nRowXMax, nUsedRowX);
			for ( size_t i=Nup;i!=MN;i++) P[i] += Nup;
			
#if DEBUG==3
			cerr<<"avant d'appliquer le pivot: P=";
			cerr<<"Nup,R,R2="<<Nup<<", "<<R<<", "<<R2<<endl;
			for ( size_t i=0; i<Nup+R2;i++)
				cerr<<P[i]<<" ";
			cerr<<endl;
#endif
			flaswp(F, Nup, A, lda, Nup, Nup+R2, P, 1);
			
#if DEBUG==3
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
