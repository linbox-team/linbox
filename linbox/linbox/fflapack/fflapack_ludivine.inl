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
LinBox::FFLAPACK::TURBO( const Field& F, const size_t M, const size_t N,		
			 typename Field::Element * NW, const size_t ld1,
			 typename Field::Element * NE, const size_t ld2,
			 typename Field::Element * SW, const size_t ld3,
			 typename Field::Element * SE, const size_t ld4	 ){
	
	if ( !(M && N) ) return 0;
	typedef typename Field::Element elt;
	static elt Mone, one, zero;
	F.init(Mone, -1);
	F.init(one,1);
	F.init(zero,0);
	// Column permutation
	size_t * P = new size_t[N];
	// Row Permutation
	size_t * rowP = new size_t[M];
	
	size_t q1,q2,q3,q3b,q4;
	q1=q2=q3=q3b=q4=0;
	
	size_t mo2 = (M>>1);
	size_t no2 = (N>>1);
	
		
	// Step 1: NW = L1.Q1.U1.P1
	size_t mloc = mo2;
	size_t nloc = no2;
// 	Timer tim;
// 	tim.clear();
// 	tim.start();
	q1 = LUdivine( F, FflasNonUnit, mloc, no2, NW, ld1, P, FflapackLQUP, rowP );
	
// 	tim.stop();
// 	cerr<<"LQUP1:"<<tim.realtime()<<endl;
// 	tim.start();
#if DEBUG
	cerr<<"NW= L1.Q1.U1.P1"<<endl;
	write_field(F,cerr,NW,mloc,no2,ld1);
#endif	
	// B1 = L^-1.NE
#if DEBUG
	cerr<<"avant B1 = L^-1.NE"<<endl;
	write_field(F,cerr,NE,mloc,N-no2,ld2);
#endif	
	solveLB( F, FflasLeft, mo2, N-no2, q1, NW, ld1, rowP, NE, ld2);
#if DEBUG
	cerr<<"B1 = L^-1.NE"<<endl;
	write_field(F,cerr,NE,mloc,N-no2,ld2);
#endif	

	// NE = Q^-1.NE
	
	applyP( F, FflasLeft, FflasNoTrans, N-no2, 0, mo2, NE, ld2, rowP );		
#if DEBUG
	cerr<<"NE=Q^-1.NE"<<endl;
	write_field(F,cerr,NE,mloc,N-no2,ld2);
#endif	

	// SW = SW.P1
	applyP( F, FflasRight, FflasTrans, M-mo2, 0, q1, SW, ld3, P );		
#if DEBUG
	cerr<<"SW = SW.P1"<<endl;
	write_field(F,cerr,SW,M-mo2,no2,ld3);
#endif	
	//flaswp(F,M-mo2,SW,ld3,0,q1,P,1);

// 	tim.stop();
// 	cerr<<"L^-1:"<<tim.realtime()<<endl;
// 	tim.start();
	
	// N1 = SW_{1,q1} . U1^-1
	ftrsm( F, FflasRight, FflasUpper, FflasNoTrans, FflasNonUnit, M-mo2, q1, one, NW, ld1 , SW, ld3 );
#if DEBUG
	cerr<<" N1 = SW_{1,q1} . U1^-1"<<endl;
	write_field(F,cerr,SW,M-mo2,no2,ld3);
#endif	

// 	tim.stop();
// 	cerr<<"trsm:"<<tim.realtime()<<endl;
// 	tim.start();
	
	// I1 = SW_{q1+1,n} - N1.G1  
	fgemm(F, FflasNoTrans, FflasNoTrans, M-mo2,  no2-q1, q1, Mone, SW, ld3, NW+q1, ld1, one, SW+q1, ld3);
#if DEBUG
	cerr<<" I1 = SW_{q1+1,n} - N1.G1"<<endl;
	write_field(F,cerr,SW,M-mo2,no2,ld3);
#endif	
// 	tim.stop();
// 	cerr<<"fgemm1:"<<tim.realtime()<<endl;
// 	tim.start();
			
	// E1 = SE - N1.B1_{1,q1}
	fgemm( F, FflasNoTrans, FflasNoTrans, M-mo2, N-no2, q1, Mone, SW, ld3, NE, ld2, one, SE, ld4);
#if DEBUG
	cerr<<"  E1 = SE - N1.B1_{1,q1}"<<endl;
	write_field(F,cerr,SE,M-mo2,N-no2,ld4);
#endif	
// 	tim.stop();
// 	cerr<<"fgemm2:"<<tim.realtime()<<endl;
// 	tim.start();

	//Step 2: E1 = L2.Q2.U2.P2
	mloc = M-mo2;
	nloc = N-no2;
	q2 = LUdivine( F, FflasNonUnit, mloc, nloc, SE, ld4, P+no2, FflapackLQUP, rowP+mo2 );
#if DEBUG
	cerr<<"  E1 = L2.Q2.U2.P2"<<endl;
	write_field(F,cerr,SE,M-mo2,N-no2,ld4);	
#endif	
// 	tim.stop();
// 	cerr<<"LQUP2:"<<tim.realtime()<<endl;
// 	tim.start();

	// [I2;F2] = L2^-1.I1
	solveLB( F, FflasLeft, mloc, no2-q1, q2, SE, ld4, rowP+mo2, SW+q1, ld3);
#if DEBUG
	cerr<<"  [I2;F2] = L2^-1.I1"<<endl;
	write_field(F,cerr,SW,M-mo2,no2,ld3);	
#endif	

	// I1 = Q2^-1.I1
	applyP( F, FflasLeft, FflasNoTrans, no2-q1, 0, mloc, SW+q1, ld3, rowP+mo2 );
#if DEBUG
	cerr<<"I1 = Q2^-1.I1"<<endl;
	write_field(F,cerr,SW,mloc,no2,ld3);
 #endif	

	// B1 = B1.P2
	applyP( F, FflasRight, FflasTrans, mo2, 0, q2, NE, ld2, P+no2 );
#if DEBUG
	cerr<<"I1 = Q2^-1.I1"<<endl;
	write_field(F,cerr,NE,mo2,N-no2,ld2);
#endif	
 	//flaswp(F,mo2,NE,ld2,0,q2,P+no2,1); 
	// Updating P
	//	for (size_t i=no2;i<N;++i)
	//	P[i] += no2;
// 	tim.stop();
// 	cerr<<"L2^-1:"<<tim.realtime()<<endl;
// 	tim.start();

	//alternative: de 0 a q2 avant
	// N2 = B1_{q1+1,mo2} . V2^-1
	ftrsm(F, FflasRight, FflasUpper,FflasNoTrans,FflasNonUnit, mo2-q1, q2, one, SE, ld4, NE+q1*ld2,ld2);
// 	tim.stop();
// 	cerr<<"trsm2:"<<tim.realtime()<<endl;
// 	tim.start();
		
	// H2 = B1_{q1+1,mo2;q2,N-no2} - N2.E2  
	fgemm(F, FflasNoTrans, FflasNoTrans, mo2-q1, N-no2-q2, q2, Mone, NE+q1*ld2, ld2, SE+q2, ld4, one, NE+q1*ld2+q2, ld2);

// 	tim.stop();
// 	cerr<<"fgemm12:"<<tim.realtime()<<endl;
// 	tim.start();
	// O2 = NW_{q1+1,mo2;q1+1,N-no2} = - N2.I2  
	fgemm(F, FflasNoTrans, FflasNoTrans, mo2-q1, no2-q1, q2, Mone, NE+q1*ld2, ld2, SW+q1, ld3, zero,
	      NW+q1*(ld1+1), ld1);
// 	tim.stop();
// 	cerr<<"fgemm22:"<<tim.realtime()<<endl;
// 	tim.start();

	//Step 3: F2 = L3.Q3.U3.P3
	mloc = M-mo2-q2;
	nloc = no2-q1;
	q3 = LUdivine( F, FflasNonUnit, mloc, nloc, SW+q2*ld3+q1, ld3, P+q1, FflapackLQUP, rowP+mo2+q2 );
	
	// Updating P
	for (size_t i=q1;i<no2;++i)
		P[i] += q1;
		
	//Step 3bis: H2 = L3b.Q3b.U3b.P3b
	mloc = mo2-q1;
	nloc = N-no2-q2;
	size_t * rP3b = new size_t[mo2-q1];
	for (size_t j=0;j<mo2-q1;++j)
		rP3b[j]=0;
	q3b = LUdivine( F, FflasNonUnit, mloc, nloc, NE+q1*ld2+q2, ld2, P+no2+q2, FflapackLQUP, rP3b );
	
// 	tim.stop();
// 	cerr<<"LQUP3et3bis:"<<tim.realtime()<<endl;
// 	tim.start();
		
	if (( q3 < no2-q1) && (q3b<mo2-q1)){
			
		// [O3;_] = L3b^-1.O2
		if (q3b>0){
			if ( mo2-q1 < N-no2-q2+q1) 
				// L is expanded to a Lower triangular matrix
				solveLB( F, FflasLeft,mloc, no2-q1, q3b, NE+q1*ld2+q2 , ld2, rP3b, NW+q1*(ld1+1), ld1);
			else{
				cerr<<"USING SOLVELB2"<<endl;
				//no modification of L
				solveLB2( F, FflasLeft,mloc, no2-q1, q3b, NE+q1*ld2+q2 , ld2, rP3b, NW+q1*(ld1+1), ld1);
			}
#if DEBUG
			cerr<<"O2 avant="<<endl;
			write_field(F,cerr,NW+q1*(ld1+1),mloc,no2-q1,ld1);
#endif	
	
			// O2 = Q3b^-1.O2
			applyP( F, FflasLeft, FflasNoTrans, no2-q1, 0, mloc, NW+q1*(ld1+1), ld1, rP3b );
#if DEBUG
			cerr<<"O2 apres="<<endl;
			write_field(F,cerr,NW+q1*(ld1+1),mloc,no2-q1,ld1);
#endif	
	
			//updating rowP
			size_t tmp;
			for (size_t j=0;j<mo2-q1;++j)
				if (rP3b[j]!=j){
					//	cerr<<"(rP3b["<<j<<"]="<<rP3b[j]<<endl;
					tmp = rowP[j+q1];
					rowP[j+q1] = rowP[rP3b[j]+q1];
					rowP[rP3b[j]+q1] = tmp;
				}
				
			// X2 = X2.P3
			// Si plusieurs niveaux rec, remplacer X2 par [NW;I2]
			applyP( F, FflasRight, FflasTrans, mo2-q1-q3b, q1, q1+q3, NW+(q1+q3b)*ld1, ld1, P );
			//flaswp(F,mo2-q1-q3b,NW+(q1+q3b)*ld1,ld1,q1,q1+q3,P,1); 
				
			// Updating P
			//for (size_t i=no2+q2;i<N;++i)
			//P[i] += no2+q2;
	
			// A faire si plusieurs niveaux recursifs
			// B2 = B2.P3b
			//flaswp(F,q1,NE,lda,no2+q2,no2+q2+q3b,P,1); 
			// E2 = E2.P3b
			//flaswp(F,q2,SE+q2,lda,no2+q2,no2+q2+q3b,P,1); 
		}
					
		// N3 = X2 . D3^-1
		ftrsm( F, FflasRight, FflasUpper, FflasNoTrans, FflasNonUnit, mo2-q1-q3b, q3, one, SW+q2*ld3+q1, ld3 ,NW+(q1+q3b)*ld1+q1,ld1);

		// T2 = T2 - N3.F3
		fgemm( F, FflasNoTrans, FflasNoTrans, mo2-q1-q3b, no2-q1-q3,q3, Mone, NW+(q1+q3b)*ld1+q1, ld1, SW+q2*ld3+q3+q1, ld3, one, NW+(q1+q3b)*ld1+q1+q3, ld1 );
				
		//Step 4: T2 = L4.Q4.U4.P4
		mloc = mo2-q1-q3b;
		nloc = no2-q1-q3;
		q4 = LUdivine( F, FflasNonUnit, mloc, nloc, NW+(q1+q3b)*ld1+q1+q3, ld1, P+q1+q3, FflapackLQUP, rowP+q1+q3b );

		// Updating P
		//for (size_t i=q1+q3;i<no2;++i)
		//	P[i] += q1+q3;

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
// LUdivine: LQUP factorisation of A. In-place computation:
// A is replaced by (LQ) and U.
// Q,P are permutation matrices stored in arrays of indexes in the
// lapack permutation storage style
//---------------------------------------------------------------------
template <class Field>
inline size_t 
LinBox::FFLAPACK::LUdivine( const Field& F, const enum FFLAS_DIAG Diag,
		    const size_t M, const size_t N,		
		    typename Field::Element * A, const size_t lda, size_t*P, 
		    const enum FFLAPACK_LUDIVINE_TAG LuTag, size_t *Q){
	
	if ( !(M && N) ) return 0;
	typedef typename Field::Element elt;
	static elt Mone, one, zero;
	F.init(Mone, -1);
	F.init(one,1);
	F.init(zero,0);

#if DEBUG==3
	cerr<<"Entering LUdivine with LUtag, M, N ="<<LuTag<<" "
	    <<M<<" "<<N<<endl;
#endif
	size_t MN = MIN(M,N);

	if (MN == 1){
		size_t ip=0;
		while (ip<N && F.isZero(*(A+ip))){ip++;}
		if (LuTag == FflapackLQUP)
			*Q=0;
		if (ip==N){ // current row is zero
			//*P=0;
			if (N==1){
				while (ip<M && F.isZero(*(A+ip*lda))){
					if (LuTag == FflapackLQUP)
						Q[ip]=ip;
					ip++;
						
				}
				if (ip==M){
					return 0;
				}
				else{
					if (LuTag == FflapackLQUP){
						//Q[ip]=0;
						*Q=ip;
					}
					if ( Diag == FflasNonUnit ){
						elt piv = *(A+ip*lda);
						while(++ip<M){
							F.divin(*(A+lda*ip), piv);
						}
					}
					return 1;
				}
			}
			else{
				if (LuTag == FflapackLQUP)
					*Q=0;
				return 0;
			}
		}
		*P=ip;
		if (ip!=0){
			// swap the pivot
			//size_t tmpi = P[0];
			//P[0] = ip;
			//P[ip] = tmpi;
			typename Field::Element tmp=*A;
			*A = *(A+ip);
			*(A+ip) = tmp;
		}
		if ( Diag == FflasUnit ){
			// Normalisation of the row
			for (size_t k=1; k<N; k++)
				F.divin(*(A+k), *A);
		}
		else if ( N==1 )
			while(++ip<M)
				F.divin(*(A+lda*ip), *A);
		return 1;
	}
		
	else{ // MN>1
		size_t Nup = MN>>1;
		size_t Ndown =  M - Nup;

		// Recursive call on NW
		size_t R = LUdivine/*_base*/(F, Diag, Nup, N, A, lda, P, LuTag, Q);

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
			applyP( F, FflasRight, FflasTrans, Ndown, 0, R, Ar, lda, P ); 
			//flaswp(F,Ndown,Ar,lda,0,R,P,1);
		      
			if ( LuTag == FflapackLSP ){
				// The triangle is not contiguous
				// => Copy into a temporary
				size_t ldt=MAX(N-R,R);
				elt * T = new elt[R*ldt]; // contiguous matrix
#ifdef CHECK_MEMORY
				CUR_MEMORY+= (R*ldt)*sizeof(elt);
				if (CUR_MEMORY > MAX_MEMORY)
					MAX_MEMORY=CUR_MEMORY;
#endif
				TriangleCopy( F, FflasUpper, Diag, R, 
					      T, ldt, A, lda); 
				// Triangular block inversion of NW and apply to SW
				// Ar <- Ar.U1^-1
				ftrsm( F, FflasRight, FflasUpper, 
				       FflasNoTrans, Diag, Ndown, R, 
				       one, T, ldt, Ar, lda);
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
				       Ar, lda, T, ldt, one, An, lda);
#if DEBUG==3
				cerr<<"Apres le FFFMMBLAS"<<endl;
				write_field(F,cerr,A,M,N,lda);
#endif
				delete[] T;
#ifdef CHECK_MEMORY
				CUR_MEMORY-=(R*ldt)*sizeof(elt);
#endif
			}
			else{
				// The triangle is contiguous
				
				// Triangular block inversion of NW and apply to SW
				// Ar <- Ar.U1^-1
				ftrsm( F, FflasRight, FflasUpper, 
				       FflasNoTrans, Diag, Ndown, R, 
				       one, A, lda, Ar, lda);
#if DEBUG==3
				cerr<<"Apres le Ftrsm"<<endl;
				write_field(F,cerr,A,M,N,lda);
#endif
				// Updating SE
				// An <- An - Ar*Ac
				fgemm( F, FflasNoTrans, FflasNoTrans, Ndown, N-R, R,
				       Mone, Ar, lda, Ac, lda, one, An, lda);
#if DEBUG==3
				cerr<<"Apres le FFFMMBLAS"<<endl;
				write_field(F,cerr,A,M,N,lda);
#endif
			}
		}
		// Recursive call on SE
		size_t R2=LUdivine/*_base*/( F, Diag, Ndown, N-R, An, lda,P+R,LuTag, 
				    (LuTag == FflapackLQUP)?Q+Nup:Q);
		for (size_t i=R;i<R+R2;++i)
			P[i] += R;

		if (R2){
			//if (LuTag == FflapackLQUP)
			// Apply P2 on An
			// An <- An.P2
			applyP( F, FflasRight, FflasTrans, Nup, R, R+R2, A, lda, P ); 
			//flaswp(F, Nup, A, lda, R, R+R2, P, 1);
		}
		else if( LuTag == FflapackSingular )
			return 0;
#if DEBUG==3
		cerr<<"Apres le deuxieme LUdivine rec et flaswp"<<endl;
		write_field(F,cerr,A,M,N,lda);
#endif
		
		// Non zero row permutations
		if ( LuTag == FflapackLQUP ){
			for (size_t i=Nup;i<M;i++)
				Q[i] += Nup;
			if (R<Nup){
				// Permutation of the 0 rows
				
				for ( size_t i=Nup, j=R ; i<Nup+R2;++i, ++j){
					fcopy( F, N-j,A+j*(lda+1),1, A+i*lda+j,1);
					for (typename Field::Element *Ai = A+i*lda+j;
					     Ai!=A+i*lda+N; ++Ai)
						F.assign(*Ai, zero);
					// Essai de lapack_style pour les lignes
					size_t t = Q[j];
					Q[j]=Q[i];
					Q[i] = t;
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
// LUP factorisation of the Krylov base matrix of A^t and v.
// When all rows have been factorized in A, and rank is full,
// then new krylov vectors are computed and then triangularized
// P is the permutation matrix stored in the lapack style
// nRowX is the number of Krylov vectors already computed,
// nUsedRowX is the number of Krylov vectors already triangularized
//---------------------------------------------------------------------

template <class Field>
size_t
LinBox::FFLAPACK::LUdivine_construct( const Field& F, const enum FFLAS_DIAG Diag,
				      const size_t M, const size_t N,
				      const typename Field::Element * A, const size_t lda,
				      typename Field::Element * X, const size_t ldx,
				      typename Field::Element * u, size_t* P,
				      bool computeX, const enum FFLAPACK_MINPOLY_TAG MinTag = FflasDense,
				      const size_t kg_mc =0, const size_t kg_mb=0, const size_t kg_j=0){

	static typename Field::Element Mone, one, zero;
	F.init(Mone, -1);
	F.init(one, 1);
	F.init(zero,0);
	size_t MN = MIN(M,N);

	if (MN == 1){ 
		size_t ip=0;
		while (ip<N && F.isZero(*(X+ip))){ip++;}
		if (ip==N){ // current row is zero
			*P=0;
			return 0;
		}
		*P=ip;
		if (ip!=0){
			// swap the pivot
			typename Field::Element tmp=*X;
			*X = *(X+ip);
			*(X+ip) = tmp;
		}
		if ( Diag == FflasUnit ){
			// Normalisation of the row
			for (size_t k=1; k<N; k++)
				F.divin(*(X+k), *X);
		}
 		if (N==1 && M>1 && computeX)// Only appends when A is 1 by 1
			F.mul(*(X+ldx),*X, *A);
				
		return 1;
	}
	else{ // MN>1
		size_t Nup = MN>>1;
		size_t Ndown =  M - Nup;
		
		// Recursive call on NW
		size_t R = LUdivine_construct(F, Diag, Nup, N, A, lda, X, ldx, u, 
					       P, computeX, MinTag, kg_mc, kg_mb, kg_j );
		if (R==Nup){
			typename Field::Element * Xr = X + Nup*ldx; //  SW
			typename Field::Element * Xc = X + Nup;     //  NE
			typename Field::Element * Xn = Xr + Nup;    //  SE
			typename Field::Element * Xi = Xr;
			if ( computeX )
				for (size_t i=0; i< Ndown; ++i, Xi+=ldx){
					if (MinTag == FflapackDense)
						fgemv(F, FflasNoTrans, N, N, one, 
						      A, lda, u, 1, zero, Xi,1);
  					else // Keller-Gehrig Fast algorithm's matrix
 						fgemv_kgf( F, N, A, lda, u, 1, Xi, 1, 
  							   kg_mc, kg_mb, kg_j );
					// cerr<<"u="<<endl;
// 					write_field(F,cerr,u,1,N,N);
//  					cerr<<"X="<<endl;
//  					write_field(F,cerr,Xi,1,N,N);
					
					fcopy(F, N, u,1,Xi, 1);
					// cerr<<"Apres fgemv X="<<endl;
					// write_field( F, cerr, X, 1, N, N );
				} 
 // 			cerr<<"X="<<endl;
// 			write_field( F, cerr, X, M, N, ldx );
			//cerr<<"B="<<endl;
			//write_field( F, cerr, A, N, N, lda );
			
			// Apply the permutation on SW
			applyP( F, FflasRight, FflasTrans, Ndown, 0, R, Xr, ldx, P); 
			//flaswp(F,Ndown,Ar,ldx,0,R,P,1);
#if DEBUG==3
			cerr<<"Apres le premier LUdivinerec et le laswp"<<endl;
			write_field(F,cerr,X,M,N,ldx);
#endif
			// Triangular block inversion of NW and apply to SW
			// Xr <- Xr.U1^-1
			ftrsm( F, FflasRight, FflasUpper, FflasNoTrans, Diag,
			       Ndown, R, one, X, ldx, Xr, ldx);
#if DEBUG==3
			cerr<<"Apres le Ftrsm"<<endl;
			write_field(F,cerr,X,M,N,ldx);
#endif
			
			// Update of SE
			// Xn <- Xn - Xr*Xc
			fgemm( F, FflasNoTrans, FflasNoTrans, Ndown, N-Nup, Nup,
			       Mone, Xr, ldx, Xc, ldx, one, Xn, ldx);
#if DEBUG==3
			cerr<<"Apres le FFFMMBLAS"<<endl;
			write_field(F,cerr,X,M,N,ldx);
#endif
			// Recursive call on SE
			
			size_t R2 = LUdivine_construct(F, Diag, Ndown, N-Nup, A, lda,
							Xn, ldx, u, P + Nup, 
							false, MinTag, kg_mc, kg_mb, kg_j);
			for ( size_t i=R;i<R+R2;++i) P[i] += R;
			
#if DEBUG==3
			cerr<<"Nup,R,R2="<<Nup<<", "<<R<<", "<<R2<<endl;
			cerr<<"avant d'appliquer le pivot: P=";
			for ( size_t i=0; i<N;i++)
				cerr<<P[i]<<" ";
			cerr<<endl;
#endif
			applyP( F, FflasRight, FflasTrans, Nup, R, R+R2, X, ldx, P); 
			//flaswp(F, Nup, X, ldx, Nup, Nup+R2, P, 1);
			
#if DEBUG==3
			cerr<<"Apres le deuxieme LSP rec et flaswp"<<endl;
			write_field(F,cerr,X,M,N,ldx);
#endif
			
			return R+=R2;
		}
		else 
			return R;
		// Rank deficient matrices can only be factorized 
		// under the condition: the first R rows are linearly independent
		// If not, the lower block is never factorized as soon as the
		// upper block is rank defficient
	}
}
