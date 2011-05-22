/* ffpack/ffpack_ludivine.inl
 * Copyright (C) 2005 Clement Pernet
 *
 * Written by Clement Pernet <Clement.Pernet@imag.fr>
 *
 * See COPYING for license information.
 */

#ifndef __LINBOX_ffpack_ludivine_INL
#define __LINBOX_ffpack_ludivine_INL

#ifndef MIN
#define MIN(a,b) (a<b)?a:b
#endif
#ifndef MAX
#define MAX(a,b) (a<b)?b:a
#endif

//#define LB_DEBUG

template<class Field>
inline size_t
FFPACK::LUdivine_gauss( const Field& F, const FFLAS_DIAG Diag,
			const size_t M, const size_t N,		
			typename Field::Element * A, const size_t lda, size_t*P, 
			size_t *Q, const FFPACK_LUDIVINE_TAG LuTag)
{
	typename Field::Element mone,one,zero;
	F.init(one,1.0);
	F.init(zero,0.0);
	F.neg (mone, one);
	size_t MN = MIN(M,N);
	typename Field::Element * Acurr = A;
	size_t r = 0;
	
	for (size_t k = 0; k < MN; ++k){
		size_t p = r;
		Acurr = A+k*lda+r;
		while ((p < N) && F.isZero (*(Acurr++)))
			p++;
		if (p < N){
			P[r] = p;
			if (r < k){
				fcopy (F, N-r, (A + r*(lda+1)), 1, (A+k*lda+r),1);
				Acurr = A+r+k*lda;
				for (size_t i=r; i<N; ++i)
					F.assign(*(Acurr++),zero);
			}
			
			fswap (F, M, A+r, lda, A+p, lda);
			Q[r] = k;
			r++;
		}
		if (k+1<M){
			ftrsv (F, FflasUpper, FflasTrans, FflasNonUnit, r, A, lda, A+(k+1)*lda, 1);
			fgemv (F, FflasTrans, r, N-r, mone, A+r, lda, A+(k+1)*lda, 1, one, A+(k+1)*lda+r, 1);
		} else
			return r;
	}
			
	return r;
}


template<class Field>
inline size_t
FFPACK::LUdivine_small( const Field& F, const FFLAS_DIAG Diag, const FFLAS_TRANSPOSE trans,
			const size_t M, const size_t N,		
			typename Field::Element * A, const size_t lda, size_t*P, 
			size_t *Q, const FFPACK_LUDIVINE_TAG LuTag)
{
	return callLUdivine_small <typename Field::Element> ()
		(F, Diag, trans, M, N, A, lda, P, Q, LuTag);
}

template<class Element>
class FFPACK::callLUdivine_small 
{
public:
	template <class Field>
	inline size_t 
	operator()( const Field& F, const FFLAS_DIAG Diag, const FFLAS_TRANSPOSE trans,
		    const size_t M, const size_t N,		
		    typename Field::Element * A, const size_t lda, size_t*P, 
		    size_t *Q, const FFPACK_LUDIVINE_TAG LuTag){

		if ( !(M && N) ) return 0;
		typedef typename Field::Element elt;
		elt mone,zero,one;;
		F.init (one, 1.0);
		F.neg(mone, one);
		F.init (zero, 0.0);
		elt * Aini = A;
		elt * Acurr;
		size_t rowp = 0;
		size_t colp;
		size_t R = 0;
		size_t k = 0;
		//size_t kmax = DotProdBound (F, 0, one) -1; // the max number of delayed operations
		while ((rowp<M) && (k<N)){

			//Find non zero pivot
			colp = k;
			Acurr = Aini;
			while ((F.isZero(*Acurr)) || (F.isZero (F.init (*Acurr, *Acurr))))
				if (++colp == N){
					if (rowp==M-1)
						break;
					colp=k; ++rowp;
					Acurr = Aini += lda;
				}
				else
					++Acurr;
		
			if ((rowp == M-1)&&(colp == N))
				break;
			R++;
			P[k] = colp;
			Q[k] = rowp;

			// Permutation of the pivot column
			fswap (F, M, A+k, lda, A + colp , lda);

			//Normalization
			elt invpiv;
			F.init(*Aini,*Aini);
			F.inv (invpiv,*Aini);

			for (size_t j=1; j<N-k; ++j)
				if (!F.isZero(*(Aini+j)))
					F.init(*(Aini+j), *(Aini+j));
			for (size_t i=lda; i<(M-rowp)*lda; i+=lda)
				if (!F.isZero(*(Aini+i)))
					F.init(*(Aini+i), *(Aini+i));


			if (Diag == FflasUnit) {
				for (size_t j=1; j<N-k; ++j)
					if (!F.isZero(*(Aini+j)))
						F.mulin (*(Aini+j),invpiv);
			} else 
				for (size_t i=lda; i<(M-rowp)*lda; i+=lda)
					if (!F.isZero(*(Aini+i)))
						F.mulin (*(Aini+i),invpiv);
		
			//Elimination
			//Or equivalently, but without delayed ops :
			fger (F, M-rowp-1, N-k-1, mone, Aini+lda, lda, Aini+1, 1, Aini+(lda+1), lda);
		
			Aini += lda+1; ++rowp; ++k;
		}

		// Compression the U matrix
		size_t l;
		if (Diag == FflasNonUnit){
			Aini = A;
			l = N;
		} else {
			Aini = A+1;
			l=N-1;
		}
		for (size_t i=0; i<R; ++i, Aini += lda+1) {
			if (Q[i] > i){
				fcopy (F, l-i, Aini, 1, Aini+(Q[i]-i)*lda, 1);
				for (size_t j=0; j<l-i; ++j)
					F.assign (*(Aini+(Q[i]-i)*lda+j), zero);
			}
		}
		return R;
	}
};

template<>
class FFPACK::callLUdivine_small<double>
{
public:
	template <class Field>
	inline size_t 
	operator()( const Field& F, const FFLAS_DIAG Diag,  const FFLAS_TRANSPOSE ,
		    const size_t M, const size_t N,		
		    typename Field::Element * A, const size_t lda, size_t*P, 
		    size_t *Q, const FFPACK_LUDIVINE_TAG ){

		if ( !(M && N) ) return 0;
		typedef typename Field::Element elt;
		elt mone,zero,one;;
		F.init (one, 1.0);
		F.neg(mone, one);
		F.init (zero, 0.0);
		elt * Aini = A;
		elt * Acurr;
		size_t rowp = 0;
		size_t colp;
		size_t R = 0;
		size_t k = 0;
		size_t delay =0;
		size_t kmax = DotProdBound (F, 0, one, FflasDouble) -1; // the max number of delayed operations
		while ((rowp<M) && (k<N)){

			//Find non zero pivot
			colp = k;
			Acurr = Aini;
			while ((F.isZero(*Acurr)) || (F.isZero (F.init (*Acurr, *Acurr))))
				if (++colp == N){
					if (rowp==M-1)
						break;
					colp=k; ++rowp;
					Acurr = Aini += lda;
				}
				else
					++Acurr;
		
			if ((rowp == M-1)&&(colp == N))
				break;
			R++;
			P[k] = colp;
			Q[k] = rowp;

			// Permutation of the pivot column
			fswap (F, M, A+k, lda, A + colp , lda);

			//Normalization
			elt invpiv;
			F.init(*Aini,*Aini);
			F.inv (invpiv,*Aini);

			for (size_t j=1; j<N-k; ++j)
				if (!F.isZero(*(Aini+j)))
					F.init(*(Aini+j), *(Aini+j));
			for (size_t i=lda; i<(M-rowp)*lda; i+=lda)
				if (!F.isZero(*(Aini+i)))
					F.init(*(Aini+i), *(Aini+i));


			if (Diag == FflasUnit) {
				for (size_t j=1; j<N-k; ++j)
					if (!F.isZero(*(Aini+j)))
						F.mulin (*(Aini+j),invpiv);
			} else 
				for (size_t i=lda; i<(M-rowp)*lda; i+=lda)
					if (!F.isZero(*(Aini+i)))
						F.mulin (*(Aini+i),invpiv);
		
			if (delay++ >= kmax){ // Reduction has to be done
				delay = 0;
				for (size_t i=1; i<M-rowp; ++i)
					for (size_t j=1; j<N-k; ++j)
						F.init(	*(Aini+i*lda+j),*(Aini+i*lda+j));
			}
			//Elimination
			for (size_t i=1; i<M-rowp; ++i)
				for (size_t j=1; j<N-k; ++j)
					*(Aini+i*lda+j) -= *(Aini+i*lda) * *(Aini+j);
			//Or equivalently, but without delayed ops :
			//fger (F, M-rowp-1, N-k-1, mone, Aini+lda, lda, Aini+1, 1, Aini+(lda+1), lda);
		
			Aini += lda+1; ++rowp; ++k;
		}

		// Compression the U matrix
		size_t l;
		if (Diag == FflasNonUnit){
			Aini = A;
			l = N;
		} else {
			Aini = A+1;
			l=N-1;
		}
		for (size_t i=0; i<R; ++i, Aini += lda+1) {
			if (Q[i] > i){
				fcopy (F, l-i, Aini, 1, Aini+(Q[i]-i)*lda, 1);
				for (size_t j=0; j<l-i; ++j)
					F.assign (*(Aini+(Q[i]-i)*lda+j), zero);
			}
		}
		return R;
	}
};

template<>
class FFPACK::callLUdivine_small<float>
{
public:
	template <class Field>
	inline size_t 
	operator()( const Field& F, const FFLAS_DIAG Diag, const FFLAS_TRANSPOSE trans,
		    const size_t M, const size_t N,		
		    typename Field::Element * A, const size_t lda, size_t*P, 
		    size_t *Q, const FFPACK_LUDIVINE_TAG LuTag){

		if ( !(M && N) ) return 0;
		typedef typename Field::Element elt;
		elt mone,zero,one;;
		F.init (one, 1.0);
		F.neg(mone, one);
		F.init (zero, 0.0);
		elt * Aini = A;
		elt * Acurr;
		size_t rowp = 0;
		size_t colp;
		size_t R = 0;
		size_t k = 0;
		size_t delay =0;
		size_t kmax = DotProdBound (F, 0, one, FflasFloat) -1; // the max number of delayed operations
		while ((rowp<M) && (k<N)){

			//Find non zero pivot
			colp = k;
			Acurr = Aini;
			while ((F.isZero(*Acurr)) || (F.isZero (F.init (*Acurr, *Acurr))))
				if (++colp == N){
					if (rowp==M-1)
						break;
					colp=k; ++rowp;
					Acurr = Aini += lda;
				}
				else
					++Acurr;
		
			if ((rowp == M-1)&&(colp == N))
				break;
			R++;
			P[k] = colp;
			Q[k] = rowp;

			// Permutation of the pivot column
			fswap (F, M, A+k, lda, A + colp , lda);

			//Normalization
			elt invpiv;
			F.init(*Aini,*Aini);
			F.inv (invpiv,*Aini);

			for (size_t j=1; j<N-k; ++j)
				if (!F.isZero(*(Aini+j)))
					F.init(*(Aini+j), *(Aini+j));
			for (size_t i=lda; i<(M-rowp)*lda; i+=lda)
				if (!F.isZero(*(Aini+i)))
					F.init(*(Aini+i), *(Aini+i));


			if (Diag == FflasUnit) {
				for (size_t j=1; j<N-k; ++j)
					if (!F.isZero(*(Aini+j)))
						F.mulin (*(Aini+j),invpiv);
			} else 
				for (size_t i=lda; i<(M-rowp)*lda; i+=lda)
					if (!F.isZero(*(Aini+i)))
						F.mulin (*(Aini+i),invpiv);
		
			if (delay++ >= kmax){ // Reduction has to be done
				delay = 0;
				for (size_t i=1; i<M-rowp; ++i)
					for (size_t j=1; j<N-k; ++j)
						F.init(	*(Aini+i*lda+j),*(Aini+i*lda+j));
			}
			//Elimination
			for (size_t i=1; i<M-rowp; ++i)
				for (size_t j=1; j<N-k; ++j)
					*(Aini+i*lda+j) -= *(Aini+i*lda) * *(Aini+j);
			//Or equivalently, but without delayed ops :
			//fger (F, M-rowp-1, N-k-1, mone, Aini+lda, lda, Aini+1, 1, Aini+(lda+1), lda);
		
			Aini += lda+1; ++rowp; ++k;
		}

		// Compression the U matrix
		size_t l;
		if (Diag == FflasNonUnit){
			Aini = A;
			l = N;
		} else {
			Aini = A+1;
			l=N-1;
		}
		for (size_t i=0; i<R; ++i, Aini += lda+1) {
			if (Q[i] > i){
				fcopy (F, l-i, Aini, 1, Aini+(Q[i]-i)*lda, 1);
				for (size_t j=0; j<l-i; ++j)
					F.assign (*(Aini+(Q[i]-i)*lda+j), zero);
			}
		}
		return R;
	}
};

template <class Field>
inline size_t 
FFPACK::LUdivine (const Field& F, const FFLAS_DIAG Diag, const FFLAS_TRANSPOSE trans,
		  const size_t M, const size_t N,		
		  typename Field::Element * A, const size_t lda, size_t*P, 
		  size_t *Q, const FFPACK_LUDIVINE_TAG LuTag, const size_t cutoff)
{
	
	if ( !(M && N) ) return 0;
	typedef typename Field::Element elt;
	elt Mone, one, zero;
	F.init(Mone, -1.0);
	F.init(one,1.0);
	F.init(zero,0.0);
	size_t MN = MIN(M,N);

	size_t incRow, incCol, rowDim, colDim;
	if (trans == FflasTrans){
		incRow = 1;
		incCol = lda;
		colDim = M;
		rowDim = N; 
	} else {
		incRow = lda;
		incCol = 1;
		colDim = N;
		rowDim = M;
	}
	
	if ((rowDim < cutoff) && (colDim < 2*cutoff)) // the coeff 2 is experimentally determined!
		return LUdivine_small (F, Diag, trans, M, N, A, lda, P, Q, LuTag);
	else if (MN == 1){
		size_t ip=0;
		//while (ip<N && !F.isUnit(*(A+ip)))ip++;
		while (F.isZero (*(A+ip*incCol)))
			if (++ip == colDim)
				break;
		*Q=0;
		if (ip == colDim){ // current row is zero
			*P=0;
			if (colDim == 1){
				//while (ip<M && !F.isUnit(*(A+ip*lda))){
				while (ip<rowDim && F.isZero(*(A + ip*incRow))){
					Q[ip]=ip;
					ip++;
				}
				if (ip == rowDim) {
					return 0;
				} else {
					size_t oldip = ip;
					if ( Diag == FflasNonUnit ){
						elt invpiv;
						F.inv(invpiv,*(A+ip*incRow));
						while(++ip<rowDim) 
							F.mulin(*(A + ip*incRow), invpiv);
						elt tmp;
						F.assign(tmp, *(A+oldip*incRow));
						F.assign( *(A+oldip*incRow), *A);
						F.assign( *A, tmp);
					}
					*Q=oldip; 
					
					return 1;
				}
			}
			else{ *Q=0; return 0;}
		}
		*P=ip;
		if (ip!=0){
			// swap the pivot
			typename Field::Element tmp=*A;
			*A = *(A + ip*incCol);
			*(A + ip*incCol) = tmp;
		}
		elt invpiv;
		F.inv(invpiv, *A);
		if ( Diag == FflasUnit ){
			// Normalisation of the row
			for (size_t k=1; k<colDim; k++)
				F.mulin(*(A+k*incCol), invpiv);
		}
		else if ( colDim==1 )
			while(++ip<rowDim) 
				F.mulin(*(A + ip*incRow), invpiv);
		return 1;
	} else { // MN>1
		size_t Nup = rowDim >> 1;
		size_t Ndown =  rowDim - Nup;
		// Recursive call on NW
		size_t R, R2;
		if (trans == FflasTrans){
			R = LUdivine (F, Diag, trans, colDim, Nup, A, lda, P, Q, LuTag, cutoff);
			typename Field::Element *Ar = A + Nup*incRow; // SW
			typename Field::Element *Ac = A + R*incCol;     // NE
			typename Field::Element *An = Ar + R*incCol;    // SE
			if (!R){
				if (LuTag == FfpackSingular ) 
					return 0;
			} else {			
				applyP (F, FflasLeft, FflasNoTrans, Ndown, 0, R, Ar, lda, P); 
				// Ar <- L1^-1 Ar
				ftrsm( F, FflasLeft, FflasLower, 
				       FflasNoTrans, Diag, R, Ndown,  
				       one, A, lda, Ar, lda);
				// An <- An - Ac*Ar
				fgemm( F, FflasNoTrans, FflasNoTrans, colDim-R, Ndown, R,
				       Mone, Ac, lda, Ar, lda, one, An, lda);
			}
			// Recursive call on SE
			R2 = LUdivine (F, Diag, trans, colDim-R, Ndown, An, lda, P + R, Q + Nup, LuTag, cutoff);
			for (size_t i = R; i < R + R2; ++i)
				P[i] += R;
			if (R2)
				// An <- An.P2
				applyP (F, FflasLeft, FflasNoTrans, Nup, R, R+R2, A, lda, P); 
			else if (LuTag == FfpackSingular)
				return 0;
			
		} else{
			R = LUdivine (F, Diag, trans, Nup, colDim, A, lda, P, Q, LuTag, cutoff);
			typename Field::Element *Ar = A + Nup*incRow; // SW
			typename Field::Element *Ac = A + R*incCol;     // NE
			typename Field::Element *An = Ar + R*incCol;    // SE
			if (!R){
				if (LuTag == FfpackSingular ) 
					return 0;
			} else {			
				// Ar <- Ar.P
				applyP (F, FflasRight, FflasTrans, Ndown, 0, R, Ar, lda, P); 
				// Ar <- Ar.U1^-1
				ftrsm( F, FflasRight, FflasUpper, 
				       FflasNoTrans, Diag, Ndown, R, 
				       one, A, lda, Ar, lda);
				// An <- An - Ar*Ac
				fgemm( F, FflasNoTrans, FflasNoTrans, Ndown, colDim-R, R,
				       Mone, Ar, lda, Ac, lda, one, An, lda);
				
			}
			// Recursive call on SE
			R2=LUdivine (F, Diag, trans, Ndown, N-R, An, lda,P+R, Q+Nup, LuTag, cutoff);
			for (size_t i = R; i < R + R2; ++i)
				P[i] += R;
			if (R2)
				// An <- An.P2
				applyP (F, FflasRight, FflasTrans, Nup, R, R+R2, A, lda, P); 
			else if (LuTag == FfpackSingular)
				return 0;
			
		}
		// Non zero row permutations
		for (size_t i = Nup; i < Nup + R2; i++)
			Q[i] += Nup;
		if (R < Nup){
			// Permutation of the 0 rows
			if (Diag == FflasNonUnit){
				for ( size_t i = Nup, j = R ; i < Nup + R2; ++i, ++j){
					fcopy( F, colDim - j, A + j * (lda + 1), incCol, A + i*incRow + j*incCol, incCol);
 					for (typename Field::Element *Ai = A + i*incRow + j*incCol;
 					     Ai != A + i*incRow + colDim*incCol; Ai+=incCol)
 						F.assign (*Ai, zero);
					size_t t = Q[j];
					Q[j]=Q[i];
					Q[i] = t;
				}
			} else {
				for ( size_t i = Nup, j = R+1 ; i < Nup + R2; ++i, ++j){
					fcopy( F, colDim - j,
					       A + (j-1)*incRow + j*incCol, incCol,
					       A + i*incRow + j*incCol, incCol);
					for (typename Field::Element *Ai = A + i*incRow + j*incCol;
					     Ai != A + i*incRow + colDim*incCol; Ai+=incCol)
						F.assign (*Ai, zero);
					size_t t = Q[j-1];
					Q[j-1]=Q[i];
					Q[i] = t;
				}
			}
		}
		return R + R2;
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
FFPACK::LUdivine_construct( const Field& F, const FFLAS_DIAG Diag,
				      const size_t M, const size_t N,
				      const typename Field::Element * A, const size_t lda,
				      typename Field::Element * X, const size_t ldx,
				      typename Field::Element * u, size_t* P,
				      bool computeX, const FFPACK_MINPOLY_TAG MinTag = FfpackDense,
				      const size_t kg_mc =0, const size_t kg_mb=0, const size_t kg_j=0)
{

            typename Field::Element Mone, one, zero;
	F.init(Mone, -1.0);
	F.init(one, 1.0);
	F.init(zero,0.0);
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
			typename Field::Element invpiv;
			F.inv(invpiv, *X);
			
			// Normalisation of the row
			for (size_t k=1; k<N; k++)
				F.mulin(*(X+k), invpiv);
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
			if ( computeX ){
				if (MinTag == FfpackDense)
					for (size_t i=0; i< Ndown; ++i, Xi+=ldx){
						fgemv(F, FflasNoTrans, N, N, one, 
						      A, lda, u, 1, zero, Xi,1);
						fcopy(F, N, u,1,Xi, 1);
					}
				else // Keller-Gehrig Fast algorithm's matrix
					for (size_t i=0; i< Ndown; ++i, Xi+=ldx){
						fgemv_kgf( F, N, A, lda, u, 1, Xi, 1, 
  							   kg_mc, kg_mb, kg_j );
						fcopy(F, N, u,1,Xi, 1);
					}
			}
			// Apply the permutation on SW
			applyP( F, FflasRight, FflasTrans, Ndown, 0, R, Xr, ldx, P); 
			// Triangular block inversion of NW and apply to SW
			// Xr <- Xr.U1^-1
			ftrsm( F, FflasRight, FflasUpper, FflasNoTrans, Diag,
			       Ndown, R, one, X, ldx, Xr, ldx);
			
			// Update of SE
			// Xn <- Xn - Xr*Xc
			fgemm( F, FflasNoTrans, FflasNoTrans, Ndown, N-Nup, Nup,
			       Mone, Xr, ldx, Xc, ldx, one, Xn, ldx);

			// Recursive call on SE
			
			size_t R2 = LUdivine_construct(F, Diag, Ndown, N-Nup, A, lda,
							Xn, ldx, u, P + Nup, 
							false, MinTag, kg_mc, kg_mb, kg_j);
			for ( size_t i=R;i<R+R2;++i) P[i] += R;
			
			applyP( F, FflasRight, FflasTrans, Nup, R, R+R2, X, ldx, P); 
			
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

//---------------------------------------------------------------------
// TURBO: rank computation algorithm 
//---------------------------------------------------------------------

template <class Field>
inline size_t 
FFPACK::TURBO (const Field& F, const size_t M, const size_t N,
	       typename Field::Element* A, const size_t lda, size_t * P, size_t * Q, const size_t cutoff) 
{

	size_t mo2 = (M>>1);
	size_t no2 = (N>>1);

	typename Field::Element * NW = A;
	typename Field::Element * NE = A + no2;
	typename Field::Element * SW = A + mo2*lda;
	typename Field::Element * SE = SW + no2;
	
	size_t ld1, ld2, ld3, ld4;
	ld1 = ld2 = ld3 = ld4 = lda;

	if ( !(M && N) ) return 0;
	typedef typename Field::Element elt;
	elt Mone, one, zero;
	F.init(Mone, -1.0);
	F.init(one,1.0);
	F.init(zero,0.0);

	// Column permutation
	size_t * P1 = new size_t[no2];
	size_t * P2 = new size_t[N-no2];
	// Row Permutation
	size_t * Q1 = new size_t[mo2];
	size_t * Q2 = new size_t[M-mo2];
	for (size_t i=0; i<mo2; ++i)
		Q1[i] = 0;
	for (size_t i=0; i<M-mo2; ++i)
		Q2[i] = 0;
	size_t q1,q2,q3,q3b,q4;
	q1=q2=q3=q3b=q4=0;
		
	// Step 1: NW = L1.Q1.U1.P1
	size_t mloc = mo2;
	size_t nloc = no2;
// 	Timer tim;
// 	tim.clear();
// 	tim.start();
	q1 = LUdivine( F, FflasNonUnit, FflasNoTrans, mloc, no2, NW, ld1, P1, Q1, FfpackLQUP, cutoff);
	
// 	tim.stop();
// 	cerr<<"LQUP1:"<<tim.realtime()<<std::endl;
// 	tim.start();
#ifdef LB_DEBUG
	std::cerr<<"NW= L1.Q1.U1.P1"<<std::endl;
	write_field(F,std::cerr,NW,M,N,lda);
#endif	
	// B1 = L^-1.NE
#ifdef LB_DEBUG
	std::cerr<<"avant B1 = L^-1.NE"<<std::endl;
	write_field(F,std::cerr,NE,mloc,N-no2,ld2);
#endif	
	solveLB( F, FflasLeft, mo2, N-no2, q1, NW, ld1, Q1, NE, ld2);
#ifdef LB_DEBUG
	std::cerr<<"B1 = L^-1.NE"<<std::endl;
	write_field(F,std::cerr,NE,mloc,N-no2,ld2);
#endif	

	// NE = Q^-1.NE
	
	applyP( F, FflasLeft, FflasNoTrans, N-no2, 0, mo2, NE, ld2, Q1);		
#ifdef LB_DEBUG
	std::cerr<<"NE=Q^-1.NE"<<std::endl;
	write_field(F,std::cerr,NE,mloc,N-no2,ld2);
#endif	

	// SW = SW.P1
	applyP( F, FflasRight, FflasTrans, M-mo2, 0, q1, SW, ld3, P1 );
#ifdef LB_DEBUG
	std::cerr<<"SW = SW.P1"<<std::endl;
	write_field(F,std::cerr,SW,M-mo2,no2,ld3);
#endif	

// 	tim.stop();
// 	std::cerr<<"L^-1:"<<tim.realtime()<<std::endl;
// 	tim.start();
	
	// N1 = SW_{1,q1} . U1^-1
	ftrsm( F, FflasRight, FflasUpper, FflasNoTrans, FflasNonUnit, M-mo2, q1, one, NW, ld1 , SW, ld3 );
#ifdef LB_DEBUG
	std::cerr<<" N1 = SW_{1,q1} . U1^-1"<<std::endl;
	write_field(F,std::cerr,SW,M-mo2,no2,ld3);
#endif	

// 	tim.stop();
// 	std::cerr<<"trsm:"<<tim.realtime()<<std::endl;
// 	tim.start();
	
	// I1 = SW_{q1+1,n} - N1.G1  
	fgemm(F, FflasNoTrans, FflasNoTrans, M-mo2,  no2-q1, q1, Mone, SW, ld3, NW+q1, ld1, one, SW+q1, ld3);
#ifdef LB_DEBUG
	std::cerr<<" I1 = SW_{q1+1,n} - N1.G1"<<std::endl;
	write_field(F,std::cerr,SW,M-mo2,no2,ld3);
#endif	
// 	tim.stop();
// 	std::cerr<<"fgemm1:"<<tim.realtime()<<std::endl;
// 	tim.start();
			
	// E1 = SE - N1.B1_{1,q1}
	fgemm( F, FflasNoTrans, FflasNoTrans, M-mo2, N-no2, q1, Mone, SW, ld3, NE, ld2, one, SE, ld4);
#ifdef LB_DEBUG
	std::cerr<<"  E1 = SE - N1.B1_{1,q1}"<<std::endl;
	write_field(F,std::cerr,SE,M-mo2,N-no2,ld4);
#endif	
// 	tim.stop();
// 	std::cerr<<"fgemm2:"<<tim.realtime()<<std::endl;
// 	tim.start();


	//Step 2: E1 = L2.Q2.U2.P2
	mloc = M-mo2;
	nloc = N-no2;
	q2 = LUdivine( F, FflasNonUnit, FflasNoTrans, mloc, nloc, SE, ld4, P2, Q2, FfpackLQUP, cutoff);
#ifdef LB_DEBUG
	std::cerr<<"  E1 = L2.Q2.U2.P2"<<std::endl;
	write_field(F,std::cerr,SE,M-mo2,N-no2,ld4);	
#endif	
// 	tim.stop();
// 	std::cerr<<"LQUP2:"<<tim.realtime()<<std::endl;
// 	tim.start();

	// [I2;F2] = L2^-1.I1
	solveLB( F, FflasLeft, mloc, no2-q1, q2, SE, ld4, Q2, SW+q1, ld3);
#ifdef LB_DEBUG
	std::cerr<<"  [I2;F2] = L2^-1.I1"<<std::endl;
	write_field(F,std::cerr,SW,M-mo2,no2,ld3);	
#endif	
	// I1 = Q2^-1.I1
	applyP( F, FflasLeft, FflasNoTrans, no2-q1, 0, mloc, SW+q1, ld3, Q2 );
#ifdef LB_DEBUG
	std::cerr<<"I1 = Q2^-1.I1"<<std::endl;
	write_field(F,std::cerr,SW,mloc,no2,ld3);
 #endif	

	// B1 = B1.P2
	applyP( F, FflasRight, FflasTrans, mo2, 0, q2, NE, ld2, P2 );
#ifdef LB_DEBUG
	std::cerr<<"B1 = B1.P2"<<std::endl;
	write_field(F,std::cerr,NE,mo2,N-no2,ld2);
#endif	
	// Updating P
	//	for (size_t i=no2;i<N;++i)
	//	P[i] += no2;
// 	tim.stop();
// 	std::cerr<<"L2^-1:"<<tim.realtime()<<std::endl;
// 	tim.start();

	//alternative: de 0 a q2 avant
	// N2 = B1_{q1+1,mo2} . V2^-1
	ftrsm(F, FflasRight, FflasUpper,FflasNoTrans,FflasNonUnit, mo2-q1, q2, one, SE, ld4, NE+q1*ld2,ld2);
// 	tim.stop();
// 	std::cerr<<"trsm2:"<<tim.realtime()<<std::endl;
// 	tim.start();
		
	// H2 = B1_{q1+1,mo2;q2,N-no2} - N2.E2  
	fgemm(F, FflasNoTrans, FflasNoTrans, mo2-q1, N-no2-q2, q2, Mone, NE+q1*ld2, ld2, SE+q2, ld4, one, NE+q1*ld2+q2, ld2);

// 	tim.stop();
// 	std::cerr<<"fgemm12:"<<tim.realtime()<<std::endl;
// 	tim.start();
	// O2 = NW_{q1+1,mo2;q1+1,N-no2} = - N2.I2  
	//write_field (F,cerr<<"avant O2"<<endl, A, M, N, lda);

	fgemm(F, FflasNoTrans, FflasNoTrans, mo2-q1, no2-q1, q2, Mone, NE+q1*ld2, ld2, SW+q1, ld3, zero,
	      NW+q1*(ld1+1), ld1);
	//	write_field (F,cerr<<"apres O2"<<endl, A, M, N, lda);
// 	tim.stop();
// 	std::cerr<<"fgemm22:"<<tim.realtime()<<std::endl;
// 	tim.start();

	
	//Step 3: F2 = L3.Q3.U3.P3
	mloc = M-mo2-q2;
	nloc = no2-q1;
	q3 = LUdivine( F, FflasNonUnit, FflasNoTrans, mloc, nloc, SW+q2*ld3+q1, ld3, P1+q1, Q2+q2, FfpackLQUP, cutoff);
	
	// Updating P1,Q2
	for (size_t i=q1;i<no2;++i)
		P1[i] += q1;
	for (size_t i=q2;i<q2+q3;++i)
		Q2[i] += q2;
		
	//Step 3bis: H2 = L3b.Q3b.U3b.P3b
	mloc = mo2-q1;
	nloc = N-no2-q2;

	q3b = LUdivine( F, FflasNonUnit, FflasNoTrans, mloc, nloc, NE+q1*ld2+q2, ld2, P2+q2, Q1+q1, FfpackLQUP, cutoff);
	
	// Updating P2, Q1
	for (size_t i = q2; i < q2+q3b; ++i)
		P2[i] += q2;

// 	tim.stop();
// 	std::cerr<<"LQUP3et3bis:"<<tim.realtime()<<std::endl;
// 	tim.start();
		
	if (( q3 < no2-q1) && (q3b<mo2-q1)){
			
		// [O3;_] = L3b^-1.O2
		if (q3b>0){
// 			if ( mo2-q1 < N-no2-q2+q1) 
// 				// L is expanded to a Lower triangular matrix
// 				solveLB( F, FflasLeft,mloc, no2-q1, q3b, NE+q1*ld2+q2 , ld2, rP3b, NW+q1*(ld1+1), ld1);
//			else{
			//std::cerr<<"USING SOLVELB2"<<std::endl;
			//no modification of L
			solveLB2( F, FflasLeft,mloc, no2-q1, q3b, NE+q1*ld2+q2 , ld2, Q1+q1, NW+q1*(ld1+1), ld1);
//			}
#ifdef LB_DEBUG
			std::cerr<<"O2 avant="<<std::endl;
			write_field(F,std::cerr,NW+q1*(ld1+1),mloc,no2-q1,ld1);
#endif	
	
			// O2 = Q3b^-1.O2
			applyP( F, FflasLeft, FflasNoTrans, no2-q1, 0, mloc, NW+q1*(ld1+1), ld1, Q1+q1 );
#ifdef LB_DEBUG
			std::cerr<<"O2 apres="<<std::endl;
			write_field(F,std::cerr,NW+q1*(ld1+1),mloc,no2-q1,ld1);
#endif	
	
			//updating Q
#if 0
 			size_t tmp;
 			for (size_t j=0;j<mo2-q1;++j)
 				if (rP3b[j]!=j){
 					//	std::cerr<<"(rP3b["<<j<<"]="<<rP3b[j]<<std::endl;
 					tmp = Q[j+q1];
 					Q[j+q1] = Q[rP3b[j]+q1];
 					Q[rP3b[j]+q1] = tmp;
 				}
#endif
				
			// X2 = X2.P3
			// Si plusieurs niveaux rec, remplacer X2 par [NW;I2]
			applyP( F, FflasRight, FflasTrans, mo2-q1-q3b, q1, q1+q3,
				NW/*+(q1+q3b)*ld1*/, ld1, P1);
			applyP( F, FflasRight, FflasTrans, q2, q1, q1+q3,
				SW/*+(q1+q3b)*ld1*/, ld3, P1);
			
	
			// A faire si plusieurs niveaux recursifs
			// B2 = B2.P3b
			applyP (F, FflasRight, FflasTrans, q1, q2, q2+q3b,
				NW, ld2, P2);
			//flaswp(F,q1,NE,lda,no2+q2,no2+q2+q3b,P,1); 
			// E2 = E2.P3b
			applyP (F, FflasRight, FflasTrans, q2, q2, q2+q3b,
				SE, ld4, P2);
			//flaswp(F,q2,SE+q2,lda,no2+q2,no2+q2+q3b,P,1); 
		}
					
		// N3 = X2 . D3^-1
		ftrsm( F, FflasRight, FflasUpper, FflasNoTrans, FflasNonUnit, mo2-q1-q3b, q3, one, SW+q2*ld3+q1, ld3 ,NW+(q1+q3b)*ld1+q1,ld1);

		// T2 = T2 - N3.F3
		fgemm( F, FflasNoTrans, FflasNoTrans, mo2-q1-q3b, no2-q1-q3,q3, Mone, NW+(q1+q3b)*ld1+q1, ld1, SW+q2*ld3+q3+q1, ld3, one, NW+(q1+q3b)*ld1+q1+q3, ld1 );


		//Step 4: T2 = L4.Q4.U4.P4
		mloc = mo2-q1-q3b;
		nloc = no2-q1-q3;
		
		// size_t * rP4 = new size_t[mloc];
// 		for (size_t j=0;j<mo2-q1;++j)
// 			rP4[j]=0;
		q4 = LUdivine( F, FflasNonUnit, FflasNoTrans, mloc, nloc, NW+(q1+q3b)*ld1+q1+q3, ld1, P1+q1+q3, Q1+q1+q3b, FfpackLQUP, cutoff);

		// Updating P
		for (size_t i=q1+q3;i<q1+q3+q4;++i)
			P1[i] += q3;

//		size_t tmp;
// 			if (rP4[j]!=j){
// 				//	std::cerr<<"(rP3b["<<j<<"]="<<rP3b[j]<<std::endl;
// 				tmp = Q[j+q1+q3b];
// 				Q[j+q1+q3b] = Q[rP3b[j]+q1+q3b];
// 				Q[rP3b[j]+q1+q3b] = tmp;
// 			}
		
		// A faire si plusieurs niveaux recursifs
		// [G1;O3] = [G1;O3].P4
		applyP (F, FflasRight, FflasTrans, q1+q3b, q1+q3, q1+q3+q4,
			NW, ld1, P1);
		//flaswp(F,q1+q3b,NE,lda,no2+q2,no2+q2+q3b,P,1); 
		// [I2;F3] = [I2;F3].P4
		applyP (F, FflasRight, FflasTrans, q2+q3, q1+q3, q1+q3+q4,
			SW, ld3, P1);
		//flaswp(F,q2,SE+q2,lda,no2+q2,no2+q2+q3b,P,1); 
	}
	//!!!!!! Attention a appliquer Q4, Q2, Q3, Q3b a gauche !!!!!!!
	
	//updating Q1
	for (size_t i = q1; i < q1+q3b; ++i)
		Q1[i] += q1;
	for (size_t i=q1+q3b;i<q1+q3b+q4;++i)
		Q1[i] += q1 + q3b;

	for (size_t i=0; i<q1; ++i)
		P[i] = P1[i];
	for (size_t i=q1; i<q1+q2; ++i)
		P[i] = P2[i-q1] + no2;
	for (size_t i=q1+q2; i<q1+q2+q3; ++i)
		P[i] = P1[i-q2];
	for (size_t i=q1+q2+q3; i<q1+q2+q3+q3b; ++i)
		P[i] = P2[i-q1-q3]+no2;
	for (size_t i=q1+q2+q3+q3b; i<q1+q2+q3+q3b+q4; ++i)
		P[i] = P1[i-q2-q3b];
	delete[] P1;
	delete[] P2;

	for (size_t i=0; i<q1; ++i)
		Q[i] = Q1[i];
	for (size_t i=q1; i<q1+q2; ++i)
		Q[i] = Q2[i-q1] + mo2;
	for (size_t i=q1+q2; i<q1+q2+q3; ++i)
		Q[i] = Q2[i-q1] + mo2;
	for (size_t i=q1+q2+q3; i<q1+q2+q3+q3b; ++i)
		Q[i] = Q1[i-q2-q3];
	for (size_t i=q1+q2+q3+q3b; i<q1+q2+q3+q3b+q4; ++i)
		P[i] = Q1[i-q2-q3];
	delete[] Q1;
	delete[] Q2;
	
	
	//write_field (F, cerr<<"avant reordonnancement"<<endl, A, M,N, lda)<<endl;
	typename Field::Element * R = new typename Field::Element[M*N];
	size_t ldr = N;
	// Copying first q1 cols
	for (size_t i=0; i<q1; ++i)
		fcopy (F, q1, R+i*ldr, 1, NW+i*ld1,1);
	for (size_t i=q1; i<q1+q2+q3; ++i)
		fcopy (F, q1, R+i*ldr, 1, SW+(i-q1)*ld3,1);
	for (size_t i=q1+q2+q3; i<q2+q3+mo2; ++i)
		fcopy (F, q1, R+i*ldr, 1, NW+(i-q2-q3)*ld1,1);
	for (size_t i=q2+q3+mo2; i<M; ++i)
		fcopy (F, q1, R+i*ldr, 1, SW+(i-mo2)*ld3,1);
	// Copying q1..q2 cols 
	for (size_t i=0; i<q1; ++i)
		fcopy (F, q2, R+q1+i*ldr, 1, NE+i*ld2,1);
	for (size_t i=q1; i<q1+q2+q3; ++i)
		fcopy (F, q2, R+q1+i*ldr, 1, SE+(i-q1)*ld4,1);
	for (size_t i=q1+q2+q3; i<q2+q3+mo2; ++i)
		fcopy (F, q2, R+q1+i*ldr, 1, NE+(i-q2-q3)*ld2,1);
	for (size_t i=q2+q3+mo2; i<M; ++i)
		fcopy (F, q2, R+q1+i*ldr, 1, SE+(i-mo2)*ld4,1);
	// Copying q2..q3 cols 
	for (size_t i=0; i<q1; ++i)
		fcopy (F, q3, R+q1+q2+i*ldr, 1, NW+q1+i*ld1,1);
	for (size_t i=q1; i<q1+q2+q3; ++i)
		fcopy (F, q3, R+q1+q2+i*ldr, 1, SW+q1+(i-q1)*ld3,1);
	for (size_t i=q1+q2+q3; i<q2+q3+mo2; ++i)
		fcopy (F, q3, R+q1+q2+i*ldr, 1, NW+q1+(i-q2-q3)*ld1,1);
	for (size_t i=q2+q3+mo2; i<M; ++i)
		fcopy (F, q3, R+q1+q2+i*ldr, 1, SW+q1+(i-mo2)*ld3,1);
	// Copying q3..q3b cols 
	for (size_t i=0; i<q1; ++i)
		fcopy (F, q3b, R+q1+q2+q3+i*ldr, 1, NE+q2+i*ld2,1);
	for (size_t i=q1; i<q1+q2+q3; ++i)
		fcopy (F, q3b, R+q1+q2+q3+i*ldr, 1, SE+q2+(i-q1)*ld4,1);
	for (size_t i=q1+q2+q3; i<q2+q3+mo2; ++i)
		fcopy (F, q3b, R+q1+q2+q3+i*ldr, 1, NE+q2+(i-q2-q3)*ld2,1);
	for (size_t i=q2+q3+mo2; i<M; ++i)
		fcopy (F, q3b, R+q1+q2+q3+i*ldr, 1, SE+q2+(i-mo2)*ld4,1);
	// Copying q3b..q4 cols 
	for (size_t i=0; i<q1; ++i)
		fcopy (F, q4, R+q1+q2+q3+q3b+i*ldr, 1, NW+q1+q3+i*ld1,1);
	for (size_t i=q1; i<q1+q2+q3; ++i)
		fcopy (F, q4, R+q1+q2+q3+q3b+i*ldr, 1, SW+q1+q3+(i-q1)*ld3,1);
	for (size_t i=q1+q2+q3; i<q2+q3+mo2; ++i)
		fcopy (F, q4, R+q1+q2+q3+q3b+i*ldr, 1, NW+q1+q3+(i-q2-q3)*ld1,1);
	for (size_t i=q2+q3+mo2; i<M; ++i)
		fcopy (F, q4, R+q1+q2+q3+q3b+i*ldr, 1, SW+q1+q3+(i-mo2)*ld3,1);
	// Copying the last cols 
	for (size_t i=0; i<q1; ++i)
		fcopy (F, no2-q1-q3-q4, R+q1+q2+q3+q3b+q4+i*ldr, 1, NW+q1+q3+q4+i*ld1,1);
	for (size_t i=q1; i<q1+q2+q3; ++i)
		fcopy (F, no2-q1-q3-q4, R+q1+q2+q3+q3b+q4+i*ldr, 1, SW+q1+q3+q4+(i-q1)*ld3,1);
	for (size_t i=q1+q2+q3; i<q2+q3+mo2; ++i)
		fcopy (F, no2-q1-q3-q4, R+q1+q2+q3+q3b+q4+i*ldr, 1, NW+q1+q3+q4+(i-q2-q3)*ld1,1);
	for (size_t i=q2+q3+mo2; i<M; ++i)
		fcopy (F, no2-q1-q3-q4, R+q1+q2+q3+q3b+q4+i*ldr, 1, SW+q1+q3+q4+(i-mo2)*ld3,1);
	// Copying the last cols 
	for (size_t i=0; i<q1; ++i)
		fcopy (F, N-no2-q2-q3b, R+no2+q2+q3b+i*ldr, 1, NE+q2+q3b+i*ld2,1);
	for (size_t i=q1; i<q1+q2+q3; ++i)
		fcopy (F, N-no2-q2-q3b, R+no2+q2+q3b+i*ldr, 1, SE+q2+q3b+(i-q1)*ld4,1);
	for (size_t i=q1+q2+q3; i<q2+q3+mo2; ++i)
		fcopy (F, N-no2-q2-q3b, R+no2+q2+q3b+i*ldr, 1, NE+q2+q3b+(i-q2-q3)*ld2,1);
	for (size_t i=q2+q3+mo2; i<M; ++i)
		fcopy (F, N-no2-q2-q3b, R+no2+q2+q3b+i*ldr, 1, SE+q2+q3b+(i-mo2)*ld4,1);

	// A=R : to be improved (avoid allocation of R). To be changed if rec data structure are used
	for (size_t i=0; i<M; ++i)
		fcopy (F, N, A+i*lda, 1, R+i*ldr,1);
	
	delete[] R;
	//delete[] Q;
	// Necessaire:
	// 1 traiter les flaswp manquants
	// Facultatif:
	// 2 permutations de lignes doivent etre coherentes
	// 3 effectuer les dernieres permutations lignes et colonnes
	//std::cerr<<q1<<" "<<q2<<" "<<q3<<" "<<q3b<<" "<<q4<<std::endl;
	return q1+q2+q3+q3b+q4;
}

#undef LB_DEBUG
#endif //__LINBOX_ffpack_ludivine_INL

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
