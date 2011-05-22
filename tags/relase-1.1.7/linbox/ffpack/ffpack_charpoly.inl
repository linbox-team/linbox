
/* ffpack/ffpack_charpoly.inl
 * Copyright (C) 2005 Clement Pernet
 *
 * Written by Clement Pernet <Clement.Pernet@imag.fr>
 *
 * See COPYING for license information.
 */

template <class Field, class Polynomial>
std::list<Polynomial>&
FFPACK::CharPoly (const Field& F, std::list<Polynomial>& charp, const size_t N,
		  typename Field::Element * A, const size_t lda,
		  const FFPACK_CHARPOLY_TAG CharpTag){
	switch (CharpTag) {
	case FfpackLUK:{
		typename Field::Element * X = new typename Field::Element[N*(N+1)];
		LUKrylov (F, charp, N, A, lda, X, N);
		delete[] X;
		return charp;
	}
	case FfpackKG:{
		return KellerGehrig (F, charp, N, A, lda);
		//break;
	}
	case FfpackDanilevski:{
		return Danilevski (F, charp, N, A, lda);
		//break;
	}
	case FfpackKGFast:{
		size_t mc, mb, j;
		if (KGFast (F, charp, N, A, lda, &mc, &mb, &j)){
			std::cerr<<"NON GENERIC MATRIX PROVIDED TO KELLER-GEHRIG-FAST"<<std::endl;
		}
		return charp;
		//break;
	}
 	case FfpackKGFastG:{
 		return KGFast_generalized (F, charp, N, A, lda);
	}
	case FfpackHybrid:{
		typename Field::Element * X = new typename Field::Element[N*(N+1)];
		LUKrylov_KGFast (F, charp, N, A, lda, X, N);
		delete[] X;
		return charp;
	}
	case FfpackArithProg:{
		size_t attempts=0;
		bool cont = false;
		FFLAS_INT_TYPE p;
		F.characteristic(p);
		// Heuristic condition (the pessimistic theoretical one being p<2n^2.
		if ((unsigned long) (p) < N)
			return CharPoly (F, charp, N, A, lda, FfpackLUK);

		do{
			try {
				CharpolyArithProg (F, charp, N, A, lda, __FFPACK_CHARPOLY_THRESHOLD);
			}
			catch (CharpolyFailed){
				if (attempts++ < 2)
					cont = true;
				else
					return CharPoly(F, charp, N, A, lda, FfpackLUK);
				
			}
		} while (cont);
		return charp;
	}
	default:{
		typename Field::Element * X = new typename Field::Element[N*(N+1)];
		LUKrylov (F, charp, N, A, lda, X, N);
		delete[] X;
		return charp;
	}
	}
}

template <class Field, class Polynomial>
std::list<Polynomial>&
FFPACK::LUKrylov (const Field& F, std::list<Polynomial>& charp, const size_t N,
		 typename Field::Element * A, const size_t lda,
		 typename Field::Element * X, const size_t ldx){
	
	typedef typename Field::Element elt;
	elt* Ai, *Xi, *X2=X;
	static elt Mone, one, zero;
	F.init(zero,0.0);
	F.init(one, 1.0);
	F.neg(Mone,one);
	int Ncurr=N;
	charp.clear();
	int nbfac = 0;
	while (Ncurr > 0){
		size_t P[Ncurr];
		Polynomial minP;//=new Polynomial();
		MinPoly (F, minP, Ncurr, A, lda, X2, ldx, P);
		int k = minP.size()-1; // degre of minpoly
		if ((k==1) && F.isZero ((minP)[0])){ // minpoly is X
			Ai = A;
			int j = Ncurr*Ncurr;
			while (j-- && F.isZero(*(Ai++))) ;
			if (!j){ // A is 0, CharPoly=X^n
				minP.resize(Ncurr+1);
				(minP)[1] = zero;
				(minP)[Ncurr] = one;
				k=Ncurr;
			}
		}
		nbfac++;
		charp.push_front (minP);
		if (k==Ncurr){
			return charp;
		}
		size_t Nrest = Ncurr-k;
		elt * X21 = X2 + k*ldx;
		elt * X22 = X21 + k;
		// Compute the n-k last rows of A' = PA^tP^t in X2_
		// A = A . P^t
		applyP (F, FflasRight, FflasTrans, Ncurr, 0, k, A, lda, P);
		// Copy X2_ = (A'_2)^t
		for (Xi = X21, Ai = A+k; Xi != X21 + Nrest*ldx; Ai++, Xi+=ldx-Ncurr)
			for (size_t jj=0; jj<Ncurr*lda; jj+=lda)
				*(Xi++) = *(Ai+jj);
		// A = A . P : Undo the permutation on A
		applyP (F, FflasRight, FflasNoTrans, Ncurr, 0, k, A, lda, P);
		// X2_ = X2_ . P^t (=  (P A^t P^t)2_) 
		applyP (F, FflasRight, FflasTrans, Nrest, 0, k, X21, ldx, P);
		// X21 = X21 . S1^-1
		ftrsm(F, FflasRight, FflasUpper, FflasNoTrans, FflasUnit, Nrest, k,
		      one, X2, ldx, X21, ldx); 
		// Creation of the matrix A2 for recurise call 
		for (Xi = X22, Ai = A;
		     Xi != X22 + Nrest*ldx;
		     Xi += (ldx-Nrest), Ai += (lda-Nrest))
			for (size_t jj=0; jj<Nrest; ++jj)
				*(Ai++) = *(Xi++);
		fgemm (F, FflasNoTrans, FflasNoTrans, Nrest, Nrest, k, Mone,
		       X21, ldx, X2+k, ldx, one, A, lda);
		X2 = X22;
		Ncurr = Nrest;
	}
	return charp;
}

template <class Field, class Polynomial>
std::list<Polynomial>&
FFPACK::LUKrylov_KGFast (const Field& F, std::list<Polynomial>& charp, const size_t N,
			 typename Field::Element * A, const size_t lda,
			 typename Field::Element * X, const size_t ldx){
	
	typedef typename Field::Element elt;
	
	static elt Mone, one, zero;
	F.init(zero,0.0);
	F.init(one, 1.0);
	F.neg(Mone,one);
	size_t kg_mc, kg_mb, kg_j;
	
	if (!KGFast (F, charp, N, A, lda, &kg_mc, &kg_mb, &kg_j))
		return charp;
	else{// Matrix A is not generic
		Polynomial *minP = new Polynomial();
		const elt* Ai;
		elt* A2i, *Xi;
		size_t *P = new size_t[N];

		MinPoly (F, *minP, N, A, lda, X, ldx, P, FfpackKGF, kg_mc, kg_mb, kg_j);
		size_t k = minP->size()-1; // degre of minpoly
		if ((k==1) && F.isZero ((*minP)[0])){ // minpoly is X
			Ai = A;
			int j = N*N;
			while (j-- && F.isZero(*(Ai++))) ;
			if (!j){ // A is 0, CharPoly=X^n
				minP->resize(N+1);
				(*minP)[1] = zero;
				(*minP)[N] = one;
				k=N;
			}
		}
		
		if (k==N){
			charp.clear();
			charp.push_front(*minP); // CharPoly = MinPoly
			delete[] P;
			return charp;
		}

		size_t Nrest = N-k;
		elt * X21 = X + k*ldx;
		elt * X22 = X21 + k;
		
		// Creates the matrix A
		//size_t lambda = MAX(0,N - kg_mc*(kg_j+1) - kg_mb);  // uint >= 0 !!!
		size_t lambda =   kg_mc*(kg_j+1) + kg_mb;
		if (lambda > N) 
			lambda = 0 ;	
		else
			lambda = N - lambda ;

		size_t imax = kg_mc+kg_mb;
		// First Id
		for (size_t j = 0; j < lambda; ++j){
			for (size_t i=0; i<imax; ++i)
				F.assign (*(A+j+i*lda), zero);
			F.assign (*(A+j+imax*lda), one);
			for (size_t i=imax+1; i<N; ++i)
				F.assign (*(A+j+i*lda), zero);
			++imax;
		}
		// Column block B
		for (typename Field::Element* Ai=A; Ai<A+N*lda; Ai+=lda)
			fcopy (F, kg_mb, Ai+lambda, 1, Ai+N-kg_mc-kg_mb, 1);

		// Second Id block
		imax = N- kg_j*kg_mc;
		for (size_t j = 0; j< kg_j*kg_mc; ++j){
			for (size_t i = 0; i<imax; ++i)
				F.assign (*(A+lambda+kg_mb+j+i*lda), zero);
			F.assign (*(A+lambda+kg_mb+j+imax*lda), one);
			for (size_t i = imax+1; i<N; ++i)
				F.assign (*(A+lambda+kg_mb+j+i*lda), zero);
			++imax;
		}
		
		// Compute the n-k last rows of A' = PA^tP^t in X2_
		
		// A = P . A 
		applyP (F, FflasLeft, FflasNoTrans, N, 0, k, 
			const_cast<typename Field::Element* &>(A), lda, P);
		
		// Copy X2_ = (A'2_)
		for (Xi = X21, Ai = A+k*lda; Xi != X21 + Nrest*ldx; Ai+=lda-N, Xi+=ldx-N){
			for (size_t jj=0; jj<N; ++jj){
				*(Xi++) = *(Ai++);
			}
		}

		// A = P^t . A : Undo the permutation on A
		applyP (F, FflasLeft, FflasTrans, N, 0, k, 
			const_cast<typename Field::Element* &>(A), lda, P);
	
		// X2_ = X2_ . P^t (=  (P A P^t)2_) 
		applyP (F, FflasRight, FflasTrans, Nrest, 0, k, X21, ldx, P);

		// X21 = X21 . S1^-1
		ftrsm(F, FflasRight, FflasUpper, FflasNoTrans, FflasUnit, Nrest, k,
		      one, X, ldx, X21, ldx);  
	
		// Creation of the matrix A2 for recurise call 
		elt * A2 = new elt[Nrest*Nrest];
	
		for (Xi = X22, A2i = A2;
		     Xi != X22 + Nrest*ldx;
		     Xi += (ldx-Nrest)){
			for (size_t jj=0; jj<Nrest; ++jj){
				*(A2i++) = *(Xi++);
			}
		}
		fgemm (F, FflasNoTrans, FflasNoTrans, Nrest, Nrest, k, Mone,
		       X21, ldx, X+k, ldx, one, A2, Nrest);
	
		// Recursive call on X22
		LUKrylov_KGFast (F, charp, Nrest, A2, Nrest, X22, ldx);
		charp.push_front (*minP);
		delete[] P;
		delete[] A2;
		return charp;
	}
}
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
