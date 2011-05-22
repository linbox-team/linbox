
/* linbox/ffpack/ffpack_frobenius.inl
 * Copyright (C) 2007 Clement Pernet
 *
 * Written by Clement Pernet <cpernet@uwaterloo.ca>
 *
 * See COPYING for license information.
 */

#ifndef MIN
#define MIN(a,b) (a<b)?a:b
#endif

//---------------------------------------------------------------------
// CharpolyArithProg: Las Vegas algorithm to compute the Charpoly 
// over a large field (Z/pZ, s.t.  p > 2n^2)
//---------------------------------------------------------------------
template <class Field, class Polynomial>
std::list<Polynomial>&
FFPACK::CharpolyArithProg (const Field& F, std::list<Polynomial>& frobeniusForm, 
			   const size_t N, typename Field::Element * A, const size_t lda,
			   const size_t c){
	
	static typename Field::Element one, zero, mone;
	F.init(one, 1UL);
	F.neg(mone, one);
	F.init(zero, 0UL);

	size_t * rp = new size_t[2*N];
	size_t noc = static_cast<size_t>(ceil(double(N)/double(c)));

	// Building the workplace matrix 
	typename Field::Element *K = new typename Field::Element[N*(noc*c)];
	typename Field::Element *K2 = new typename Field::Element[N*(noc*c)];
	size_t ldk = N;

	size_t *dA = new size_t[N]; //PA
	size_t *dK = new size_t[noc*c];
	for (size_t i=0; i<noc; ++i)
		dK[i]=0;

	// Picking a random noc x N block vector U^T
	typename Field::RandIter g (F);
	typename Field::NonZeroRandIter nzg (F,g);
	for (size_t i = 0; i < noc; ++i)
 		for (size_t j = 0; j < N; ++j)
 			g.random( *(K + i*ldk +j) );
	for (size_t i = 0; i < noc; ++i)
		nzg.random (*(K + i*ldk +i));

	// Computing the bloc Krylov matrix [U AU .. A^(c-1) U]^T
	for (size_t i = 1; i<c; ++i){
		fgemm( F, FflasNoTrans, FflasTrans,  noc, N, N, one,
		       K+(i-1)*noc*ldk, ldk, A, lda, zero, K+i*noc*ldk, ldk);
	}
	// K2 <- K (re-ordering)
	size_t w_idx = 0;
	for (size_t i=0; i<noc; ++i)
		for (size_t j=0; j<c; ++j, w_idx++)
			fcopy(F, N, (K2+(w_idx)*ldk), 1, (K+(i+j*noc)*ldk), 1);

	// Copying K <- K2
	for (size_t i=0; i<noc*c; ++i)
		fcopy (F, N, (K+i*ldk), 1, K2+i*ldk, 1);

	size_t * Pk = new size_t[N];
	size_t * Qk = new size_t[N];
	for (size_t i=0; i<N; ++i)
		Qk[i] = 0;
	for (size_t i=0; i<N; ++i)
		Pk[i] = 0;
	
	size_t R = LUdivine(F, FflasNonUnit, FflasNoTrans, N, N, K, ldk, Pk, Qk, FfpackLQUP);
			
	size_t row_idx = 0;
	size_t i=0;
	size_t dold = c;
	size_t nb_full_blocks = 0;
	size_t Mk = 0;
	// Determining the degree sequence dK
	for (size_t k = 0; k<noc; ++k){
		size_t d = 0;
		while ( (d<c) && (row_idx<R) && (Qk[row_idx] == i)) {i++; row_idx++; d++;}
		if (d > dold){
			// std::cerr << "FAIL in preconditionning phase:"
			//           << " degree sequence is not monotonically not increasing"
			// 	     << std::endl;
			delete[] rp; delete[] K;
			delete[] Pk; delete[] Qk; delete[] dA; delete[] dK;
			throw CharpolyFailed();
		}
		dK[k] = dold = d;
		Mk++;
		if (d == c)
			nb_full_blocks++;
		if (row_idx < N)
			i = Qk[row_idx];
	}

	// Selection of the last iterate of each block

	typename Field::Element * K3 = new typename Field::Element[Mk*N];
	typename Field::Element * K4 = new typename Field::Element[Mk*N];
	size_t bk_idx = 0;
	for (size_t i = 0; i < Mk; ++i){
		fcopy (F, N, (K3+i*ldk), 1, (K2 + (bk_idx + dK[i]-1)*ldk), 1);
		bk_idx += c;
	}
	delete[] K2;

	// K <- K A^T 
	fgemm( F, FflasNoTrans, FflasTrans, Mk, N, N, one,  K3, ldk, A, lda, zero, K4, ldk);

	// K <- K P^T
	applyP (F, FflasRight, FflasTrans, Mk, 0, R, K4, ldk, Pk);

	// K <- K U^-1
	ftrsm (F, FflasRight, FflasUpper, FflasNoTrans, FflasNonUnit, Mk, R, one, K, ldk, K4, ldk);

	// L <-  Q^T L
	applyP(F, FflasLeft, FflasNoTrans, N, 0, R, K, ldk, Qk);

	// K <- K L^-1
	ftrsm (F, FflasRight, FflasLower, FflasNoTrans, FflasUnit, Mk, R, one, K, ldk, K4, ldk);

	//undoing permutation on L
	applyP(F, FflasLeft, FflasTrans, N, 0, R, K, ldk, Qk);

	// Recovery of the completed invariant factors
	size_t Ma = Mk;
	size_t Ncurr = R;
	size_t offset = Ncurr-1;
	for (size_t i=Mk-1; i>=nb_full_blocks+1;  --i){
		if (dK[i] >= 1){ 
			for (size_t j = offset+1; j<R; ++j)
				if (!F.isZero(*(K4 + i*ldk + j))){
					//std::cerr<<"FAIL C != 0 in preconditionning"<<std::endl;
					delete[] K3; delete[] K4; delete[] K;
					delete[] Pk; delete[] Qk; delete[] rp;
					delete[] dA; delete[] dK;
					throw CharpolyFailed();
				}
			Polynomial P (dK [i]+1);
			F.assign(P[dK[i]], one);
			for (size_t j=0; j < dK [i]; ++j)
				F.neg (P [dK [i]-j-1], *(K4 + i*ldk + (offset-j)));
			frobeniusForm.push_front(P);
			offset -= dK [i];
			Ncurr -= dK [i];
			Ma--;
		}
	}
	Mk = Ma;

	if (R<N){
		for (size_t i=0; i<nb_full_blocks + 1; ++i)
			for (size_t j=R; j<N; ++j){
				if (!F.isZero( *(K4+i*ldk+j) )){
					delete[] K3; delete[] K4; delete[] K;
					delete[] Pk; delete[] Qk; delete[] rp;
					delete[] dA; delete[] dK;
					throw CharpolyFailed();
				}
			}
		
		//std::cerr<<"Preconditionning failed; missing rank = "<<N-R
		//	 <<" completing the Krylov matrix"
		//	 <<std::endl;
		size_t Nrest = N-R;
		typename Field::Element * K21 = K + R*ldk;
		typename Field::Element * K22 = K21 + R;
		typename Field::Element * Ki, *Ai;

		//  Compute the n-k last rows of A' = P A^T P^T in K2_
		// A = A . P^t
		applyP( F, FflasRight, FflasTrans, N, 0, R, A, lda, Pk);

		// Copy K2_ = (A'_2)^t
		for (Ki = K21, Ai = A+R; Ki != K21 + Nrest*ldk; Ai++, Ki+=ldk-N)
			for ( size_t j=0; j<N*lda; j+=lda )
				*(Ki++) = *(Ai+j);

		// A = A . P : Undo the permutation on A
		applyP( F, FflasRight, FflasNoTrans, N, 0, R, A, lda, Pk);

		// K2_ = K2_ . P^t (=  ( P A^t P^t )2_ ) 
		applyP( F, FflasRight, FflasTrans, Nrest, 0, R, K21, ldk, Pk);

		// K21 = K21 . S1^-1
		ftrsm (F, FflasRight, FflasUpper, FflasNoTrans, FflasNonUnit, Nrest, R,
		       one, K, ldk, K21, ldk);  

		typename Field::Element * Arec = new typename Field::Element[Nrest*Nrest];
		size_t ldarec = Nrest;
		
		// Creation of the matrix A2 for recursive call 
		for (Ki = K22,  Ai = Arec;
		     Ki != K22 + Nrest*ldk;
		     Ki += (ldk-Nrest) )
			for ( size_t j=0; j<Nrest; ++j )
				*(Ai++) = *(Ki++);
		fgemm (F, FflasNoTrans, FflasNoTrans, Nrest, Nrest, R, mone,
		       K21, ldk, K+R, ldk, one, Arec, ldarec);

		std::list<Polynomial> polyList;
		polyList.clear();

		// Recursive call on the complementary subspace
		CharPoly(F, polyList, Nrest, Arec, ldarec);
		delete[] Arec;
		frobeniusForm.merge(polyList);
	}

	delete[] Pk;
	delete[] Qk;
	size_t deg = c+1;
	for (size_t i=0; i<Mk; ++i)
 		dA[i] = dK[i];
	bk_idx = 0;
	
	typename Field::Element *Arp = new typename Field::Element[Ncurr*Ma];
	typename Field::Element *Ac = new typename Field::Element[Ncurr*Ma];
	size_t ldac = Ma;
	size_t ldarp = Ncurr;
	
	for (size_t i=0; i < Ncurr; ++i)
 		for (size_t j=0; j<Ma; ++j)
			*(K+i*ldk+j) = *(Ac + i*Ma +j) = *(K4 + i + (j)*ldk);
	delete[] K4;

	size_t block_idx, it_idx, rp_val;

	// Main loop of the arithmetic progession
	while ((nb_full_blocks >= 1) && (Mk > 1)) {
		delete[] K;
		delete[] K3;
		K = new typename Field::Element[Ncurr*Ma];
		K3 = new typename Field::Element[Ncurr*Ma];
		ldk = Ma;
 
		// Computation of the rank profile
		for (size_t i=0; i < Ncurr; ++i)
			for (size_t j=0; j < Ma; ++j)
				*(Arp + j*ldarp + Ncurr-i-1) = *(Ac + i*ldac + j);
		for (size_t i=0; i<2*Ncurr; ++i)
			rp[i] = 0;
		size_t R;
		try{
			R = SpecRankProfile (F, Ma, Ncurr, Arp, ldarp, deg-1, rp);
		} catch (CharpolyFailed){
			delete[] Arp; delete[] Ac; delete[] K; delete[] K3;
			delete[] rp; delete[] dA; delete[] dK;
			throw CharpolyFailed();
		}
		if (R < Ncurr){
			//std::cerr<<"FAIL R<Ncurr"<<std::endl;
			delete[] Arp; delete[] Ac; delete[] K; delete[] K3;
			delete[] rp; delete[] dA; delete[] dK;
			throw CharpolyFailed();
		}

		// Computation of the degree vector dK
		it_idx = 0;
		rp_val = 0;
		size_t g = 0;
		size_t dtot=0;
		block_idx = 0;
		nb_full_blocks = 0;
		while (dtot<Ncurr){
			do {g++; rp_val++; it_idx++;}
			while ( /*(g<Ncurr ) &&*/ (rp[g] == rp_val) && (it_idx < deg ));
			if ((block_idx)&&(it_idx > dK[block_idx-1])){
				delete[] Arp; delete[] Ac;delete[] K; delete[] K3;
				delete[] rp; delete[] dA; delete[] dK;
				throw CharpolyFailed();
				//std::cerr<<"FAIL d non decroissant"<<std::endl;
				//exit(-1);
			}
			dK[block_idx++] = it_idx;
			dtot += it_idx;
			if (it_idx == deg)
				nb_full_blocks ++;
			it_idx=0;
			rp_val = rp[g];
		}

		Mk = block_idx;
				
		// Selection of dense colums of K 
		for (size_t i=0; i < nb_full_blocks; ++i){
			fcopy (F, Ncurr, K+i, ldk, Ac+i, ldac);
		}
		
		// K <- QK K
		size_t pos = nb_full_blocks*(deg-1);
		for (size_t i = nb_full_blocks; i < Mk; ++i){
			for (size_t j=0; j<Ncurr; ++j)
				F.assign (*(K + i + j*ldk), zero);
			F.assign (*(K + i + (pos + dK[i]-1)*ldk), one);
			pos += dA[i];
		}

		// Copying K3 <- K
		for (size_t i=0; i<Mk; ++i)
			fcopy (F, Ncurr, K3+i, ldk, K+i, ldk);
		CompressRowsQK (F, Mk, K3 + nb_full_blocks*(deg-1)*ldk, ldk,
				Arp, ldarp, dK+nb_full_blocks, deg, Mk-nb_full_blocks);

		// K <- PA K
		CompressRows (F, nb_full_blocks, K, ldk, Arp, ldarp, dA, Ma);
		
		// A <- newQA^T K (compress)
		CompressRowsQA (F, Ma, Ac, ldac, Arp, ldarp, dA, Ma);
		
		// K <- A K
		fgemm (F, FflasNoTrans, FflasNoTrans, Ncurr-Ma, nb_full_blocks, Ma, one,
		       Ac, ldac, K+(Ncurr-Ma)*ldk, ldk, one, K, ldk);
		fgemm (F, FflasNoTrans, FflasNoTrans, Ma, nb_full_blocks, Ma, one,
		       Ac+(Ncurr-Ma)*ldac, ldac, K+(Ncurr-Ma)*ldk, ldk, zero, Arp, ldarp);
		for (size_t i=0; i< Ma; ++i)
			fcopy(F, nb_full_blocks, K+(Ncurr-Ma+i)*ldk, 1, Arp+i*ldarp, 1);
		
		// Copying the last rows of A times K
		offset = (deg-2)*nb_full_blocks;
		for (size_t i = nb_full_blocks; i < Mk; ++i) {
			for (size_t j=0; j<Ncurr; ++j)
				F.assign(*(K+i+j*ldk), zero);
			if (dK[i] == dA[i]) // copy the column of A
				fcopy (F, Ncurr, K+i, ldk, Ac+i, ldac);
			else{
				F.assign (*(K + i + (offset+dK[i]-1)*ldk), one);
			}
			offset += dA[i]-1;
		}
				
		// K <- QA K
		DeCompressRowsQA (F, Mk, Ncurr, K, ldk, Arp, ldarp, dA, Ma);

		// K <- QK^T K
		CompressRowsQK (F, Mk, K + nb_full_blocks*(deg-1)*ldk, ldk, Arp, ldarp,
				dK+nb_full_blocks, deg, Mk-nb_full_blocks);
		
		// K <- K^-1 K
		size_t *P=new size_t[Mk];
		size_t *Q=new size_t[Mk];
		if (LUdivine (F, FflasNonUnit, FflasNoTrans, Mk, Mk , K3 + (Ncurr-Mk)*ldk, ldk, P, Q, FfpackLQUP) < Mk){
			// should never happen (not a LAS VEGAS check)
			//std::cerr<<"FAIL R2 < MK"<<std::endl;
			//			exit(-1);
		}
		ftrsm (F, FflasLeft, FflasLower, FflasNoTrans, FflasUnit, Mk, Mk, one,
		       K3 + (Ncurr-Mk)*ldk, ldk, K+(Ncurr-Mk)*ldk, ldk);
		ftrsm (F, FflasLeft, FflasUpper, FflasNoTrans, FflasNonUnit, Mk, Mk, one,
		       K3+(Ncurr-Mk)*ldk, ldk, K+(Ncurr-Mk)*ldk, ldk);
		applyP (F, FflasLeft, FflasTrans, Mk, 0, Mk, K+(Ncurr-Mk)*ldk,ldk, P);
		fgemm (F, FflasNoTrans, FflasNoTrans, Ncurr-Mk, Mk, Mk, mone,
		       K3, ldk, K+(Ncurr-Mk)*ldk,ldk, one, K, ldk);
		delete[] P;
		delete[] Q;
		
		// K <- PK^T K
		DeCompressRows (F, Mk, Ncurr, K, ldk, Arp, ldarp, dK, Mk);
		
		// K <- K PK (dA <- dK)
		if (nb_full_blocks*deg < Ncurr)
			Ma = nb_full_blocks+1;
		else
			Ma = nb_full_blocks;
		
		for (size_t i=0; i< Ma; ++i)
			dA[i] = dK[i];

		// Recovery of the completed invariant factors
		offset = Ncurr-1;
		size_t oldNcurr = Ncurr;
		for (size_t i=Mk-1; i>=nb_full_blocks+1;  --i)
			if (dK[i] >= 1){ 
				Polynomial  P (dK [i]+1);
				F.assign(P[dK[i]], one);
				for (size_t j=0; j < dK[i]; ++j)
					F.neg( P[dK[i]-j-1], *(K + i + (offset-j)*ldk));
				frobeniusForm.push_front(P);
				offset -= dK[i];
				Ncurr -= dK[i];
			}
		for (size_t i= offset+1; i<oldNcurr; ++i)
			for (size_t j=0; j<nb_full_blocks+1; ++j){
				if (!F.isZero( *(K+i*ldk+j) )){
					//std::cerr<<"FAIL C != 0"<<std::endl;
					delete[] rp; delete[] Arp; delete[] Ac;
					delete[] K; delete[] K3;
					delete[] dA; delete[] dK;
					throw CharpolyFailed();
				}
			}
		
		// A <- K
		delete[] Ac; delete[] Arp;
		Ac = new typename Field::Element[Ncurr*Mk];
		ldac = Mk;
		Arp = new typename Field::Element[Ncurr*Mk];
		ldarp=Ncurr;
		for (size_t i=0; i < Ncurr; ++i )
			fcopy (F, Mk, Ac + i*ldac, 1, K + i*ldk, 1);

		deg++;
			
	}

	// Recovery of the first invariant factor
	Polynomial Pl(dK [0]+1);
	F.assign(Pl[dK[0]], one);
	for (size_t j=0; j < dK[0]; ++j)
		F.neg( Pl[j], *(K  + j*ldk));
	frobeniusForm.push_front(Pl);
	delete[] rp; delete[] Arp; delete[] Ac; delete[] K; delete[] K3;
	delete[] dA; delete[] dK;
	return frobeniusForm;
}

template <class Field>
void FFPACK::CompressRowsQK (Field& F, const size_t M,
			   typename Field::Element * A, const size_t lda,
			   typename Field::Element * tmp, const size_t ldtmp,
			   const size_t * d, const size_t deg,const size_t nb_blocs){

	int currtmp = 0;
	size_t currw = d[0]-1;
	size_t currr = d[0]-1;
	for (int i = 0; i< int(nb_blocs)-1; ++i){
		for (int j = d[i]-1; j<int(deg)-1; ++j, currr++, currtmp++)
			fcopy(F, M, tmp + currtmp*ldtmp, 1,  A + currr*lda, 1);
		for (int j=0; j < int(d[i+1]) -1; ++j, currr++, currw++){
			fcopy(F, M, A + (currw)*lda, 1, A+(currr)*lda, 1);
		}
	}
	for (int i=0; i < currtmp; ++i, currw++){
		fcopy (F, M, A + (currw)*lda, 1, tmp + i*ldtmp, 1);
	}
}

template <class Field>
void FFPACK::CompressRows (Field& F, const size_t M,
			     typename Field::Element * A, const size_t lda,
			     typename Field::Element * tmp, const size_t ldtmp,
			     const size_t * d, const size_t nb_blocs){

	size_t currd = d[0]-1;
	size_t curri = d[0]-1;
	for (int i = 0; i< int(nb_blocs)-1; ++i){
		fcopy(F, M, tmp + i*ldtmp, 1,  A + currd*lda, 1);
		for (int j=0; j < int(d[i+1]) -1; ++j){
			fcopy(F, M, A + (curri++)*lda, 1, A+(currd+j+1)*lda, 1);
		}
		currd += d[i+1];
	}
	for (int i=0; i < int(nb_blocs)-1; ++i){
		fcopy (F, M, A + (curri++)*lda, 1, tmp + i*ldtmp, 1);
	}
}

template <class Field>
void FFPACK::DeCompressRows (Field& F, const size_t M, const size_t N,
			     typename Field::Element * A, const size_t lda,
			     typename Field::Element * tmp, const size_t ldtmp,
			     const size_t * d, const size_t nb_blocs){
	
	for (int i=0; i<int(nb_blocs)-1; ++i)
		fcopy(F, M, tmp + i*ldtmp, 1, A + (N-nb_blocs+i)*lda, 1);
	
	size_t w_idx = N - 2;
	size_t r_idx = N - nb_blocs - 1;
	for (int i = int(nb_blocs)-2; i>=0; --i){
		for (size_t j = 0; j<d[i+1]-1; ++j)
			fcopy (F, M, A + (w_idx--)*lda, 1, A + (r_idx--)*lda, 1);
		fcopy (F, M, A + (w_idx--)*lda, 1, tmp + i*ldtmp, 1);
	}
}

template <class Field>
void FFPACK::DeCompressRowsQK (Field& F, const size_t M, const size_t N,
			       typename Field::Element * A, const size_t lda,
			       typename Field::Element * tmp, const size_t ldtmp,
			       const size_t * d, const size_t deg,const size_t nb_blocs){
	
	size_t zeroblockdim = 1; // the last block contributes with 1
	size_t currtmp = 0;
	for (int i=0; i<int(nb_blocs)-1; ++i)
		zeroblockdim += deg - d[i];
	for (int i=0; i < zeroblockdim - 1; ++i, ++currtmp)
		fcopy(F, M, tmp + currtmp*ldtmp, 1,  A + (N - zeroblockdim +i)*lda, 1);
	currtmp--;
	size_t w_idx = N - 2;
	size_t r_idx = N - zeroblockdim - 1;

	for (int i = int(nb_blocs)-2; i>=0; --i){
		for (size_t j = 0; j < d [i+1] - 1; ++j)
			fcopy (F, M, A + (w_idx--)*lda, 1, A + (r_idx--)*lda, 1);
		for (size_t j = 0; j < deg - d[i]; ++j)
			fcopy (F, M, A + (w_idx--)*lda, 1, tmp + (currtmp--)*ldtmp, 1);
	}
}

template <class Field>
void FFPACK::CompressRowsQA (Field& F, const size_t M,
			     typename Field::Element * A, const size_t lda,
			     typename Field::Element * tmp, const size_t ldtmp,
			     const size_t * d, const size_t nb_blocs){

	size_t currd = 0;
	size_t curri = 0;
	for (size_t i = 0; i< nb_blocs; ++i){
		fcopy(F, M, tmp + i*ldtmp, 1,  A + currd*lda, 1);
		for (size_t j=0; j < d[i] -1; ++j)
			fcopy(F, M, A + (curri++)*lda, 1, A+(currd+j+1)*lda, 1);
		currd += d[i];
	}
	for (size_t i=0; i < nb_blocs; ++i)
		fcopy (F, M, A + (curri++)*lda, 1, tmp + i*ldtmp, 1);
}

template <class Field>
void FFPACK::DeCompressRowsQA (Field& F, const size_t M, const size_t N,
			       typename Field::Element * A, const size_t lda,
			       typename Field::Element * tmp, const size_t ldtmp,
			       const size_t * d, const size_t nb_blocs){
	
	for (size_t i=0; i<nb_blocs; ++i)
		fcopy(F, M, tmp + i*ldtmp, 1, A + (N-nb_blocs+i)*lda, 1);

	size_t w_idx = N - 1;
	size_t r_idx = N - nb_blocs - 1;
	for (int i = int(nb_blocs)-1; i>=0; --i){
		for (size_t j = 0; j<d[i]-1; ++j)
			fcopy (F, M, A + (w_idx--)*lda, 1, A + (r_idx--)*lda, 1);
		fcopy (F, M, A + (w_idx--)*lda, 1, tmp + i*ldtmp, 1);
	}
}
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
