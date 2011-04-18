/* ffpack/ffpack_krylovelim.inl
 * Copyright (C) 2007 Clement Pernet
 *
 * Written by Clement Pernet <Clement.Pernet@imag.fr>
 *
 * See COPYING for license information.
 */

#ifndef __LINBOX_ffpack_krylovelim_INL
#define __LINBOX_ffpack_krylovelim_INL

#ifndef MIN
#define MIN(a,b) (a<b)?a:b
#endif
#ifndef MAX
#define MAX(a,b) (a<b)?b:a
#endif

//#define LB_DEBUG

// A is m x n with m <= n
// Ensures : rankprof is the row rankprofil of the matrix k x n matrix B formed as follows (k = sum d_i):
// for d_i <  j < d_{i+1} the jth row B_j = is e_j
// B_{d_i} = A_i
// iterates must be initialized by [1, 2, ...]
// inviterates is the inverse finction of iterates
template <class Field>
inline size_t 
FFPACK::KrylovElim( const Field& F, const size_t M, const size_t N,		
		    typename Field::Element * A, const size_t lda, size_t*P, 
		    size_t *Q, const size_t deg, size_t *iterates,  size_t * inviterates,size_t  maxit,
		    size_t virt)
{

	if ( !(M && N) ) return 0;
	typedef typename Field::Element elt;
	static elt Mone, one, zero;
	F.init(Mone, -1.0);
	F.init(one,1.0);
	F.init(zero,0.0);

	if (M == 1){
		virt += deg;
		for (size_t i=0; i<virt; ++i)
			if (iterates[i]){
				//	cerr<<"A["<<N-i-1<<"]=0"<<endl;
				F.assign (A [N-iterates[i]], zero);
			}
		size_t ip=0;
		//while (ip<N && !F.isUnit(*(A+ip)))ip++;
		while (ip<N && F.isZero (*(A+ip))) {
			//			cout<<(*(A+ip))<<" ";
			ip++;
		}
		//cout<<endl;
		*Q=0;
		//cerr<<"ip = "<<ip<<endl;
		if (ip==N){ // current row is zero
			*P=0;
			//cerr<<"return 0"<<endl;
			return 0;
		}
		*P = ip;
			// 	cerr<<"iterates = ";
// 		cerr<<"iterates avant"<<endl;
// 		for (size_t i=0; i<N; ++i)
//  			cerr<<iterates[i]<<" ";
//  		cerr<<endl;
//  		cerr<<"inviterates["<<N<<"-"<<ip<<"] ="<<inviterates[N-ip]<<endl;
//  		cerr<<"iterates[inviterates[N-ip]] ="<<iterates[inviterates[N-ip]]<<endl;
		iterates [inviterates[N-ip] -1 ] = 0;
		if (ip >0){
			iterates [inviterates[N] -1 ] = N-ip;
			inviterates[N-ip] = inviterates[N];
		}
// 		cerr<<"iterates apres"<<endl;
// 		for (size_t i=0; i<N; ++i)
//  			cerr<<iterates[i]<<" ";
//  		cerr<<endl;
		//cout<<"iterates ["<<N<<"-1-"<<ip<<"] = 0"<<endl;

		if (ip!=0){
			// swap the pivot
			typename Field::Element tmp=*A;
			*A = *(A+ip);
			*(A+ip) = tmp;
		}
		return 1;
	} else { // MN>1
		size_t Nup = M>>1;
		size_t Ndown =  M - Nup;

		// Recursive call on NW
		size_t R = KrylovElim (F,  Nup, N, A, lda, P, Q, deg, iterates, inviterates, maxit, virt);
		
		typename Field::Element *Ar = A + Nup*lda; // SW
		typename Field::Element *Ac = A + R;     // NE
		typename Field::Element *An = Ar + R;    // SE

		if (R){
			// Ar <- Ar.P
			applyP (F, FflasRight, FflasTrans, Ndown, 0, R, Ar, lda, P); 
			// Ar <- Ar.U1^-1
			ftrsm( F, FflasRight, FflasUpper, 
			       FflasNoTrans, FFLAS::FflasNonUnit, Ndown, R, 
			       one, A, lda, Ar, lda);
			// An <- An - Ar*Ac
			fgemm( F, FflasNoTrans, FflasNoTrans, Ndown, N-R, R,
			       Mone, Ar, lda, Ac, lda, one, An, lda);
		}
		// Recursive call on SE
		size_t R2 = KrylovElim (F, Ndown, N-R, An, lda,P+R, Q+Nup, deg, iterates, inviterates, maxit, MIN(maxit-deg,(virt+Nup*deg)));

		for (size_t i = R; i < R + R2; ++i)
			P[i] += R;
		if (R2)
			// An <- An.P2
			applyP (F, FflasRight, FflasTrans, Nup, R, R+R2, A, lda, P); 
		
		// Non zero row permutations
		for (size_t i = Nup; i < M; i++)
			Q[i] += Nup;
		if (R < Nup){
			// Permutation of the 0 rows
			for ( size_t i = Nup, j = R ; i < Nup + R2; ++i, ++j){
				fcopy( F, N - j, A + j*(lda + 1), 1, A + i*lda + j, 1);
				for (typename Field::Element *Ai = A + i*lda + j;
				     Ai != A + i*lda + N; ++Ai)
					F.assign (*Ai, zero);
				size_t t = Q[j];
				Q[j]=Q[i];
				Q[i] = t;
			}
		}
		return R + R2;
	}
}

template <class Field>
size_t 
FFPACK::SpecRankProfile (const Field& F, const size_t M, const size_t N,
			 typename Field::Element * A, const size_t lda, const size_t deg,
			 size_t *rankProfile)
{

	//size_t deg = (N-1)/M+1; // Number of trivial iterates per blocs
	size_t * Q = new size_t[M];
	size_t * P = new size_t[N];
	size_t * iterates = new size_t[N];
	size_t * inviterates = new size_t[N+1];
	for (size_t i=0; i < N; ++i)
		inviterates[i+1] = iterates[i] = i+1;

	size_t R = KrylovElim (F, M, N, A, lda, P, Q, deg, iterates, inviterates, N,0);
	//cerr<<"Apres tout iterates = "<<endl;

// 	for (size_t i=0; i<N; ++i)
// 		cerr<<iterates[i]<<" ";
// 	cerr<<endl;

	size_t curr_row = 0;
	size_t it_idx = 0;
	size_t bk_idx = 0;
	size_t rp_idx = 0;
	bool dependent = false;
	for (size_t i=0; i<M; ++i){
		//	cerr<<"Block "<<i<<endl;
		for (size_t j=0; j<deg; ++j){
			if (curr_row < N+M -1){
				if (iterates[it_idx++]){
					rankProfile [rp_idx++] = curr_row;
					if (dependent){
#ifdef LB_DEBUG
						std::cerr<<"FAIL itere dependant intercale"<<std::endl;
#endif
						delete[] P;
						delete[] Q;
						delete[] iterates;
						delete[] inviterates;
						throw CharpolyFailed();
					}
#ifdef LB_DEBUG
					std::cerr<<"X";
#endif					
				}
				else{
					dependent = true;
#ifdef LB_DEBUG
					std::cerr<<"O";
#endif
				}
				curr_row++;
			}
		}
		dependent = false;
		if ((Q [bk_idx] == i)&&(i<R)){
			rankProfile [rp_idx++] = curr_row;
#ifdef LB_DEBUG
 			std::cerr<<"V"<<std::endl;
#endif
 			bk_idx++;
 		}
#ifdef LB_DEBUG
 		else
 			std::cerr<<"W"<<std::endl;
#endif
		curr_row++;
	}
	delete[] P;
	delete[] Q;
	delete[] inviterates;
	delete[] iterates;
	
	return rp_idx;
}

#endif //__LINBOX_ffpack_krylovelim_INL

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
