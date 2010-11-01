/* linbox/ffpack/ffpack_charpoly_kgfast.inl
 * Copyright (C) 2004 Clement Pernet
 *
 * Written by Clement Pernet <Clement.Pernet@imag.fr>
 *
 * See COPYING for license information.
 */

#ifndef __LINBOX_ffpack_charpoly_kgfastgeneralized_INL
#define __LINBOX_ffpack_charpoly_kgfastgeneralized_INL

#ifndef MIN
#define MIN(a,b) (a<b)?a:b
#endif

//---------------------------------------------------------------------
// CharPoly: Compute the characteristic polynomial of A using 
// Keller-Gehrig's fast algorithm. 
//---------------------------------------------------------------------

//#define LB_DEBUG

#ifdef LB_DEBUG
#include "tests/Matio.h"

template <class Field>
void printA(const Field& F,
	    std::ostream& os,
	    const typename Field::Element * E,
	    const typename Field::Element * C, 
	    const size_t lda,
	    const size_t*B, 
	    const size_t*T, 
	    const size_t me,const size_t mc, const size_t lambda, const size_t mu){
	
	typename Field::Element * A = buildMatrix(F,E,C,lda,B,T,me,mc,lambda,mu);
	size_t N = mc+me+lambda+mu;
	write_field(F,os,A,N,N,N);
	delete[] A;
}
#endif

template <class Field>
typename Field::Element * buildMatrix (const Field& F,
				       const typename Field::Element * E,
				       const typename Field::Element * C, 
				       const size_t lda,
				       const size_t*B, 
				       const size_t*T, 
				       const size_t me,
				       const size_t mc, 
				       const size_t lambda, 
				       const size_t mu)
{
	
	typename Field::Element zero,one;
	F.init (zero,0UL);
	F.init (one,1UL);
	size_t N = mc+me+lambda+mu;
	typename Field::Element * A = new typename Field::Element[N*N];
	for (size_t j=0; j<lambda+me;++j)
		if (B[j] < N){
			for (size_t i=0;i<N;++i)
				F.assign( *(A+i*N+j),zero);
			F.assign( *(A+B[j]*lda+j), one);
		} else {
			FFLAS::fcopy (F, N, A+j, N, E+B[j]-N, lda);
		}
	for (size_t j=lambda+me; j<lambda+me+mu; ++j)
		for (size_t i=0;i<N;++i)
			F.assign( *(A+i*N+j),zero);
	for (size_t i=0; i<mu; ++i)
		F.assign( *(A+(lambda+me+mc+i)*lda+lambda+me+T[i]), one);
	for (size_t j=0; j<mc; ++j)
		FFLAS::fcopy(F,N,A+N-mc+j,N,C+j,lda);
	return A;
}

template <class Field, class Polynomial>
std::list<Polynomial>&
FFPACK::KGFast_generalized (const Field& F, std::list<Polynomial>& charp, 
			    const size_t N,
			    typename Field::Element * A, const size_t lda)
{
	
	//std::cerr<<"Dans KGFast"<<std::endl;
	static typename Field::Element one, zero, mone;
	F.init(one, 1UL);
	F.neg(mone, one);
	F.init(zero, 0UL);
	size_t mc=N>>1; // Matrix A is transformed into a mc_Frobenius form
	size_t me=N-mc;
	// B[i] = j, the row of the 1 if the col Ai is sparse; 
	// B[i] = n+k, if the col Ai is the kth col of E
	size_t * B = new size_t[N];
	bool * allowedRows = new bool[N];
	for (size_t i=0;i<(N+1)/2;++i) 
		allowedRows[i]=true;
	// T[i] = j si T_i,j = 1
	size_t * T = new size_t[N];
	for (size_t i=0;i<N;++i) 
		T[i]=i;
	size_t lambda=0;
	
	typename Field::Element * C, *E = A;
#ifdef LB_DEBUG
	std::cerr<<"Debut KGFG"<<std::endl
	  <<" ----------------------------"<<std::endl;
#endif
	while (mc > 0) {
#ifdef LB_DEBUG
 		std::cerr<<"Boucle1: mc,me,lambda="<<mc<<" "<<me<<" "<<lambda<<std::endl;
		// 		write_field (F, std::cerr, A, N, N, lda);
#endif
		size_t mu=0;
		C = A + (N-mc);
		for (size_t i = 0; i<me;++i)
			B[lambda+i] = N+i;
#ifdef LB_DEBUG
		for (size_t i=0;i<lambda+me;++i)
			std::cerr<<"B["<<i<<"] = "<<B[i]<<std::endl;
		//std::cerr<<std::endl<<"mc="<<mc<<":";
#endif
		while (mu < N-mc) {
#ifdef LB_DEBUG
 			std::cerr<<"Boucle2: mu,me,lambda="<<mu<<" "<<me<<" "<<lambda<<std::endl;
			printA(F,std::cerr<<"A="<<std::endl,E,C,lda,B,T,me,mc,lambda,mu);

			// B1 <- C1^-1.B1
			std::cerr<<"Forming LUP";
#endif
			size_t ncols = ((mu==0)||(mc<=mu))?mc:mc-mu;
			typename Field::Element * LUP = new typename Field::Element[(lambda+me)*ncols];
			for (size_t i=0;i < lambda + me; ++i)
				if (allowedRows[i])
					fcopy (F, ncols, LUP+i*ncols, 1, C+i*lda, 1);
				else
					for (size_t j = 0; j < ncols; ++j)
						F.assign (*(LUP+i*ncols+j), zero);
#ifdef LB_DEBUG
			std::cerr<<"..done"<<std::endl;
			write_field (F,std::cerr<<"LUP="<<std::endl,LUP,lambda+me,ncols,ncols);
			std::cerr<<"LQUP(C1)";
#endif
			size_t * P = new size_t[ncols];
			size_t * Q = new size_t[lambda+me];
			for (size_t i=0; i<ncols;++i)
				P[i]=0;
			for (size_t i=0; i<lambda+me;++i)
				Q[i]=0;

			size_t r = LUdivine (F, FflasNonUnit, FflasNoTrans, lambda + me, ncols, LUP, ncols, 
					   P, Q, FfpackLQUP);
#ifdef LB_DEBUG
			std::cerr<<"..done"<<std::endl;
#endif

			if (r==0){
				if ((lambda == 0) && (ncols == mc)){
					std::cerr<<"BLOCAGE lambda=0!!!"<<std::endl;
					//Rec call on the leading block
					KGFast_generalized (F, charp, me, A, lda);
					
					//Rec call on the trailing block
					typename Field::Element * At = buildMatrix(F,E,C,lda,B,T,me,mc,lambda,mu);
					KGFast_generalized (F, charp, N-me, At+me*(lda+1), lda);
					delete[] At;
					exit(-1);

				} else if (me != 0) {
					std::cerr<<"BLOCAGE me!=0!!!"<<std::endl;
					exit(-1);
								
				}
				else {
					for (int i=mu; i>=0; --i)
						T[i+lambda] = T[i]+lambda;
					for (size_t i=0; i< lambda; ++i)
						T[B[i]-mc-1] = i;
					mu += lambda;
					lambda = 0;
					break;
				}
				//std::cerr<<"BLOCAGE !!!"<<std::endl;
				//exit(-1);
			}

#ifdef LB_DEBUG
			std::cerr<<"Forming genreric rank profil C1";
			// form the generic rank profil block C1 Q^TPAP^TQ
			for (size_t i=0;i<r;++i)
				std::cerr<<"P["<<i<<"] = "<<P[i]<<std::endl;
#endif
			applyP (F, FflasRight, FflasTrans, N, 0, r, C, lda, P);
#ifdef LB_DEBUG
			std::cerr<<".";
#endif
			//printA(F,cerr<<"A="<<std::endl,E,C,lda,B,T,me,mc,lambda,mu);
			// (E, C) <- P(E, C)
			applyP (F, FflasLeft, FflasNoTrans, me, 0, r, E+(N-mc)*lda, lda, P);
#ifdef LB_DEBUG
			std::cerr<<".";
			//printA(F,std::cerr<<"A="<<std::endl,E,C,lda,B,T,me,mc,lambda,mu);
#endif
			applyP (F, FflasLeft, FflasNoTrans, mc, 0, r, C+(N-mc)*lda, lda, P);
#ifdef LB_DEBUG
			std::cerr<<".";
#endif
			//printA(F,std::cerr<<"A="<<std::endl,E,C,lda,B,T,me,mc,lambda,mu);
			// T <- P T

			// !!!!!!!Attention -> ajouter le traitement du cas 0<mu<mc
			for (size_t k = 0; k<r; ++k)
				if (P[k] > (size_t) k){
					if ((mu>=mc-k)){
#ifdef LB_DEBUG
						std::cerr<<"// on permute LN-mc+k et L_N-mc+P[k]"<<std::endl;
#endif
						size_t tmp = T[mu-mc+k];
						T[mu-mc+k] = T[mu-mc+P[k]];
						T[mu-mc+P[k]] = tmp;
					}
					else if (mu){
						std::cerr<<"CAS MU < MC - k"<<std::endl;
						exit(-1);
					}
					// Updating B to be improved (tabulated B^-1)
					for (size_t i=0; i<lambda+me; ++i){
						if (B[i] == N-mc+k)
							B[i] = N-mc+P[k];
						else if (B[i] == N-mc+P[k])
							B[i] = N-mc+k;
					}
					
				}
#ifdef LB_DEBUG
			std::cerr<<".";
			//printA(F,std::cerr<<"A="<<std::endl,E,C,lda,B,T,me,mc,lambda,mu);
#endif

			// (E, C) <- Q^T(E, C) 
			applyP (F, FflasLeft, FflasTrans, me, 0, r, E, lda, Q);
#ifdef LB_DEBUG
			std::cerr<<".";
			//printA(F,std::cerr<<"A="<<std::endl,E,C,lda,B,T,me,mc,lambda,mu);
#endif
			applyP (F, FflasLeft, FflasTrans, mc, 0, r, C, lda, Q);
#ifdef LB_DEBUG
			std::cerr<<".";
#endif
			// F <- Q^T F
			size_t * tempP = new size_t[lambda+me+mc];
			for (size_t i=0; i< lambda+me+mc; ++i)
				tempP[i] = i;
			for (int i = r-1; i>=0; --i)
				if (Q[i] > (size_t) i){
#ifdef LB_DEBUG
					std::cerr<<"Permutation de tempP["<<i
					    <<"] et tempP["<<Q[i]<<"]"<<std::endl;
#endif
					// on permute LN-mc+k et L_N-mc+P[k]
					size_t tmp = tempP[i];
					tempP[i] = tempP[Q[i]];
					tempP[Q[i]] = tmp;
				}
				
#ifdef LB_DEBUG
			std::cerr<<".";
#endif
			for (size_t i=0; i < lambda+me; ++i)
				if (B[i] < N)
					B[i] = tempP[B[i]];
#ifdef LB_DEBUG
			std::cerr<<".";
#endif
			delete[] tempP;

#ifdef LB_DEBUG
			std::cerr<<std::endl<<"Avant B<-BQ"<<std::endl;
			for (size_t i=0; i<lambda+me;++i)
				std::cerr<<"B["<<i<<"] = "<<B[i]<<std::endl;
#endif
			// B <- B Q
			for (int k = r-1; k>=0; --k)
				if (Q[k] > (size_t) k){
					// on permute Ck et C_Q[k]
					size_t tmp = B[k];
					B[k] = B[Q[k]];
					B[Q[k]] = tmp;
				}
#ifdef LB_DEBUG
			std::cerr<<"Apres"<<std::endl;
			for (size_t i=0; i<lambda+me;++i)
				std::cerr<<"B["<<i<<"] = "<<B[i]<<std::endl;

			std::cerr<<".";
#endif

			// grouping the bloc L in LUP
			for (size_t i=0; i<r; ++i)
				if (Q[i]>i)
					fcopy(F, i, LUP+i*mc, 1, LUP+Q[i]*mc,1);
#ifdef LB_DEBUG
			std::cerr<<"..done"<<std::endl;
			
			printA(F,std::cerr<<"A="<<std::endl,E,C,lda,B,T,me,mc,lambda,mu);


// 			std::cerr<<"LUP="<<std::endl;
// 			write_field (F, std::cerr, LUP, mc, mc, mc);
            //std::cerr<<" "<<r;

			// E'1 <- C11^-1 E1
			std::cerr<<"// E'1 <- C11^-1 E1";
#endif

			ftrsm(F, FflasLeft, FflasLower, FflasNoTrans, FflasUnit,
			   r, me, one, LUP, mc , E, lda);
			ftrsm(F, FflasLeft, FflasUpper, FflasNoTrans, FflasNonUnit, 
			   r, me, one, LUP, mc , E, lda);
#ifdef LB_DEBUG
			std::cerr<<"..done"<<std::endl;
			printA(F,std::cerr<<"A="<<std::endl,E,C,lda,B,T,me,mc,lambda,mu);

			// C'12 <- C11^-1 C12 
			std::cerr<<"// C'12 <- C11^-1 C12"; 
#endif
			ftrsm(F, FflasLeft, FflasLower, FflasNoTrans, FflasUnit,
			   r, mc-r, one, LUP, mc , C+r, lda);
			ftrsm(F, FflasLeft, FflasUpper, FflasNoTrans, FflasNonUnit, 
			   r, mc-r, one, LUP, mc , C+r, lda);
			delete[] LUP;
			delete[] P;
			delete[] Q;
#ifdef LB_DEBUG
			std::cerr<<"..done"<<std::endl;
			printA(F,std::cerr<<"A="<<std::endl,E,C,lda,B,T,me,mc,lambda,mu);
#endif

// 			std::cerr<<"Apres B1<-C1^-1"<<std::endl;
// 			write_field (F, std::cerr, A, N, N, lda);
            
			// E'2 <- E2 - C21.E'1
#ifdef LB_DEBUG
			std::cerr<<"// E'2 <- E2 - C21.E'1";
#endif
			fgemm(F, FflasNoTrans, FflasNoTrans, N-r, me, r, 
			   mone, C+r*lda, lda, E, lda, 
			   one, E+r*lda, lda);
#ifdef LB_DEBUG
			std::cerr<<"..done"<<std::endl;			
			printA(F,std::cerr<<"A="<<std::endl,E,C,lda,B,T,me,mc,lambda,mu);
			// C'22 <- C22 - C21.C'12
			std::cerr<<"// C'22 <- C22 - C21.C'12";
#endif
			fgemm(F, FflasNoTrans, FflasNoTrans, N-r, mc-r, r, 
			   mone, C+r*lda, lda, C+r, lda, 
			   one, C+r*(lda+1), lda);
#ifdef LB_DEBUG
			std::cerr<<"..done"<<std::endl;
			printA(F,std::cerr<<"A="<<std::endl,E,C,lda,B,T,me,mc,lambda,mu);
			
// 			std::cerr<<"Apres B2<-B2-C2.B1"<<std::endl;
//             write_field (F, std::cerr, A, N, N, lda);

			// Shifting E: E1;E2 -> E2;E1
			std::cerr<<"// Shifting E: E1;E2 -> E2;E1";
#endif
			typename Field::Element * tmp = new typename Field::Element[r*me];
			for (size_t i=0; i<r; ++i)
				fcopy (F, me, tmp+i*me, 1, E+i*lda, 1);
			for (size_t i=r; i< N; ++i)		
				fcopy (F, me, E+(i-r)*lda, 1, E+i*lda, 1);
			for (size_t i=0; i<r; ++i)
				fcopy (F, me, E+(i+N-r)*lda, 1, tmp+i*me, 1);
			delete[] tmp;
#ifdef LB_DEBUG
			std::cerr<<"..done"<<std::endl;

			// Shifting C_{*,2}: C_{1,2};C_{2,2} -> C_{2,2};C_{1,2}
			std::cerr<<"// Shifting C_{*,2}: C_{1,2};C_{2,2} -> C_{2,2};C_{1,2}";
#endif
			tmp = new typename Field::Element[r*(mc-r)];
			for (size_t i=0; i<r; ++i)
				fcopy (F, mc-r, tmp+i*(mc-r), 1, C+r+i*lda, 1);
			for (size_t i=r; i< N; ++i)		
				fcopy (F, mc-r, C+r+(i-r)*lda, 1, C+r+i*lda, 1);
			for (size_t i=0; i<r; ++i)
				fcopy (F, mc-r, C+r+(i+N-r)*lda, 1, tmp+i*(mc-r), 1);
			delete[] tmp;
#ifdef LB_DEBUG
			std::cerr<<"..done"<<std::endl;
			
			printA(F,std::cerr<<"A="<<std::endl,E,C,lda,B,T,me,mc,lambda,mu);

// 			std::cerr<<"Apres shift de B"<<std::endl;
//             write_field (F, std::cerr, A, N, N, lda);

			// C'2 <- T C2
			std::cerr<<"// C'2 <- T C2";
#endif
			// To be improved!!!
			tmp = new typename Field::Element[mu*r];
			typename Field::Element * C2 = C+(N-mu-mc)*lda;
			for (size_t i=0; i<mu; ++i)
				fcopy (F, r, tmp+i*r, 1, C2+T[i]*lda, 1); 
			for (size_t i=0; i<mu; ++i)
				fcopy (F, r, C2+i*lda, 1, tmp+i*r, 1); 
			delete[] tmp;
#ifdef LB_DEBUG
			std::cerr<<"..done"<<std::endl;
			
			// [C'2;C'3] += [E2;E3].C 
			std::cerr<<"// [C'2;C'3] += [E2;E3].C";
#endif
			tmp = new typename Field::Element[me*r];
			for (size_t i=0; i<lambda+me; ++i)
				if (B[i] >= N){
					fcopy (F, r, tmp+(B[i]-N)*r, 1, C+i*lda, 1);
				}
			fgemm (F, FflasNoTrans, FflasNoTrans, mu + r, r, me, 
			    one, E+(N-mu-r)*lda, lda, tmp, r,
			    one, C+(N-mu-mc)*lda, lda);

			delete[] tmp;
#ifdef LB_DEBUG
			std::cerr<<"..done"<<std::endl;

			// shifting [C'2;C'3]
			std::cerr<<"// shifting [C'2;C'3]";
#endif
			tmp = new typename Field::Element[(mc-r)*r];
			typename Field::Element * C4 = C + (N-mc+r)*lda;
			for (size_t i=0; i < (mc-r); ++i){
				fcopy (F, r, tmp+i*r, 1, C4 + i*lda, 1);
			}
			for (int i = N-1; i >= (int) (N -mu-r); --i)
				fcopy (F, r, C+i*lda, 1, C+(i-mc+r)*lda, 1);
#ifdef LB_DEBUG
			std::cerr<<"..done"<<std::endl;

			printA(F,std::cerr<<"A="<<std::endl,E,C,lda,B,T,me,mc,lambda,mu);
			

			// tmp2 <- C'1 (the rows corresponding to E)
			std::cerr<<"// tmp2 <- C'1 (the rows corresponding to E)";
#endif
			typename Field::Element * tmp2 = new typename Field::Element[me*r];
			for (size_t i = 0; i < lambda+me; ++i)
				if (B[i] >= N){
#ifdef LB_DEBUG
					std::cerr<<"saving in row "<<B[i]-N<<std::endl;
#endif
					fcopy (F, r, tmp2+(B[i]-N)*r, 1, C+i*lda, 1);
				}
#ifdef LB_DEBUG
			std::cerr<<"..done"<<std::endl;

			// C'_F[i] <- C_i
			std::cerr<<"// C'_F[i] <- C_i";
			std::cerr<<"lambda,r,me = "<<lambda<<" "<<r<<" "<<me<<std::endl;
#endif
			typename Field::Element * tmp3 = new typename Field::Element[(lambda+me)*r];

			for (size_t i = 0; i < lambda+me; ++i)
				if (B[i] < N){
#ifdef LB_DEBUG
					std::cerr<<"copie de la ligne "<<i<<std::endl;
#endif
					fcopy (F, r, tmp3 + i*r, 1, C + i*lda, 1);
				}
#ifdef LB_DEBUG
			std::cerr<<"1"<<std::endl;
#endif
			for (size_t i = 0; i < N-mu-r; ++i)
				for (size_t j = 0; j < r; ++j)
					F.assign (*(C+i*lda+j), zero);
#ifdef LB_DEBUG
			std::cerr<<"2"<<std::endl;
#endif
			for (size_t i = 0; i < lambda+me; ++i){
#ifdef LB_DEBUG
				std::cerr<<"B["<<i<<"] = "<<B[i]<<std::endl;
#endif
				if (B[i] < N)
					fcopy (F, r, C+(B[i]-r)*lda, 1, tmp3+i*r, 1);
			}
#ifdef LB_DEBUG
			std::cerr<<"3"<<std::endl;
#endif
			delete[] tmp3;
#ifdef LB_DEBUG
			std::cerr<<"..done"<<std::endl;
 
			printA(F,std::cerr<<"A="<<std::endl,E,C,lda,B,T,me,mc,lambda,mu);

			// C'1 += E1 tmp2
			std::cerr<<"// C'1 += E1 tmp2";
#endif
			fgemm(F, FflasNoTrans, FflasNoTrans, N-mu-r, r, me, 
			   one, E, lda, tmp2, r, one, C, lda);
			delete[] tmp2;
#ifdef LB_DEBUG
			std::cerr<<"..done"<<std::endl;
 
			printA(F,std::cerr<<"A="<<std::endl,E,C,lda,B,T,me,mc,lambda,mu);

			// C'_1 += C_2 C4
			std::cerr<<"// C'_1 += C_2 C4";
#endif
			fgemm(F, FflasNoTrans, FflasNoTrans, N, r, mc-r, 
			   one, C+r, lda, tmp, r, one, C, lda);
			delete[] tmp;
#ifdef LB_DEBUG
			std::cerr<<"..done"<<std::endl;

			printA(F,std::cerr<<"A="<<std::endl,E,C,lda,B,T,me,mc,lambda,mu);

			// switching C_1 <-> C_2
			std::cerr<<"// switching C_1 <-> C_2";
#endif
			tmp = new typename Field::Element[N*r];
			for (size_t j = 0; j<r; ++j)
				fcopy (F, N, tmp+j, r, C+j, lda);
			for (size_t j = r; j<mc; ++j)
				fcopy (F, N, C+j-r, lda, C+j, lda);
			for (size_t j = 0; j<r; ++j)
				fcopy (F, N, C+mc-r+j, lda, tmp+j, r);
			delete[] tmp;
#ifdef LB_DEBUG
			std::cerr<<"..done"<<std::endl;

			
			printA(F,std::cerr<<"A="<<std::endl,E,C,lda,B,T,me,mc,lambda,mu);

			
			// update the datastructure:
			std::cerr<<"// update the datastructure:";
#endif
			mu += r;
			tmp2 = new typename Field::Element[N*me];
			size_t nlambda= 0, nme=0;
			for (size_t i=0;i<lambda+me;++i)
				allowedRows[i]=true;
			for (size_t j=r; j < lambda + me; ++j){
				if (B[j] >= N){
#ifdef LB_DEBUG
					std::cerr<<"B["<<j-r<<"] = "<<N+nme<<std::endl;
#endif
					fcopy (F, N, tmp2+nme, me, E+(B[j]-N), lda);
					B[j-r] = N + nme;
					nme++;
				} else {
#ifdef LB_DEBUG
					std::cerr<<"B["<<j-r<<"] = "<<B[j]<<std::endl;
#endif
					B[j-r] = B[j]-r;
					allowedRows[B[j]-r] = false;
					nlambda++;
				}
			}
			for (size_t j=0; j<nme; ++j)
				fcopy (F, N, E+j, lda, tmp2+j, me);
			lambda = nlambda;
			me = nme;
#ifdef LB_DEBUG
			std::cerr<<"..done"<<std::endl;
#endif
			delete[] tmp2;
		}
		// update the datastructure: F <- T
		for (size_t i=0; i<mu; ++i){
#ifdef LB_DEBUG
			std::cerr<<"B[T["<<i<<"]] = "<<"B["<<T[i]<<"] = "<<mc+i<<std::endl;
#endif

			B[T[i]] = mc+i;
			T[i]=i;
		}
		E=C;
		me = mc;
		mc>>=1;
		me -= mc;
		lambda = mu;
		for (size_t i=0;i<me+mc;++i)
			allowedRows[i]=true;
		for (size_t i=me+mc;i<lambda+me+mc;++i)
			allowedRows[i]=false;
		
	}

	Polynomial *minP = new Polynomial();
	minP->resize(N+1);
	minP->operator[](N) = one;
	typename Polynomial::iterator it = minP->begin();
	for (size_t j=0; j<N; ++j, it++){
		F.neg(*it, *(A+N-1+j*lda));
	}
	charp.clear();
	charp.push_front(*minP);
	return charp;
}

#undef LB_DEBUG

#endif // __LINBOX_ffpack_charpoly_kgfastgeneralized_INL

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
