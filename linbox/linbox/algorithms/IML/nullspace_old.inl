/* ---------------------------------------------------------------------
 *
 * -- Integer Matrix Library (IML)
 *    (C) Copyright 2004, 2006 All Rights Reserved
 *
 * -- IML routines -- Version 1.0.1 -- November, 2004
 *
 * Author         : Zhuliang Chen
 * Contributor(s) : Arne Storjohann
 * University of Waterloo -- School of Computer Science
 * Waterloo, Ontario, N2L3G1 Canada
 *
 * ---------------------------------------------------------------------
 *
 * -- Copyright notice and Licensing terms:
 *
 *  Redistribution  and  use in  source and binary forms, with or without
 *  modification, are  permitted provided  that the following  conditions
 *  are met:
 *
 * 1. Redistributions  of  source  code  must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce  the above copyright
 *    notice,  this list of conditions, and the  following disclaimer in
 *    the documentation and/or other materials provided with the distri-
 *    bution.
 * 3. The name of the University,  the IML group,  or the names of its
 *    contributors  may not be used to endorse or promote products deri-
 *    ved from this software without specific written permission.
 *
 * -- Disclaimer:
 *
 * THIS  SOFTWARE  IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES,  INCLUDING,  BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE UNIVERSITY
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,  INDIRECT, INCIDENTAL, SPE-
 * CIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
 * TO,  PROCUREMENT  OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
 * OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEO-
 * RY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT  (IN-
 * CLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
 * THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */

#ifndef __LINBOX_algorithms_iml_nullspace_INL
#define __LINBOX_algorithms_iml_nullspace_INL


#include "linbox/matrix/factorized-matrix.h"
#include "linbox/algorithms/matrix-hom.h"
#include "linbox/algorithms/rational-solver.h"
#include "linbox/randiter/random-prime.h"
#include "linbox/integer.h"
#include "linbox/util/debug.h"
// #include "linbox/algorithms/integer-tools.h"



namespace LinBox { namespace iml {
namespace Protected {
#define _LB_MAX_NS_TRY 5

	template<class Matrix>
	struct StackMat {
		Matrix & first ;
		Matrix & second;
		StackMat(Matrix & A, Matrix & B) :
			first(A),second(B)
		{}
	} ;

	/*! Extraction of submatrices according to permuations.
	 * Let \c P and \c Q be two permutations and \c A be a dense matrix.
	 * Let \c QAP be the resulting permuted matrix. For a given \c r,
	 * we want to extract the left upper \c r square matrix \c C in \c QAP and
	 * the rectangular matrix \c B formed by the remaining elements in the first \c r
	 * rows.
	 * We use the formula \f$\left(QAP\right)_{i,j} = A_{Q^{-1}(i),P(j)}\f$.
	 * @param[out] QAP such that \c QAP=pair(B,C).
	 * @param A input matrix
	 * @param Qt row permutation
	 * @param P column permutation
	 */
	template <class DenseMatrix, class DenseMatrix2, class Permutation>
	StackMat<DenseMatrix> & extractMatrix(  StackMat<DenseMatrix> & QAP,
						const DenseMatrix2 & A,
						const Permutation & Qt, const Permutation & P)
	{
		// std::cout                  << "Extracting matrix..." << std::endl;
		// A.write(std::cout          << "A is : "              << std::endl);
		// P.write(std::cout          << "P^{-1} : " , false )          << std::endl;
		// Qt.write(std::cout          << "Q^{-1} : " , false )          << std::endl;
		DenseMatrix & D = QAP.first ;
		DenseMatrix & G = QAP.second ;
		// Permutation P(Pt) ;
		// P.write(std::cout<< 'a')<<std::endl;
		// P.Invert();
		// P.Compress();
		// P.write(std::cout<< 'b')<<std::endl;
		// Pt.write(std::cout<< 'a')<<std::endl;
		size_t r = D.rowdim();
		if (!P.isIdentity()) {
			// std::cout << "P != ID" << std::endl;
			if (!Qt.isIdentity()) {
				// std::cout << "Qt != ID" << std::endl;
				for (size_t i = 0 ; i < r ; ++i){
					size_t i_r = Qt[i] ;
					size_t j = 0 ;
					for ( ; j < r ; ++j){
						D.setEntry(i,j,A.getEntry(i_r,P[j]));
					}
					for ( ; j < A.coldim() ; ++j){
						G.setEntry(i,j-r,A.getEntry(i_r,P[j]));
					}

				}
			}
			else { // Qt ==id
				// std::cout << "Qt == ID" << std::endl;
				for (size_t i = 0 ; i < r ; ++i){
					size_t j = 0 ;
					for ( ; j < r ; ++j){
						D.setEntry(i,j,A.getEntry(i,P[j]));
					}
					for ( ; j < A.coldim() ; ++j){
						G.setEntry(i,j-r,A.getEntry(i,P[j]));
					}

				}
			}
		}
		else { // P == ID
			// std::cout << "P == ID" << std::endl;
			if (!Qt.isIdentity()) {
				// std::cout << "Qt != ID" << std::endl;
				for (size_t i = 0 ; i < r ; ++i){
					size_t i_r = Qt[i] ;
					size_t j = 0 ;
					for ( ; j < r ; ++j){
						D.setEntry(i,j,A.getEntry(i_r,j));
					}
					for ( ; j < A.coldim() ; ++j){
						G.setEntry(i,j-r,A.getEntry(i_r,j));
					}

				}
			}
			else {
				// std::cout << "Qt == ID" << std::endl;
				std::cerr << "you should not be here" << std::endl;
				for (size_t i = 0 ; i < r ; ++i){
					size_t j = 0 ;
					for ( ; j < r ; ++j){
						D.setEntry(i,j,A.getEntry(i,j));
					}
					for ( ; j < A.coldim() ; ++j){
						G.setEntry(i,j-r,A.getEntry(i,j));
					}

				}
			}
		}

		// std::cout                  << "... matrix extracted." << std::endl;
		// QAP.first.write(std::cout  << "C is : "              << std::endl,false) << ';' << std::endl;
		// QAP.second.write(std::cout << "B is : "              << std::endl,false) << ';' << std::endl;

		return QAP ;
	}



		size_t
		nullspaceRight(
			  BlasMatrix<PID_integer>         & N
			  , const BlasMatrix<PID_integer> & A
			 )
	{
		size_t lA ;
		MaxElementSize(lA,A) ;
		std::vector<Integer> d ;

		size_t m = A.rowdim();
		size_t n = A.coldim();
		size_t r = 0  ; // rank
		size_t s ; // ker dim, the result.

		//typedef double Element ;
		typedef Modular<double> Field;     // dans lequel on réduit A
		typedef Field::Element  Element ;

		PID_integer ZZ ;
		MatrixDomain<PID_integer>        ZMD(ZZ);
		BlasMatrixDomain< PID_integer > BZMD(ZZ) ;

		size_t essai = _LB_MAX_NS_TRY ;

		while (1) { /* main loop */
			BlasPermutation<size_t> P (n);
			BlasPermutation<size_t> Qt(m);
			// ReduceForm RR(A);
			// RR.getDims(r,s);
			// RR.getPermuts(P,Qt);
			// RR.getCB(C,B);

			//!@todo use MatrixHom::map ?
			{ /* LU mod p  */
				RandomPrimeIter PrimeGen(19) ;
				Element p = (Element) PrimeGen.random_between(15);
				Field F(p);
				BlasMatrix<Field> Ap(F,m,n);
				MatrixHom::map(Ap,A,F);
				// Ap.write(std::cout << "A mod " << p << " :=",true) << std::endl;
				LQUPMatrix<Field> LUp(Ap,P,Qt);
				r = LUp.getRank();

			}

			s = n - r; // taille du noyal.
			d.resize(s);


			if (s == 0) { // A était inversible !
				N= BlasMatrix<PID_integer>(ZZ,0,0);
				break;
			}
			else if (r == 0) { // A est nulle ?
				bool nulle = ZMD.isZero(A) ;
				if (!nulle) {
					--essai;
					std::cout << "Bad prime (very bad)" << std::endl;
					if (!essai)
						throw(LinBoxFailure(__func__,__FILE__,__LINE__,
								    "could not find a nullspace..."));
					continue; // A non nulle, pas de bol hein ?
				}
				N = BlasMatrix<PID_integer>(ZZ,n,n);
				Integer one(1);
				for (size_t i = 0 ; i < n ; ++i){
					N.setEntry(i,i,one);
				}
				break ;
			}
			else { // A a un noyau non trivial !

				BlasMatrix<PID_integer> * C = NULL ;
				BlasMatrix<PID_integer> * B = NULL ;

				P.Compress();
				Qt.Compress();

				if (!Qt.isIdentity()) {
					Qt.Invert();
					Qt.Compress();
				}

				if (P.isIdentity() && Qt.isIdentity()) {
					C = (new BlasMatrix<PID_integer>(A,0,0,r,r));
					B = (new BlasMatrix<PID_integer>(A,0,r,r,s));
				}
				else {
					C = (new BlasMatrix<PID_integer>(ZZ,r,r)) ;
					B = (new BlasMatrix<PID_integer>(ZZ,r,s)) ;
				}

				BlasMatrix<PID_integer> NN(ZZ,n,s) ; // le noyal !!

				StackMat<BlasMatrix<PID_integer> >  QAP(*C,*B) ; // A = (C|B\?|?)
				if (!P.isIdentity() || !Qt.isIdentity()) {
					extractMatrix(QAP,A,Qt,P) ;
				}


				Timer Solver ; Solver.clear();
				Solver.start();

#ifdef __LINBOX_USE_IML
				mpz_t * C_iml = reinterpret_cast<mpz_t*>(C->getWritePointer());
				mpz_t * B_iml,*N_iml, D_iml;
				mpz_init(D_iml);
#ifndef __LINBOX_CUTS_HAIR_IN_FOUR
				B_iml = reinterpret_cast<mpz_t*>(B->getWritePointer());
				N_iml = reinterpret_cast<mpz_t*>(NN.getWritePointer()) ;
				nonsingSolvLlhsMM(IML::RightSolu, r, s, C_iml, B_iml, N_iml, D_iml);
				Solver.stop(); std::cout << "total solver : " << Solver << std::endl;
				Integer di ;
				mpz_set(di.get_mpz(),D_iml);
				Integer::negin(di) ;

				for (size_t i = 0 ; i < s ; ++i) {
					NN.setEntry( r+i , i , di) ;
				}
				mpz_clear(D_iml);
#else
				/*  B^t */
				BlasMatrix<PID_integer> Bt(ZZ,s,r) ;
				for (size_t i = 0 ; i < s ; ++i)
					for (size_t j = 0 ; j < r ; ++j)
						Bt.setEntry( i,j,B->getEntry(j,i) );
				BlasMatrix<PID_integer> Nt(ZZ,s,r) ;
				B_iml = reinterpret_cast<mpz_t*>(Bt.getWritePointer());
				N_iml = reinterpret_cast<mpz_t*>(Nt.getWritePointer()) ;
				for (size_t j = 0 ; j < s ; ++j) {
					nonsingSolvLlhsMM(RightSolu, r, 1, C_iml, B_iml, N_iml, D_iml);
					B_iml += r;
					N_iml += r;
					mpz_set(d[j].get_mpz(),D_iml);
				}
				for (size_t i = 0 ; i < r ; ++i)
					for (size_t j = 0 ; j < s ; ++j)
						NN.setEntry( i,j,Nt.getEntry(j,i) );

				for (size_t i = 0 ; i < s ; ++i) {
					Integer::negin(d[i]) ;
					NN.setEntry( r+i , i , d[i]) ;
				}

				Solver.stop(); std::cout << "total solver : " << Solver << std::endl;
#endif
#else
				// d.CxN=B
				RationalSolver<PID_integer,Modular<double>,RandomPrimeIterator,DixonTraits> QD ;
				RationalSolver<PID_integer,Modular<double>,RandomPrimeIterator,WanTraits> QS ;
				// if (lA<50)
					// QS.hasSmallCoeffs(); //!@bug B two...

				BlasMatrix<PID_integer> CC (QAP.first);
#if 0
				Timer Solver_in ; Solver_in.clear();
				Timer Chrono ;
#endif

				bool recontructed = true ;
				for (size_t j = 0 ; j < s ; ++j) {
					std::vector<Integer> v_n(r);
					std::vector<Integer> v_b(r);
					for (size_t i = 0 ; i < r ; ++i)
						v_b[i] = B->getEntry(i,j) ;
					// QQ.solveNonsingular(v_n,d,QAP.first,v_b) ; //! @todo
					enum SolverReturnStatus SRS ;
#if 0
					Chrono.clear(); Chrono.start();
#endif
					if (lA<50) {
						// std::cout << "trying fast solve" << std::endl;
						SRS = QS.solve(v_n,d[j],CC,v_b) ; //! @todo
						if (SRS != SS_OK){
							// std::cout << "failed, falling back to standard" << std::endl;
							SRS = QD.solveNonsingular(v_n,d[j],CC,v_b,false,5) ; //! @todo
						}
					}
					else {
						SRS = QD.solveNonsingular(v_n,d[j],CC,v_b,false,5) ; //! @todo
					}
#if 0
					Chrono.stop(); Solver_in += Chrono ;
#endif
#if 0 /*  debug */
					std::vector<Integer> v_r(r);
					CC.apply(v_r,v_n);
					for (size_t i = 0 ; i < r ; ++i)
						if (v_b[i]*d[j] != v_r[i])
							std::cout << 'X' ;
						else
							std::cout << 'O' ;
					std::cout << '\n' ;
#endif
					if (SRS != SS_OK){
						--essai;
						std::cout << "Recoonstruction failed" << std::endl;
						if (!essai)
							throw(LinBoxFailure(__func__,__FILE__,__LINE__,"could not find a nullspace..."));

						recontructed = false ;
						break;

					}
					for (size_t i = 0 ; i < r ; ++i)
						NN.setEntry(i,j,v_n[i])  ;

				}
#if 0
				std::cout << "only solver : " << Solver_in << std::endl;
				std::cout << "per iter : " << Solver_in.realtime()/s << std::endl;
				Solver.stop(); std::cout << "total solver : " << Solver << std::endl;
#endif
				if (!recontructed){
					continue ;
				}
				delete B ;
				delete C ;


				// on rajoute l'identité
				for (size_t i = 0 ; i < s ; ++i) {
					Integer::negin(d[i]) ;
					NN.setEntry( r+i , i , d[i]) ;
				}
#endif


				/* on applique P^{-1} */
				TransposedBlasMatrix<BlasPermutation<size_t > > PP(P);
				BZMD.mulin_right (PP,NN);


				{/*  Checking with other rows */
					size_t reste = m-r;
					const size_t t = (size_t) std::log((double)reste)+1;
					BlasMatrix<PID_integer>Ax(ZZ,t,n);

					for (int i=(int)t ; i-- ;) {
						//!@todo comment faire sans copie ? créer deux vecteurs 'statiques' et jouer aux pointeurs ?
						int i_q = (int)Qt[(size_t)(r+Integer::random_lessthan(Integer(reste)))] ;
						linbox_check(!(i_q<0) && i_q<(int)m);
						for ( size_t j = 0 ; j < A.coldim() ; ++j)
							Ax.setEntry(i,j, A.getEntry(i_q,j));
					}

					RandomPrimeIter PrimeGen(19) ;
					Element p = (Element) PrimeGen.random_between(15);
					Field F(p);
					BlasMatrix<Field>     Axp(Ax,F);
					BlasMatrix<Field>      Np(NN,F);
					BlasMatrix<Field>       Rp(F,t,s);
					BlasMatrixDomain<Field>   BFMD(F);
					MatrixDomain<Field>        FMD(F);

					BFMD.mul(Rp,Axp,Np);
					if(!FMD.isZero(Rp)){
#if 0 /*  debug */
						std::cout << "p := " << p << ";" << std::endl;
						Axp.write(std::cout << "Axp :=  ",true) << ';' << std::endl;
						Np.write(std::cout << "Np :=  ",true) << ';' << std::endl;
						Rp.write(std::cout << "Rp :=  ",true) << ';' << std::endl;
#endif
						--essai ;
						std::cout << "trying again..." << std::endl;
						if (!essai)
							throw(LinBoxFailure(__func__,__FILE__,__LINE__,"could not find a nullspace..."));

						continue;
					}
				}

				N = NN ;
				break;
			}
		}
		return s;

	}
} // Protected

template<class ZZ>
size_t nullspace (
		BlasMatrix<ZZ>        & N
		,const BlasMatrix<ZZ> & A
		, const Tag::Side Side
		)
{
	if (Side == Tag::Right)
		return Protected::nullspaceRight(N,A);
	else {
		//! XXX throw error
		return -1 ;
	}
}


#if 0
		/*
		 * Calling Sequence:
		 *   nullspaceLong(n, m, A, mp_N_pass)
		 *
		 * Summary:
		 *   Compute the right nullspace of A.
		 *
		 * Input:  n: long, row dimension of A
		 *         m: long, column dimension of A
		 *         A: 1-dim signed long array length n*m, representing n x m matrix
		 *            in row major order
		 *
		 * Output:
		 *   - *mp_N_pass: points to a 1-dim mpz_t array of length m*s, where s is the
		 *                dimension of the right nullspace of A
		 *   - the dimension s of the nullspace is returned
		 *
		 * Notes:
		 *   - The matrix A is represented by one-dimension array in row major order.
		 *   - Space for what mp_N_points to is allocated by this procedure: if the
		 *     nullspace is empty, mp_N_pass is set to NULL.
		 */
		long
		nullspaceLong(const long n, const long m, const long *A, mpz_t * *mp_N_pass)
		{
			long i, j, k, r, s, *P, *rp, *Pt, *rpt, *C, flag, temp;
			double *DA;
			FiniteField p, d;
			mpz_t *mp_B, *mp_N, mp_D, mp_t1, mp_t2;

			P = XCALLOC(long, n + 1);
			rp = XCALLOC(long, n + 1);
			while (1) {
				p = RandPrime(15, 19);
				DA = XCALLOC(double, n * m);
				for (i = 0; i < n * m; i++)
					DA[i] =
					(double) ((temp =
						   (A[i] % ((long) p))) >=
						  0 ? temp : ((long) p) + temp);
				for (i = 0; i < n + 1; i++) {
					P[i] = i;
					rp[i] = 0;
				}
				d = 1;
				RowEchelonTransform(p, DA, n, m, 1, 1, 0, 0, P, rp, &d);
				XFREE(DA);
				r = rp[0];
				s = m - r;
				if (s == 0) {
					*mp_N_pass = NULL;
				} else if (r == 0) {
					flag = 1;
					for (i = 0; i < n * m; i++)
						if (A[i] != 0)
							flag = 0;
					if (!flag)
						continue;
					mp_N = XCALLOC(mpz_t, m * m);
					for (i = 0; i < m; i++) {
						for (j = 0; j < m; j++)
							mpz_init_set_ui(mp_N[i * m + j], 0);
						mpz_init_set_ui(mp_N[i * m + i], 1);
					}
					*mp_N_pass = mp_N;
				} else {		/* r>0 and s>0 */

					Pt = revseq(r, n, P);
					rpt = revseq(r, m, rp);

					C = XCALLOC(long, r * r);
					for (i = 0; i < r; i++)
						for (j = 0; j < r; j++)
							C[i * r + j] = A[Pt[i] * m + rpt[j]];

					mp_B = XCALLOC(mpz_t, r * s);
					for (i = 0; i < r; i++)
						for (j = 0; j < s; j++)
							mpz_init_set_si(mp_B[i * s + j],
									A[Pt[i] * m + rpt[r + j]]);

					mpz_init(mp_D);
					mp_N = XCALLOC(mpz_t, m * s);
					for (i = 0; i < m * s; i++)
						mpz_init(mp_N[i]);

					nonsingSolvMM(RightSolu, r, s, C, mp_B, mp_N, mp_D);
					mpz_neg(mp_D, mp_D);
					for (i = 0; i < s; i++)
						mpz_set(mp_N[(r + i) * s + i], mp_D);

					XFREE(C);
					for (i = 0; i < r * s; i++)
						mpz_clear(mp_B[i]);
					XFREE(mp_B);
					mpz_clear(mp_D);

					for (i = r; i >= 1; i--)
						for (j = 0; j < s; j++)
							mpz_swap(mp_N[(i - 1) * s + j],
								 mp_N[(rp[i] - 1) * s + j]);

					*mp_N_pass = mp_N;

					flag = 1;
					mpz_init(mp_t1);
					mpz_init(mp_t2);
					for (i = r; i < n && flag; i++) {
						for (j = 0; j < s && flag; j++) {
							mpz_set_ui(mp_t2, 0);
							for (k = 0; k < m; k++) {
								mpz_mul_si(mp_t1, mp_N[k * s + j],
									   A[Pt[i] * m + k]);
								mpz_add(mp_t2, mp_t2, mp_t1);
							}
							if (mpz_sgn(mp_t2) != 0)
								flag = 0;
						}
					}
					mpz_clear(mp_t1);
					mpz_clear(mp_t2);

					XFREE(Pt);
					XFREE(rpt);

					if (!flag) {
						for (i = 0; i < m * s; i++)
							mpz_clear(mp_N[i]);
						XFREE(mp_N);
						continue;
					}
				}
				break;
			}
			XFREE(P);
			XFREE(rp);

			return s;

		}


		/*
		 * Calling Sequence:
		 *   nullspaceMP(n, m, A, mp_N_pass)
		 *
		 * Summary: Compute the right nullspace of A. In this function A is a
		 *          1-dimensional mpz_t array.
		 *
		 * Input:  n: long, row dimension of A
		 *         m: long, column dimension of A
		 *         A: 1-dim mpz_t array length n*m, representing n x m matrix
		 *            in row major order
		 *
		 * Output:
		 *   - *mp_N_pass: points to a 1-dim mpz_t array of length m*s, where s is the
		 *                dimension of the right nullspace of A
		 *   - the dimension s of the nullspace is returned
		 *
		 * Notes:
		 *   - The matrix A is represented by one-dimension array in row major order.
		 *   - Space for what mp_N_points to is allocated by this procedure: if the
		 *     nullspace is empty, mp_N_pass is set to NULL.
		 */

		long
		nullspaceMP(const long n, const long m, const mpz_t *A, mpz_t * *mp_N_pass)
		{
			long i, j, k, r, s, *P, *rp, *Pt, *rpt, flag, temp;
			double *DA;
			FiniteField p, d = 1;
			mpz_t *mp_B, *mp_N, mp_D, mp_t1, mp_t2, *C, mp_r;

			mpz_init(mp_r);

			P = XCALLOC(long, n + 1);
			rp = XCALLOC(long, n + 1);
			while (1) {
				p = RandPrime(15, 19);
				DA = XCALLOC(double, n * m);
				for (i = 0; i < n * m; i++) {
					mpz_mod_ui (mp_r, A[i], p);
					DA[i] = mpz_get_d(mp_r);
				}
				for (i = 0; i < n + 1; i++) {
					P[i] = i;
					rp[i] = 0;
				}
				d = 1;
				RowEchelonTransform(p, DA, n, m, 1, 1, 0, 0, P, rp, &d);
				XFREE(DA);
				r = rp[0];
				s = m - r;
				if (s == 0) {
					*mp_N_pass = NULL;
				} else if (r == 0) {
					flag = 1;
					for (i = 0; i < n * m; i++)
						if ( mpz_cmp_si(A[i],0) )
							flag = 0;
					if (!flag)
						continue;
					mp_N = XCALLOC(mpz_t, m * m);
					for (i = 0; i < m; i++) {
						for (j = 0; j < m; j++)
							mpz_init_set_ui(mp_N[i * m + j], 0);
						mpz_init_ui(mp_N[i * m + i], 1);
					}
					*mp_N_pass = mp_N;
				} else {		/* r>0 and s>0 */

					Pt = revseq(r, n, P);
					rpt = revseq(r, m, rp);

					C = XCALLOC(mpz_t, r * r);
					for (i = 0; i < r; i++)
						for (j = 0; j < r; j++)
							mpz_init_set(C[i * r + j], A[Pt[i] * m + rpt[j]]);

					mp_B = XCALLOC(mpz_t, r * s);
					for (i = 0; i < r; i++)
						for (j = 0; j < s; j++)
							mpz_init_set(mp_B[i * s + j], A[Pt[i] * m + rpt[r + j]]);

					mpz_init(mp_D);
					mp_N = XCALLOC(mpz_t, m * s);
					for (i = 0; i < m * s; i++)
						mpz_init(mp_N[i]);


					nonsingSolvLlhsMM(RightSolu, r, s, C, mp_B, mp_N, mp_D);

					mpz_neg(mp_D, mp_D);
					for (i = 0; i < s; i++)
						mpz_set(mp_N[(r + i) * s + i], mp_D);

					for (i = 0; i < r*r; i++)
						mpz_clear(C[i]);
					XFREE(C);

					for (i = 0; i < r * s; i++)
						mpz_clear(mp_B[i]);
					XFREE(mp_B);
					mpz_clear(mp_D);

					for (i = r; i >= 1; i--)
						for (j = 0; j < s; j++)
							mpz_swap(mp_N[(i - 1) * s + j],
								 mp_N[(rp[i] - 1) * s + j]);

					*mp_N_pass = mp_N;

					flag = 1;
					mpz_init(mp_t1);
					mpz_init(mp_t2);
					for (i = r; i < n && flag; i++) {
						for (j = 0; j < s && flag; j++) {
							mpz_set_ui(mp_t2, 0);
							for (k = 0; k < m; k++) {
								mpz_mul(mp_t1, mp_N[k * s + j],  A[Pt[i] * m + k]);
								mpz_add(mp_t2, mp_t2, mp_t1);
							}
							if (mpz_sgn(mp_t2) != 0)
								flag = 0;
						}
					}
					mpz_clear(mp_t1);
					mpz_clear(mp_t2);

					XFREE(Pt);
					XFREE(rpt);

					if (!flag) {
						for (i = 0; i < m * s; i++)
							mpz_clear(mp_N[i]);
						XFREE(mp_N);
						continue;
					}
				}
				break;
			}
			XFREE(P);
			XFREE(rp);

			mpz_clear(mp_r);
			return s;

		}


		/*
		 * Calling Sequence:
		 *   kernelLong(n, m, A, mp_N_pass, reduce)
		 *
		 * Summary:
		 *   Compute the right kernel of A.
		 *   If reduce is true the kernel will be reduced using LLL.
		 *
		 * Input:  n: long, row dimension of A
		 *         m: long, column dimension of A
		 *         A: 1-dim signed long array length n*m, representing n x m matrix
		 *            in row major order
		 *
		 * Output:
		 *   - *mp_N_pass: points to a 1-dim mpz_t array of length m*s, where s is the
		 *                dimension of the right kernel of A
		 *   - the dimension s of the kernel is returned
		 *
		 * Notes:
		 *   - The matrix A is represented by one-dimension array in row major order.
		 *   - Space for what mp_N_points to is allocated by this procedure: if the
		 *     nullspace is empty, mp_N_pass is set to NULL.
		 */
		long
		kernelLong(const long n, const long m, const long *A, mpz_t * *mp_N_pass, const long reduce)
		{
			long i, j, k, r, s, *P, *rp, *Pt, *rpt, *C, flag, temp;
			double *DA;
			FiniteField p, d;
			mpz_t *mp_B, *mp_N, *mp_M, mp_D, mp_t1, mp_t2;

			P = XCALLOC(long, n + 1);
			rp = XCALLOC(long, n + 1);
			while (1) {
				p = RandPrime(15, 19);
				DA = XCALLOC(double, n * m);
				for (i = 0; i < n * m; i++)
					DA[i] =
					(double) ((temp =
						   (A[i] % ((long) p))) >=
						  0 ? temp : ((long) p) + temp);
				for (i = 0; i < n + 1; i++) {
					P[i] = i;
					rp[i] = 0;
				}
				d = 1;
				RowEchelonTransform(p, DA, n, m, 1, 1, 0, 0, P, rp, &d);
				XFREE(DA);
				r = rp[0];
				s = m - r;
				if (s == 0) {
					*mp_N_pass = NULL;
				} else if (r == 0) {
					flag = 1;
					for (i = 0; i < n * m; i++)
						if (A[i] != 0)
							flag = 0;
					if (!flag)
						continue;
					mp_N = XCALLOC(mpz_t, m * m);
					for (i = 0; i < m; i++) {
						for (j = 0; j < m; j++)
							mpz_init_set_ui(mp_N[i * m + j], 0);
						mpz_init_set_ui(mp_N[i * m + i], 1);
					}
					*mp_N_pass = mp_N;
				} else {		/* r>0 and s>0 */

					Pt = revseq(r, n, P);
					rpt = revseq(r, m, rp);

					C = XCALLOC(long, r * r);
					for (i = 0; i < r; i++)
						for (j = 0; j < r; j++)
							C[i * r + j] = A[Pt[i] * m + rpt[j]];

					mp_B = XCALLOC(mpz_t, r * s);
					for (i = 0; i < r; i++)
						for (j = 0; j < s; j++)
							mpz_init_set_si(mp_B[i * s + j],
									A[Pt[i] * m + rpt[r + j]]);

					mpz_init(mp_D);
					mp_N = XCALLOC(mpz_t, m * s);
					for (i = 0; i < m * s; i++)
						mpz_init(mp_N[i]);

					nonsingSolvMM(RightSolu, r, s, C, mp_B, mp_N, mp_D);
					mpz_neg(mp_D, mp_D);
					for (i = 0; i < s; i++)
						mpz_set(mp_N[(r + i) * s + i], mp_D);

					XFREE(C);
					for (i = 0; i < r * s; i++)
						mpz_clear(mp_B[i]);
					XFREE(mp_B);
					mpz_clear(mp_D);

					// Convert nullspace to kernel
					mp_M = XCALLOC(mpz_t, m * s);
					for (i = 0; i < m * s; i++)
						mpz_init_set(mp_M[i], mp_N[i]);

					specialHermite(0, m-s, s, 0, mp_M, NULL);

					kernelBasis(m-s, s, 0, mp_M, mp_N);

					if (reduce==1)  {
						for (i=0; i<m; i++)
							for (j=0; j<s; j++)
								mpz_set(mp_M[j*m+i], mp_N[i*s+j]);

						ired(mp_M, s, m, s);

						for (i=0; i<m; i++)
							for (j=0; j<s; j++)
								mpz_set(mp_N[i*s+j], mp_M[j*m+i]);
					}

					for (i = 0; i < m * s; i++)
						mpz_clear(mp_M[i]);
					XFREE(mp_M);

					for (i = r; i >= 1; i--)
						for (j = 0; j < s; j++)
							mpz_swap(mp_N[(i - 1) * s + j],
								 mp_N[(rp[i] - 1) * s + j]);

					*mp_N_pass = mp_N;

					flag = 1;
					mpz_init(mp_t1);
					mpz_init(mp_t2);
					for (i = r; i < n && flag; i++) {
						for (j = 0; j < s && flag; j++) {
							mpz_set_ui(mp_t2, 0);
							for (k = 0; k < m; k++) {
								mpz_mul_si(mp_t1, mp_N[k * s + j],
									   A[Pt[i] * m + k]);
								mpz_add(mp_t2, mp_t2, mp_t1);
							}
							if (mpz_sgn(mp_t2) != 0)
								flag = 0;
						}
					}
					mpz_clear(mp_t1);
					mpz_clear(mp_t2);

					XFREE(Pt);
					XFREE(rpt);

					if (!flag) {
						for (i = 0; i < m * s; i++)
							mpz_clear(mp_N[i]);
						XFREE(mp_N);
						continue;
					}
				}
				break;
			}
			XFREE(P);
			XFREE(rp);

			return s;

		}
#endif

	} // IML
} // LinBox

#endif // __LINBOX_algorithms_iml_nullspace_INL

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
