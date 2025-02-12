/* Copyright (C)  LinBox
 * Written by Zhendong Wan
 *			  Jean-Guillaume Dumas
 *
 * ========LICENCE========
 * This file is part of the library LinBox.
 *
  * LinBox is free software: you can redistribute it and/or modify
 * it under the terms of the  GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 */

#ifndef __LINBOX_signature_H
#define __LINBOX_signature_H
/* Function related to the signature computation of symmetric matrices */

#include "linbox/ring/modular.h"
#include "linbox/algorithms/cra-domain.h"
#include "linbox/algorithms/cra-builder-early-multip.h"
#include <fflas-ffpack/ffpack/ffpack.h>
#include "linbox/randiter/random-prime.h"
#include "linbox/matrix/dense-matrix.h"
#include "linbox/matrix/matrix-domain.h"
#include "linbox/polynomial/dense-polynomial.h"
#include "linbox/solutions/minpoly.h"

namespace LinBox
{

	/*! @brief CRA iteration to get a diagonal with the same signature.
	 */
	template <class Matrix, class Field = Givaro::Modular<double>>
	class SemiDIteration {
	private:
		const Matrix& _M;
		size_t _n;

	public:
		using Element = typename Field::Element;

		SemiDIteration(const Matrix& M) :
			_M(M),
			_n(M.rowdim())
		{ }

		template <typename Vector>
		IterationResult operator() (Vector& v, const Field& K) {
			// XXX it would be nice to just create this matrix once
			// and re-use it, but then the iteration wouldn't be thread safe.
			std::vector<Element> FA (_n*_n);
			{
				auto raw_p = _M.Begin();
				for (auto p = FA.begin(); p != FA.end(); ++ p, ++ raw_p)
				{
					K. init (*p, *raw_p);
				}
			}

			if (! FFPACK::fsytrf(K, FFLAS::FflasUpper, _n, FA.data(), _n)) {
				return IterationResult::SKIP;
			}

			v.resize(_n);
			Element tmp;
			K.init(tmp, K.one);
			{
				auto vp = v.begin();
				for (size_t j = 0; vp != v.end(); ++j, ++vp) {
					K.mulin(tmp, FA[j * _n + j]);
					K.assign(*vp, tmp);
				}
			}

			return IterationResult::CONTINUE;
		}
	};

	class Signature {
	public:

		class BLAS_LPM_Method {};
		class Minpoly_Method {};

		template <class Matrix>
		bool isPosDef (const Matrix& M);

		template <class Matrix>
		static bool isPosDef (const Matrix& M, const BLAS_LPM_Method& meth)
            {
                PrimeIterator<IteratorCategories::HeuristicTag>::setSeed(static_cast<uint64_t>(std::time(nullptr)));
                size_t n = M. rowdim();
                std::vector<int> P;
                symmetricLU (P, M);
                if (P. size () < n)
                    return false;

                DenseVector<typename Matrix::Field> D(M.field(),n);
                semiD(D, M);

                    //std::clog << "All principal minors are: [";
                    //for (int i = 0; i < n; ++ i)
                    //	std::clog << D[(size_t)i] << ", ";
                    //std::clog << "]\n";

                if (allPos(D)) return true;
                else return false;
            }


		template <class Matrix>
		static bool isPosDef (const Matrix& M, const Minpoly_Method& meth)
            {

                typedef DensePolynomial<typename Matrix::Field> Poly;
                Poly p(M.field());
                minpoly (p, M);
                typename Poly::reverse_iterator p_p;
                typename Matrix::Field R = M. field();
                bool flip = false;
                for (p_p = p .rbegin(); p_p != p. rend(); ++ p_p) {
                    if (flip)
                        R. negin(*p_p);
                    flip = 1 - flip;
                }

                if(allPos(p)) return true;
                else return false;
            }


		template <class Matrix>
		static bool isPosSemiDef (const Matrix& M, const BLAS_LPM_Method& meth)
            {
                // FIXME this doesn't actually affect the seeding for the algorithms
                PrimeIterator<IteratorCategories::HeuristicTag>::setSeed(static_cast<uint64_t>(std::time(nullptr)));
                size_t n = M. rowdim();
                std::vector<int> P;
                size_t r = (size_t)rank_random (M);
                    //std::clog << "Rank:= " << r << std::endl;
                if (r == 0)
                    return true;
                symmetricLU (P, M);
                if (P. size () < r)
                    return false;

                typename Matrix::Field R = M. field();
                DenseVector<typename Matrix::Field> D(R, P.size());

                    //std::clog << "Begin semiD:\n";
                if(P. size() == n)
                    semiD(D, M);
                else {
                    Matrix PM (R, P.size(), P.size());
                    typename Matrix::RowIterator cur_r; int j = 0;
                    for (cur_r = PM. rowBegin(); cur_r != PM. rowEnd(); ++ cur_r, ++j) {
                        typename Matrix::ConstRowIterator m_r = M. rowBegin() + P[(size_t)j];
                        for (size_t k = 0; k < P.size(); ++ k)
                            R. assign ((*cur_r)[k],
                                       (*m_r)[(size_t)P[(size_t)k]]);
                        // R. assign (cur_r -> operator[] (k),
                        //                m_r -> operator[] ((size_t)P[(size_t)k]));
                    }
                    semiD (D, PM);
                }
                    //std::clog << "End semiD:\n";

                return allPos(D);
            }

		template <class Matrix>
		static bool isPosSemiDef (const Matrix& M, const Minpoly_Method& meth)
            {

                typedef DensePolynomial<typename Matrix::Field> Poly;
                Poly p;
                minpoly (p, M);
                typename Poly::reverse_iterator p_p;
                typename Matrix::Field R = M. field();
                bool flip = false;
                for (p_p = p .rbegin(); p_p != p. rend(); ++ p_p) {
                    if (flip)
                        R. negin(*p_p);
                    flip = 1 - flip;
                }

                if(allNonNeg(p)) return true;
                else return false;
            }

	protected:

		template <class Vector>
		static bool allPos (const Vector& v)
            {

                for (auto p = v. begin(); p != v. end(); ++ p) {
                    if (*p <= 0)
                        return false;
				}

                return true;
            }

		template <class Vector>
		static bool allNonNeg (const Vector& v)
            {

                typename Vector::const_iterator p;
                for (p = v. begin(); p != v. end(); ++ p)
                    if (*p < 0)
                        return false;

                return true;
            }

            /* Compute the equivalent diagonal matrix
             * ie. with the same signature
             * Assume M is non-singular and symmetric with generic rank profile
             */
		template <class Matrix, class Vector>
		static Vector& semiD (Vector& out, const Matrix& M)
            {
                typedef Givaro::Modular<double> Field;

				ChineseRemainder<CRABuilderEarlyMultip<Field>> cra(LINBOX_DEFAULT_EARLY_TERMINATION_THRESHOLD);
				SemiDIteration<Matrix,Field> iter(M);
				PrimeIterator<IteratorCategories::HeuristicTag> primes(FieldTraits<Field>::bestBitSize(M.coldim()));
				return cra(out, iter, primes);
			}

            //only works with symmetric integer matrix
            // return a permutation matrix which is represented as a vector, such that
            // all principal of PAP^T are non-zero, up to a maximal.
		template <class Vector, class Matrix>
		static Vector& symmetricLU (Vector& v, const Matrix& IM)
            {

                typedef Givaro::Modular<int32_t> Field;
                    // typedef Givaro::Modular<double> Field;
                typedef Field::Element Element;
                typedef DenseMatrix<Field> FMatrix;
                PrimeIterator<IteratorCategories::HeuristicTag> primeg(FieldTraits<Field>::bestBitSize(IM.coldim()));
                Field F (*primeg);
                FMatrix FM(F, IM.rowdim(), IM.coldim());
                    //std::clog << "Random prime " << p << "\n";

                MatrixHom::map (FM, IM);
                VectorDomain<Field> VD(F);
                FMatrix& M = FM;

                    //typename FMatrix::RowIterator cur_r, tmp_r;
                    // typedef FMatrix::Row Row;
                    //the index is 0-based.
                int i = 0;
                int n = (int) M. rowdim();
                std::vector<int> P((size_t)n);

                for (i = 0; i < n; ++ i)
                    P[(size_t)i] = i;

                    //M. write(std::clog);
                for (i = 0; i < n; ++ i) {
                        //std::clog << "i= " << i << "\n";
                    int j;
                        //find a pivot
                    for (j = i; j < n; ++ j) {
                        if (!F. isZero(M[(size_t)j][(size_t)j])) break;
                    }

                        //no piviot
                    if (j == n) break;

                        // a pivot
                    if (j != i) {
                        VD. swap (*(M. colBegin() + j), *(M. colBegin() + i));
                        VD. swap (*(M. rowBegin() + j), *(M. rowBegin() + i));
                    }
                        //std::clog << "Pivot= " << j << '\n';
                        //M. write(std::clog);

                    P[(size_t)i] = j;
                    Element tmp;
                    F. inv (tmp, M[(size_t)i][(size_t)i]);
                    F. negin(tmp);
                    VD. mulin(*(M. rowBegin() + i), tmp);
                        //M. write(std::clog);

                    for (j = i + 1; j < n; ++ j) {
                        F. assign (tmp,  M[(size_t)j][(size_t)i]);
                        VD. axpyin (*(M. rowBegin() + j), tmp,
                                    *(M. rowBegin() + i));
                    }

                        //not necessary
                        //M. write(std::clog);
                    for (j = i + 1; j < n; ++ j)
                        F. assign (M[(size_t)i][(size_t)j], F.zero);
                }

                v. resize ((size_t)n);
                std::vector<int>::iterator i_p; int j;
                for (i_p = v. begin(), j = 0; i_p != v. end(); ++ i_p, ++ j)
                    *i_p = j;

                for (j = 0; j < i; ++ j) {
                    if (j != P[(size_t)j])
                        std::swap (v[(size_t)j], v[(size_t)P[(size_t)j]]);
                }

                v. resize ((size_t)i);

                    //std::clog << "Pseud-rank: " << i << "\n[";
                    //for (i_p = v. begin(); i_p != v. end(); ++ i_p)
                    //	std::clog << *i_p << ", ";
                    //std::clog << "]\n";

                return v;
            }

            // This assumes Matrix is BlasMatrix
            // (that it's rawiterator will go thru n^2 values row by row.)
		template <class Matrix>
		static long rank_random (const Matrix& M)
            {

                    // typedef typename Matrix::Field Ring;
                    // typedef typename Ring::Element Integer_t;
                typedef Givaro::Modular<double> Field;
                    // typedef Field::Element Element;

                int n = (int)M. rowdim();

                Field::Element* FA = new Field::Element[n*n], *p;

                    // get a prime.
                    // Compute the rank mod that prime. Accumulate into v with CRA.
                PrimeIterator<IteratorCategories::HeuristicTag> primeg(FieldTraits<Field>::bestBitSize(n));


                Field K(*primeg);

                typename Matrix::ConstIterator raw_p;
                for (p = FA, raw_p = M. Begin(); p != FA + (n*n); ++ p, ++ raw_p)
                    K. init (*p, *raw_p);

                long r = (long)FFPACK::Rank( K, (size_t)n, (size_t)n, FA, (size_t)n);

                delete[] FA;
                return r;

            }
	}; // end of class Signature

} //end of namespace LinBox

#endif //__LINBOX_signature_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
