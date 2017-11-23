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
#include "linbox/algorithms/cra-early-multip.h"
#include <fflas-ffpack/ffpack/ffpack.h>
#include "linbox/randiter/random-prime.h"
#include "linbox/matrix/dense-matrix.h"
#include "linbox/matrix/matrix-domain.h"
#include "linbox/polynomial/dense-polynomial.h"
#include "linbox/solutions/minpoly.h"

namespace LinBox
{

	class Signature {
	public:

		class BLAS_LPM_Method {};
		class Minpoly_Method {};

		template <class Matrix>
		bool isPosDef (const Matrix& M);

		template <class Matrix>
		static bool isPosDef (const Matrix& M, const BLAS_LPM_Method& meth)
            {
                RandomPrimeIterator::setSeed((size_t)time(0));
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
                RandomPrimeIterator::setSeed((size_t)time(0));
                size_t n = M. rowdim();
                std::vector<int> P;
                size_t r = (size_t)rank_random (M);
                    //std::clog << "Rank:= " << r << std::endl;
                if (r == 0)
                    return true;
                symmetricLU (P, M);
                if (P. size () < r)
                    return false;

                typedef typename Matrix::Field::Element Int;
                std::vector<Int> D(P.size());

                typename Matrix::Field R = M. field();

                    //std::clog << "Begin semiD:\n";
                if(P. size() == n)
                    semiD(D, M);
                else {
                    Matrix PM (R, P.size(), P.size());
                    typename Matrix::RowIterator cur_r; int j = 0;
                    for (cur_r = PM. rowBegin(); cur_r != PM. rowEnd(); ++ cur_r, ++j) {
                        typename Matrix::ConstRowIterator m_r = M. rowBegin() + P[(size_t)j];
                        for (size_t k = 0; k < P.size(); ++ k)
                            R. assign (cur_r -> operator[] (k),
                                       m_r -> operator[] ((size_t)P[(size_t)k]));
                    }
                    semiD (D, PM);
                }

                    //std::clog << "End semiD:\n";

                if (allPos(D)) return true;
                else return false;
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

                typename Vector::const_iterator p;
                for (p = v. begin(); p != v. end(); ++ p)
                    if (*p <= 0)
                        return false;

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

#ifdef _LB_DEBUG
                std::clog << "Debug begin with input matrix:\n";
                M. write (std::clog);
#endif
                typedef typename Matrix::Field Ring;
                typedef typename Ring::Element Integer_t;
                typedef Givaro::Modular<double> Field;
                typedef Field::Element Element;

                size_t n = M. rowdim();
                RandomPrimeIterator primeg;
                if( ! primeg.template setBitsDelayedField<Field>(n) )
                    primeg.template setBitsField<Field>();



                Element* FA = new Element[n*n];
                size_t* P= new size_t[n], *PQ = new size_t[n];

                Element* p; Element tmp;
                EarlyMultipCRA< Field > cra(3UL);

                Integer_t m = 1;
                    // std::vector<Element> v(n);
                typedef Givaro::ZRing<Element> NoField;
                NoField unF ;
                DenseVector<NoField> v(unF,n);
                size_t j = 0;
                Field K2;
                bool faithful = true;
                typename Matrix::ConstIterator raw_p;

                do {
                        // get a prime.
                        // Compute mod that prime. Accumulate into v with CRA.
                    ++primeg ; while(cra.noncoprime(*primeg)) ++primeg;
                    Field K1(*primeg);
                    K2 = K1;

#ifdef _LB_DEBUG
                    std::clog << "Computing blackbox matrix mod " << *primeg;
#endif
                    for (p = FA, raw_p = M. Begin(); p != FA + (n*n); ++ p, ++ raw_p)
                        K1. init (*p, *raw_p);

                    faithful = FFPACK::fsytrf(K1, FFLAS::FflasUpper, n, FA, n);

#ifdef _LB_DEBUG
                    if (!faithful) {
                        std::clog << "Not a faithful prime\n";
                    }
#endif

                } while(! faithful);
                K2. assign(tmp, K2.one);

                typename DenseVector<NoField>::iterator vp;
                for (j = 0, vp = v.begin(); vp != v.end(); ++j, ++vp) {
                    K2.mulin(tmp, *(FA + (j * n + j)));
                    K2.assign(*vp, tmp);
                }
                cra. initialize(K2, v );

                while (! cra.terminated() ){
                        // get a prime.
                    ++primeg; while(cra.noncoprime(*primeg)) ++primeg;
                    Field K3(*primeg);
#ifdef _LB_DEBUG
                    std::clog << "Computing blackbox matrix mod " << *primeg;
#endif
                    for (p = FA, raw_p = M. Begin(); p != FA + (n*n); ++ p, ++ raw_p)
                        K3. init (*p, *raw_p);

#ifdef _LB_DEBUG
                    std::clog << "\rComputing lup mod " << *primeg << ". ";
#endif
                    faithful = FFPACK::fsytrf(K3, FFLAS::FflasUpper, n, FA, n);
                    if (!faithful) {
#ifdef _LB_DEBUG
                        std::clog << "Not a faithful prime\n";
#endif
                        continue;
                    }

                    K3. assign(tmp, K3.one);

                    for (j = 0, vp = v.begin(); vp != v.end(); ++j, ++vp) {
                        K3.mulin(tmp, *(FA + (j * n + j)));
                        K3.assign(*vp, tmp);
                    }
#ifdef _LB_DEBUG
                    std::clog << "Faithful image:[";
                    for (size_t l = 0; l < v. size(); ++ l)
                        std::clog << v[l] << ", ";
                    std::clog << "]\n";
#endif
                    cra. progress(K3, v);
                }

                delete[] FA;
                delete[] P;
                delete[] PQ;
                    //std::clog << "Compute the final answer.\n";
                cra.result(out);

#ifdef _LB_CRATIMING
                cra.reportTimes(std::clog);
#endif
                return out;
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
                RandomPrimeIterator primeg; primeg.template setBitsField<Field>();
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
                RandomPrimeIterator primeg;
                if( ! primeg.template setBitsDelayedField<Field>(n) )
                    primeg.template setBitsField<Field>();

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

// vim:sts=4:sw=4:ts=4:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:

