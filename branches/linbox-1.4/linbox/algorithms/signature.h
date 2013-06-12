/* Copyright (C)  LinBox
 * Written by Zhendong Wan
 *
 *
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

#include "linbox/field/modular.h"
#include "linbox/algorithms/cra-early-multip.h"
#include <fflas-ffpack/ffpack/ffpack.h>
#include "linbox/randiter/random-prime.h"
#include "linbox/matrix/blas-matrix.h"
#include "linbox/algorithms/blas-domain.h"
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

			typedef typename Matrix::Field::Element Int;
			// std::vector<Int> D(n);
			BlasVector<typename Matrix::Field> D(M.field(),n);
			semiD(D, M);

			//std::cout << "All principal minors are: [";
			//for (int i = 0; i < n; ++ i)
			//	std::cout << D[(size_t)i] << ", ";
			//std::cout << "]\n";

			if (allPos(D)) return true;
			else return false;
		}


		template <class Matrix>
		static bool isPosDef (const Matrix& M, const Minpoly_Method& meth)
		{

			typedef typename Matrix::Field::Element Int;
			typedef std::vector<Int> Poly;
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

			//std::cout << "Begin semiD:\n";
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

			//std::cout << "End semiD:\n";

			if (allPos(D)) return true;
			else return false;
		}

		template <class Matrix>
		static bool isPosSemiDef (const Matrix& M, const Minpoly_Method& meth)
		{

			typedef typename Matrix::Field::Element Int;
			typedef std::vector<Int> Poly;
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

	private:

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

			//std::cout << "Debug begin with input matrix:\n";
			//M. write (std::cout);
			typedef typename Matrix::Field Ring;
			typedef typename Ring::Element Integer_t;
			typedef Modular<double> Field;
			typedef Field::Element Element;

			size_t n = M. rowdim();
	                RandomPrimeIterator primeg;
                	if( ! primeg.template setBitsDelayedField<Field>(n) )
                       		primeg.template setBitsField<Field>();



			Element* FA = new Element[n*n];
			size_t* P= new size_t[n], *PQ = new size_t[n];
			size_t* P_p, * PQ_p;

			Element* p; Element tmp;
			EarlyMultipCRA< Field > cra(3UL);

			Integer_t m = 1;
			// std::vector<Element> v(n);
			typedef UnparametricField<Element> NoField;
			NoField unF ;
			BlasVector<NoField> v(unF,n);
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

				//clog << "Computing blackbox matrix mod " << prime;
				for (p = FA, raw_p = M. Begin(); p != FA + (n*n); ++ p, ++ raw_p)
					K1. init (*p, *raw_p);

				//clog << "\rComputing lup mod " << prime << ". ";
				FFPACK::LUdivine((typename Field::Father_t)K1, FFLAS::FflasNonUnit, FFLAS::FflasNoTrans, n, n, FA, n, P, PQ, FFPACK::FfpackLQUP);

				faithful = true;
				for ( j = 0, P_p = P, PQ_p = PQ; j < n; ++ j, ++ P_p, ++ PQ_p)
					if ((*P_p != j) || (*PQ_p != j))  {
						faithful = false;
						break;
					}

			} while(! faithful);
			K2. init (tmp, 1UL);

			// typename std::vector<Element>::iterator vp;
			typename BlasVector<NoField>::iterator vp;
			for (j = 0, vp = v.begin(); vp != v.end(); ++j, ++vp) {
				K2.mulin(tmp, *(FA + (j * n + j)));
				K2.assign(*vp, tmp);
			}
			// BlasVector<Field> v2(K2,v);//!@bug should not do it like that...
			cra. initialize(K2, v );

			while (! cra.terminated() ){
				// get a prime.
				++primeg; while(cra.noncoprime(*primeg)) ++primeg;
				Field K3(*primeg);
				//clog << "Computing blackbox matrix mod " << prime;
				for (p = FA, raw_p = M. Begin(); p != FA + (n*n); ++ p, ++ raw_p)
					K3. init (*p, *raw_p);

				//clog << "\rComputing lup mod " << prime << ". ";
				FFPACK::LUdivine((typename Field::Father_t)K3, FFLAS::FflasNonUnit, FFLAS::FflasNoTrans, n, n, FA, n, P, PQ, FFPACK::FfpackLQUP);

				faithful = true;
				for ( j = 0, P_p = P, PQ_p = PQ; j < n; ++ j, ++ P_p, ++ PQ_p)
					if ((*P_p != j) || (*PQ_p != j))  {
						faithful = false;
						break;
					}

				if (!faithful) {
					//std::cout << "Not a faithful prime\n";
					continue;
				}

				K3. init (tmp, 1UL);

				for (j = 0, vp = v.begin(); vp != v.end(); ++j, ++vp) {
					K3.mulin(tmp, *(FA + (j * n + j)));
					K3.assign(*vp, tmp);
				}
#if 0
				std::cout << "Faithful image:[";
				for (int l = 0; l < v. size(); ++ l)
					std::cout << v[l] << ", ";
				std::cout << "]\n";
#endif
				// BlasVector<Field> v3(K3,v); //!@bug should not be doing that...
				cra. progress(K3, v);
			}

			delete[] FA;
			delete[] P;
			delete[] PQ;
			//std::cout << "Compute the final answer.\n";
			cra.result(out);
			return out;
		}

		//only works with symmetric integer matrix
		// return a permutation matrix which is represented as a vector, such that
		// all principal of PAP^T are non-zero, up to a maximal.
		template <class Vector, class Matrix>
		static Vector& symmetricLU (Vector& v, const Matrix& IM)
		{

			typedef Modular<int32_t> Field;
			// typedef Modular<double> Field;
			typedef Field::Element Element;
			typedef BlasMatrix<Field> FMatrix;
			RandomPrimeIterator primeg; primeg.template setBitsField<Field>();
			Field F ((unsigned long)*primeg);
			FMatrix FM(F, IM.rowdim(), IM.coldim());
			//std::cout << "Random prime " << p << "\n";

			Element zero; F. init (zero, 0);
			MatrixHom::map (FM, IM, F);
			VectorDomain<Field> VD(F);
			FMatrix& M = FM;

			//typename FMatrix::RowIterator cur_r, tmp_r;
			typedef FMatrix::Row Row;
			//the index is 0-based.
			int i = 0;
			int n = (int) M. rowdim();
			std::vector<int> P((size_t)n);

			for (i = 0; i < n; ++ i)
				P[(size_t)i] = i;

			//M. write(std::cout);
			for (i = 0; i < n; ++ i) {
				//std::cout << "i= " << i << "\n";
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
				//std::cout << "Pivot= " << j << '\n';
				//M. write(std::cout);

				P[(size_t)i] = j;
				Element tmp;
				F. inv (tmp, M[(size_t)i][(size_t)i]);
				F. negin(tmp);
				VD. mulin(*(M. rowBegin() + i), tmp);
				//M. write(std::cout);

				for (j = i + 1; j < n; ++ j) {
					F. assign (tmp,  M[(size_t)j][(size_t)i]);
					VD. axpyin (*(M. rowBegin() + j), tmp,
						    *(M. rowBegin() + i));
				}

				//not necessary
				//M. write(std::cout);
				for (j = i + 1; j < n; ++ j)
					F. assign (M[(size_t)i][(size_t)j], zero);
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

			//std::cout << "Pseud-rank: " << i << "\n[";
			//for (i_p = v. begin(); i_p != v. end(); ++ i_p)
			//	std::cout << *i_p << ", ";
			//std::cout << "]\n";

			return v;
		}

		// This assumes Matrix is BlasMatrix
		// (that it's rawiterator will go thru n^2 values row by row.)
		template <class Matrix>
		static long rank_random (const Matrix& M)
		{

			typedef typename Matrix::Field Ring;
			// typedef typename Ring::Element Integer_t;
			typedef Modular<double> Field;
			typedef Field::Element Element;

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

			long r = (long)FFPACK::Rank((typename Field::Father_t) K, (size_t)n, (size_t)n, FA, (size_t)n);

			delete[] FA;
			return r;

		}
	}; // end of class Signature

} //end of namespace LinBox

#endif //__LINBOX_signature_H

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

