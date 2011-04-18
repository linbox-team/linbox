/* Copyright (C)  LinBox
 * Written by Zhendong Wan 
 *
 *
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#ifndef __LINBOX_signature_H
#define __LINBOX_signature_H
/* Function related to the signature computation of symmetric matrices */

#include <linbox/field/modular-double.h>
#include <linbox/field/modular-int32.h>
#include <linbox/algorithms/cra-early-multip.h>
#include <linbox/ffpack/ffpack.h>
#include <linbox/randiter/random-prime.h>
#include <linbox/matrix/blas-matrix.h>
#include <linbox/algorithms/blas-domain.h>
#include <linbox/solutions/minpoly.h>

namespace LinBox 
{

class Signature {
public:

	class BLAS_LPM_Method {};
	class Minpoly_Method {};
	
	template <class Matrix>
	bool isPosDef (const Matrix& M);

	template <class Matrix>
	static bool isPosDef (const Matrix& M, const BLAS_LPM_Method& meth) {
		RandomPrimeIterator::setSeed(time(0));
		size_t n = M. rowdim();
		std::vector<int> P;
		symmetricLU (P, M);
		if (P. size () < n) 
			return false;

		typedef typename Matrix::Field::Element Int;
		std::vector<Int> D(n);
		semiD(D, M);

		//std::cout << "All principal minors are: [";
		//for (int i = 0; i < n; ++ i)
		//	std::cout << D[i] << ", ";
		//std::cout << "]\n";

		if (allPos(D)) return true;
		else return false;
	}	
	

	template <class Matrix>
	static bool isPosDef (const Matrix& M, const Minpoly_Method& meth) {

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
	static bool isPosSemiDef (const Matrix& M, const BLAS_LPM_Method& meth) {
		RandomPrimeIterator::setSeed(time(0));
		size_t n = M. rowdim();
		std::vector<int> P;
		size_t r = rank_random (M);
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
				typename Matrix::ConstRowIterator m_r = M. rowBegin() + P[j];
				for (size_t k = 0; k < P.size(); ++ k) 
					R. assign (cur_r -> operator[] (k),
								m_r -> operator[] (P[k]));
			}
			semiD (D, PM);
		}

		//std::cout << "End semiD:\n";

		if (allPos(D)) return true;
		else return false;
	}	

	template <class Matrix>
	static bool isPosSemiDef (const Matrix& M, const Minpoly_Method& meth) {

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
	static bool allPos (const Vector& v) {
		
		typename Vector::const_iterator p;
		for (p = v. begin(); p != v. end(); ++ p)
			if (*p <= 0)
				return false;

		return true;
	}

	template <class Vector>
	static bool allNonNeg (const Vector& v) {
		
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
	static Vector& semiD (Vector& out, const Matrix& M) {
		
		//std::cout << "Debug begin with input matrix:\n";
		//M. write (std::cout);
		typedef typename Matrix::Field Ring;
		typedef typename Ring::Element Integer;
		typedef Modular<double> Field;
		typedef Field::Element Element;

		size_t n = M. rowdim();
			
		integer mmodulus;
		FieldTraits<Field>::maxModulus(mmodulus);
		long bit1 = (long) floor (log((double)mmodulus)/M_LN2);
		long bit2 = (long) floor (log(sqrt(double(4503599627370496LL/n)))/M_LN2);
		RandomPrimeIterator primeg(bit1 < bit2 ? bit1 : bit2); 

		Field::Element* FA = new Field::Element[n*n];
		size_t* P= new size_t[n], *PQ = new size_t[n];
		size_t* P_p, * PQ_p;

		Field::Element* p; Field::Element tmp;
		EarlyMultipCRA< Field > cra(3UL);

		Integer m = 1;
		std::vector<Field::Element> v(n);
		size_t j = 0;
		Field K2;
		bool faithful = true;
		typename Matrix::ConstRawIterator raw_p;

		do {
		// get a prime. 
		// Compute mod that prime. Accumulate into v with CRA. 
		++primeg ; while(cra.noncoprime(*primeg)) ++primeg; 
		Field K1(*primeg); 
		K2 = K1;

		//clog << "Computing blackbox matrix mod " << prime;
		for (p = FA, raw_p = M. rawBegin(); p != FA + (n*n); ++ p, ++ raw_p)
		  K1. init (*p, *raw_p);

		//clog << "\rComputing lup mod " << prime << ". ";
		FFPACK::LUdivine(K1, FFLAS::FflasNonUnit, FFLAS::FflasNoTrans, n, n, FA, n, P, PQ, FFPACK::FfpackLQUP);

		faithful = true;
		for ( j = 0, P_p = P, PQ_p = PQ; j < n; ++ j, ++ P_p, ++ PQ_p) 
		  if ((*P_p != j) || (*PQ_p != j))  {
		    faithful = false;
		    break;
		  }

		} while(! faithful);
		K2. init (tmp, 1UL);

		typename std::vector<Field::Element>::iterator vp;
		for (j = 0, vp = v.begin(); vp != v.end(); ++j, ++vp) {
		  K2.mulin(tmp, *(FA + (j * n + j)));
		  K2.assign(*vp, tmp);
		}
		cra. initialize(K2, v); 

		while (! cra.terminated() ){
			// get a prime. 
		   ++primeg; while(cra.noncoprime(*primeg)) ++primeg; 
		  Field K(*primeg); 
		  //clog << "Computing blackbox matrix mod " << prime;
		  for (p = FA, raw_p = M. rawBegin(); p != FA + (n*n); ++ p, ++ raw_p)
		    K. init (*p, *raw_p);
		  
		  //clog << "\rComputing lup mod " << prime << ". ";
		  FFPACK::LUdivine(K, FFLAS::FflasNonUnit, FFLAS::FflasNoTrans, n, n, FA, n, P, PQ, FFPACK::FfpackLQUP);
		  
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
		  
		  K. init (tmp, 1UL);

		  for (j = 0, vp = v.begin(); vp != v.end(); ++j, ++vp) {
		    K.mulin(tmp, *(FA + (j * n + j)));
		    K.assign(*vp, tmp);
		  }
// 		  std::cout << "Faithful image:[";
// 		  for (int l = 0; l < v. size(); ++ l)
// 		  std::cout << v[l] << ", ";
// 		  std::cout << "]\n";
		  cra. progress(K, v); 
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
	static Vector& symmetricLU (Vector& v, const Matrix& IM) {

		typedef Modular<int32> Field;
		typedef Field::Element Element;
		typedef DenseMatrix<Field> FMatrix;
		RandomPrimeIterator primeg(20);
		Field F (*primeg);
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
		int n = M. rowdim();
		std::vector<int> P(n);

		for (i = 0; i < n; ++ i)
			P[i] = i;
		
		//M. write(std::cout);
		for (i = 0; i < n; ++ i) {
			//std::cout << "i= " << i << "\n";
			int j;
			//find a pivot
			for (j = i; j < n; ++ j) {
				if (!F. isZero(M[j][j])) break;
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

			P[i] = j;
			Element tmp;
			F. inv (tmp, M[i][i]);
			F. negin(tmp);
			VD. mulin(*(M. rowBegin() + i), tmp);
			//M. write(std::cout);
			
			for (j = i + 1; j < n; ++ j) {
				F. assign (tmp,  M[j][i]);
				VD. axpyin (*(M. rowBegin() + j), tmp,
							*(M. rowBegin() + i));
			}

			//not necessary
			//M. write(std::cout);
			for (j = i + 1; j < n; ++ j)
				F. assign (M[i][j], zero);
		}

		v. resize (n);
		std::vector<int>::iterator i_p; int j;
		for (i_p = v. begin(), j = 0; i_p != v. end(); ++ i_p, ++ j)
			*i_p = j;

		for (j = 0; j < i; ++ j) {
			if (j != P[j])
				std::swap (v[j], v[P[j]]);
		}

		v. resize (i);

		//std::cout << "Pseud-rank: " << i << "\n[";
		//for (i_p = v. begin(); i_p != v. end(); ++ i_p)
		//	std::cout << *i_p << ", ";
		//std::cout << "]\n";

		return v;
	}

	// This assumes Matrix is DenseMatrix 
	// (that it's rawiterator will go thru n^2 values row by row.)
	template <class Matrix>
	static long rank_random (const Matrix& M) {

		typedef typename Matrix::Field Ring;
		typedef typename Ring::Element Integer;
		typedef Modular<double> Field;
		typedef Field::Element Element;

		int n = M. rowdim();
			
		integer mmodulus;
		FieldTraits<Field>::maxModulus(mmodulus);
		long bit1 = (long) floor (log((double)mmodulus)/M_LN2);
		long bit2 = (long) floor (log(sqrt(double(4503599627370496LL/n)))/M_LN2);
		RandomPrimeIterator primeg(bit1 < bit2 ? bit1 : bit2); 

		Field::Element* FA = new Field::Element[n*n], *p;

		// get a prime. 
		// Compute the rank mod that prime. Accumulate into v with CRA. 
		Field K(*primeg); 

		typename Matrix::ConstRawIterator raw_p;
		for (p = FA, raw_p = M. rawBegin(); p != FA + (n*n); ++ p, ++ raw_p)
			K. init (*p, *raw_p);

		long r = FFPACK::Rank( K, n, n, FA, n);

		delete[] FA;
		return r;

	}
}; // end of class Signature

} //end of namespace LinBox

#endif //__LINBOX_signature_H
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
