/* Copyright (C) LinBox
 *
 * author: Zhendong Wan
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



#ifndef __LINBOX_minpoly_integer_H
#define __LINBOX_minpoly_integer_H

#include <iostream>
#include <math.h>

/*! \file algorithms/minpoly-integer.h
 * Compute the minpoly of a matrix over an integer ring using modular arithmetic
 * @todo better filter out repeated primes
 */



#include "linbox/field/field-traits.h"
#include "linbox/algorithms/matrix-hom.h"
#include "linbox/vector/vector-domain.h"
#include "linbox/randiter/random-prime.h"
//#include "linbox/solutions/minpoly.h"
#include "linbox/util/commentator.h"
#include <fflas-ffpack/ffpack/ffpack.h>
#include "linbox/algorithms/cra-early-multip.h"

namespace LinBox
{

	/* compute the minpoly of a matrix over the Integer ring
	 * via modular method over Field.
	 */
	template <class _Integer, class _Field>
	class MinPoly {
	public:
		typedef _Field Field;
		//typedef _Integer Integer;
		typedef typename Field::Element Element;

		/* Blackbox case */
		template<class Poly, class IMatrix>
		static Poly& minPoly(Poly& y, const IMatrix& M);

		template <class IMatrix>
		static int minPolyDegree (const IMatrix& M, int n_try = 1);

		template<class Poly, class IMatrix>
		static Poly& minPoly(Poly& y, const IMatrix& M, int degree);

		template<class Poly, class IMatrix>
		static Poly& minPolyNonSymmetric(Poly& y, const IMatrix& M, int degree);

		template<class Poly, class IMatrix>
		static Poly& minPolySymmetric(Poly& y, const IMatrix& M, int degree);

		template <class IMatrix>
		static bool isSymmetric(const IMatrix& M, int n_try = 1);
	};

	template <class _Integer, class _Field>
	class MinPolyBlas {
	public:
		typedef _Field Field;
		typedef _Integer Integer;
		typedef typename Field::Element Element;

		template <class Poly, class Ring>
		static Poly& minPolyBlas (Poly& y, const BlasMatrix<Ring>& M);

		template <class Poly, class Ring>
		static Poly& minPolyBlas (Poly& y, const BlasMatrix<Ring>& M, int degree);

		template <class Ring>
		static int minPolyDegreeBlas (const BlasMatrix<Ring>& M, int n_try = 1);
	};

	template<class _Integer, class _Field>
	template<class Poly, class IMatrix>
	Poly& MinPoly<_Integer, _Field>::minPoly(Poly& y, const IMatrix& M)
	{
		int degree = minPolyDegree (M);
		minPoly(y, M, degree);
		return y;
	}

	template <class _Integer, class _Field>
	template <class IMatrix>
	int MinPoly<_Integer, _Field>::minPolyDegree (const IMatrix& M, int n_try)
	{
		int degree = 0;
		typedef typename IMatrix::template rebind<Field>::other FBlackbox;
		typedef std::vector<Element> FPoly;
		FPoly fp;
		RandomPrimeIterator primeg; primeg.template setBitsField<Field>();
		for (int i = 0; i < n_try; ++ i) {
			++primeg;
			Field F(*primeg);
			FBlackbox  fbb(M, F);
			minpoly (fp, fbb);
			if (degree < ((int) fp.size() - 1)) degree = fp.size() -1;
		}
		return degree;
	}

	template <class _Integer, class _Field>
	template<class Poly, class IMatrix>
	Poly& MinPoly<_Integer, _Field>::minPoly(Poly& y, const IMatrix& M, int degree)
	{
		if (isSymmetric(M))  {
			//std::cout << "Symmetric:\n";
			minPolySymmetric(y, M, degree);
		}
		else {
			//std::cout << "NonSymmetric:\n";
			minPolyNonSymmetric(y, M, degree);
		}
		return y;
	}

	template <class _Integer, class _Field>
	template<class Poly, class IMatrix>
	Poly& MinPoly<_Integer, _Field>::minPolyNonSymmetric(Poly& y, const IMatrix& M, int degree)
	{

		typedef typename IMatrix::template rebind<Field>::other FBlackbox;
		typedef std::vector<Element> FPoly;

		RandomPrimeIterator primeg; primeg.template setBitsField<Field>();

		FPoly fp (degree + 1);
		typename FPoly::iterator fp_p;
		y.resize (degree + 1);

		EarlyMultipCRA< _Field > cra(3UL);
		do {
			++primeg;
			Field F(*primeg);
			FBlackbox fbb(M, F);
			minpoly (fp, fbb);
			cra.initialize(F, fp);
		} while( (int)fp.size() - 1 != degree); // Test for Bad primes

		while(! cra.terminated()) {
			++primeg; while(cra.noncoprime(*primeg)) ++primeg;
			Field F(*primeg);
			FBlackbox fbb(M, F);
			minpoly (fp, fbb);
			if ((int)fp.size() - 1 != degree) {
				commentator().report (Commentator::LEVEL_IMPORTANT,
						    INTERNAL_DESCRIPTION) << "Bad prime.\n";
				continue;
			}
			cra.progress(F, fp);
		}

		cra. result (y);
		// commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION) <<  "Number of primes needed: " << cra. steps() << std::endl;
		return y;
	}

	template <class _Integer, class _Field>
	template<class Poly, class IMatrix>
	Poly& MinPoly<_Integer, _Field>::minPolySymmetric(Poly& y, const IMatrix& M, int degree)
	{

		typedef typename IMatrix::template rebind<Field>::other FBlackbox;
		typedef std::vector<Element> FPoly;

		RandomPrimeIterator primeg; primeg.template setBitsField<Field>();

		FPoly fp (degree + 1);
		typename FPoly::iterator fp_p;
		y.resize (degree + 1);

		EarlyMultipCRA< _Field > cra(3UL);
		do {
			++primeg;
			Field F(*primeg);
			FBlackbox fbb(M,F);
			minpolySymmetric (fp, fbb);
			cra.initialize(F, fp);
		} while( (int)fp.size() - 1 != degree); // Test for Bad primes

		while(! cra.terminated()) {
			++primeg; while(cra.noncoprime(*primeg)) ++primeg;
			Field F(*primeg);
			FBlackbox fbb(M,F);
			minpolySymmetric (fp, fbb);
			if ((int)fp.size() - 1 != degree) {
				commentator().report (Commentator::LEVEL_IMPORTANT,
						    INTERNAL_DESCRIPTION) << "Bad prime.\n";
				continue;
			}
			cra.progress(F, fp);
		}
		cra. result (y);
		//std::cout << "Number of primes needed: " << cra. steps() << std::endl;
		return y;
	}


	template <class _Integer, class _Field>
	template <class IMatrix>
	bool MinPoly<_Integer, _Field>::isSymmetric(const IMatrix& M, int n_try)
	{
		typedef typename IMatrix::Field Ring;
		typedef typename Ring::Element Element;
		Ring R(M. field()); int order = M. rowdim();
		std::vector<Element> x(order), mx (order), xm (order);
		typename std::vector<Element>::iterator x_p;
		VectorDomain<Ring> RVD (R);
		if (M. rowdim() != M. coldim()) return false;

		for (int i = 0; i < n_try; ++ i) {
			for (x_p = x. begin(); x_p != x. end(); ++ x_p)
				R. init (*x_p, rand());
			M. apply (mx, x);
			M. applyTranspose (xm, x);
			if (!RVD.areEqual(mx, xm)) return false;
		}
		return true;
	}

	template <class _Integer, class _Field>
	template <class Poly, class Ring>
	Poly& MinPolyBlas<_Integer, _Field>::minPolyBlas (Poly& y, const BlasMatrix<Ring>& M)
	{
		int degree = minPolyDegreeBlas (M);
		minPolyBlas (y, M, degree);
		return y;
	}

	template <class _Integer, class _Field>
	template <class Poly, class Ring>
	Poly& MinPolyBlas<_Integer, _Field>::minPolyBlas (Poly& y, const BlasMatrix<Ring>& M, int degree)
	{

		y. resize (degree + 1);
		size_t n = M. rowdim();
		RandomPrimeIterator primeg;
		if( ! primeg.template setBitsDelayedField<Field>(n) )
			primeg.template setBitsField<Field>();
		Element* FA = new Element [n*n];
		Element* X = new Element [n*(n+1)];
		size_t* Perm = new size_t[n];
		Element* p;
		typename BlasMatrix<Ring>::ConstIterator raw_p;
		std::vector<Element> poly (degree + 1);
		typename std::vector<Element>::iterator poly_ptr;

		EarlyMultipCRA< _Field > cra(3UL);
		do {
			++primeg; while(cra.noncoprime(*primeg)) ++primeg;
			Field F(*primeg);
			for (p = FA, raw_p = M. Begin();
			     p != FA + (n*n); ++ p, ++ raw_p)

				F. init (*p, *raw_p);

			FFPACK::MinPoly((typename _Field::Father_t) F, poly, n, FA, n, X, n, Perm);

			cra.initialize(F, poly);
		} while( poly. size() != degree + 1) ; // Test for Bad primes

		while (! cra. terminated()) {
			++primeg; while(cra.noncoprime(*primeg)) ++primeg;
			Field F(*primeg);
			for (p = FA, raw_p = M. Begin();
			     p != FA + (n*n); ++ p, ++ raw_p)

				F. init (*p, *raw_p);

			FFPACK::MinPoly((typename _Field::Father_t) F, poly, n, FA, n, X, n, Perm);

			if(poly. size() != degree + 1) {
				commentator().report (Commentator::LEVEL_IMPORTANT,
						    INTERNAL_DESCRIPTION) << "Bad prime.\n";
				continue;
			}
			cra.progress(F, poly);
		}
		cra. result(y);
		//std::cout << "Number of primes needed: " << cra. steps() << std::endl;
		delete FA; delete X; delete Perm;

		return y;
	}


	template <class _Integer, class _Field>
	template <class Ring>
	int MinPolyBlas<_Integer, _Field>::minPolyDegreeBlas (const BlasMatrix<Ring>& M, int n_try)
	{
		size_t n = M. rowdim();
		int degree = 0;
		Element* FA = new Element [n*n];
		Element* X = new Element [n*(n+1)];
		size_t* Perm = new size_t[n];
		Element* p;
		std::vector<Element> Poly;

                RandomPrimeIterator primeg;
                if( ! primeg.template setBitsDelayedField<Field>(n) )
                        primeg.template setBitsField<Field>();

		typename BlasMatrix<Ring>::ConstIterator raw_p;
		for (int i = 0; i < n_try; ++ i) {
			++primeg;
			Field F(*primeg);
			for (p = FA, raw_p = M. Begin();
			     p!= FA + (n*n); ++ p, ++ raw_p)
				F. init (*p, *raw_p);

			FFPACK::MinPoly((typename _Field::Father_t) F, Poly, n, FA, n, X, n, Perm);

			if (degree < Poly. size() - 1)
				degree = Poly. size() -1;
		}
		delete FA; delete X; delete Perm;

		return degree;
	}
} // LinBox

#endif //__LINBOX_minpoly_integer_H

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

