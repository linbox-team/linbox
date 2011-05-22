/* Copyright (C) LinBox
 *
 * author: Zhendong Wan
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



#ifndef __LINBOX_minpoly_integer_H
#define __LINBOX_minpoly_integer_H

#include <iostream>
#include <math.h>

/*! \file algorithms/minpoly-integer.h
 * Compute the minpoly of a matrix over an integer ring using modular arithmetic 
 * @todo better filter out repeated primes
 */



#include <linbox/field/field-traits.h>
#include <linbox/algorithms/matrix-hom.h>
#include <linbox/vector/vector-domain.h>
#include <linbox/randiter/random-prime.h>
//#include <linbox/solutions/minpoly.h>
#include <linbox/util/commentator.h>
#include <linbox/ffpack/ffpack.h>
#include <linbox/algorithms/cra-early-multip.h>

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
		static Poly& minPolyBlas (Poly& y, const DenseMatrix<Ring>& M);

		template <class Poly, class Ring>
		static Poly& minPolyBlas (Poly& y, const DenseMatrix<Ring>& M, int degree);

		template <class Ring>
		static int minPolyDegreeBlas (const DenseMatrix<Ring>& M, int n_try = 1);
	};

	template<class _Integer, class _Field>
	template<class Poly, class IMatrix>
	Poly& MinPoly<_Integer, _Field>::minPoly(Poly& y, const IMatrix& M) {
		int degree = minPolyDegree (M);
		minPoly(y, M, degree);
		return y;
	}

	template <class _Integer, class _Field>
	template <class IMatrix>
	int MinPoly<_Integer, _Field>::minPolyDegree (const IMatrix& M, int n_try) {
		int degree = 0;
                typedef typename IMatrix::template rebind<Field>::other FBlackbox;
		typedef std::vector<Element> FPoly;
		FPoly fp;
		integer mmodulus; 
		FieldTraits<Field>::maxModulus(mmodulus);
		long bits = (long) floor (log((double)mmodulus)/M_LN2);
		RandomPrimeIterator primeg(bits); 
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
	Poly& MinPoly<_Integer, _Field>::minPoly(Poly& y, const IMatrix& M, int degree) {
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
	Poly& MinPoly<_Integer, _Field>::minPolyNonSymmetric(Poly& y, const IMatrix& M, int degree) {

                typedef typename IMatrix::template rebind<Field>::other FBlackbox;
		typedef std::vector<Element> FPoly;

		integer mmodulus; 
		FieldTraits<Field>::maxModulus(mmodulus);
		long bits = (long) floor (log((double)mmodulus)/M_LN2);

		RandomPrimeIterator primeg(bits); 
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
                        commentator.report (Commentator::LEVEL_IMPORTANT,
                                            INTERNAL_DESCRIPTION) << "Bad prime.\n";
                        continue;
                    }
                    cra.progress(F, fp);
		}
                
		cra. result (y);
                    // commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION) <<  "Number of primes needed: " << cra. steps() << std::endl;
		return y;
	}

	template <class _Integer, class _Field>
	template<class Poly, class IMatrix>
	Poly& MinPoly<_Integer, _Field>::minPolySymmetric(Poly& y, const IMatrix& M, int degree) {

                typedef typename IMatrix::template rebind<Field>::other FBlackbox;
		typedef std::vector<Element> FPoly;


		integer mmodulus; 
		FieldTraits<Field>::maxModulus(mmodulus);
		long bits = (long) floor (log((double)mmodulus)/M_LN2);

		RandomPrimeIterator primeg(bits); 
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
                        commentator.report (Commentator::LEVEL_IMPORTANT,
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
	bool MinPoly<_Integer, _Field>::isSymmetric(const IMatrix& M, int n_try) {
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
	Poly& MinPolyBlas<_Integer, _Field>::minPolyBlas (Poly& y, const DenseMatrix<Ring>& M) {
		int degree = minPolyDegreeBlas (M);
		minPolyBlas (y, M, degree);
		return y;
	}

	template <class _Integer, class _Field>
	template <class Poly, class Ring>
	Poly& MinPolyBlas<_Integer, _Field>::minPolyBlas (Poly& y, const DenseMatrix<Ring>& M, int degree) {

		y. resize (degree + 1);
		size_t n = M. rowdim();
		integer mmodulus; 
		FieldTraits<Field>::maxModulus(mmodulus);
		long bit1 = (long) floor (log((double)mmodulus)/M_LN2);
		long bit2 = (long) floor (log(sqrt(double(4503599627370496LL/n)))/M_LN2);
		RandomPrimeIterator primeg(bit1 < bit2 ? bit1 : bit2);
		Element* FA = new Element [n*n];
		Element* X = new Element [n*(n+1)];
		size_t* Perm = new size_t[n];
		Element* p;
		typename DenseMatrix<Ring>::ConstRawIterator raw_p;
		std::vector<Element> poly (degree + 1);
		typename std::vector<Element>::iterator poly_ptr;

                EarlyMultipCRA< _Field > cra(3UL);
                do {
                    ++primeg; while(cra.noncoprime(*primeg)) ++primeg;   
                    Field F(*primeg);
                    for (p = FA, raw_p = M. rawBegin(); 
                         p != FA + (n*n); ++ p, ++ raw_p)
                        
                        F. init (*p, *raw_p);
                    
                    FFPACK::MinPoly( F, poly, n, FA, n, X, n, Perm);

                    cra.initialize(F, poly);
                } while( poly. size() != degree + 1) ; // Test for Bad primes

		while (! cra. terminated()) {
                    ++primeg; while(cra.noncoprime(*primeg)) ++primeg;   
                    Field F(*primeg);
                    for (p = FA, raw_p = M. rawBegin(); 
                         p != FA + (n*n); ++ p, ++ raw_p)
                        
                        F. init (*p, *raw_p);
                    
                    FFPACK::MinPoly( F, poly, n, FA, n, X, n, Perm);
                    
                    if(poly. size() != degree + 1) {
                        commentator.report (Commentator::LEVEL_IMPORTANT, 
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
	int MinPolyBlas<_Integer, _Field>::minPolyDegreeBlas (const DenseMatrix<Ring>& M, int n_try) {
		size_t n = M. rowdim();
		int degree = 0;
		Element* FA = new Element [n*n];
		Element* X = new Element [n*(n+1)];
		size_t* Perm = new size_t[n];
		Element* p;
		std::vector<Element> Poly;

		integer mmodulus; 
		FieldTraits<Field>::maxModulus(mmodulus);
		long bit1 = (long) floor (log((double)mmodulus)/M_LN2);
		long bit2 = (long) floor (log(sqrt(double(4503599627370496LL/n)))/M_LN2);
		RandomPrimeIterator primeg(bit1 < bit2 ? bit1 : bit2); 
		
		typename DenseMatrix<Ring>::ConstRawIterator raw_p;
		for (int i = 0; i < n_try; ++ i) {
			++primeg;
			Field F(*primeg);
			for (p = FA, raw_p = M. rawBegin(); 
				 p!= FA + (n*n); ++ p, ++ raw_p)
				F. init (*p, *raw_p);

			FFPACK::MinPoly( F, Poly, n, FA, n, X, n, Perm);

			if (degree < Poly. size() - 1) 
				degree = Poly. size() -1;
		}
		delete FA; delete X; delete Perm;

		return degree;
	}
} // LinBox

#endif //__LINBOX_minpoly_integer_H
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
