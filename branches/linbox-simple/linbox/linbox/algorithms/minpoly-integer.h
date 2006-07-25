/* -*- mode: C++; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*- */
/*better filter out repeated primes*/
/* Compute the minpoly of a matrix over an integer ring using modular arithmetic */
/* author: Zhendong Wan*/

#ifndef _LINBOX_MINPOLY_INTEGER_H__
#define _LINBOX_MINPOLY_INTEGER_H__

#include <iostream>
#include <math.h>

#include <linbox/field/field-traits.h>
#include <linbox/algorithms/matrix-hom.h>
#include <linbox/vector/vector-domain.h>
#include <linbox/randiter/random-prime.h>
//#include <linbox/solutions/minpoly.h>
#include <linbox/util/commentator.h>
#include <linbox/ffpack/ffpack.h>
#include <linbox/algorithms/cra-domain.h>

namespace LinBox {

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
		FBlackbox* fbb; FPoly fp;
		integer mmodulus; 
		FieldTraits<Field>::maxModulus(mmodulus);
		long bits = (long) floor (log((double)mmodulus)/M_LN2);
		RandomPrime primeg(bits); integer  prime;
		for (int i = 0; i < n_try; ++ i) {
			primeg. randomPrime (prime);
			Field F(prime);
			MatrixHom::map (fbb, M, F);
			//LinBox::minpoly (fp, *fbb); delete fbb;
			minpoly (fp, *fbb); delete fbb;
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

		RandomPrime primeg(bits); integer prime;
		FBlackbox* fbb; 
		FPoly fp (degree + 1);
		typename FPoly::iterator fp_p;
		y.resize (degree + 1);

                ChineseRemainder< _Field > cra(3UL, degree+1);
		while(! cra.terminated()) {
                    primeg. randomPrime(prime);
                    while(cra.noncoprime(prime))
                        primeg. randomPrime(prime);   
                    Field F(prime);
                    MatrixHom::map (fbb, M, F);
                    minpoly (fp, *fbb); delete fbb;  
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

		RandomPrime primeg(bits); integer prime;
		FBlackbox* fbb; 
		FPoly fp (degree + 1);
		typename FPoly::iterator fp_p;
		y.resize (degree + 1);
                ChineseRemainder< _Field > cra(3UL, degree+1);
		while(! cra.terminated()) {
                    primeg. randomPrime(prime); 
                    while(cra.noncoprime(prime))
                        primeg. randomPrime(prime);   
                    Field F(prime); 
                    MatrixHom::map (fbb, M, F); 
                    minpolySymmetric (fp, *fbb); delete fbb;
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
		RandomPrime primeg(bit1 < bit2 ? bit1 : bit2); integer prime;
		Element* FA = new Element [n*n];
		Element* X = new Element [n*(n+1)];
		size_t* Perm = new size_t[n];
		Element* p;
		typename DenseMatrix<Ring>::ConstRawIterator raw_p;
		std::vector<Element> poly (degree + 1);
		typename std::vector<Element>::iterator poly_ptr;
                ChineseRemainder< _Field > cra(3UL, degree+1);
		while (! cra. terminated()) {
                    primeg. randomPrime(prime);
                    while(cra.noncoprime(prime))
                        primeg. randomPrime(prime);   
                    Field F(prime);
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
		RandomPrime primeg(bit1 < bit2 ? bit1 : bit2); integer prime;
		
		typename DenseMatrix<Ring>::ConstRawIterator raw_p;
		for (int i = 0; i < n_try; ++ i) {
			primeg. randomPrime(prime);
			Field F(prime);
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

#endif
