/* -*- mode: C++; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*- */
/*better filter out repeated primes*/
#include <iostream>
#include <math.h>

#include <linbox/field/field-traits.h>
#include <linbox/algorithms/matrix-mod.h>
#include <linbox/vector/vector-domain.h>
#include <linbox/randiter/random-prime.h>
#include <linbox/solutions/minpoly.h>
#include <linbox/util/commentator.h>
#include "cra.h"
#include <linbox/fflapack/fflapack.h>

namespace LinBox {

	/* compute the minpoly of a matrix over the Integer ring
	 * via modular method over Fieled.
	 */
	template <class _Integer, class _Field>
	class MinPoly {
	public:
		typedef _Field Field;
		typedef _Integer Integer;
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
		typedef typename MatrixModTrait<IMatrix, Field>::value_type FBlackbox;
		typedef std::vector<Element> FPoly;
		FBlackbox* fbb; FPoly fp;
		integer mmodulus; 
		FieldTraits<Field>::maxModulus(mmodulus);
		long bits = (long) floor (log((double)mmodulus)/M_LN2);
		RandomPrime primeg(bits); long prime;
		for (int i = 0; i < n_try; ++ i) {
			primeg. randomPrime (prime);
			Field F(prime);
			MatrixMod::mod (fbb, M, F);
			LinBox::minpoly (fp, *fbb, F); delete fbb;
			if (degree < (fp.size() - 1)) degree = fp.size() -1;
		}
		return degree;
	}

	template <class _Integer, class _Field>
	template<class Poly, class IMatrix>
	Poly& MinPoly<_Integer, _Field>::minPoly(Poly& y, const IMatrix& M, int degree) {
		if (isSymmetric(M))  {
			std::cout << "Symmetric:\n";
			minPolySymmetric(y, M, degree);
		}
		else {
			std::cout << "NonSymmetric:\n";
			minPolyNonSymmetric(y, M, degree);
		}
		return y;
	}

	template <class _Integer, class _Field>
	template<class Poly, class IMatrix>
	Poly& MinPoly<_Integer, _Field>::minPolyNonSymmetric(Poly& y, const IMatrix& M, int degree) {

		typedef typename MatrixModTrait<IMatrix, Field>::value_type FBlackbox;
		typedef std::vector<Element> FPoly;

		integer mmodulus; 
		FieldTraits<Field>::maxModulus(mmodulus);
		long bits = (long) floor (log((double)mmodulus)/M_LN2);

		RandomPrime primeg(bits); long prime;
		FBlackbox* fbb; 
		FPoly fp (degree + 1);
		typename FPoly::iterator fp_p;
		std::vector<Integer> v(degree + 1);
		typename std::vector<Integer>::iterator v_p;
		y.resize (degree + 1);

		CRA<Integer> cra;
		Integer m = 1;
		while(! cra.terminated()) {
			primeg. randomPrime(prime);
			Field F(prime);
			if ((m % prime) != 0) {
				m *= prime;
				MatrixMod::mod (fbb, M, F);
				LinBox::minpoly (fp, *fbb, F); delete fbb;
				if (fp.size() - 1 != degree) {
					commentator.report (Commentator::LEVEL_IMPORTANT,
										INTERNAL_DESCRIPTION) << "Bad prime.\n";
					continue;
				}
			}
			else {
				commentator.report (Commentator::LEVEL_IMPORTANT,
									INTERNAL_DESCRIPTION) << "Repeated prime.\n";
				continue;
			}
			for (fp_p = fp.begin(), v_p = v.begin(); fp_p != fp.end(); ++fp_p, ++v_p)
				F. convert (*v_p, *fp_p);

			cra.step(prime, v);
		}
			
		cra. result (y);
		std::cout << "Number of primes needed: " << cra. steps() << std::endl;
		return y;
	}

	template <class _Integer, class _Field>
	template<class Poly, class IMatrix>
	Poly& MinPoly<_Integer, _Field>::minPolySymmetric(Poly& y, const IMatrix& M, int degree) {

		typedef typename MatrixModTrait<IMatrix, Field>::value_type FBlackbox;
		typedef std::vector<Element> FPoly;

		CRA<Integer> cra; 
		integer mmodulus; 
		FieldTraits<Field>::maxModulus(mmodulus);
		long bits = (long) floor (log((double)mmodulus)/M_LN2);

		RandomPrime primeg(bits); long prime;
		FBlackbox* fbb; 
		FPoly fp (degree + 1);
		typename FPoly::iterator fp_p;
		std::vector<Integer> v(degree + 1);
		typename std::vector<Integer>::iterator v_p;
		y.resize (degree + 1);

		Integer m = 1;
		while(! cra.terminated()) {
			primeg. randomPrime(prime); 
			Field F(prime); 
			if ((m % prime) != 0) {
				m *= prime;
				MatrixMod::mod (fbb, M, F); 
				LinBox::minpolySymmetric (fp, *fbb, F); delete fbb;
				if (fp.size() - 1 != degree) {
					commentator.report (Commentator::LEVEL_IMPORTANT,
										INTERNAL_DESCRIPTION) << "Bad prime.\n";
					continue;
				}
			}
			else {
				commentator.report (Commentator::LEVEL_IMPORTANT,
									INTERNAL_DESCRIPTION) << "Repeated prime.\n";
				continue;
			}

			for (fp_p = fp.begin(), v_p = v.begin(); fp_p != fp.end(); ++fp_p, ++v_p)
				F. convert (*v_p, *fp_p);

			cra.step(prime, v);
		}
			
		cra. result (y);
		std::cout << "Number of primes needed: " << cra. steps() << std::endl;
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
		RandomPrime primeg(bit1 < bit2 ? bit1 : bit2); long prime;
		Element* FA = new Element [n*n];
		Element* X = new Element [n*(n+1)];
		size_t* Perm = new size_t[n];
		Element* p;
		typename DenseMatrix<Ring>::ConstRawIterator raw_p;
		std::vector<Element> poly (degree + 1);
		typename std::vector<Element>::iterator poly_ptr;
		std::vector<Integer> v(degree + 1);
		typename std::vector<Integer>::iterator v_p;
		CRA<Integer> cra;
		
		Integer m = 1; 
		while (! cra. terminated()) {
			primeg. randomPrime(prime);
			Field F(prime);
			if ((m % prime) != 0) {
				m *= prime;
				for (p = FA, raw_p = M. rawBegin(); 
					 p != FA + (n*n); ++ p, ++ raw_p)

					F. init (*p, *raw_p);

				FFLAPACK::MinPoly( F, poly, n, FA, n, X, n, Perm);

				if(poly. size() != degree + 1) {
					commentator.report (Commentator::LEVEL_IMPORTANT, 
									    INTERNAL_DESCRIPTION) << "Bad prime.\n";
					continue;
				}
			}
			else {
				commentator.report (Commentator::LEVEL_IMPORTANT, 
								    INTERNAL_DESCRIPTION) << "Repeated prime.\n";
				continue;
			}
			for (poly_ptr = poly.begin(), v_p = v.begin(); 
				 poly_ptr != poly.end(); ++poly_ptr, ++v_p)
				F.convert(*v_p, *poly_ptr);

			cra.step(prime, v);
		}
		cra. result(y);
		std::cout << "Number of primes needed: " << cra. steps() << std::endl;
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

		CRA<Integer> cra; 
		integer mmodulus; 
		FieldTraits<Field>::maxModulus(mmodulus);
		long bit1 = (long) floor (log((double)mmodulus)/M_LN2);
		long bit2 = (long) floor (log(sqrt(double(4503599627370496LL/n)))/M_LN2);
		RandomPrime primeg(bit1 < bit2 ? bit1 : bit2); long prime;
		
		typename DenseMatrix<Ring>::ConstRawIterator raw_p;
		for (int i = 0; i < n_try; ++ i) {
			primeg. randomPrime(prime);
			Field F(prime);
			for (p = FA, raw_p = M. rawBegin(); 
				 p!= FA + (n*n); ++ p, ++ raw_p)
				F. init (*p, *raw_p);

			FFLAPACK::MinPoly( F, Poly, n, FA, n, X, n, Perm);

			if (degree < Poly. size() - 1) 
				degree = Poly. size() -1;
		}
		delete FA; delete X; delete Perm;

		return degree;
	}
}
