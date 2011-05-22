/* linbox/algorithms/rational-solver2.h
 * Copyright (C) 2010 LinBox
 * Author Z. Wan
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





#ifndef __LINBOX_rational_solver2__H
#define __LINBOX_rational_solver2__H

/* 
 * Implementation the algorithm in manuscript, available at http://www.cis.udel.edu/~wan/jsc_wan.ps
 */

#include <memory.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <linbox/integer.h>
#include <linbox/algorithms/rational-reconstruction2.h>

namespace LinBox 
{

#if __LINBOX_HAVE_DGETRF && __LINBOX_HAVE_DGETRI
	template <class Ring, class Field, class RandomPrime>
	inline int RationalSolver<Ring, Field, RandomPrime, NumericalTraits>::cblas_dgeinv(double* M, int n) {
		enum CBLAS_ORDER order = CblasRowMajor;
		int lda = n;
		int P[n];
		int ierr = clapack_dgetrf (order, n, n, M, lda, P);
		if (ierr != 0) {
                    std::cerr << "In RationalSolver::cblas_dgeinv Matrix is not full rank" << std::endl;
                    return -1;
		}
		clapack_dgetri (order, n, M, lda, P);
		return 0;
	}

	template <class Ring, class Field, class RandomPrime>
	inline int RationalSolver<Ring, Field, RandomPrime, NumericalTraits>::cblas_rsol (int n, const double* M, integer* numx, integer& denx, double* b) {
		if (n < 1) return 0;
		double* IM = new double[n * n];
		memcpy ((void*)IM, (const void*)M, sizeof(double)*n*n);
		int ret;
		//compute the inverse by flops
		ret = cblas_dgeinv (IM, n); 
		if (ret != 0) {delete[] IM; return 1;}

		double mnorm = cblas_dOOnorm(M, n, n);
		// residual
		double* r = new double [n];
		// A^{-1}r
		double* x = new double [n];
		//ax  = A x
		double* ax = new double [n];
		// a digit, d \approx \alpha x
		double* d = new double [n];

		const double* p2;
		double* pd;
		const double T = 1 << 30;
	
		integer* num = new integer [n];
		integer* p_mpz;
		integer tmp_mpz, den, denB, B;

		den = 1;
		// compute the hadamard bound
		cblas_hbound (denB, n, n, M);      
		B = denB * denB;
		// shouble be a check for tmp_mpz
		tmp_mpz = 2 * mnorm + cblas_dmax (n, b, 1);
		B <<= 1; B *= tmp_mpz; //B *= tmp_mpz;
	
		//double log2 = log (2.0);
		double log2 = M_LN2;
		// r = b
		memcpy ((void*) r, (const void*) b, sizeof(double)*n);

		do  {
			cblas_dapply (n, n, IM, r, x);
			// compute ax
			cblas_dapply (n, n, M, x, ax);
			// compute ax = ax -r, the negative of residual
			cblas_daxpy (n, -1, r, 1, ax, 1);
			// compute possible shift
			double normr1, normr2, normr3, shift1, shift2;
			normr1 = cblas_dmax(n, r, 1);
			normr2 = cblas_dmax(n, ax, 1);
			normr3 = cblas_dmax(n, x, 1);
			//try to find a good scalar
			int shift = 30;
			if (normr2 <.0000000001) 
				shift = 30;
			else {
				shift1 = floor(log (normr1 / normr2) / log2) - 2;
				shift = (int)(30 < shift1 ? 30 : shift1);
			}

			normr3 = normr3 > 2 ? normr3 : 2;
			shift2 = floor(53. * log2 / log (normr3));
			shift = (int)(shift < shift2 ? shift : shift2);

			if (shift <= 0) {
#ifdef DEBUGRC
				printf ("%s", "Bad scalar \n");
				printf("%f, %f\n", normr1, normr2);
				printf ("%d, shift = ", shift);
				printf ("OO-norm of matrix: %f\n", cblas_dOOnorm(M, n, n));
				printf ("OO-norm of inverse: %f\n", cblas_dOOnorm(IM, n, n));
				printf ("Error, abort\n"); 
#endif
				delete[] IM; delete[] r; delete[] x; delete[] ax; delete[] d; delete[] num; 
				return 2;
			}

			int scalar = ((long long int)1 << shift);
			for (pd = d, p2 = x; pd != d + n; ++ pd, ++ p2) 
				//better use round, but sun sparc machine doesnot supprot it
				*pd = floor (*p2 * scalar);

			// update den
			den <<= shift; 
			//update num
			update_num (num, n, d, shift);
		
#ifdef DEBUGRC
			printf ("in iteration\n");
			printf ("residual=\n");
			printvec (r,  n);
			printf ("A^(-1) r\n");
			printvec (x,  n);
			printf ("scalar= ");
			printf ("%d \n", scalar);
			printf ("One digit=\n");
			printvec (d, n);
			printf ("Current bound= \n");
			std::cout << B;
			printf ("den= \n");
			std::cout << den;
			printf ("accumulate numerator=\n");
			printvec (num, n);
#endif
			// update r = r * shift - M d
			double tmp = 2 * mnorm + cblas_dmax (n, r, 1);
			if (tmp < T) update_r_int (r, n, M, d, shift);
			else update_r_ll (r, n, M, d, shift);
			//update_r_ll (r, n, M, d, shift);
		} while (den < B);
	
		integer q, rem, den_lcm, tmp_den;
		integer* p_x, * p_x1;
		p_mpz = num;
		p_x = numx;
		// construct first answer
		rational_reconstruction (*p_x, denx, *p_mpz, den, denB);
		++ p_mpz;
		++ p_x;

		int sgn;
		for (; p_mpz != num + n; ++ p_mpz, ++ p_x)  {
			sgn = sign (*p_mpz);
			tmp_mpz = denx * (*p_mpz);
			tmp_mpz = abs (tmp_mpz);
			integer::divmod (q, rem, tmp_mpz, den);
		
			if ( rem < denx)  {
				if (sgn >= 0)
					*p_x = q;
				else 
					*p_x = -q;
			}
			else {
				rem = den - rem;
				q += 1;
				if (rem < denx) {
					if (sgn >= 0)
						*p_x = q;
					else
						*p_x = -q;
				}
				else {
					rational_reconstruction (*p_x, tmp_den, *p_mpz, den, denB);
					lcm (den_lcm, tmp_den, denx);
					integer::divexact (tmp_mpz, den_lcm, tmp_den);
					integer::mul (*p_x, *p_x, tmp_mpz);
					integer::divexact (tmp_mpz, den_lcm, denx);
					denx = den_lcm;
					for (p_x1 = numx; p_x1 != p_x; ++ p_x1)
						integer::mul (*p_x1, *p_x1, tmp_mpz);
				}
			}	
		}
#ifdef DEBUGRC
		std::cout << "rational answer\nCommon den = ";
		std::cout << denx;
		std::cout << "\nNumerator= \n";
		printvec (numx, n);
#endif
	
		//normalize the answer
		if (denx != 0) {	
			integer g; g = denx;
			for (p_x = numx; p_x != numx + n; ++ p_x)
				g = gcd (g, *p_x);
			for (p_x = numx; p_x != numx + n; ++ p_x)
				integer::divexact (*p_x, *p_x, g);
			integer::divexact (denx, denx, g);
		}
	
		//check if the answer is correct, not necessary
		cblas_mpzapply (n, n, M, (const integer*)numx, num);
		integer* sb = new integer [n];
		double* p;
		for (p_mpz = sb, p = b; p_mpz != sb + n; ++ p_mpz, ++ p) {
			*p_mpz = *p;
			integer::mulin(*p_mpz, denx);
		}
		ret = 0;
		for (p_mpz = sb, p_x = num; p_mpz != sb + n; ++ p_mpz, ++ p_x) 
			if (*p_mpz != *p_x) {
				ret = 3;
				break;
			}
#ifdef DEBUGRC
		if (ret == 3) {
	
			std::cout << "Input matrix:\n";
			for (int i = 0; i < n; ++ i) {
				const double* p = M + (i * n);
				printvec (p, n);
			}
			std::cout << "Input rhs:\n";
			printvec (b, n);
			std::cout << "Common den: " << denx << '\n';
			std::cout << "Numerator: ";
			printvec (numx, n);
			std::cout << "A num: ";
			printvec (num, n);
			std::cout << "denx rhs: ";
			printvec (sb, n);
		}
#endif

		// garbage collector
		delete[] IM; delete[] r; delete[] x; delete[] ax; delete[] d; delete[] num; delete[] sb;

		return ret;
	}

#endif


	template <class Ring, class Field, class RandomPrime>
	/* apply  y <- Ax */
	inline int RationalSolver<Ring, Field, RandomPrime, NumericalTraits>::cblas_dapply (int m, int n, const double* A, const double* x, double* y) {
		cblas_dgemv (CblasRowMajor, CblasNoTrans, m, n, 1, A, n, x, 1, 0, y, 1);
		return 0;
	}

	template <class Ring, class Field, class RandomPrime>
	inline int RationalSolver<Ring, Field, RandomPrime, NumericalTraits>::cblas_mpzapply (int m, int n, const double* A, const integer* x, integer* y) {
		const double* p_A;
		const integer* p_x;
		integer* p_y;
		integer tmp;
		for (p_A = A, p_y = y; p_y != y + m; ++ p_y) {
			*p_y = 0;
			for (p_x = x; p_x != x + n; ++ p_x, ++ p_A) {
				//mpz_set_d (tmp, *p_A);
				//mpz_addmul_si (*p_y, *p_x, (int)(*p_A));
				tmp = *p_x  * (long long int)(*p_A);
				integer::addin (*p_y, tmp);
			}
		}
		return 0;
	}

	template <class Ring, class Field, class RandomPrime>
	template <class Elt>
	inline int RationalSolver<Ring, Field, RandomPrime, NumericalTraits>::printvec (const Elt* v, int n) {
		const Elt* p;
		std::cout << "\[";
		for (p = v; p != v + n; ++ p)
			std::cout << *p << ' ';
		std::cout << "]\n";
		return 0;
	}

	template <class Ring, class Field, class RandomPrime>
	//update num, *num <- *num * 2^shift + d
	inline int RationalSolver<Ring, Field, RandomPrime, NumericalTraits>::update_num (integer* num, int n, const double* d, int shift) {
		integer* p_mpz;
		integer tmp_mpz;
		const double* pd;
		for (p_mpz = num, pd = d; p_mpz != num + n; ++ p_mpz, ++ pd) {
			(*p_mpz) = (*p_mpz) << shift;
			tmp_mpz = *pd;
			integer::add (*p_mpz, *p_mpz, tmp_mpz);
		} 
		return 0;
	}

	template <class Ring, class Field, class RandomPrime>
	//update r = r * shift - M d
	inline int RationalSolver<Ring, Field, RandomPrime, NumericalTraits>::update_r_int (double* r, int n, const double* M, const double* d, int shift) {
		int tmp;
		double* p1;
		const double* p2;
		const double* pd;
		for (p1 = r, p2 = M; p1 != r + n; ++ p1) {
			tmp = (int)(long long int) *p1;
			tmp <<= shift;
			for (pd = d; pd != d + n; ++ pd, ++ p2) {
				tmp -= (int)(long long int)*pd * (int)(long long int)*p2;
			}
			*p1 = (double)tmp;
		}
		return 0;
	}

	template <class Ring, class Field, class RandomPrime>
	//update r = r * shift - M d
	inline int RationalSolver<Ring, Field, RandomPrime, NumericalTraits>::update_r_ll (double* r, int n, const double* M, const double* d, int shift) {
		long long int tmp;
		double* p1;
		const double* p2;
		const double* pd;
		for (p1 = r, p2 = M; p1 != r + n; ++ p1) {
			tmp = (long long int) *p1;
			tmp <<= shift;
			for (pd = d; pd != d + n; ++ pd, ++ p2) {
				tmp -= (long long int)*pd * (long long int) *p2;
			}
			*p1 = tmp;
		}
		return 0;
	}

	template <class Ring, class Field, class RandomPrime>
	inline double RationalSolver<Ring, Field, RandomPrime, NumericalTraits>::cblas_dOOnorm(const double* M, int m, int n) {
		double norm = 0;
		double old = 0;
		const double* p;
		for (p = M; p != M + (m * n); ) {
			old = norm;
			norm = cblas_dasum (n, p ,1);
			if (norm < old) norm = old;
			p += n;
		}	
		return norm;
	}

	template <class Ring, class Field, class RandomPrime>
	inline double RationalSolver<Ring, Field, RandomPrime, NumericalTraits>::cblas_dmax (const int N, const double* a, const int inc) {
		return fabs(a[cblas_idamax (N, a, inc)]);
	}

	template <class Ring, class Field, class RandomPrime>
	inline int RationalSolver<Ring, Field, RandomPrime, NumericalTraits>::cblas_hbound (integer& b, int m, int n, const double* M){
		double norm = 0;
		const  double* p;
		integer tmp;
		b = 1;
		for (p = M; p != M + (m * n); ) {
			norm = cblas_dnrm2 (n, p ,1);
			tmp =  norm;
			integer::mulin (b, tmp);
			p += n;
		}

		return 0;
	}
}//LinBox

#endif //__LINBOX_rational_solver2__H


/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
