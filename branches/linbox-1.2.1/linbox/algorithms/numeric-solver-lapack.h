/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/* Copyright (C) 2011 LinBox
 * Written Bryan Youse <>
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
 * Free Software Foundation, Inc., 51 Franklin Street - Fifth Floor,
 * Boston, MA 02110-1301, USA.
 */

/* numeric-solver-lapack.h
 *  numeric solver using lapack routines
 *  to support development of iterative numeric/symbolic solver
 */

#ifndef __LINBOX_numeric_solver_lapack_H
#define __LINBOX_numeric_solver_lapack_H

#include "linbox/blackbox/blas-blackbox.h"

namespace LinBox {

	template <class Matrix>
	struct LPS {

		LPS() : _Ap(NULL), _IM(NULL), _m(0), _n(0) {}
		~LPS() {
			if(_IM) {
				// std::cout << "delete" << std::endl;
				delete[] _IM;
			}
		}
		LPS(Matrix& A) { init(A); }
		int init(Matrix & A); // set up for solving - expect multiple subsequent calls to solve() and apply().

		template<class Vector> int solve(Vector& x, const Vector& b); // x such that Ax = b (approx)
		template<class Vector> Vector& apply(Vector& y, const Vector& x); // y = Ax (approx)

	protected:
		Matrix* _Ap; // for right now, assume this points to the input, double matrix A.
		double *_IM;

		size_t _m, _n;
	};

	template <class Matrix>
	int LPS<Matrix>::init(Matrix& A)
	{
		_Ap = &A; // would need memcpy if *_Ap to get mods.
		_m = A.rowdim();
		_n = A.coldim();

		//  kludgey pointer to beginning of double vals
		void *thedata = &*(_Ap->Begin());

		linbox_check(_n);
		_IM = new double[_n * _n];
		memcpy((void *)_IM, thedata, sizeof(double)*_m*_n);

		// time to set up inverse of matrix
		int lda = (int)_n;
		int * P = new int[_n];
		// std::cerr << "Bef getrf: M_0,0 " << *_IM << ", M_n-1,n-1 " << *(_IM+_n*_n-1) << std::endl;
		int ierr = clapack_dgetrf (CblasRowMajor, (int)_n, (int)_n, _IM, lda, P);
		// std::cerr << "Aft getrf: M_0,0 " << *_IM << ", M_n-1,n-1 " << *(_IM+_n*_n-1) << std::endl;
		if (ierr != 0) {
			// std::cerr << "In LPS::init Matrix is not full rank" << std::endl;
			delete[] P ;
			return -1;
		}
		clapack_dgetri (CblasRowMajor, (int)_n, _IM, lda, P);
		delete[] P ;

		return 0;
	}

	template <class Matrix>
	template<class Vector>
	int LPS<Matrix>::solve(Vector& x, const Vector& b)
	{
		linbox_check(typeid(typename Vector::value_type)==typeid(double));
		// std::cout << "input :" << b << std::endl;
		const double * bdata = &*(b.begin());
		double * xdata = &*(x.begin());

		cblas_dgemv(CblasRowMajor, CblasNoTrans, (int)_m, (int)_n, 1, _IM, (int)_n, bdata, 1, 0, xdata, 1);
		// std::cout << "result to solve :" << x << std::endl;

		return 0;
	}

	template <class Matrix>
	template<class Vector>
	Vector& LPS<Matrix>::apply(Vector& y, const Vector& x)
	{
		// std::cout << "input :" << x << std::endl;
		const double * xdata = &*(x.begin());
		double * ydata = &*(y.begin());
		double *thedata = &*(_Ap->Begin());

		cblas_dgemv(CblasRowMajor, CblasNoTrans, (int)_m, (int)_n, 1, thedata, (int)_n, xdata, 1, 0, ydata, 1);
		// std::cout << "result to apply :" << y << std::endl;

		return y;
	}

} // namespace LinBox

#endif // __LINBOX_numeric_solver_lapack_H
