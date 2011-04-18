/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
/* numeric-solver-lapack.h
 *  numeric solver using lapack routines
 *  to support development of iterative numeric/symbolic solver
 */

#ifndef __NUMERIC_SOLVER_LAPACK_H
#define __NUMERIC_SOLVER_LAPACK_H

#include <linbox/blackbox/blas-blackbox.h>

namespace LinBox {

template <class Matrix>
struct LPS {
	
	LPS() : _Ap(NULL), _IM(NULL), _m(0), _n(0) {}
	~LPS() { delete _IM; } 
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
int LPS<Matrix>::init(Matrix& A) { 
	_Ap = &A; // would need memcpy if *_Ap to get mods.
	_m = A.rowdim();
	_n = A.coldim();
	
   //  kludgey pointer to beginning of double vals
   void *thedata = &*(_Ap->rawBegin());

	_IM = new double[_n * _n];
   memcpy((void *)_IM, thedata, sizeof(double)*_m*_n);

	// time to set up inverse of matrix
	int lda = _n;
	int P[_n];
	//std::cerr << "Bef getrf: M_0,0 " << *_IM << ", M_n-1,n-1 " << *(_IM+_n*_n-1) << std::endl;
	int ierr = clapack_dgetrf (CblasRowMajor, _n, _n, _IM, lda, P);
	//std::cerr << "Aft getrf: M_0,0 " << *_IM << ", M_n-1,n-1 " << *(_IM+_n*_n-1) << std::endl;
	if (ierr != 0) {
		//std::cerr << "In LPS::init Matrix is not full rank" << std::endl;
		return -1;
	}
	clapack_dgetri (CblasRowMajor, _n, _IM, lda, P);

	return 0;
} 

template <class Matrix>
template<class Vector>
int LPS<Matrix>::solve(Vector& x, const Vector& b) {
	const double * bdata = &*(b.begin());
	double * xdata = &*(x.begin());
	
	cblas_dgemv(CblasRowMajor, CblasNoTrans, _m, _n, 1, _IM, _n, bdata, 1, 0, xdata, 1);

	return 0;
}

template <class Matrix>
template<class Vector>
Vector& LPS<Matrix>::apply(Vector& y, const Vector& x) { 
	const double * xdata = &*(x.begin());
	double * ydata = &*(y.begin());
   double *thedata = &*(_Ap->rawBegin());

	cblas_dgemv(CblasRowMajor, CblasNoTrans, _m, _n, 1, thedata, _n, xdata, 1, 0, ydata, 1);

	return y;
}

} // namespace LinBox

#endif // __NUMERIC_SOLVER_LAPACK_H
