#ifndef __SLICED_DOMAIN_H
#define __SLICED_DOMAIN_H

#include "dense-sliced.h"

/*
The SlicedDomain template over a Field type parameter has constructor from an instance of Field which must represent GF(3).
The SlicedDomain::Matrix subtype meets the sliced dense matrix concept
and it is the interface for working with sliced matrices.
A Blackbox is any matrix type that has apply and applyTranspose applicable to Matrix.  That is to say the signature of apply is 
  Matrix& Blackbox::apply(Matrix& Y, const Matrix& X)

SlicedDomain provides, for A a Matrix and B a Blackbox(preconditioner) 
  mulin(A, B) // A *= B
  addin(A, A2) // A += A2
*/

template<class _Field, typename _WordT = unsigned long long int>
// _Field must be able to represent GF(3).
struct SlicedDomain : public _Field
{
	typedef _Field Field;
	typedef _WordT Word_T;
	typedef typename _Field::Element Element;
	typedef Element Scalar;
	typedef Sliced<SlicedDomain> Matrix;

	SlicedDomain () : Field(3) {}
	//SlicedDomain (size_t p, size_t e = 1) : Field(p, e) {}
	SlicedDomain (size_t p, size_t e = 1) : Field(p) 
	{  if (p != 3 || e != 1) throw LinBoxError("bad field"); }

	SlicedDomain& operator=(const SlicedDomain<Field, Word_T>& R) 
	{ static_cast<Field*>(this)->operator= (R); return *this; }

// A domain provides field functions and (dense) matrix functions.

// field functions
	using Field::mOne;
	using Field::zero;
	using Field::one;
	using Field::init;
	using Field::add;	
	using Field::sub;	
	using Field::neg;	
	using Field::mul;	
	using Field::div;	
	using Field::inv;	
	using Field::axpy;	
	using Field::addin;	
	using Field::subin;	
	using Field::negin;	
	using Field::mulin;	
	using Field::divin;	
	using Field::invin;	
	using Field::axpyin;	
	using Field::isZero;	
	using Field::isOne;	
	using Field::areEqual;	
	using Field::cardinality;
	using Field::characteristic;
	using Field::write;
	using Field::read;

// matrix functions

	// Y += X, where X and Y are conformally packed (both row or both col).
	Matrix& addin (Matrix& Y, Matrix& X) { 
		return Y.addin(X);
	}	

	// X *= a
	Matrix& smulin (Matrix& X, Scalar& a) { 
		return X.smulin(a);
	}

	// X -= X
	Matrix& neg (Matrix& X) { 
		return X.smulin(mOne);
	}

	//  return C <-- A * B
	//  is assuming row sliced
	template <class Gettable> // any matrix rep with getEntry().
	Matrix & mul (Matrix& C, Gettable& A, Matrix& B){
		return C.mul(A,B);
	}

	// A += x*B
	Matrix& axpyin( Matrix& A, Scalar& x, Matrix &B) {
		typename Matrix::RawIterator Ab(A.rawBegin()), Ae(A.rawEnd()), Bb(B.rawBegin());
		return A.axpyin(Ab, Ae, x, Bb);
	}

	//  C += A * B
	Matrix& axpyin(Matrix& C, Matrix& A, Matrix &B) {
		//  temp mat to store mul
		Matrix T;
		T.init(A.rowdim(), B.coldim());
		mul(T, A, B);
		return addin(C, T);	
	}	

	Matrix& random(Matrix &A, size_t seed=0) const {
		return A.random(seed);
	}

	Matrix& clear(Matrix &A) const {
		return A.zero();
	}

	// assignment operator will not be a deep copy
	Matrix& deepcopy(Matrix &dst, Matrix& src) const {
		return dst.deepcopy(src);
	}

	// submatrix is a Matrix member function, not here.
	/*
	Matrix& submatrix(Matrix& super, size_t i, size_t j, size_t m, size_t n) const {
		return *(new Matrix(super, i, j, m, n));
	}
	*/

	bool areEqual(Matrix& A, Matrix &B) const {
		return A.isEqual(B);
	}

	std::ostream& write(std::ostream& os, Matrix& A) const { 
		return A.write(os << A.rowdim() << " " << A.coldim() << std::endl); 
	}
	std::istream& read(std::istream& is, Matrix& A) const { 
		size_t r, c;
		Element x; init(x);
		is >> r >> c;
		A.init(r, c);
		for (size_t i = 0; i < r; ++i)
			for (size_t j = 0; j < c; ++j){
				read(is, x);
				A.setEntry(i, j, x);
			}
		return is; 
	}
};
	
#endif // __SLICED_DOMAIN_H
