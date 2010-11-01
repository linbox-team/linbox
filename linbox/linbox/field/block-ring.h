/* linbox/fields/blas-ring.h
 * Copyright (C) 2007 LinBox Team
 *
 * Written by JP May, with tweaks by D. Saunders, Z. Wan
 *
 * See COPYING for license information.
 */


#ifndef __LINBOX_blockring_H
#define __LINBOX_blockring_H
#include <iostream>
#include "linbox/matrix/blas-matrix.h"
#include "linbox/field/field-interface.h"
//#include "linbox/matrix/blas-domain.h"
#include "linbox/algorithms/blas-domain.h"
#include "linbox/fflas/fflas.h"

namespace LinBox 
{
	
	
/** Elements are wrapped BlasMatrix objects.

    Operations expect conformal sizes for inputs and outputs.

    You can expect good performance due to BLAS usage,
    especially when the Field is Modular<double> or Modular<float>.
 */

template < class _Field >
class BlockRing : public FieldInterface
{
    public:
	_Field _F;
	BlasMatrixDomain<_Field> _D;        
	size_t _b;
	
	typedef BlasMatrix<typename _Field::Element> Matrix;
	typedef typename _Field::Element Scalar;
	
	
	
	/// default constructable wrapper for BlasMatrix
	class Element 
	{ 
		
	    public:
		
		typedef _Field Field;
		typedef Scalar Entry;
		
		Element() : matrix(0) {}
		
		~Element() { 
			release(); 
		}
		
		// copy constructor
		Element(const Element& e) : matrix(0){
			if (e.matrix != 0) {
				matrix = new Matrix(*(e.matrix));
				// memory leak of previous value?
			}
		}
		
		// overload assignment
		Element& operator= (const Element& e) {
			if (matrix == e.matrix) {
				return *this;
			}  else if (e.matrix == 0) {
				release();
				return *this;
			} else {
				//set(new Matrix(*(e.matrix))); // does this really copy?
				clone(e);
				return *this;
			}
		}
		
		void clone(const Element& A) {
			// make this a deep copy of A
			// BlasMatrix copy constructor is shallow!
			
			release();
			
			if(A.matrix == 0) return;
			
			size_t rows = A.matrix->rowdim();
			size_t cols = A.matrix->coldim();
			
			set(new Matrix(rows, cols));
			
			Scalar* a=A.matrix->getPointer();
			Scalar* b=  matrix->getPointer();
			
			for(size_t i=0; i < rows*cols; i++) {
				*b=*a;
				a++; b++;
			}
			
			
			}
		
		// cleanly deletes the current stored value
		// before assigning to the new value
		void set(Matrix* thematrix) {
			this -> release();
			matrix = thematrix;
		}
		
		Matrix* matrix;
		
	    private:
		
		void release() {
			if (matrix != 0 ) 
				delete matrix;
			matrix = 0;
		}  
		
	}; // class Element
	
	
	class RandIter 
	{
		typedef typename _Field::RandIter FieldRandIter;
		
		FieldRandIter r;
		size_t dim;
		
	    public:
		RandIter(const BlockRing<_Field>& BR, 
			 const integer size=0, 
			 const integer seed=0) :
			r(BR._F, size, seed), dim(BR._b) {}
		
		Element& random(Element& e) const {
			// e must be init'd
			for(size_t i=0; i < e. matrix -> rowdim(); i++)
				for(size_t j=0; j < e. matrix -> coldim(); j++)
					r.random(e.matrix->refEntry(i,j));
			return e;
		}
		
	}; //class RandIter
	
	
	BlockRing(const _Field& F, size_t d=1)
		: _F(F), _D(F), _b(d) {}
	
	Element& init(Element& B) const {
		// B is garbage from memory
		B.set(new Matrix(_b,_b));
		return B;
	}
	
	template <typename ints>
		Element& init(Element& B, ints n, size_t r = 0, size_t c = 0) const
		// n supposed to be integer, r num rows, c num cols
	{ 
		// default block dim is default dim of ring, but others are allowed.
		if (r == 0) r = _b; 
		if (c == 0) c = _b; 
		
		B.set(new Matrix(r,c));
		
		size_t k = ( (r < c) ? r : c );
		
		typename _Field::Element N; _F.init(N, n);
		
		for (size_t i = 0; i < k; ++i) (B.matrix)->setEntry(i, i, N);
		
		return B;
		
	}
	
	
	template <typename ints>
	ints& convert(ints& x) const
	{
		return _F.convert(x);
	}
	
	
	template <typename ints>
	ints& convert(ints& x, const Element &A) const 
	{
		return _F.convert(x, *(A.matrix->getPointer()));
	}
	
	Element& assign(Element &A, const Element &B) const
	{
		return A = B;
	}
	
	
	integer& cardinality(integer &c) const
	{
		// c = p^(b^2)
		
		_F.cardinality(c);
		
			if(c > 1) // _F is a finite field
			{
				integer tmp, n;
				n = _b*_b;
				c = expt(tmp, c, n);
				} // else c  = -1
			
			return c;
	}
	
	
	integer& characteristic(integer &c) const
	{
		return _F.characteristic(c);
	}
	
	size_t dim() const { return _b; }
	
	
	//Operations
	// 
	// All operations will work for matrices of dim != _b
	// but assume that the dimensions of all the given matrices
	// are compatible with the dimensions of the A operand.
	
	//Operations from the Matrix Domain:
	
	Element& mul(Element& C, const Element& A, const Element& B) const
	{
		_D.mul(*(C.matrix), *(A.matrix), *(B.matrix)); 
		return C;
	}
	
	
	//non-commutative: use mulin_left: A = A*B
	Element& mulin(Element& A, const Element& B) const
	{
		_D.mulin_left(*(A.matrix), *(B.matrix)); 
		return A;
	}
	
	// D = A*X+Y
	Element& axpy(Element& D, const Element& A, const Element& X, const Element& Y) const
		{
			_D.axpy(*(D.matrix), *(A.matrix), *(X.matrix), *(Y.matrix)); 
			return D;
		}
	
	// R = A*X+R
	Element& axpyin(Element& R, const Element& A, const Element& X) const
	{
			_D.axpyin(*(R.matrix), *(A.matrix), *(X.matrix)); 
			return R;
	}
	
	
	// These operations will not work for all elements
	// and no checks are provided!
	
	// B = A^{-1}
	Element& inv(Element& B, const Element& A) const {
		
		int nullflag = 0;
		
		_D.inv(*(B.matrix), *(A.matrix), nullflag);

                if (nullflag)
                  throw PreconditionFailed(__FUNCTION__,__LINE__,"InvMatrix: inverse undefined");
			
		return B;
	}
	
	// A=A^{-1} not really inplace!
	Element& invin(Element& A) const {
		
		int nullflag = 0;
		
		//_D.invin(A, A, nullflag);
		
		Element B;  init(B, A.matrix->rowdim(), A.matrix->coldim());
                _D.inv(*(B.matrix), *(A.matrix), nullflag);
                
                if (nullflag)
                     throw PreconditionFailed(__FUNCTION__,__LINE__,"InvMatrix: inverse undefined");
                
                A=B;
			
		return A;
	}
	
	
		// C = A*B^{-1}
	Element& div(Element& C, const Element& A, const Element& B) const {
		
		_D.right_solve(*(C.matrix),*(B.matrix),*(A.matrix));
		return C;
	}
		
	
	//A = A*B^{-1};
	Element& divin( Element& A, const Element& B) const {
		_D.right_solve(*(B.matrix),*(A.matrix)); 
		return A;
	}
	
		
	
	
	// Unwrapped operations using simple loops:
	
	// C = A + B
	Element& add(Element& C, const Element& A, const Element& B) const
	{
		size_t rows = A.matrix->rowdim();
		size_t cols = A.matrix->coldim();
		
		Scalar* a=A.matrix->getPointer();
		Scalar* b=B.matrix->getPointer();
		Scalar* c=C.matrix->getPointer();
		
		//FFLAS::fcopy(_F, rows*cols, b, 1, c, 1); // C = B
		
		
		for(size_t i=0; i < rows*cols; i++) {
			_F.add(*c,*a,*b);
				a++; b++; c++;
		}
		
		//Scalar alpha; _F.init(alpha, 1);
		//FFLAS::faxpy(_F, rows*cols, alpha, a, 1, c, 1);
		
		return C;
	}
	
	// A = A + B
	Element& addin(Element& A, const Element& B) const
	{
		size_t r = A.matrix->rowdim();
		size_t c = A.matrix->coldim();
		
		Scalar* a=A.matrix->getPointer();
		Scalar* b=B.matrix->getPointer();
		
		for(size_t i=0; i < r*c; i++) {
			_F.addin(*a,*b);
			a++; b++;
		}
		
		return A;
		
	}
	
	
	// C = A - B
	Element& sub(Element& C, const Element& A, const Element& B) const
	{
		size_t rows = A.matrix->rowdim();
		size_t cols = A.matrix->coldim();
		
		Scalar* a=A.matrix->getPointer();
		Scalar* b=B.matrix->getPointer();
		Scalar* c=C.matrix->getPointer();
		
		
		for(size_t i=0; i < rows*cols; i++) {
			_F.sub(*c,*a,*b);
			a++; b++; c++;
		}
		
		return C;
	}
	
	
	// A = A - B
	Element& subin(Element& A, const Element& B) const
	{
		size_t r = A.matrix->rowdim();
		size_t c = A.matrix->coldim();
		
		Scalar* a=A.matrix->getPointer();
		Scalar* b=B.matrix->getPointer();
		
		for(size_t i=0; i < r*c; i++) {
			_F.subin(*a,*b);
			a++; b++;
		}
		
		return A;
	}
	
	
	//B = -1*A
	Element& neg(Element& B, const Element& A) const
	{
		size_t r = A.matrix->rowdim();
		size_t c = A.matrix->coldim();
		
		Scalar* a=A.matrix->getPointer();
		Scalar* b=B.matrix->getPointer();
		
		for(size_t i=0; i < r*c; i++) {
			_F.neg(*b,*a);
			a++; b++;
		}
		
		return B;
	}
	
	
	//A = -1*A
	Element& negin(Element& A) const
	{
		size_t r = A.matrix->rowdim();
		size_t c = A.matrix->coldim();
		
		Scalar* a=A.matrix->getPointer();
		
		for(size_t i=0; i < r*c; i++) {
			_F.negin(*a);
			a++;
		}
		
		return A;
	}
	
	bool areEqual(const Element& A, const Element& B) const
	{
		
		size_t r = A.matrix->rowdim();
		size_t c = A.matrix->coldim();
		
		Scalar* a=A.matrix->getPointer();
		Scalar* b=B.matrix->getPointer();
		
		for(size_t i=0; i < r*c; i++) {
			
			if(!_F.areEqual(*a,*b)) {
				return false;
			}
			
			a++; b++;
		}
		return true;
	}
	
	
	bool isOne(const Element& X) const
	{
		size_t n = X.matrix->rowdim();
		
		if(n != X.matrix->coldim()) {
			return false;
		}
		
		Scalar* x=X.matrix->getPointer();
		
		for(size_t i=1; i <= n; i++) 
			for(size_t j=1; j <= n; j++) 
			{
				if(i==j) { // on the diagonal
					if(!_F.isOne(*x)) {
						return false;
					}
				}
				else {
					if(!_F.isZero(*x)) {
						return false;
					}
				}
				
				x++;
			}
		
		
		return true;
	}
	
	
	bool isZero(const Element& X) const
	{
			
		size_t r = X.matrix->rowdim();
		size_t c = X.matrix->coldim();
		
		Scalar* x=X.matrix->getPointer();
		
		for(size_t i=0; i < r*c; i++) 
		{
			if(!_F.isZero(*x)) {
				return false;
			}
			
			x++;  
		}
		
		return true;
	}
	
	
	//stubs for read and write field
	std::ostream& write(std::ostream& os) const 
	{
		return _F.write(os << "Dimension " << _b << " square matrices over ");
	}
	
	
	
	std::istream& read(std::istream& is) 
	{
		return is;
	}
	
	// wrapped read and write element
	std::ostream& write(std::ostream& os, const Element& A) const 
	{
		return (A.matrix)->write(os << std::endl, _F);
	}
	
	
	std::istream& read(std::istream& is, const Element& A) const 
	{
		
		return (A.matrix)->read(is, _F);
	}
	
	
    private:
	
	
	// recursive helper to compute exponentiation of integers
	static	  integer& expt (integer& res, integer& a, integer& n)
	{
		if (n == 0) {
			res=1;
		}
		else if (n == 1) {
			res=a;
		}
		else if (n[0] & 1) {
			n -= 1;
			expt(res, a, n);
			res*=a;
		} else {
			n /= 2;
			expt(res, a, n);
			res*=res;
		}
		
		return res;
	}
	
	
	
  }; // BlockRing
	
} // LinBox

#endif // __LINBOX_blockring_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
