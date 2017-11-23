#ifndef __LINBOX_plain_domain_h
#define __LINBOX_plain_domain_h

#include "linbox/util/error.h"
#include "linbox/matrix/plain-matrix.h"

namespace LinBox {
/*
   PlainDomain is a reference implementation of a matrix domain.
   It makes a domain from any field by contributing very basic matrix functions using
   the getEntry/setEntry functions of the DenseMatrix interface.

PlainDomain has a constructor from field,
it's PlainDomain::Matrix type meets the dense submatrix concept
and it is the interface for working with dense matrices.  When you need to allocate a new matrix use
PlainDomain::NewMatrix.

A Blackbox over PlainDomain<Field> is any matrix type that has apply and applyTranspose applicable to PlainDomain<Field>::Matrix.  That is to say the signature of apply is
  Matrix& Blackbox::apply(Matrix& Y, const Matrix& X)

PlainDomain provides, for dense Submatrices A,B,C
  mul(C, A, B) // C = A * B
  axpyin(C, A, B) // C += A * B
  addin(A, B) // A += B

  functions on Submatrix
  add sub neg mul div inv axpy (both scalar and matrix
  inplace forms of those
  areEqual isZero

  alternative functions
  gemm (rplace axpy) (trans)
  additional functions
  trsm (trans, uplo) trmm(trans, uplo) symm (trans, uplo)
  (in particular triangular solving must be included
*/

template<class Field_>
struct PlainDomain : public Field_
{
	typedef PlainDomain<Field_> Self_t;
	typedef Field_ Father_t;
	typedef size_t Index;
	// A domain provides distinct types: Scalar, Matrix, and Blackbox.
	typedef typename Father_t::Element Scalar;
	typedef typename Father_t::Element Element;
	typedef PlainSubmatrix<Self_t> Submatrix;
	typedef PlainMatrix<Self_t> Matrix;
	// but matrix should admit triangular forms

	// constructors and assignment
	//PlainDomain (size_t p = 0, size_t e = 1): Field(p) {}
	PlainDomain()
	{}
	PlainDomain(const Element& p)
	: Father_t(p) {}
	PlainDomain(const Father_t& F)
	: Father_t(F) {}
	PlainDomain(const PlainDomain<Father_t>& D)
	: Father_t(D) {}
	//~PlainDomain()
	// {} // default dstor is fine.
	using Father_t::operator=;
	const Father_t& field() const // transitional
	{ return *this; }

// A domain provides field functions and (dense) matrix functions.

// field functions
	using Father_t::init;
	using Father_t::add;
	using Father_t::sub;
	using Father_t::neg;
	using Father_t::mul;
	using Father_t::div;
	using Father_t::inv;
	using Father_t::axpy;
	using Father_t::addin;
	using Father_t::subin;
	using Father_t::negin;
	using Father_t::mulin;
	using Father_t::divin;
	using Father_t::invin;
	using Father_t::axpyin;
	using Father_t::isZero;
	using Father_t::isOne;
	using Father_t::areEqual;
	using Father_t::cardinality;
	using Father_t::characteristic;
	using Father_t::write;
	using Father_t::read;

// matrix arithmetic functions

	Submatrix& add (Submatrix& C, const Submatrix& A, const Submatrix& B) const // C = A + B, where A, B, C have the same shape.
	{	C = B; return addin(C,A);	}
	Submatrix& neg (Submatrix& B, const Submatrix& A) const  // B = -A, where B and A have the same shape.
	{	B = A; return negin(B);		}
	Submatrix& sub (Submatrix& C, const Submatrix& A, const Submatrix& B) const // C = A - B, where A, B, C have the same shape.
	{	C = B; return subin(C,A);	}
	Submatrix& smul (Submatrix& B, const Scalar& a, const Submatrix& A) const // B = a*A.
	{	B = A; return smulin(B,a);	}
	Submatrix& saxpy (Submatrix& C, const Scalar& a, const Submatrix& A, const Submatrix& B) const // C = a*A + B.
	{	C = B; return saxpyin(C,a,A);	}
	Submatrix& mul (Submatrix& C, const Submatrix& A, const Submatrix& B) const // C = A*B, conformal shapes required.
	{	for (Index i = 0; i < C.rowdim(); ++i)
		for (Index j = 0; j < C.coldim(); ++j)
		{	Scalar x,y,z; field().assign(x, field().zero); field().assign(y,field().zero); field().assign(z,field().zero);
			for (Index k = 0; k < B.coldim(); ++k)
				field().axpyin(x, A.getEntry(y,i,k), B.getEntry(z,k,j));
			C.setEntry(i,j,x);
		}
		return C;
	}
	Submatrix& inv (Submatrix& B, const Submatrix& A) const // B = A^{-1}
	{	B = A; return invin(B);		}
	/*?*/Submatrix& div (Submatrix& C, const Submatrix& A, const Submatrix& B) const // C = A/B, B nonsingular required.
	{	throw(LinboxError("PlainDomain - what is div?")); return C; }
	Submatrix& axpy( Submatrix& D, const Submatrix& C, const Submatrix& A, const Submatrix &B) const // D = C + A*B, conformal shapes required.
	{	D = C; return axpyin(D,A,B);		}
	Submatrix& addin (Submatrix& B, const Submatrix& A) const // B += A, where B and A have the same shape.
	{	for (Index i = 0; i < B.rowdim(); ++i)
		for (Index j = 0; j < B.coldim(); ++j)
		{	Scalar x = field().zero, y = field().zero;
			A.getEntry(x,i,j); B.getEntry(y,i,j);
			B.setEntry(i,j, field().addin(x, y));
		}
		return B;
	}
	Submatrix& negin (Submatrix& A) const // A = -A.
	{	for (Index i = 0; i < A.rowdim(); ++i)
		for (Index j = 0; j < A.coldim(); ++j)
		{	Scalar x = field().zero;
			A.getEntry(x,i,j);
			A.setEntry(i,j,field().negin(x));
		}
		return A;
	}
	Submatrix& subin (Submatrix& B, const Submatrix& A) const // B -= A, where B and A have the same shape.
	{	for (Index i = 0; i < B.rowdim(); ++i)
		for (Index j = 0; j < B.coldim(); ++j)
		{	Scalar x = field().zero, y = field().zero;
			A.getEntry(x,i,j);
			B.setEntry(i,j, field().subin(x, y));
		}
		return B;
	}
	Submatrix& smulin (Submatrix& B, const Scalar& a) const // B = aB
	{	for (Index i = 0; i < B.rowdim(); ++i)
		for (Index j = 0; j < B.coldim(); ++j)
		{	Scalar x; init(x);
			B.getEntry(x,i,j);
			B.setEntry(i,j, mulin(x, a));
		}
		// this could use a _private_ refEntry...
		return B;
	}
	Submatrix& saxpyin( Submatrix& A, const Scalar& a, const Submatrix &B) const // A += a*B, shapes must conform.
	{	for (Index i = 0; i < B.rowdim(); ++i)
		for (Index j = 0; j < B.coldim(); ++j)
		{	Scalar x = field().zero, y = field().zero;
			A.getEntry(x,i,j); B.getEntry(y,i,j);
			A.setEntry(i,j, field().axpyin(x, a, y));
		}
		// this could use a _private_ refEntry...
		return A;
	}
	Submatrix& mulin_left (Submatrix& A, const Submatrix& B) const // A = A*B, for square, same dim A and B.
	{	Matrix C(*this, B.rowdim(), A.coldim());
		mul(C,A,B); return A = C;
	}
	Submatrix& mulin_right (const Submatrix& A, Submatrix& B) const // B = A*B, for square, same dim A and B.
	{	Matrix C(*this, B.rowdim(), A.coldim());
		mul(C,A,B); return B = C;
	}
	/*?*/Submatrix& invin (Submatrix& A) const // A = A^{-1}
	{	throw(LinboxError("PlainDomain invin for nonsing Submatrix not yet impl.")); return A; }
	/*?*/Submatrix& divin (Submatrix& A, const Submatrix& B) const // A /= B, B nonsingular required.
	{	throw(LinboxError("PlainDomain what is divin?")); return A; }
	Submatrix& axpyin( Submatrix& C, const Submatrix& A, const Submatrix &B) const // C += A*B, shapes must conform.
	{	for (Index i = 0; i < C.rowdim(); ++i)
		for (Index j = 0; j < C.coldim(); ++j)
		{	Scalar x,y,z; field().assign(x, field().zero); field().assign(y,field().zero); field().assign(z,field().zero);
			C.getEntry(x,i,j);
			for (Index k = 0; k < B.coldim(); ++k)
				field().axpyin(x, A.getEntry(y,i,k), B.getEntry(z,k,j));
			C.setEntry(i,j,x);
		}
		return C;
	}
	//Submatrix& copy(Submatrix &dst, const Submatrix& src) // deep copy, same as assignment


	bool areEqual(const Submatrix& A, const Submatrix &B) const // A == B, same shape not required
	{
		Scalar a = field().zero, b = field().zero;
		for (size_t i = 0; i < A.rowdim(); ++i)
			for (size_t j = 0; j < A.coldim(); ++j)
			{	A.getEntry(a, i, j);
				B.getEntry(b, i, j);
				if (not field().areEqual(a, b) )
					return false;
			}
		return true;
	}

// simple write
	std::ostream& write(std::ostream& out, const Submatrix& A) const
	{
		Scalar a;
		//out << "%%MatrixMarketExtended array" << std::endl;
		out << A.rowdim() << " " << A.coldim() << std::endl;
		for (size_t i = 0; i < A.rowdim(); ++i){
			for (size_t j = 0; j < A.coldim(); ++j)
				write(out, A.getEntry(a, i, j)) << " ";
			out << std::endl;
		}
		return out << std::endl;
	}

// simple read.  should use matrix reader
	std::istream& read(std::istream& in, Submatrix& A) const
	{
		Scalar a;
		size_t r, c;
		in >> r >> c;
		A.init(r, c);
		for (size_t i = 0; i < A.rowdim(); ++i)
			for (size_t j = 0; j < A.coldim(); ++j){
				read(in, a);
				A.setEntry(i, j, a);
			}
		return in;
	}

}; // PlainDomain

} // LinBox
#endif // __LINBOX_plain_domain_h

