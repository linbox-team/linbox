/* iinbox/blackbox/fibb-product.h
 * Copyright (C) 2015 bds for LinBox Team.  See linbox/COPYING.LESSER for License info.
 *
 * Written by bds
 */
// todo: possiby merge with compose blackbox
#ifndef LB_FIBBProduct_H
#define LB_FIBBProduct_H

#include "linbox/blackbox/fibb.h"

namespace LinBox{

template<class Field_> 
struct FIBBProduct : public FIBB<Field_> { // Fast Inverse BlackBox
	typedef Field_ Field;
	typedef FIBB<Field> Father_t;
	typedef FIBBProduct<Field> Self_t;
  protected: 
	const FIBB<Field>** factors_;
	size_t n_;
	bool alloc_; // true only if new used within construction of factors_
	bool alloc_members_; // true only if new used to construct members.
	// (all or nothing ownership of members)
	const FIBB<Field>& head() const { return *(factors_[0]); }
	FIBBProduct<Field>& tail (FIBBProduct<Field>& A) const
	{ A.factors_ = factors_ + 1; A.n_ = n_-1; 
	  A.alloc_ = A.alloc_members_ = false; 
	  return A;
	}
  public:
    using Element = typename Father_t::Element;
    using Matrix = typename Father_t::Matrix;
    using MotherMatrix = typename Father_t::MotherMatrix;

	/* Blackbox functions */
	BBType bbTag() const { return FIBBProductTag; }
	size_t rowdim() const;
	size_t coldim() const;
	const Field& field() const;
	Matrix& applyRight(Matrix& Y, const Matrix& X) const; // Y = AX
	Matrix& applyLeft(Matrix& Y, const Matrix& X) const; // Y = XA
	// todo: rebind - hom support
	std::istream& read(std::istream& is);
	std::ostream& write(std::ostream& os) const;

	/* fibb functions */
	size_t& rank(size_t& r) const;

	Element& det(Element& d) const;

	Matrix& solveRight(Matrix& Y, const Matrix& X) const; 
	/// Y: YA = X, for this A
	Matrix& solveLeft(Matrix& Y, const Matrix& X) const; 

	/// N: AN = 0, each col random.
	Matrix& nullspaceRandomRight(Matrix& N) const; 
	/// N: NA = 0, each row random.
	Matrix& nullspaceRandomLeft(Matrix& N) const; 

	/* nullspaceBasisRight and nullspaceBasisLeft */
	/** B: columns are a right nullspace basis for this A.
		
		B is resized and filled so that:
		(1) AB = 0, (2) Ax = 0 => exists y: x = By, and (3) B has full rank.
	*/
	MotherMatrix& nullspaceBasisRight(MotherMatrix& B) const; 
	/// BA= 0 and xA = 0 => exists y: x = yB and B full rank.
	MotherMatrix& nullspaceBasisLeft(MotherMatrix& B) const; 

	/* cstors, dstor, initializers */
	FIBBProduct();
	FIBBProduct(const Field& F);
//	FIBBProduct(const FIBBProduct<Field>& A);
	FIBBProduct(const FIBB<Field>& A1, const FIBB<Field>& A2);
	FIBBProduct(const FIBB<Field>& A1, const FIBB<Field>& A2, 
				const FIBB<Field>& A3);
	FIBBProduct(const FIBB<Field>& A1, const FIBB<Field>& A2, 
				const FIBB<Field>& A3, const FIBB<Field>& A4);
	FIBBProduct(const FIBB<Field>& A1, const FIBB<Field>& A2, 
				const FIBB<Field>& A3, const FIBB<Field>& A4, 
				const FIBB<Field>& A5);
	FIBBProduct(const FIBB<Field>& A1, const FIBB<Field>& A2, 
				const FIBB<Field>& A3, const FIBB<Field>& A4, 
				const FIBB<Field>& A5, const FIBB<Field>& A6);
	~FIBBProduct();
	FIBBProduct& init();
	FIBBProduct& init(const FIBB<Field>& A1);
	FIBBProduct& init(const FIBB<Field>& A1, const FIBB<Field>& A2);
	FIBBProduct& init(const FIBB<Field>& A1, const FIBB<Field>& A2, 
					  const FIBB<Field>& A3);
	FIBBProduct& init(const FIBB<Field>& A1, const FIBB<Field>& A2, 
					  const FIBB<Field>& A3, const FIBB<Field>& A4);
	FIBBProduct& init(const FIBB<Field>& A1, const FIBB<Field>& A2, 
					  const FIBB<Field>& A3, const FIBB<Field>& A4, 
					  const FIBB<Field>& A5);
	FIBBProduct& init(const FIBB<Field>& A1, const FIBB<Field>& A2, 
					  const FIBB<Field>& A3, const FIBB<Field>& A4, 
					  const FIBB<Field>& A5, const FIBB<Field>& A6);
    // product and owner of heap allocated FIBBs
	FIBBProduct& incorporate(
			const FIBB<Field>* A1 = NULL, const FIBB<Field>* A2 = NULL, 
			const FIBB<Field>* A3 = NULL, const FIBB<Field>* A4 = NULL, 
			const FIBB<Field>* A5 = NULL, const FIBB<Field>* A6 = NULL);
	const FIBB<Field>& operator[](size_t i) const;

}; // class FIBBProduct
}// namespace LinBox
#endif // LB_FIBBProduct_H

// blackbox/FIBBProduct.inl
#ifndef LB_FIBBProduct_INL
#define LB_FIBBProduct_INL
namespace LinBox {

// Blackbox interface
template<class Field> size_t FIBBProduct<Field>:: 
rowdim() const { return n_ > 0 ? factors_[0]->rowdim() : 0; }

template<class Field> size_t FIBBProduct<Field>:: 
coldim() const { return n_ > 0 ? factors_[n_-1]->coldim() : 0; }

template<class Field> const Field& FIBBProduct<Field>:: 
field() const { return head().field(); }

template<class Field> typename FIBBProduct<Field>::Matrix& FIBBProduct<Field>:: 
applyRight(typename FIBBProduct<Field>::Matrix& Y, const typename FIBBProduct<Field>::Matrix& X) const
{	if (n_==0) return Y;
	if (n_==1) return head().applyRight(Y,X);
	const FIBB<Field>& A = head();
	FIBBProduct<Field> B; tail(B);
	MotherMatrix X1b(field(), B.rowdim(), X.coldim());
	Matrix X1(X1b);
	B.applyRight(X1, X);
	A.applyRight(Y, X1);
	return Y;
}

template<class Field> typename FIBBProduct<Field>::Matrix& FIBBProduct<Field>:: 
applyLeft(typename FIBBProduct<Field>::Matrix& Y, const typename FIBBProduct<Field>::Matrix& X) const
{	if (n_==0) return Y;
	if (n_==1) return head().applyLeft(Y,X);
	const FIBB<Field>& A = head();
	FIBBProduct<Field> B; tail(B);
	MotherMatrix X1b(field(), X.rowdim(), B.coldim());
	Matrix X1(X1b);
	A.applyLeft(X1, X);
	B.applyLeft(Y, X1);
  	return Y;
}

/* rebind - hom support
*/
// todo: read
template<class Field> std::istream& FIBBProduct<Field>:: 
read(std::istream& is) 
{ return is; }

template<class Field> std::ostream& FIBBProduct<Field>:: 
write(std::ostream& os) const
{	os << std::endl << "%%MatrixMarket matrix composite integer general" << std::endl;
	field().write(os << "% written by LinBox::FIBBProduct< ") << " >" << std::endl;
	os << n_ << " factors:";
	for (size_t i = 0; i < n_; ++i) 
		factors_[i]->write(os << std::endl);
	return os;
}

template<class Field> size_t& FIBBProduct<Field>:: 
rank( size_t& r ) const
{	if (n_ == 0) return r = 0;
	factors_[0]->rank(r);
	size_t s;
	for (size_t i = 1; i < n_; ++i) 
		r = std::min(r, factors_[i]->rank(s));
	return r;
}

template<class Field> typename FIBBProduct<Field>::Element& FIBBProduct<Field>:: 
det( typename FIBBProduct<Field>::Element& d ) const
{	if (n_==0) return field().assign(d, field().one);
	factors_[0]->det(d); 
	typename Field::Element e; field().init(e);
	for (size_t i = 1; i < n_; ++i) 
		field().mulin(d, factors_[i]->det(e));
	return d;
}

template<class Field> typename FIBBProduct<Field>::Matrix& FIBBProduct<Field>:: 
solveRight( typename FIBBProduct<Field>::Matrix& Y, const typename FIBBProduct<Field>::Matrix& X ) const
{	if (n_==0) return Y;
	if (n_==1) return head().solveRight(Y,X);
	const FIBB<Field>& A = head();
	FIBBProduct<Field> B(field()); tail(B);
	MotherMatrix Zb(field(), A.coldim(), X.coldim());
	Matrix Z(Zb);
	A.solveRight(Z,X); // A1*Z = X
	return B.solveRight(Y,Z); // A2*Y = Z
}

template<class Field> typename FIBBProduct<Field>::Matrix& FIBBProduct<Field>:: 
solveLeft( typename FIBBProduct<Field>::Matrix& Y, const typename FIBBProduct<Field>::Matrix& X ) const
{	if (n_==0) return Y;
	if (n_==1) return head().solveLeft(Y,X);
	const FIBB<Field>& A = head();
	FIBBProduct<Field> B(field()); tail(B);
	MotherMatrix Zb(field(), X.rowdim(), A.coldim()); 
	Matrix Z(Zb);
	B.solveLeft(Z,X); // Z*A2 = X
	return A.solveLeft(Y,Z); // Y*A1 = Z
}

// Randomly fill N such that ABN = 0, where this is AB.
template<class Field> typename FIBBProduct<Field>::Matrix& FIBBProduct<Field>:: 
nullspaceRandomRight( typename FIBBProduct<Field>::Matrix& N ) const 
{
 	if (n_==0) return N;
	if (n_==1) return head().nullspaceRandomRight(N);
	size_t r;
	const FIBB<Field>& A = head();
	FIBBProduct<Field> B(field()); tail(B);
	if (A.rowdim() == A.coldim() and A.rank(r) == A.coldim())
		return B.nullspaceRandomRight(N);
	else
	{	MotherMatrix N1b(field(), A.coldim(), N.coldim());
		Matrix N1(N1b);
		A.nullspaceRandomRight(N1);
		B.solveRight(N,N1);
		return N;
		// a solveRightin would be good if B is a perm.
	}
}

// Randomly fill N such that NAB = 0, where this is AB.
template<class Field> typename FIBBProduct<Field>::Matrix& FIBBProduct<Field>:: 
nullspaceRandomLeft( typename FIBBProduct<Field>::Matrix& N ) const
{	if (n_==0) return N;
	if (n_==1) return head().nullspaceRandomLeft(N);
	size_t r;
	const FIBB<Field>& A = head();
	FIBBProduct<Field> B(field()); tail(B);
	if (B.rowdim() == B.coldim() and B.rank(r) == B.coldim())
		return A.nullspaceRandomLeft(N);
	else
	{	MotherMatrix N1b(field(), N.rowdim(), B.rowdim());
		Matrix N1(N1b);
		B.nullspaceRandomLeft(N1);
		return A.solveLeft(N,N1);
		// a solveLeftin would be good if A is a perm.
	}
}

template<class Field> typename FIBBProduct<Field>::MotherMatrix& FIBBProduct<Field>:: 
nullspaceBasisRight( typename FIBBProduct<Field>::MotherMatrix& N ) const
{	if (n_==0) { N.resize(0,0); return N; }
	if (n_==1) return head().nullspaceBasisRight(N);
	size_t r;
	const FIBB<Field>& A = head();
	FIBBProduct<Field> B(field()); tail(B);
	if (A.rowdim() == A.coldim() and A.rank(r) == A.rowdim())
	 	B.nullspaceBasisRight(N);
	else 
	{	MotherMatrix N1(field());
		A.nullspaceBasisRight(N1);
		N.resize(N1.rowdim(), N1.coldim());
		Matrix Ns(N), N1s(N1);
		B.solveRight(Ns, N1s);
	}
	return N;
}

template<class Field> typename FIBBProduct<Field>::MotherMatrix& FIBBProduct<Field>:: 
nullspaceBasisLeft( typename FIBBProduct<Field>::MotherMatrix& N ) const
{	if (n_==0) { N.resize(0,0); return N; }
	if (n_==1) return head().nullspaceBasisLeft(N);
	size_t r;
	const FIBB<Field>& A = head();

	FIBBProduct<Field> B(field()); tail(B);
	if (B.rowdim() == B.coldim() and B.rank(r) == B.rowdim())
	 	A.nullspaceBasisLeft(N);
	else 
	{	MotherMatrix N1(field());
		B.nullspaceBasisLeft(N1);
		N.resize(N1.rowdim(), N1.coldim());
		Matrix Ns(N), N1s(N1);
		A.solveLeft(Ns, N1s);
	}
	return N;
}

/* cstors, dstor */
template<class Field> FIBBProduct<Field>::
FIBBProduct() 
: factors_(0), n_(0), alloc_(false), alloc_members_(false) 
{}

template<class Field> FIBBProduct<Field>::
FIBBProduct(const Field& F) 
: factors_(0), n_(0), alloc_(false), alloc_members_(false) 
{}

//template<class Field> FIBBProduct<Field>::
//FIBBProduct(const FIBBProduct<Field>& A): Ap(A.Ap), Bp(A.Bp), allocA(A.allocA), allocB(A.allocB)
//{}

template<class Field> FIBBProduct<Field>:: 
FIBBProduct( const FIBB<Field>& A1, const FIBB<Field>& A2 ) 
: factors_(0), n_(0), alloc_(false), alloc_members_(false) 
{ init(A1, A2); } 

template<class Field> FIBBProduct<Field>:: 
FIBBProduct( const FIBB<Field>& A1, const FIBB<Field>& A2, 
			 const FIBB<Field>& A3 ) 
: factors_(0), n_(0), alloc_(false), alloc_members_(false) 
{ init(A1, A2, A3); }

template<class Field> FIBBProduct<Field>:: 
FIBBProduct( const FIBB<Field>& A1, const FIBB<Field>& A2, 
			 const FIBB<Field>& A3, const FIBB<Field>& A4 ) 
: factors_(0), n_(0), alloc_(false), alloc_members_(false) 
{ init(A1, A2, A3, A4); }

template<class Field> FIBBProduct<Field>:: 
FIBBProduct( const FIBB<Field>& A1, const FIBB<Field>& A2, 
			 const FIBB<Field>& A3, const FIBB<Field>& A4, 
			 const FIBB<Field>& A5 ) 
: factors_(0), n_(0), alloc_(false), alloc_members_(false) 
{ init(A1, A2, A3, A4, A5); }

template<class Field> FIBBProduct<Field>:: 
FIBBProduct( const FIBB<Field>& A1, const FIBB<Field>& A2, 
			 const FIBB<Field>& A3, const FIBB<Field>& A4, 
			 const FIBB<Field>& A5, const FIBB<Field>& A6 ) 
: factors_(0), n_(0), alloc_(false), alloc_members_(false) 
{ init(A1, A2, A3, A4, A5, A6); }

template<class Field> FIBBProduct<Field>:: 
~FIBBProduct()
{	if (alloc_members_) 
		for(size_t i = 0; i < n_; ++i) delete factors_[i];
	if (alloc_) delete[] factors_; 
}

/* initializers */

template<class Field> FIBBProduct<Field>& FIBBProduct<Field>:: 
init( )
{	 if (alloc_ and factors_) delete[] factors_;
	factors_ = NULL; n_ = 0; 
	alloc_ = false; alloc_members_ = false;
	return *this;
}

template<class Field> FIBBProduct<Field>& FIBBProduct<Field>:: 
init( const FIBB<Field>& A )
{	if (alloc_ and factors_) delete[] factors_;
	factors_ = new const FIBB<Field>*[1]; n_ = 1; 
	alloc_ = true; alloc_members_ = false;
	factors_[0] = &A;
	return *this;
}

template<class Field> FIBBProduct<Field>& FIBBProduct<Field>:: 
init( const FIBB<Field>& A1, const FIBB<Field>& A2 )
{	if (alloc_ and factors_) delete[] factors_;
	factors_ = new const FIBB<Field>*[2]; n_ = 2; 
	alloc_ = true; alloc_members_ = false;
	factors_[0] = &A1; factors_[1] = &A2; 
	return *this;
}

template<class Field> FIBBProduct<Field>& FIBBProduct<Field>:: 
init( const FIBB<Field>& A1, const FIBB<Field>& A2, 
	  const FIBB<Field>& A3 )
{	if (alloc_ and factors_) delete[] factors_;
	factors_ = new const FIBB<Field>*[3]; n_ = 3; 
	alloc_ = true; alloc_members_ = false;
	factors_[0] = &A1; factors_[1] = &A2; 
	factors_[2] = &A3; 
	return *this;
}

template<class Field> FIBBProduct<Field>& FIBBProduct<Field>:: 
init( const FIBB<Field>& A1, const FIBB<Field>& A2, 
	  const FIBB<Field>& A3, const FIBB<Field>& A4 )
{	if (alloc_ and factors_) delete[] factors_;
	factors_ = new const FIBB<Field>*[4]; n_ = 4; 
	alloc_ = true; alloc_members_ = false;
	factors_[0] = &A1; factors_[1] = &A2; 
	factors_[2] = &A3; factors_[3] = &A4; 
	return *this;
}

template<class Field> FIBBProduct<Field>& FIBBProduct<Field>:: 
init( const FIBB<Field>& A1, const FIBB<Field>& A2, 
	  const FIBB<Field>& A3, const FIBB<Field>& A4, 
  	  const FIBB<Field>& A5 )
{	if (alloc_ and factors_) delete[] factors_;
	factors_ = new const FIBB<Field>*[5]; n_ = 5; 
	alloc_ = true; alloc_members_ = false;
	factors_[0] = &A1; factors_[1] = &A2; 
	factors_[2] = &A3; factors_[3] = &A4; 
	factors_[4] = &A5;
	return *this;
}

template<class Field> FIBBProduct<Field>& FIBBProduct<Field>:: 
init( const FIBB<Field>& A1, const FIBB<Field>& A2, 
	  const FIBB<Field>& A3, const FIBB<Field>& A4, 
  	  const FIBB<Field>& A5, const FIBB<Field>& A6 )
{	if (alloc_ and factors_) delete[] factors_;
	factors_ = new const FIBB<Field>*[6]; n_ = 6; 
	alloc_ = true; alloc_members_ = false;
	factors_[0] = &A1; factors_[1] = &A2; 
	factors_[2] = &A3; factors_[3] = &A4; 
	factors_[4] = &A5; factors_[5] = &A6;
	return *this;
}

template<class Field> FIBBProduct<Field>& FIBBProduct<Field>:: 
incorporate(const FIBB<Field>* A1, const FIBB<Field>* A2, 
			const FIBB<Field>* A3, const FIBB<Field>* A4, 
			const FIBB<Field>* A5, const FIBB<Field>* A6)
{
	if (A1 == NULL) init();
	else if (A2 == NULL) init(*A1);
	else if (A3 == NULL) init(*A1,*A2);
	else if (A4 == NULL) init(*A1,*A2,*A3);
	else if (A5 == NULL) init(*A1,*A2,*A3,*A4);
	else if (A6 == NULL) init(*A1,*A2,*A3,*A4,*A5);
	else 			init(*A1,*A2,*A3,*A4,*A5,*A6);
	alloc_members_ = true;
	return *this;
}

template<class Field> const FIBB<Field>& FIBBProduct<Field>::
operator[](size_t i) const
{ return *(factors_[i]); }

}// namespace LinBox

#endif // LB_FIBBProduct_INL
