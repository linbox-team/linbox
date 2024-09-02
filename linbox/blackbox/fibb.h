/* linbox/blackbox/fibb.h
 * Copyright (C) 2015 bds for LinBox Team.  See linbox/COPYING.LESSER for License info.
 *
 * Written by bds
 */
#ifndef LB_FIBB_H
#define LB_FIBB_H

/*
FIBB: Fast Inverse BlackBox

The FIBB functions are those of a blackbox plus rank, det, and the solvers:
Solve, NullSpaceRandom, NullSpaceBasis.
The solvers have left and right forms.

THe FIBBs are Diagonal, Permutation, Triangular, and products of FIBBs in which one or both are nonsingular.
*/
#include "linbox/blackbox/bb.h"
namespace LinBox{

template <class Ring>
struct FIBB : public BB<Ring>
{
	using Field = Ring;
	using Element = typename Ring::Element;
	using ResizableMatrix = DenseMatrix<Field>;
	using Matrix = DenseMatrix<Field>;

//	virtual const Field& field() const = 0;

	virtual ~FIBB(){}

	virtual size_t& rank(size_t& r) const
	= 0;

	virtual Element& det( Element& d) const
	= 0;

	// solveRight and solveLeft
	/** @brief Y: AY = X, for this A.
		Solve nonsingular or consistent singular system.
		If it is consistent singular, an arbitrary solution is provided.
		X and Y must have identical shape.

		Note that Y+Z is a random sample of the solution space after
		{solveRight(Y, X); nullspaceRandomRight(Z);}.

		Behaviour is unspecified for inconsistent systems (see solveMP).
	*/
	virtual Matrix& solveRight(Matrix& Y, const Matrix& X) const
	= 0;
	/// Y: YA = X, for this A
	virtual Matrix& solveLeft(Matrix& Y, const Matrix& X) const
	= 0;

	/// N: AN = 0, each col random.
	virtual Matrix& nullspaceRandomRight(Matrix& N) const
	= 0;
	/// N: NA = 0, each row random.
	virtual Matrix& nullspaceRandomLeft(Matrix& N) const
	= 0;
	// this generic is virtual so that it may be specialized for performance

	// nullspaceBasisRight and nullspaceBasisLeft

	/** B: columns are a right nullspace basis for this A.

		B is resized and filled so that:
		(1) AB = 0, (2) Ax = 0 => exists y: x = By, and (3) B has full rank.
	*/
	virtual ResizableMatrix& nullspaceBasisRight(ResizableMatrix& B) const
	= 0;
	/// BA= 0 and xA = 0 => exists y: x = yB and B full rank.
	virtual ResizableMatrix& nullspaceBasisLeft(ResizableMatrix& B) const
	= 0;
}; // class FIBB

/// N: AN = 0, each col random.
template<class Field>
DenseMatrix<Field>& genericNullspaceRandomRight(DenseMatrix<Field>& N, const FIBB<Field>& A)
//DenseMatrix<Field, std::vector<typename Field::Element> >& genericNullspaceRandomRight(DenseMatrix<Field, std::vector<typename Field::Element> >& N, const FIBB<Field>& A)
{	typedef DenseMatrix<Field> ResizableMatrix;
	typedef DenseMatrix<Field> Matrix;
	ResizableMatrix Xb(A.field(), N.rowdim(), N.coldim());
	ResizableMatrix Yb(A.field(), A.rowdim(), N.coldim());
	Matrix X(Xb); X.random();
	Matrix Y(Yb);
	A.applyRight(Y, X); // Y = AX
	A.solveRight(N, Y); // AN = AX
	BlasMatrixDomain<Field> MD(A.field());
	return MD.subin(N, X);
}

/// N: NA = 0, each row random.
template<class Field>
DenseMatrix<Field>& genericNullspaceRandomLeft(DenseMatrix<Field>& N, const FIBB<Field>& A)
//DenseMatrix<Field, std::vector<typename Field::Element> >& genericNullspaceRandomLeft(DenseMatrix<Field, std::vector<typename Field::Element> >& N, const FIBB<Field>& A)
{	typedef DenseMatrix<Field> ResizableMatrix;
	typedef DenseMatrix<Field> Matrix;
	ResizableMatrix Xb(A.field(), N.rowdim(), N.coldim());
	ResizableMatrix Yb(A.field(), N.rowdim(), A.coldim());
	Matrix X(Xb); X.random();
	Matrix Y(Yb);
	A.applyLeft(Y, X); // Y = XA
	A.solveLeft(N, Y); // NA = XA
	BlasMatrixDomain<Field> MD(A.field());
	return MD.subin(N, X);
}

} // namespace LinBox
#endif // LB_FIBB_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
