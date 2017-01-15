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
	using MotherMatrix = DenseMatrix<Ring>;
	using Matrix = typename BB<Ring>::Matrix;
	using BB<Ring>::applyLeft;
	using BB<Ring>::applyRight;

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

	/// N: Each col random subject only to AN = 0..
	virtual Matrix& nullspaceRandomRight(Matrix& N) const 
	= 0;
	/// N: Each row random subject only to NA = 0..
	virtual Matrix& nullspaceRandomLeft(Matrix& N) const 
	= 0;
	// this generic is virtual so that it may be specialized for performance

	// nullspaceBasisRight and nullspaceBasisLeft

	/** B: columns are a right nullspace basis for this A.
		
		B is resized and filled so that:
		(1) AB = 0, (2) Ax = 0 => exists y: x = By, and (3) B has full rank.
	*/
	virtual MotherMatrix& nullspaceBasisRight(MotherMatrix& B) const 
	= 0;
	/// BA= 0 and xA = 0 => exists y: x = yB and B full rank.
	virtual MotherMatrix& nullspaceBasisLeft(MotherMatrix& B) const 
	= 0;
}; // class FIBB

/** N: AN = 0, each col random.
  This is a tool Fibb classes can use unless a specialized method is desired.
*/
template<class Field>
typename FIBB<Field>::Matrix& genericNullspaceRandomRight(
		typename FIBB<Field>::Matrix& N, 
		const FIBB<Field>& A)
{	typedef typename FIBB<Field>::Matrix Matrix;
 	typedef typename FIBB<Field>::MotherMatrix MotherMatrix;
	MotherMatrix Xb(A.field(), N.rowdim(), N.coldim());
	MotherMatrix Yb(A.field(), A.rowdim(), N.coldim());
	Matrix X(Xb); X.random();
	Matrix Y(Yb); 
	A.applyRight(Y, X); // Y = AX
	A.solveRight(N, Y); // AN = AX
	BlasMatrixDomain<Field> MD(A.field());
	return MD.subin(N, X);
}

/** N: AN = 0, each row random.
  This is a tool Fibb classes can use unless a specialized method is desired.
*/
template<class Field>
typename FIBB<Field>::Matrix& genericNullspaceRandomLeft(
		typename FIBB<Field>::Matrix& N, 
		const FIBB<Field>& A)
{	typedef DenseMatrix<Field> MotherMatrix;
	typedef DenseSubmatrix<Field> Matrix;
	MotherMatrix Xb(A.field(), N.rowdim(), N.coldim());
	MotherMatrix Yb(A.field(), N.rowdim(), A.coldim());
	Matrix X(Xb); X.random();
	Matrix Y(Yb); 
	A.applyLeft(Y, X); // Y = XA
	A.solveLeft(N, Y); // NA = XA
	BlasMatrixDomain<Field> MD(A.field());
	return MD.subin(N, X);
}

} // namespace LinBox
#endif // LB_FIBB_H
