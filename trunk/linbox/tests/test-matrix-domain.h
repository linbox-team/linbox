/*
   Examples of domains
BlasMatrixDomain<Field>
PlainDomain<Field>
Sliced3Domain<F3_rep>
M4riDomain
...

   classes: Index Scalar Vector Matrix Fibb Field Polynomial
Index i,j,k
Scalar a,b,c
Vector u,v,w // Subvector
Matrix A,B,C // Submatrix
Fibb F,G
Field K
MatrixDomain MD
Polynomial P

   Functions to test:
MatrixDomain()
MatrixDomain(K)
MatrixDomain(MD)
MatrixDomain(p,e) // future
field()
   vector arith:
add(z,x,y) z = x + y
addin(x,y) x += y
sub(z,x,y) z = x - y
subin(x,y) x -= y
neg(z,x) z = -x
negin(x) x = -x
axpy(z,a,x,y) z = a*x + y
axpyin(a,x,y) y += a*x
mulladd(z,a,x,b,y) z = ax + by
mulladdin(a,x,b,y) x = ax + by
copy(x,y) x <-- y
swap(x,y) x <--> y
   matrix arith
add(C,A,B) C = A + B
addin(A,B) A += B
sub(C,A,B) C = A - B
subin(A,B) A -= B
neg(C,A) C = -A
negin(A) A = -A
mul(C,A,B) C = A*B
mulin_left(A,B) A = A*B (square)
mulin_right(A,B) B = A*B (square)
axpy(D,A,B,C) D = A*B + C
axpyin(C,A,B) C += A*B
axpy(D,a,B,C) D = a*B + C
axpyin(a,B,C) C += a*B
mulladd(D,b,C,a,A,B) D = bC + aAB
mulladdin(b,C,a,A,B) C = bC + aAB
inv(F,A) F = A^{dagger}
invin(F,A) F = A^{dagger}, A modified
inv(B,A) B = A^{dagger}
invin(B,A) B = A^{dagger}, A modified
 nullity?
div_left(C,A,B) C = A^{dagger}*B (if consistent, AC = B)
divin_left(A,B) B = A^{dagger}*B, and A modified
div_right(C,A,B) C = A*B^{dagger} (if consistent, CB = A)
divin_right(A,B) A = A*B^{dagger}, and B modified

areEqual(A, B)
copy(B,A) B <-- A
swap(B,A) B <--> A

left_solve(X,A,B) = div_left(X,A,B) X = A^{dagger}*B
right_solve(X,A,B) = div_right(X,A,B) = X = A*B^{dagger}
		Operand& left_solve (const Matrix& A, Operand& B) const
		Operand& right_solve (const Matrix& A, Operand& B) const

Matrix<MatrixDomain> could do the rest
rank(i,A)
rankin(i,A)

det(d,A)
detin(d,A)

minpoly(P,A)
charpoly(P,A)
charpoly(Plist,A) // factored form(?)

isZero(A)
	(? isIdentity, isUnit, isZeroDivisor)


	matrices set themselves, read and write themselves
		void setIdentity(Matrix & I)
		random: void setRandom(Matrix & I, opt args)
		void setZero(Matrix & I)
		inline std::ostream &write (std::ostream &os, const Matrix &A) const
		inline std::istream &read (std::istream &is, Matrix &A) const


*/
/*

	class VectorDomain : public virtual DotProductDomain<Field> {//, public virtual VectorDomainBase<Field> {
	propose types
		VectorDomain()
		void init(const Field& F) { this->_field = &F; }
		VectorDomain (const VectorDomain &VD) :
		VectorDomain &operator = (const VectorDomain &VD)
		using VectorDomainBase<Field>::field;
		inline std::ostream &write (std::ostream &os, const Vector &x) const
		inline std::istream &read (std::istream &is, Vector &x) const
		inline Vector1 &copy (Vector1 &res, const Vector2 &v) const
	x	inline Vector1 &copy (Vector1 &res, const Vector2 &v, size_t i, size_t len = 0) const
	propose subvector()

		inline bool areEqual (const Vector1 &v1, const Vector2 &v2) const
		inline bool isZero (const Vector &v) const
		inline Element &dot (Element &res, const Vector1 &v1, const Vector2 &v2) const
		inline Vector1 &add (Vector1 &res, const Vector2 &y, const Vector3 &x) const
		inline Vector1 &addin (Vector1 &y, const Vector2 &x) const
		inline Vector1 &sub (Vector1 &res, const Vector2 &y, const Vector3 &x) const
		inline Vector1 &subin (Vector1 &y, const Vector2 &x) const
		inline Vector1 &neg (Vector1 &res, const Vector2 &x) const
		inline Vector &negin (Vector &y) const
	?	inline Vector1 &mul (Vector1 &res, const Vector2 &x, const Element &a) const
	?	inline Vector &mulin (Vector &x, const Element &a) const
	propose smul, smulin
		inline Vector1 &axpy (Vector1 &res, const Element &a, const Vector2 &x, const Vector3 &y) const
		inline Vector1 &axpyin (Vector1 &y, const Element &a, const Vector2 &x) const
	propose muladd, muladdin

		inline void swap (Vector &v1, Vector &v2) const
	?	inline Vector &permute (Vector   &v, Iterator  P_start, Iterator  P_end) const
		VectorDomain (const Field &F) :
		Vector& random(Vector& v)
*/
/*
proposed categories
Field/Ring
VectorDomain
  class Vector
  class NewVector
MatrixDomain inherits from VectorDomain, thus contains vector functions
  class Matrix
  class NewMatrix
BlackBox<MatrixDomain>
  uses MatrixDomain::Matrix and MatrixDomain::Vector
Matrix<MatrixDomain>
  domain() returns the MatrixDomain
    domain().mul
  field() returns MD().field()
SparseMatrix<MatrixDomain>

 Those matrix classes know basic matrix functions rank, det, minpoly, etc
 as well as
*/

#ifndef LinBox_test_matrix_domain_h
#define LinBox_test_matrix_domain_h
#include "test-field.h"

template <class MDom>
bool testMatrixDomain(const MDom& MD, int n) {
	bool pass = true;

	typename MDom::Matrix A(MD,n,n), B(MD,n,n), C(MD,n,n), D(MD, n, n);
	A.random(); B.random();
	C.zero();
	typename MDom::Scalar a,b,c,d;
	MD.init(a); MD.init(b); MD.init(c); MD.init(d);

	MD.add(C, A, B);
	MD.sub(C, A, B);
	MD.neg(C, A);
	MD.smul(C, a, B);
	MD.saxpy(C, a, A, B);
	MD.mul(C, A, B);
	MD.axpy(D, A, B, C);

	MD.addin(A, B);
	MD.subin(A, B);
	MD.negin(A);
	MD.smulin(B, a);
	MD.saxpyin(B, a, A);
	MD.mulin_left(A, B);
	MD.mulin_right(A, B);
	MD.axpyin(C, A, B);
	std::cout << MD.cardinality() << std::endl;
	std::cout << MD.characteristic() << std::endl;

	pass = pass and runBasicRingTests(MD, "matrix domain", 1, false);
	pass = pass and runFieldTests<MDom>(MD, "matrix domain", 1,0, false);
	return pass;
}
#endif //LinBox_test_matrix_domain_h
