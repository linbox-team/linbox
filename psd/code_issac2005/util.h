#ifndef __UTIL_H__
#define __UTIL_H__

#include <linbox/integer.h>
#include <linbox/blackbox/sparse.h>
#include <linbox/algorithms/matrix-mod.h>
#include <linbox/field/gmp-rational.h>
#include "lie-matrix.h"
// conversion to C++ using LinBox rationals by bds and zw

namespace LinBox {
// d <= 2(v1 . v2)/ (v2. v2)
template <class Field, class V1, class V2>
void chdot (typename Field::Element& d, const V1& v1, const V2& v2, const Field& f);

// r list of chdot[lambda,root[j]], j_l corresponding index j
template <class Vect, class JV, class Field, class NV, class LRV, class RV>
void compute_r (Vect& r, JV& j_l,
				const NV& nu, const LRV& root, const RV& rho, const Field& f);

/* compute the lcm of denom of all entries of a rational matrix*/ 
/* error probability \approx 2^{-n_try + 1} */
template <class IMatrix>
integer& lcmDenEntries (integer& l, const IMatrix& A, int n_try = 20);

/* scalar <- lcm (denominator of all entries of M) */
template <class Ring>
integer& lcmDenEntries(integer& scalar, const SparseMatrix<Ring>& M);

/* gcd of numerators of all entries */
template <class GMatrix>
void gcdNumEntries (integer& g, const GMatrix& M, int n_try = 20);

/*
template <class Ring>
void gcdNumEntries (integer& g, const SparseMatrix<Ring>& M, int n_try);
*/

/* RM <- scalar * QM */
/*
template <class RMatrix, class QMatrix>
void scalarMul (RMatrix& RM, const QMatrix& QM, const integer& scalar);
*/

/* M <- a * M + b * I */
template <class ZMatrix, class Scalar>
void addIdentity (ZMatrix& M, const Scalar& a, const Scalar& b);

/* M <- D * M, D is a diagonal matrix, represented by a vector v*/
template <class Matrix, class Vector>
void diagonalMulIn (Matrix& M, const Vector& v);

/* test if  a matrix is a diagonal matrix*/
template <class Matrix>
bool isDiagonal (const Matrix& M, int n_try = 1);

template <class Ring>
bool isDiagonal (const SparseMatrix<Ring>& M);

template <class Ops, class Ring>
void optimize(Ops& op, const Ring& R);

template <class QOPVect, class JVect, class QRVect, class QOps>
void buildLieMatrixGMP(LieMatrix<GMPRationalField, SparseMatrix<GMPRationalField> >*& M,
					   const QOPVect& Qop, const JVect& op_i, 
			  		   const QRVect& Qr, QOps* const QF);

struct Signature{size_t pos; size_t zero; size_t neg;};

std::ostream& operator<<(std::ostream& out, const Signature& s) {
	out << "p_z_n[" <<s. pos <<", " <<s. zero <<", " <<s. neg <<"]";
	return out;
}

template <class Ring>
int nonZeroEntries(const SparseMatrix<Ring>& M);

template <class Ring, class Blackbox>
int nonZeroEntries(const LieMatrix<Ring, Blackbox>& M);

template <class Poly, class Realring>
size_t alternations(const Poly& p, const Realring& R);

template <class Poly, class Realring>
Signature& signature(Signature& s, const Poly& p, Realring R);

template <class Vector>
bool isAllPositive (const Vector& v);

template <class Vector>
bool isAlternativeSign (const Vector& v);

template <class Matrix, class MinPoly>
bool check_minpoly(const Matrix& A, const MinPoly& m, int n_try = 1);
}

#include "util.inl"
#endif
