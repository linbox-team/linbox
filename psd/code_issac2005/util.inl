#ifndef __UTIL_INL__
#define __UTIL_INL__

#include "util.h"
#include <linbox/vector/vector-domain.h>

namespace LinBox{

template <class Field, class V1, class V2>
void chdot (typename Field::Element& d, const V1& v1, const V2& v2, const Field& f) {
	VectorDomain <Field> vd (f);
	typename Field::Element len, in;
	vd. dot (in, v1, v2);
	vd. dot (len, v2, v2);
	f. init (d, 2);
	f. mulin (d, in);
	f. divin (d, len);
}

template <class Vect, class JV, class Field, class NV, class LRV, class RV>
void compute_r (Vect& r, JV& j_l,
				const NV& nu, const LRV& root, const RV& rho, const Field& f) {

	typedef typename Field::Element Element;
	NV lambda = nu;
	RV weight = rho;
	int rank = root. size ();
	bool test =true;
	LinBox::VectorDomain<Field> vd(f);

	r. clear ();
	j_l. clear();

	while (test) {
		test=false;
 		int j = 0;
		while (j < rank ) {
			Element d;
			vd. dot (d, root[j], weight);
			if (f. sign(d) > 0) { //??
				test=true;
				chdot (d, lambda, root[j], f);				
				r. push_back (d);
				j_l. push_back (j);
				f. negin (d);
				vd. axpyin (lambda, d, root[j]);
				chdot (d, weight, root[j], f);				
				f. negin (d);
				vd. axpyin (weight, d, root[j]);
				j = rank;
				
			}
			++ j;
		}
	}
}

template <class GMatrix>
integer& lcmDenEntries (integer& l, const GMatrix& A, int n_try) {

	typedef typename GMatrix::Field Ring;
	typedef typename Ring::Element Rational;
	l = 1; Ring Q(A. field()); integer tmp;
	int n = A. rowdim(); int m = A. coldim();
	std::vector<Rational> x(n), y(m);
	typename std::vector<Rational>::iterator x_p, y_p;
	for (int i = 0; i < n_try; ++ i) {
		for (y_p = y. begin(); y_p != y. end(); ++ y_p)
			Q. init (*y_p, rand());
		A. apply (x, y);
		for (x_p = x. begin(); x_p != x. end(); ++ x_p) {
			Q. get_den (tmp, *x_p);
			lcm (l, l, tmp);
		}
	}
			
	return l;
}

template <class Ring>
integer& lcmDenEntries(integer& scalar, const SparseMatrix<Ring>& M) {
	typedef SparseMatrix<Ring> Matrix;
	typename Matrix::ConstRawIterator raw_p;
	typename Matrix::Field Q(M. field());
	scalar = 1; integer den;
	for (raw_p = M. rawBegin(); raw_p != M. rawEnd(); ++ raw_p) {
		Q.get_den (den, *raw_p);
		lcm (scalar, scalar, den);
	}

	return scalar;
}

template <class LMatrix>
void gcdNumEntries (integer& g, const LMatrix& M, int n_try) {
		
	typedef typename LMatrix::Field Ring;
	typedef typename Ring::Element  Element;
	typedef std::vector<Element> R_Vect;
	typedef typename R_Vect::iterator RIterator;

	Ring R(M. field());
	R_Vect x(M. coldim()), y (M. rowdim());
	RIterator x_p, y_p;

	g = 0; integer tmp;
	for (int i = 0; i < n_try; ++ i) {

	 	for (x_p = x. begin(); x_p != x. end(); ++ x_p)
			R. init (*x_p, rand());

		M. apply (y, x);
		for (y_p = y. begin(); y_p != y. end(); ++ y_p) {
			R. get_num (tmp, *y_p);
			gcd (g, g, tmp);
		}

		if (g == 1) return;
	}
}

template <class ZMatrix, class Scalar>
void addIdentity (ZMatrix& M, const Scalar& a, const Scalar& b) {

	typename ZMatrix::Field::Element ref;
	typename ZMatrix::Field F(M. field());
	typename ZMatrix::RawIterator raw_p;
	for (raw_p = M. rawBegin(); raw_p != M. rawEnd(); ++ raw_p)
		F. mulin (*raw_p, a);

	int dim = M. rowdim();
	for (int i = 0; i < dim; ++ i) 
		F. addin (M.refEntry(i, i), b);
}

template <class Matrix, class Vector>
void diagonalMulIn (Matrix& M, const Vector& v) {
	typedef typename Matrix::Field Ring;
	typedef typename Ring::Element Element;
	typedef std::vector<Element> R_Vect;
	Ring R(M. field());
	VectorDomain<Ring> RVD(R);

	typename Matrix::RowIterator row_p;
	typename Vector::const_iterator v_p;
	for (row_p = M. rowBegin(), v_p = v. begin();
		 row_p != M. rowEnd(); ++ row_p, ++ v_p)
		
		RVD. mulin (*row_p, *v_p);
}
			
template <class Matrix>
bool isDiagonal (const Matrix& M, int n_try = 1) {
	typedef typename Matrix::Field Ring;
	typedef typename Ring::Element Element;
	typedef std::vector<Element> R_Vect;
	Ring R(M. field());
	VectorDomain<Ring> RVD(R);
	Element one; R. init (one, 1);
	int dim = M. rowdim();
	if (dim != M. coldim()) return false;
	R_Vect d(dim), x(dim), Mx(dim), dx(dim);
	typename R_Vect::iterator x_p, d_p, dx_p;
	for (x_p = x. begin(); x_p != x. end(); ++ x_p)
		R. assign (*x_p, one);
	M. apply (d, x);
	
	for (int i = 0; i < n_try; ++ i) {
		for (x_p = x. begin(); x_p != x. end(); ++ x_p)
			R. init (x, rand());
		M. apply (Mx, x);
		for (x_p = x. begin(), d_p = d. begin(), dx_p = dx. begin();
			 x_p != x. end(); ++ x_p, ++ d_p, ++ dx_p)
			 R. mul (*dx_p, *x_p, *d_p);
		
		if (!areEqual (dx, Mx)) return false;
	}

	return true;
}

template <class Ring>
bool isDiagonal (const SparseMatrix<Ring>& M) {
	typename SparseMatrix<Ring>::ConstRawIndexedIterator raw_p;
	for (raw_p = M. rawIndexedBegin(); raw_p != M. rawIndexedEnd(); ++ raw_p)
	if (raw_p. rowIndex() != raw_p. colIndex())
		return false;

	return true;
}

template <class Ops, class Ring>
void optimize(Ops& op, const Ring& R) {

	int rest = op. size();
	int str = op. size() -1;
	typedef typename Ring::Element Element;
	typedef std::vector<Element> R_Vect;
	
	int cur, i;
	while (rest > 1) {
		cur = str;
		if (isDiagonal(*(op[cur]))) {
			R_Vect d(op[cur] -> rowdim());
			typename R_Vect::iterator d_p; int i;
			for (i = 0, d_p = d. begin(); 
				 i < op[cur] -> rowdim(); ++ i, ++ d_p)
				op[cur] -> getEntry (*d_p, i, i);
			diagonalMulIn(*(op[cur - 1]), d);
			delete op[cur];
			op. erase (op. begin() + cur);
		}

		-- str;
		-- rest;
	}
}

template <class QOPVect, class JVect, class QRVect, class QOps>
void buildLieMatrixGMP(LieMatrix<GMPRationalField, SparseMatrix<GMPRationalField> >*& M,
		  			   const QOPVect& Qop, const JVect& op_i, 
		  			   const QRVect& Qr, QOps* const QF) {

	typedef GMPRationalField Ring;
	typedef Ring::Element Element;
	typedef QOps QOperator;
	typedef typename QOperator::Field Rationals;
	typedef QOPVect QOList;
	typedef SparseMatrix<Ring> ROperator;
	typedef std::vector<ROperator*> ROList;
	typedef std::vector<Element> RVector;
	typedef LieMatrix<Ring, SparseMatrix<Ring> > LMatrix;
	Rationals Q(QF -> field()); Ring R;

	int length = Qr. size(), dim = Qop. front() -> rowdim();
	int nop = Qop. size();
	typename QOPVect::const_iterator Qop_p; ROList::iterator Zop_p;
	typename QRVect::const_iterator Qr_p;
	typename JVect::const_iterator op_i_ptr;

	ROList ops(length + 1); 
	ops. front() = new ROperator (*QF);

	ROList::iterator ops_p; Element tmp1, tmp2, one; R. init (one, 1);
	for (ops_p = ops. begin() + 1, op_i_ptr = op_i. begin(),
	     Qr_p = Qr. begin(); ops_p != ops. end(); 
		 ++ ops_p, ++ op_i_ptr, ++ Qr_p) {

		 R. add (tmp2, *Qr_p, one);
		 R. invin (tmp2);
		 R. mul (tmp1, *Qr_p, tmp2);

		 *ops_p = new ROperator(*(Qop[*op_i_ptr]));
		 addIdentity (**ops_p, tmp1, tmp2);
	}

	/*
	std::cout << "Operators: \n";
	for (ops_p = ops. begin(); ops_p != ops. end(); ++ ops_p) {
		std::cout << ops_p - ops. begin() << "th operator\n";
		(*ops_p) -> write (std::cout);
	}
	*/
	std::clog << "Before optimization, chain length. " << ops. size();
	std::clog << std::endl;
	optimize (ops, R);
	std::clog << "After optimization, chain length. " << ops. size();
	std::clog << std::endl;
	/*
	std::cout << "Operators: \n";
	for (ops_p = ops. begin(); ops_p != ops. end(); ++ ops_p) {
		std::cout << ops_p - ops. begin() << "th operator\n";
		(*ops_p) -> write (std::cout);
	}
	*/

	M = new LMatrix(ops, R);
}

template <class Vector>
bool isAlternativeSign (const Vector& v) {

	if (v. size() < 1) return true;

	typename Vector::const_iterator p;
	int cur_sign = sign(v. front());

	if (cur_sign == 0) return false;
	for (p = v. begin() + 1; p != v. end(); ++ p) {
		cur_sign =  -cur_sign;
		if (sign (*p) != cur_sign) return false;
	}

	return true;
}
		
template <class Poly, class Realring>
size_t alternations(const Poly& p, const Realring& R){
	size_t alts = 0;
	int curSign, prevSign = R.sign(p.back());
	for (typename Poly::const_reverse_iterator ptr = p.rbegin() + 1; 
	      ptr != p.rend(); ++ptr){
		  curSign = R.sign(*ptr);
		  if (curSign != 0 && curSign != prevSign) {
		    ++alts;
			prevSign = curSign;
		  }
	}
	return alts;
}

template <class Poly, class Realring>
Signature& signature(Signature& s, const Poly& p, Realring R){
	s.pos = alternations(p, R);

	typename Poly::const_iterator pptr;
	for(pptr = p.begin(); ! R.isZero(*pptr); ++pptr);
	s.zero = pptr - p.begin();

	s.neg = p.size() - 1 - s.pos - s.zero;
	return s;
}

template <class Vector>
bool isAllPositive (const Vector& v) {

	if (v. size() < 1) return true;

	typename Vector::const_iterator p;
	for (p = v. begin(); p != v. end(); ++ p) {
		if (sign (*p) <= 0) return false;
	}

	return true;
}

template <class Matrix, class MinPoly>
bool check_minpoly(const Matrix& A, const MinPoly& m, int n_try) {

	linbox_check (A. rowdim() == A. coldim());
	typedef typename Matrix::Field Field;
	typedef typename Field::Element Element;
	int n = A. rowdim();
	std::vector<Element> x(n), y(n), tmp(n);
	Field F(A. field());
	VectorDomain<Field> VCD(F);
	Element e;

	typename MinPoly::const_reverse_iterator m_p;
	typename std::vector<Element>::iterator x_p, y_p;

	for (int i = 0; i < n_try; ++ i) {

		for(x_p = x. begin(); x_p != x. end(); ++ x_p) 
			F. init (*x_p, rand());

		F. init (e, m. back());
		VCD. mul (y, x, e);

		for (m_p = m. rbegin() + 1; m_p != m. rend(); ++ m_p) {
			A. apply (tmp, y);
			F. init (e, *m_p);
			VCD. axpy (y, e, x, tmp);
		}

		if(!VCD.isZero(y)) return false;
	}

	return true;
}

template <class Ring>
int nonZeroEntries (const SparseMatrix<Ring>& M) {
	int E = 0;
	typename SparseMatrix<Ring>::ConstRawIterator raw_p;
	for (raw_p = M. rawBegin(); raw_p != M. rawEnd(); ++ raw_p)
		++ E;

	return E;
}

template <class Ring, class Blackbox>
int nonZeroEntries (const LieMatrix<Ring, Blackbox>& M) {
	
	typedef typename LieMatrix<Ring, Blackbox>::Ops Ops;
	typename Ops::const_iterator ops_p;

	int E = 0;
	for (ops_p = M. getOps(). begin(); ops_p != M. getOps(). end(); ++ ops_p)
		E += nonZeroEntries (**ops_p);
	
	return E;
}

}
#endif
