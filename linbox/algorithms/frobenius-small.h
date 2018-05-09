/* linbox/algorithms/coppersmith-invariant-factors.h
 * Copyright (C) 2018 Gavin Harrison
 *
 * Written by Gavin Harrison <gavin.har@gmail.com>
 *
 * ========LICENCE========
 * This file is part of the library LinBox.
 *
  * LinBox is free software: you can redistribute it and/or modify
 * it under the terms of the  GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 */

#ifndef __LINBOX_frobenius_small_H
#define __LINBOX_frobenius_small_H

#include <list>
#include <vector>
#include <math.h> 

#include "linbox/algorithms/block-coppersmith-domain.h"
#include "linbox/algorithms/blackbox-container.h"
#include "linbox/algorithms/blackbox-block-container-smmx2.h"
#include "linbox/matrix/random-matrix.h"

namespace LinBox
{

template<class _Field, class _PolynomialRing>
class FrobeniusSmall {
public:
	typedef _Field Field;
	typedef typename Field::Element Element;
	typedef typename Field::RandIter RandIter;
	typedef MatrixDomain<Field> MatrixDom;
	typedef typename MatrixDom::OwnMatrix Matrix;
	
	typedef SparseMatrix<Field, SparseMatrixFormat::CSR> Blackbox;
	typedef FFLAS::Sparse<Field, FFLAS::SparseMatrix_t::CSR> FSparseMat;
	
	typedef BlasVector<Field> Vector;
	
	typedef _PolynomialRing PolynomialRing;
	typedef typename PolynomialRing::Element Polynomial;
	typedef typename PolynomialRing::Coeff Coeff;
		
protected:
	Field _F;
	RandIter _RI;
	PolynomialRing _R;
	MatrixDom _MD;
	VectorDomain<Field> _VD;
	
public:
	FrobeniusSmall(const Field &F, const PolynomialRing &R) : _F(F), _RI(F), _R(R), _MD(F), _VD(F) {}

//protected:

	void convert(FSparseMat &A, const Blackbox &BB, bool transpose = false) const {
		std::vector<index_t> st1 = BB.getStart();
		std::vector<index_t> col1 = BB.getColid();
		std::vector<Element> data1 = BB.getData();
		
		uint64_t nnz = data1.size();
		
		index_t *row = FFLAS::fflas_new<index_t>(nnz);
		for (size_t j = 0; j < BB.rowdim(); ++j) {
			for (index_t k = st1[j] ; k < st1[j+1]; ++k) {
				row[k] = j;
			}
		}
		
		index_t *col = FFLAS::fflas_new<index_t>(nnz);
		for (size_t i = 0; i < col1.size(); i++) {
			col[i] = col1[i];
		}
		
		typename Field::Element_ptr data = FFLAS::fflas_new<Element>(nnz);
		for (size_t i = 0; i < data1.size(); i++) {
			data[i] = data1[i];
		}
		
		if (transpose) {
			FFLAS::sparse_init(_F, A, col, row, data, BB.rowdim(), BB.coldim(), nnz);
		} else {
			FFLAS::sparse_init(_F, A, row, col, data, BB.rowdim(), BB.coldim(), nnz);
		}
	}
	
	void convert(Vector &p, const Polynomial &f) {
		p.resize(_R.deg(f) + 1);
		for (size_t i = 0; i <= _R.deg(f); i++) {
			integer ei;
			Element e;
			Coeff c;
			
			_R.getCoeffField().convert(ei, _R.getCoeff(c, f, i));
			_F.init(e, i);
			
			p.setEntry(i, e);
		}
	}

	size_t trialbound() const {
		size_t p = _F.cardinality();
		
		if (p == 2) {
			return 6;
		} else if (p == 3) {
			return 5;
		} else if (p == 4) {
			return 4;
		} else if (5 <= p && p <= 9) {
			return 3;
		}
		
		return 2;
	}

	void randomVector(Vector &v) const {
		for (size_t i = 0; i < v.size(); i++) {
			Element e;
			_RI.random(e);
			
			v.setEntry(i, e);
		}
	}

	template<class Blackbox1>
	void minpolyseq(Polynomial &f, const Vector &u, const Blackbox1 &A, const Vector &v) const {
		typedef BlackboxContainer<Field, Blackbox1> Sequence;
		
		Sequence seq(&A, _F, u, v);
		MasseyDomain<Field, Sequence> WD(&seq, 20);
		
		BlasVector<Field> phi(_F);
		unsigned long deg;
		WD.minpoly(phi, deg);
		
		_R.init(f, phi);
	}
	
	template<class Blackbox1>
	void minpolyvec(Polynomial &f, const Blackbox1 &A, const Vector &v, size_t r) {
		_R.assign(f, _R.one);
		for (size_t i = 0; i < r; i++) {
			Vector u(_F, A.rowdim());
			randomVector(u);
			
			Polynomial tmp;
			minpolyseq(tmp, u, A, v);
						
			_R.lcmin(f, tmp);
		}
	}
	
	// compute the minimal polynomial of uAV with random V in F^(n-by-b)
	void minpolyvec(Polynomial &f, const Vector &u, const FSparseMat &A, size_t b) {
		typedef BlackboxBlockContainerSmmx2<Field> Sequence;
		
		RandomDenseMatrix<RandIter, Field> RDM(_F, _RI);
				
		BlasMatrix<Field> V(_F, u.size(), b);
		RDM.random(V);
		
		Sequence seq(A, _F, u, V);
		
		BlockCoppersmithDomain<MatrixDom, Sequence> coppersmith(_MD, &seq, 20);
		
		std::vector<BlasMatrix<Field>> gen;
		coppersmith.right_minpoly(gen);
		
		std::vector<Element> coeffs;
		for (size_t i = 0; i < gen.size(); i++) {
			coeffs.push_back(gen[i].getEntry(0, 0));
		}
		
		_R.init(f, coeffs);
	}
	
	void filterp(Polynomial &h, const Polynomial &f, const Polynomial &g) {
		Polynomial d;
		
		_R.assign(d, g);
		_R.assign(h, f);
		
		while (_R.deg(d) > 0) {
			_R.divin(h, d);
			_R.gcdin(d, h);
		}
	}
	
	template<class Blackbox1>
	void apply(Vector &v, const Polynomial &f, const Blackbox1 &A, const Vector &vin) {
		Vector p(_F);
		convert(p, f);
		
		PolynomialBB<Blackbox1, Vector> BB(A, p);
		v.resize(vin.size());
		BB.apply(v, vin);
	}
	
	void filterv(Vector &u, Vector &v, Polynomial &f, const Blackbox &A, const Vector &uin, const Vector &vin) {
		typedef Transpose<Blackbox> Transpose;
		Transpose T(A);
		
		Polynomial fm;
		minpolyseq(fm, uin, A, vin);
		
		Polynomial fl;
		minpolyvec(fl, T, uin, 7);
		
		Polynomial fr;
		minpolyvec(fr, A, vin, 7);
		
		//_R.write(std::cout << "fm: ", fm) << std::endl;
		//_R.write(std::cout << "fl: ", fl) << std::endl;
		//_R.write(std::cout << "fr: ", fr) << std::endl;
		
		Polynomial g, fd;
		_R.quo(fd, fl, fm);
		_R.gcd(g, fm, fd);
		
		Polynomial tmp;
		filterp(tmp, fm, g);
		_R.assign(fm, tmp);
		
		_R.quo(fd, fr, fm);
		_R.gcd(g, fm, fd);
		filterp(f, fm, g);
		
		Polynomial gl, gr;
		_R.quo(gl, fl, f);
		_R.quo(gr, fr, f);
		
		apply(u, gl, T, uin);
		apply(v, gr, A, vin);
	}
	
	void filterv(Vector &u, Vector &v, Polynomial &f, 
		const Blackbox &A, 
		const Blackbox &T,
		const FSparseMat &FA,
		const FSparseMat &FT,
		const Vector &uin, const Vector &vin) {
		
		Polynomial fm;
		minpolyseq(fm, uin, A, vin);
		
		Polynomial fl;
		minpolyvec(fl, uin, FA, 8);
		
		Polynomial fr;
		minpolyvec(fr, vin, FT, 8);
		
		//_R.write(std::cout << "fm: ", fm) << std::endl;
		//_R.write(std::cout << "fl: ", fl) << std::endl;
		//_R.write(std::cout << "fr: ", fr) << std::endl;
		
		Polynomial g, fd;
		_R.quo(fd, fl, fm);
		_R.gcd(g, fm, fd);
		
		Polynomial tmp;
		filterp(tmp, fm, g);
		_R.assign(fm, tmp);
		
		_R.quo(fd, fr, fm);
		_R.gcd(g, fm, fd);
		filterp(f, fm, g);
		
		Polynomial gl, gr;
		_R.quo(gl, fl, f);
		_R.quo(gr, fr, f);
		
		apply(u, gl, T, uin);
		apply(v, gr, A, vin);
	}
	
	void mergev(Vector &u, Vector &v, Polynomial &f, const Blackbox &A, const Vector &u1, const Vector &v1, const Polynomial &f1, const Vector &u2, const Vector &v2, const Polynomial &f2) {
		Polynomial g1, tmp;
		_R.lcm(tmp, f1, f2);
		_R.quoin(tmp, f2);
		filterp(g1, f1, tmp);
		
		Polynomial h1;
		_R.quo(h1, f1, g1);
		
		Polynomial g, h2;
		_R.gcd(g, h1, f2);
		filterp(h2, f2, g);
		
		Polynomial g2;
		_R.quo(g2, f2, h2);
		
		typedef Transpose<Blackbox> Transpose;
		Transpose T(A);
		
		Vector tu1(_F);
		Vector tu2(_F);
		apply(tu1, g1, T, u1);
		apply(tu2, g2, T, u2);
		u.resize(tu1.size());
		_VD.add(u, tu1, tu2);
		
		Vector tv1(_F);
		Vector tv2(_F);
		apply(tv1, g1, A, v1);
		apply(tv2, g2, A, v2);
		v.resize(tv1.size());
		_VD.add(v, tv1, tv2);
		
		_R.mul(f, h1, h2);
	}
	
	void minpolspace(Vector &u, Vector &v, Polynomial &f, const Blackbox &A, const std::vector<Vector> &us, const std::vector<Vector> &vs) {
		filterv(u, v, f, A, us[0], vs[0]);
		for (size_t i = 1; i < us.size(); i++) {
			Vector tmpu(_F);
			Vector tmpv(_F);
			Polynomial tmpf;
			filterv(tmpu, tmpv, tmpf, A, us[i], vs[i]);
			mergev(u, v, f, A, u, v, f, tmpu, tmpv, tmpf);
		}
	}
	
	void minpolspace(Vector &u, Vector &v, Polynomial &f,
		const Blackbox &A, 
		const Blackbox &T,
		const FSparseMat &FA,
		const FSparseMat &FT,
		const std::vector<Vector> &us, const std::vector<Vector> &vs) {
	
		filterv(u, v, f, A, us[0], vs[0]);
		for (size_t i = 1; i < us.size(); i++) {
			Vector tmpu(_F);
			Vector tmpv(_F);
			Polynomial tmpf;
			filterv(tmpu, tmpv, tmpf, A, T, FA, FT, us[i], vs[i]);
			mergev(u, v, f, A, u, v, f, tmpu, tmpv, tmpf);
		}
	}
	
	void solve(std::vector<Polynomial> &fs, const Blackbox &A) {
		std::vector<Vector> us;
		std::vector<Vector> vs;
		
		for (size_t i = 0; i < trialbound(); i++) {
			Vector u(_F, A.rowdim());
			Vector v(_F, A.coldim());
			
			randomVector(u);
			randomVector(v);
			
			us.push_back(u);
			vs.push_back(v);
		}
		
		Vector u(_F);
		Vector v(_F);
		Polynomial f;
		minpolspace(u, v, f, A, us, vs);
		
		fs.push_back(f);
	}
	
	void solve(std::vector<Polynomial> &fs, 
		const Blackbox &A, 
		const Blackbox &T, 
		const FSparseMat &FA,
		const FSparseMat &FT) {
	
		std::vector<Vector> us;
		std::vector<Vector> vs;
		
		for (size_t i = 0; i < trialbound(); i++) {
			Vector u(_F, A.rowdim());
			Vector v(_F, A.coldim());
			
			randomVector(u);
			randomVector(v);
			
			us.push_back(u);
			vs.push_back(v);
		}
		
		Vector u(_F);
		Vector v(_F);
		Polynomial f;
		minpolspace(u, v, f, A, T, FA, FT, us, vs);
		
		fs.push_back(f);
	}
};

}

#endif //__LINBOX_frobenius_small_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
