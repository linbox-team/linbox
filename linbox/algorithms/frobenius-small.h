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

#include "givaro/givtimer.h"

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
	
	Givaro::Timer TW;
public:
	FrobeniusSmall(const Field &F, const PolynomialRing &R) : _F(F), _RI(F), _R(R), _MD(F), _VD(F) {}

//protected:

	void print(const Vector &v) {
		for (size_t i = 0; i < v.size(); i++) {
			_F.write(std::cout, v.getEntry(i)) << " ";
		}
		std::cout << std::endl;
	}
	
	void convert(Vector &p, const Polynomial &f) {
		p.resize(_R.deg(f) + 1);
		for (size_t i = 0; i <= _R.deg(f); i++) {
			integer ei;
			Element e;
			Coeff c;
			
			_R.getCoeffField().convert(ei, _R.getCoeff(c, f, i));
			_F.init(e, ei);
			
			p.setEntry(i, e);
		}
	}
	
	void copy(Vector &a, const Vector &b) {
		for (size_t i = 0; i < a.size(); i++) {
			a.setEntry(i, b.getEntry(i));
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
	
	// return true if fk(A)v <> 0
	template<class Blackbox1>
	bool test(const Polynomial &fk, const Blackbox1 &A, const Vector &v) {
		Vector r(_F, v.size());
		
		apply(r, fk, A, v);
		
		for (size_t i = 0; i < r.size(); i++) {
			if (!_F.isZero(r.getEntry(i))) {
				return true;
			}
		}
		
		return false;
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
	
	void randomU(
		Vector &u, 
		const std::vector<std::vector<Vector>> &basis_us, 
		const std::vector<std::vector<Vector>> &basis_vs, 
		const std::vector<std::vector<Element>> &ss) {
		
		randomVector(u);
		Vector const_u(_F, u.size());
		copy(const_u, u);
		
		for (size_t i = 0; i < basis_us.size(); i++) {
			for (size_t j = 0; j < basis_us[i].size(); j++) {
				Vector tmp(_F, u.size());
				Element scale;
				
				_VD.dot(scale, const_u, basis_vs[i][j]);
				_F.divin(scale, ss[i][j]);
				_VD.mul(tmp, basis_us[i][j], scale);
				
				_VD.subin(u, tmp);
			}
		}
	}
	
	void randomV(
		Vector &v,
		const std::vector<std::vector<Vector>> &basis_us, 
		const std::vector<std::vector<Vector>> &basis_vs, 
		const std::vector<std::vector<Element>> &ss) {
	
		randomVector(v);
		Vector const_v(_F, v.size());
		copy(const_v, v);
		
		for (size_t i = 0; i < basis_us.size(); i++) {
			for (size_t j = 0; j < basis_us[i].size(); j++) {
				Vector tmp(_F, v.size());
				Element scale;
				
				_VD.dot(scale, const_v, basis_us[i][j]);
				_F.divin(scale, ss[i][j]);
				_VD.mul(tmp, basis_vs[i][j], scale);
				
				_VD.subin(v, tmp);
			}
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
	void minpolyvec(Polynomial &f, const Blackbox1 &A, const Vector &v) {		
		Vector u(_F, A.rowdim());
		randomVector(u);
		minpolyseq(f, u, A, v);
		
		while (test(f, A, v)) {
			Polynomial tmp;
			randomVector(u);
			
			minpolyseq(tmp, u, A, v);
			
			_R.lcmin(f, tmp);
		}
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
	
	void filterv(Vector &u, Vector &v, Polynomial &f, const Blackbox &A, const Vector &uin, const Vector &vin) {
		typedef Transpose<Blackbox> Transpose;
		Transpose T(A);
		
		Polynomial fm;
		minpolyseq(fm, uin, A, vin);
		
		Polynomial fl;
		minpolyvec(fl, T, uin);
		
		Polynomial fr;
		minpolyvec(fr, A, vin);
		
		//std::cout << "filterv" << std::endl;
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
		
		//_R.write(std::cout << "gl: ", gl) << std::endl;
		//_R.write(std::cout << "gr: ", gr) << std::endl;
		
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
		size_t n = A.rowdim();
		
		filterv(u, v, f, A, us[0], vs[0]);	
		for (size_t i = 1; i < us.size(); i++) {			
			Vector tmpu(_F, n);
			Vector tmpv(_F, n);
			Polynomial tmpf;
			filterv(tmpu, tmpv, tmpf, A, us[i], vs[i]);
			
			Vector oldu(_F, n);
			Vector oldv(_F, n);
			Polynomial oldf;
			
			copy(oldu, u);
			copy(oldv, v);
			_R.assign(oldf, f);
			
			mergev(u, v, f, A, oldu, oldv, oldf, tmpu, tmpv, tmpf);
		}
	}
	
	void dualbasis(
		std::vector<Vector> &basis_u, 
		std::vector<Vector> &basis_v,
		std::vector<Element> &basis_s,
		const Vector &u, 
		const Blackbox &A, 
		const Vector &v, 
		size_t dk) {
	
		size_t n = A.rowdim();
	
		Vector tmpu1(u);
		Vector tmpv1(v);
		Element tmps1;
		_VD.dot(tmps1, tmpu1, tmpv1);
		
		basis_u.push_back(tmpu1);
		basis_v.push_back(tmpv1);
		basis_s.push_back(tmps1);
		
		for (size_t i = 1; i < dk; i++) {
			Vector tmpu(_F, n);
			Vector tmpv(_F, n);
			Element tmps;
			
			A.applyTranspose(tmpu, basis_u[basis_u.size() - 1]);
			basis_u.push_back(tmpu);
			
			A.apply(tmpv, basis_v[basis_v.size() - 1]);
			basis_v.push_back(tmpv);
			
			_VD.dot(tmps, tmpu, tmpv);
			basis_s.push_back(tmps);
		}
	}
	
	void solve(std::vector<Polynomial> &fs, const Blackbox &A) {
		TW.start();
		
		size_t n = A.rowdim();
		size_t k = 0;
		size_t d = 0;
		
		std::vector<Vector> us;
		std::vector<Vector> vs;
		
		std::vector<std::vector<Vector>> basis_us;
		std::vector<std::vector<Vector>> basis_vs;
		std::vector<std::vector<Element>> ss;
		
		Transpose<Blackbox> T(A);
		
		while (d < n) {
			std::vector<Vector> trial_us;
			std::vector<Vector> trial_vs;
			
			Vector old_u(_F, n);
			Vector old_v(_F, n);
			Polynomial old_f;
			
			bool dropped = false;
			for (size_t i = 0; i < trialbound(); i++) {
				Vector u(_F, A.rowdim());
				Vector v(_F, A.coldim());
				
				randomU(u, basis_us, basis_vs, ss);
				randomV(v, basis_vs, basis_vs, ss);
				
				trial_us.push_back(u);
				trial_vs.push_back(v);
				
				while (k > 0 && (test(fs[fs.size() - 1], T, u) || test(fs[fs.size() - 1], A, v))) {
					dropped = true;
					
					_R.assign(old_f, fs[fs.size() - 1]);
					copy(old_u, us[us.size() - 1]);
					copy(old_v, vs[vs.size() - 1]);
					
					fs.pop_back();
					us.pop_back();
					vs.pop_back();
					
					basis_us.pop_back();
					basis_vs.pop_back();
					ss.pop_back();
					
					d -= _R.deg(old_f);
					k--;
				}
			}
			
			Vector u(_F, n);
			Vector v(_F, n);
			Polynomial f;
			if (dropped) {
				Vector tmpu(_F);
				Vector tmpv(_F);
				Polynomial tmpf;
				
				minpolspace(tmpu, tmpv, tmpf, A, trial_us, trial_vs);
				mergev(u, v, f, A, old_u, old_v, old_f, tmpu, tmpv, tmpf);
			} else {
				minpolspace(u, v, f, A, trial_us, trial_vs);
			}
			
			fs.push_back(f);
			us.push_back(u);
			vs.push_back(v);
			
			{
				Polynomial fu;
				minpolyvec(fu, T, u);
				
				Polynomial fv;
				minpolyvec(fv, A, v);
				
				//_R.write(std::cout << "fm: ", f) << std::endl;
				//_R.write(std::cout << "fu: ", fu) << std::endl;
				//_R.write(std::cout << "fv: ", fv) << std::endl;
			}
			
			k++;
			d += _R.deg(f);
			
			/*
			std::vector<Vector> basis_u;
			std::vector<Vector> basis_v;
			std::vector<Element> basis_s;
			dualbasis(basis_u, basis_v, basis_s, u, A, v, _R.deg(f));
			basis_us.push_back(basis_u);
			basis_vs.push_back(basis_v);
			ss.push_back(basis_s);
			
			bool ifBasis = true;
			for (size_t i = 0; i < basis_s.size(); i++) {
				_F.write(std::cout << "s[" << i << "]: ", basis_s[i]) << std::endl;
				ifBasis = ifBasis && !_F.isZero(basis_s[i]);
			}
			std::cout << std::endl;
			*/
			
			if (k == 1) {
				break;
			}
		}
		
		TW.stop();
		std::cout << "frobenius-small: " << TW.usertime() << std::endl;
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
