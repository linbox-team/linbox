/* Copyright (C) 2015 LinBox
 *
 *  Author: Gavin Harrison
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 */

#ifndef __LINBOX_smith_form_iliopoulos2_H
#define __LINBOX_smith_form_iliopoulos2_H

namespace LinBox
{

template<class _Field>
class IliopoulosDomain {
public:
	typedef _Field Field;
	typedef MatrixDomain<Field> Domain;
	typedef typename Domain::OwnMatrix Matrix;
	typedef typename Field::Element Element;

protected:
	Domain _MD;
	Field _F;

public:
	void xgcd(
		Element &g,
		std::vector<Element> &ts,
		const std::vector<Element> &as)
	{
		ts.resize(as.size());
		_F.assign(ts[0], _F.one);
		_F.assign(g, as[0]);
		
		for (size_t i = 1; i < as.size(); i++) {
			Element r, t1, t2;
			_F.gcd(r, t1, t2, g, as[i]);
			
			for (size_t j = 0; j < i; j++) {
				_F.mulin(ts[j], t1);
			}
			
			_F.assign(ts[i], t2);
			_F.assign(g, r);
		}
	}
	
	void swapRows(Matrix &A, size_t a, size_t b)
	{
		size_t n = A.coldim();
		
		Element tmp1, tmp2;
		for (size_t i = 0; i < n; i++) {
			A.getEntry(tmp1, a, i);
			A.getEntry(tmp2, b, i);
			
			A.setEntry(a, i, tmp2);
			A.setEntry(b, i, tmp1);
		}
	}
	
	void swapCols(Matrix &A, size_t a, size_t b)
	{
		size_t n = A.coldim();
		
		Element tmp1, tmp2;
		for (size_t i = 0; i < n; i++) {
			A.getEntry(tmp1, i, a);
			A.getEntry(tmp2, i, b);
			
			A.setEntry(i, a, tmp2);
			A.setEntry(i, b, tmp1);
		}
	}
	
	bool isRowZero(const Matrix &A, size_t p)
	{
		size_t n = A.coldim();
		
		Element tmp;
		for (size_t i = p; i < n; i++) {
			if (!_F.isZero(_MD.getEntry(tmp, p, i))) {
				return false;
			}
		}
		
		return true;
	}
	
	bool moveNonzeroRowToPivot(Matrix &A, size_t p)
	{
		size_t n = A.coldim();
		
		if (!isRowZero(A, p)) {
			return true;
		}
		
		for (size_t i = p; i < n; i++) {
			if (!isRowZero(A, i)) {
				swapRows(A, p, i);
				return true;
			}
		}
		
		return false;
	}
	
	bool moveZeroColToPivot(Matrix &A, size_t p)
	{
		size_t n = A.coldim();
		
		Element tmp;
		if (_F.isZero(A.getEntry(tmp, p, p))) {
			return true;
		}
		
		for (size_t i = p+1; i < n; i++) {
			if (_F.isZero(A.getEntry(tmp, p, i))) {
				swapCols(A, p, i);
				return true;
			}
		}
		
		return false;
	}
	
	void makePivotColZero(Matrix &A, size_t p)
	{
		Element r, x1, x2, y1, y2, a, b;
		
		A.getEntry(a, p, p);
		A.getEntry(b, p, p+1);
		
		_F.gcd(r, x1, x2, a, b);
		
		_F.div(y1, b, r);
		_F.div(y2, a, r);
		_F.negin(y2);
		
		size_t n = A.coldim();
		
		for (size_t i = 0; i < n; i++) {
			A.getEntry(a, i, p);
			A.getEntry(b, i, p+1);
			
			Element new_a, new_b, tmp;
			_F.mul(tmp, y1, a);
			_F.mul(new_a, y2, b);
			_F.addin(new_a, tmp);
			
			_F.mul(tmp, x1, a);
			_F.mul(new_b, x2, b);
			_F.addin(new_b, tmp);
			
			A.setEntry(i, p, new_a);
			A.setEntry(i, p+1, new_b);
		}
	}
	
	void makePivotRowZero(Matrix &A, size_t p)
	{
		Element r, x1, x2, y1, y2, a, b;
		
		A.getEntry(a, p, p);
		A.getEntry(b, p+1, p);
		
		_F.gcd(r, x1, x2, a, b);
		_F.div(y1, b, r);
		_F.div(y2, a, r);
		_F.negin(y2);
		
		size_t n = A.rowdim();
		
		for (size_t i = 0; i < n; i++) {
			A.getEntry(a, p, i);
			A.getEntry(b, p+1, i);
			
			Element new_a, new_b, tmp;
			_F.mul(tmp, y1, a);
			_F.mul(new_a, y2, b);
			_F.addin(new_a, tmp);
			
			_F.mul(tmp, x1, a);
			_F.mul(new_b, x2, b);
			_F.addin(new_b, tmp);
			
			A.setEntry(p, i, new_a);
			A.setEntry(p+1, i, new_b);
		}
	}
	
	void scaleAndAddColsToPivotCol(
		Matrix &A,
		const std::vector<Element> &ts,
		size_t p)
	{	
		size_t n = A.coldim();
		
		// Add cols p+1 to n into p, scaled by ts[i]
		for (size_t j = p; j < n; j++) {
			Element a;
			A.getEntry(a, j, p);
			
			for (size_t i = p+1; i < n; i++) {
				Element b, tmp;
				A.getEntry(b, j, i);
				_F.mulin(b, ts[i - p - 1]);
				_F.addin(a, b);
			}
			
			A.setEntry(j, p, a);
		}
	}
	
	void eliminateOffPivotCols(Matrix &A, Element &s, size_t p)
	{
		size_t n = A.coldim();
		
		// zero out other columns entry in pivot row
		for (size_t i = p+1; i < n; i++) {
			Element q;
			A.getEntry(q, p, i);
			_F.divin(q, s);
			
			for (size_t j = p; j < n; j++) {
				Element a, b, tmp;
				A.getEntry(a, j, p);
				A.getEntry(b, j, i);
				_F.mul(tmp, a, q);
				_F.subin(b, tmp);
				A.setEntry(j, i, b);
			}
		}
	}
	
	void reduceEntries(Matrix &A, const Element &d, size_t p)
	{
		size_t n = A.coldim();
		
		// reduce entries mod d
		for (size_t i = p; i < n; i++) {
			for (size_t j = p; j < n; j++) {
				Element tmp;
				A.getEntry(tmp, i, j);
				_F.modin(tmp, d);
				A.setEntry(i, j, tmp);
			}
		}
	}
	
	void eliminateRow(Matrix &A, const Element &d, size_t p)
	{
		size_t n = A.coldim();
		
		if (!moveZeroColToPivot(A, p)) {
			makePivotColZero(A, p);
		}
		
		std::vector<Element> ks;
		for (size_t i = p+1; i < n; i++) {
			Element tmp;
			A.getEntry(tmp, p, i);
			ks.push_back(tmp);
		}
		
		Element s;
		std::vector<Element> ts;
		xgcd(s, ts, ks);
		
		scaleAndAddColsToPivotCol(A, ts, p);
		
		eliminateOffPivotCols(A, s, p);
		
		reduceEntries(A, d, p);
	}
	
	bool moveZeroRowToPivot(Matrix &A, size_t p)
	{
		size_t n = A.rowdim();
		
		for (size_t i = p+1; i < n; i++) {
			Element tmp;
			A.getEntry(tmp, i, p);
			
			if (_F.isZero(tmp)) {
				swapRows(A, p, i);
				return true;
			}
		}
		
		return false;
	}
	
	void scaleAndAddRowsToPivotRow(
		Matrix &A,
		const std::vector<Element> &ts,
		size_t p)
	{	
		size_t n = A.coldim();
		
		// Add rows p+1 to n into p, scaled by ts[i]
		for (size_t j = p; j < n; j++) {
			Element a;
			A.getEntry(a, p, j);
			
			for (size_t i = p+1; i < n; i++) {
				Element b, tmp;
				A.getEntry(b, i, j);
				_F.mulin(b, ts[i - p - 1]);
				_F.addin(a, b);
			}
			
			A.setEntry(p, j, a);
		}
	}
	
	void eliminateOffPivotRows(Matrix &A, Element &s, size_t p)
	{
		size_t n = A.coldim();
		
		// zero out other rows entry in pivot col
		for (size_t i = p+1; i < n; i++) {
			Element q;
			A.getEntry(q, i, p);
			_F.divin(q, s);
			
			for (size_t j = p; j < n; j++) {
				Element a, b, tmp;
				A.getEntry(a, p, j);
				A.getEntry(b, i, j);
				_F.mul(tmp, a, q);
				_F.subin(b, tmp);
				A.setEntry(i, j, b);
			}
		}
	}
	
	void eliminateCol(Matrix &A, const Element &d, size_t p)
	{
		size_t n = A.coldim();
		
		if (!moveZeroRowToPivot(A, p)) {
			makePivotRowZero(A, p);
		}
		
		std::vector<Element> ks;
		for (size_t i = p+1; i < n; i++) {
			Element tmp;
			A.getEntry(tmp, i, p);
			ks.push_back(tmp);
		}
		ks.push_back(d);
		
		Element s;
		std::vector<Element> ts;
		xgcd(s, ts, ks);
		
		scaleAndAddRowsToPivotRow(A, ts, p);
		A.setEntry(p, p, s);
		
		eliminateOffPivotRows(A, s, p);
		
		reduceEntries(A, d, p);
	}
	
	bool pivotDividesRow(Matrix &A, size_t p)
	{
		size_t n = A.coldim();
		
		Element a;
		A.getEntry(a, p, p);
		
		for (size_t i = p+1; i < n; i++) {
			Element b;
			A.getEntry(b, p, i);
			
			if (!_F.isDivisor(b, a)) {
				return false;
			}
		}
		
		return true;
	}
	
	void zeroOutRow(Matrix &A, size_t p)
	{
		size_t n = A.coldim();
		
		for (size_t i = p+1; i < n; i++) {
			A.setEntry(p, i, _F.zero);
		}
	}

public:
	IliopoulosDomain(Field &F) :
		_MD(F),
		_F(F)
	{
	}
	
	void smithFormIn(Matrix &A, const Element &d)
	{
		size_t n = A.coldim();
		
		for (size_t p = 0; p < n - 1; p++) {
			do {
				eliminateRow(A, d, p);
				eliminateCol(A, d, p);
			} while(!pivotDividesRow(A, p));
			zeroOutRow(A, p);
		}
		
		for (size_t p = 0; p < n; p++) {
			for (size_t q = p+1; q < n; q++) {
				Element h;
				A.getEntry(h, p, p);
				
				Element Aqq;
				A.getEntry(Aqq, q, q);
				
				Element App;
				_F.gcd(App, Aqq, h);
				_F.divin(h, App);
				_F.mulin(Aqq, h);
				
				A.setEntry(p, p, App);
				A.setEntry(q, q, Aqq);
			}
		}
	}
	
	template <class PolyRingVector>
	void smithForm(PolyRingVector &diag, const Matrix &A, const Element &d)
	{
		size_t n = A.coldim();
		
		Matrix B(_F, n, n);
		_MD.copy(B, A);
		
		smithFormIn(B, d);
		
		for (size_t i = 0; i < n; i++) {
			B.getEntry(diag[i], i, i);
			_F.gcdin(diag[i], d);
			_F.modin(diag[i], d);
			_F.normalizeIn(diag[i]);
		}
	}
};

}

#endif //__LINBOX_smith_form_iliopoulos2_H