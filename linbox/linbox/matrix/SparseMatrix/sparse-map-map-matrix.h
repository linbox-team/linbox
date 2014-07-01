/* Copyright (C) 2013 LinBox
 * Written by AJS <stachnik@udel.edu>
 *
 *
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

/*! @file   matrix/SparseMatrix/sparse-map-map-matrix.h
 * @ingroup sparsematrix
 * @brief Supports fast elementary row AND column operations simultaneously
 */

#ifndef __LINBOX_SPARSE_MAP_MAP_MATRIX_H
#define __LINBOX_SPARSE_MAP_MAP_MATRIX_H

#include <stdlib.h>
#include <fstream>
#include <map>
#include <iostream>

#include <linbox/linbox-config.h>
#include <linbox/matrix/sparse-formats.h>
#include <linbox/util/debug.h>
#include <linbox/util/field-axpy.h>
#include <linbox/blackbox/blackbox-interface.h>
#include <linbox/field/hom.h>


namespace LinBox
{

template<class Field_>
class SparseMatrix<Field_,SparseMatrixFormat::SMM> : public BlackboxInterface {
public:
	typedef Field_ Field;
	typedef typename Field::Element Element;
	typedef size_t Index;
	typedef SparseMatrix<Field_,SparseMatrixFormat::SMM> Self_t;
	typedef std::map<Index,Element> VectorType;
	typedef typename VectorType::iterator VectorIt;
	typedef typename VectorType::const_iterator VectorConstIt;
	typedef std::map<Index,VectorType> MapType;
	typedef typename MapType::iterator MapIt;
	typedef typename MapType::const_iterator MapConstIt;

	SparseMatrix();

	SparseMatrix(const Field& F);

	SparseMatrix(const Field& F, Index r, Index c);

	SparseMatrix(const SparseMatrix& M);

	void init(const Field& F, Index r, Index c);

	void shape(Index r, Index c);

	SparseMatrix& operator=(const SparseMatrix& M);

	~SparseMatrix();

	void finalize() {}

	Index rowdim() const;

	Index coldim() const;

	const Field& field() const;

	void setEntry(Index i, Index j, const Element& e);

	const Element& getEntry(Index i, Index j) const;

	const Element& getEntry(Element& d, Index i, Index j) const;

	bool verify();

	template<class OutVector, class InVector>
	OutVector& apply(OutVector& y, const InVector& x) const;

	template<class OutVector, class InVector>
	OutVector& applyTranspose(OutVector& y, const InVector& x) const;

	//forall r: A_{i,r}<-A_{i,r}+k*A_{j,r}
	void addRow(const Element& k, Index i, Index j);

	//forall r: A_{r,i}<-A_{r,i}+k*A_{r,j}
	void addCol(const Element& k, Index i, Index j);

	//forall r: A_{i,r}<-k*A_{i,r}
	void timesRow(const Element& k, Index i);

	//forall r: A_{r,j}<-k*A_{r,j}
	void timesCol(const Element& k, Index j);

	void swapRows(Index i, Index j);

	void swapCols(Index i, Index j);

	void scaleMat(const Element& k);

	void transpose();

	Index nnz() const;

	// A -> A' = SAS^{-1}, and A' has about nnz nonzero entries.
	void randomSim(Index nnz, int seed = 0);

	// A -> A' = UAV, with U and V nonsingular, and A' has about nnz nonzero entries.
	void randomEquiv(Index nnz, int seed = 0);

	bool areEqual(Self_t& rhs) const;

	std::istream& read(std::istream& in);

	std::ostream& print(std::ostream& out) const;

	std::ostream& write(std::ostream& out) const;

	template<class Matrix>
	void copy(Matrix& mat) const;

	template<class Matrix>
	void copyFrom(Matrix& mat);

	template<class Vector>
	void toVector(Vector& vec) const;

	template<class Vector>
	void fromVector(const Vector& vec, Index r, Index c);

	static void generateCompanion(Self_t& mat,std::vector<Element>& coeffs);

	static void generateDenseRandMat(Self_t& mat, int q);

	static void generateRandMat(Self_t& mat, int nnz, int q);

	static void generateScaledIdent(Self_t& mat, int alpha);

	static void generateSparseNonSingular(Self_t& mat, int approxNNZ, int seed=0);

protected:

	static int randRange(int start, int end);

	MatrixDomain<Field> MD_;

	MapType rowMap_;

	MapType colMap_;

	Index numCols_;

	Index numRows_;

	Index nnz_;

	Element zero_;
};

}

#include <linbox/matrix/SparseMatrix/sparse-map-map-matrix.inl>

#endif // __SPARSE_MAP_MAP_MATRIX_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
