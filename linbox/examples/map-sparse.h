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

/*! @file   examples/map-sparse.h
 * @ingroup examples
 * @brief
 */

#ifndef __LINBOX_MAP_SPARSE_H
#define __LINBOX_MAP_SPARSE_H

#include <stdlib.h>
#include <fstream>
#include <map>
#include <iostream>

#include <linbox/linbox-config.h>
#include <linbox/util/debug.h>
#include <linbox/util/field-axpy.h>
#include <linbox/blackbox/blackbox-interface.h>
#include <linbox/field/hom.h>


namespace LinBox
{

template<class Field_>
class MapSparseIF {
public:
	typedef Field_ Field;
	typedef typename Field::Element Element;
	typedef size_t Index;

        virtual ~MapSparseIF() {};

	virtual void setEntry(Index i, Index j, const Element& e) =0;

        virtual const Element& getEntry(Index i, Index j) const =0;

protected:
        typedef std::map<Index,Element> VectorType;
	typedef typename VectorType::iterator VectorIt;
	typedef typename VectorType::const_iterator VectorConstIt;
	typedef std::map<Index,VectorType> MapType;
	typedef typename MapType::iterator MapIt;
	typedef typename MapType::const_iterator MapConstIt;


};

template<class Field_>
class MapSparse : public MapSparseIF<Field_> {
public:
	typedef Field_ Field;
	typedef typename Field::Element Element;
	typedef size_t Index;
	typedef std::map<Index,Element> VectorType;
	typedef typename VectorType::iterator VectorIt;
	typedef typename VectorType::const_iterator VectorConstIt;
	typedef std::map<Index,VectorType> MapType;
	typedef typename MapType::iterator MapIt;
	typedef typename MapType::const_iterator MapConstIt;

	MapSparse();

	MapSparse(const Field& F, Index r, Index c);

	MapSparse(const MapSparse& M);

        void init(const Field& F, Index r, Index c);

        void shape(Index r, Index c);

	MapSparse& operator=(const MapSparse& M);

	~MapSparse();

	Index rowdim() const;

	Index coldim() const;

	const Field& field() const;

	void setEntry(Index i, Index j, const Element& e);

	const Element& getEntry(Index i, Index j) const;

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

        bool areEqual(MapSparse<Field_> rhs) const;

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

        static void generateDenseRandMat(MapSparse<Field_>& mat, int q);

        static void generateRandMat(MapSparse<Field>& mat, int nnz, int q);

        static void generateScaledIdent(MapSparse<Field>& mat, int alpha);

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

#include <examples/map-sparse.inl>

#endif // __LINBOX_MAP_SPARSE_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
