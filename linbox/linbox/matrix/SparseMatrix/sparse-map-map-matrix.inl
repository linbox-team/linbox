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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 */

/*! @file   matrix/SparseMatrix/sparse-map-map-matrix.inl
 * @ingroup sparsematrix
 * @brief Supports fast elementary row AND column operations simultaneously
 */

#ifndef __LINBOX_SPARSE_MAP_MAP_MATRIX_INL
#define __LINBOX_SPARSE_MAP_MAP_MATRIX_INL

#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <set>
#include <utility>
#include <linbox/util/matrix-stream.h>

namespace LinBox
{

template<class Field_>
SparseMatrix<Field_,SparseMatrixFormat::SMM>::SparseMatrix() : MD_(),numCols_(0), numRows_(0), nnz_(0) {}

template<class Field_>
SparseMatrix<Field_,SparseMatrixFormat::SMM>::SparseMatrix(const Field& F) :
	MD_(F),numCols_(0), numRows_(0), nnz_(0) {F.assign(zero_,F.zero);}

template<class Field_>
SparseMatrix<Field_,SparseMatrixFormat::SMM>::SparseMatrix(const Field& F, Index r, Index c):
	MD_(F), numCols_(c), numRows_(r), nnz_(0) {F.assign(zero_,F.zero);}

template<class Field_>
void SparseMatrix<Field_,SparseMatrixFormat::SMM>::init(const Field& F, Index r, Index c) {
        MD_=F;
        F.assign(zero_,F.zero);
        shape(r,c);
}

template<class Field_>
void SparseMatrix<Field_,SparseMatrixFormat::SMM>::shape(Index r, Index c) {
        rowMap_.clear();
        colMap_.clear();
        numRows_=r;
        numCols_=c;
        nnz_=0;
}

template<class Field_>
template<class Vector>
void SparseMatrix<Field_,SparseMatrixFormat::SMM>::fromVector(const Vector& vec, Index r, Index c) {
        shape(r,c);
        if (numCols_==1) {
                for (Index i=0;i<numRows_;++i) {
                        setEntry(i,0,vec[i]);
                }
        } else {
                for (Index j=0;j<numCols_;++j) {
                        setEntry(0,j,vec[j]);
                }
        }
}

template<class Field_>
SparseMatrix<Field_,SparseMatrixFormat::SMM>::SparseMatrix(const SparseMatrix<Field_,SparseMatrixFormat::SMM>& M):
        MD_(M.field()),
	rowMap_(M.rowMap_), colMap_(M.colMap_),
        numCols_(M.numCols_), numRows_(M.numRows_),
        nnz_(M.nnz_), zero_(M.zero_) {}

template<class Field_>
SparseMatrix<Field_,SparseMatrixFormat::SMM>& SparseMatrix<Field_,SparseMatrixFormat::SMM>::operator=(const SparseMatrix<Field_,SparseMatrixFormat::SMM>& rhs)
{
	if (rhs==this) return;
	MD_.init(rhs.MD_);
	numCols_=rhs.numCols_;
	numRows_=rhs.numRows_;
	rowMap_=rhs.rowMap_;
	colMap_=rhs.colMap_;
	nnz_=rhs.nnz_;
	zero_=rhs.zero_;

	return *this;
}

template<class Field_>
SparseMatrix<Field_,SparseMatrixFormat::SMM>::~SparseMatrix() {}

template<class Field_>
const Field_& SparseMatrix<Field_,SparseMatrixFormat::SMM>::field() const
{
	return MD_.field();
}

template<class Field_>
bool SparseMatrix<Field_,SparseMatrixFormat::SMM>::verify()
{
	for (int i=0;i<rowdim();++i) {
		for (int j=0;j<coldim();++j) {
			Element d=zero_;
			MapConstIt row=rowMap_.find(i);
			if (row != rowMap_.end()) {
				VectorConstIt entry=(row->second).find(j);
				if (entry != (row->second.end())) {
					d=entry->second;
				}
			}
			Element e=zero_;
			MapConstIt col=colMap_.find(j);
			if (col != colMap_.end()) {
				VectorConstIt entry=(col->second).find(i);
				if (entry != (col->second.end())) {
					e=entry->second;
				}
			}
			if (d!=e) {
				return false;
			}
		}
	}
        return true;
}

template<class Field_>
void SparseMatrix<Field_,SparseMatrixFormat::SMM>::setEntry(Index i, Index j, const Element& e)
{
	VectorIt it=rowMap_[i].find(j);
	if (it != rowMap_[i].end()) {
		--nnz_;
		rowMap_[i].erase(it);
		VectorIt colIt=colMap_[j].find(i);
		colMap_[j].erase(colIt);
	}

	if (!field().isZero(e)) {
		++nnz_;
		field().assign(rowMap_[i][j],e);
		field().assign(colMap_[j][i],e);
	}
}

template<class Field_> const typename SparseMatrix<Field_,SparseMatrixFormat::SMM>::Element&
SparseMatrix<Field_,SparseMatrixFormat::SMM>::getEntry(Index i, Index j) const
{
	MapConstIt row=rowMap_.find(i);
	if (row != rowMap_.end()) {
		VectorConstIt entry=(row->second).find(j);
		if (entry != (row->second.end())) {
			return (entry->second);
		}
	}
	return zero_;
}

template<class Field_> const typename SparseMatrix<Field_,SparseMatrixFormat::SMM>::Element&
SparseMatrix<Field_,SparseMatrixFormat::SMM>::getEntry(Element& d, Index i, Index j) const
{
	MapConstIt row=rowMap_.find(i);
	if (row != rowMap_.end()) {
		VectorConstIt entry=(row->second).find(j);
		if (entry != (row->second.end())) {
			d=entry->second;
			return (entry->second);
		}
	}
	d=zero_;
	return zero_;
}

template<class Field_>
template<class OutVector, class InVector>
OutVector& SparseMatrix<Field_,SparseMatrixFormat::SMM>::apply(OutVector& y, const InVector& x) const
{
	linbox_check(rowdim()==y.size());
	linbox_check(coldim()==x.size());
	for (int i=0;i<rowdim();++i) {
		Element d,e;
		field().init(d,0);
		MapConstIt rowI=rowMap_.find(i);
		if (rowI != rowMap_.end()) {
			for (VectorConstIt it=rowI->second.begin();
			     it!=rowI->second.end();++it) {
				field().mul(e,it->second,x[it->first]);
				field().addin(d,e);
			}
		}
		y[i]=d;
	}
	return y;
}

template<class Field_>
template<class OutVector, class InVector>
OutVector& SparseMatrix<Field_,SparseMatrixFormat::SMM>::applyTranspose(OutVector& y, const InVector& x) const
{
	linbox_check(coldim()==y.size());
	linbox_check(rowdim()==x.size());
	for (int i=0;i<coldim();++i) {
		Element d,e;
		field().init(d,0);
		MapConstIt colI=colMap_.find(i);
		if (colI != colMap_.end()) {
			for (VectorConstIt it=colI->second.begin();
			     it!=colI->second.end();++it) {
				field().mul(e,it->second,x[it->first]);
				field().addin(d,e);
			}
		}
		y[i]=d;
	}
	return y;
}

//forall r: A_{i,r}<-A_{i,r}+k*A_{j,r}
template<class Field_>
void SparseMatrix<Field_,SparseMatrixFormat::SMM>::addRow(const Element& k, Index i, Index j)
{
	MapIt rowJ=rowMap_.find(j);
	if (rowJ != rowMap_.end()) {
		for (VectorIt it=rowJ->second.begin();it!=rowJ->second.end();++it) {
			Index col=it->first;
			Element d=getEntry(i,col);
			field().axpyin(d,k,it->second);
			setEntry(i,col,d);
		}
	}
}

//forall r: A_{r,i}<-A_{r,i}+k*A_{r,j}
template<class Field_>
void SparseMatrix<Field_,SparseMatrixFormat::SMM>::addCol(const Element& k, Index i, Index j)
{
	MapIt colJ=colMap_.find(j);
	if (colJ != colMap_.end()) {
		for (VectorIt it=colJ->second.begin();it!=colJ->second.end();++it) {
			Index row=it->first;
			Element d = getEntry(row,i);
			field().axpyin(d,k,it->second);
			setEntry(row,i,d);
		}
	}
}

//forall r: A_{i,r}<-k*A_{i,r}
template<class Field_>
void SparseMatrix<Field_,SparseMatrixFormat::SMM>::timesRow(const Element& k, Index i)
{
	linbox_check(!(field().isZero(k)));

	MapIt row=rowMap_.find(i);
	if (row != rowMap_.end()) {
		for (VectorIt it=row->second.begin();it!=row->second.end();++it) {
			Index col=it->first;
			field().mulin(rowMap_[i][col],k);
			field().mulin(colMap_[col][i],k);
		}
	}
}

//forall r: A_{r,j}<-k*A_{r,j}
template<class Field_>
void SparseMatrix<Field_,SparseMatrixFormat::SMM>::timesCol(const Element& k, Index j)
{
	linbox_check(!(field().isZero(k)));

	MapIt col=colMap_.find(j);
	if (col != colMap_.end()) {
		for (VectorIt it=col->second.begin();it!=col->second.end();++it) {
			Index row=it->first;
			field().mulin(rowMap_[row][j],k);
			field().mulin(colMap_[j][row],k);
		}
	}
}

template<class Field_>
void SparseMatrix<Field_,SparseMatrixFormat::SMM>::swapRows(Index i, Index j)
{
	VectorType oldRowI=rowMap_[i];
	VectorType oldRowJ=rowMap_[j];

	rowMap_[i].clear();
	rowMap_[j].clear();

	for (VectorIt it=oldRowI.begin();it!=oldRowI.end();++it) {
		Index col=it->first;
		colMap_[col].erase(colMap_[col].find(i));
	}

	for (VectorIt it=oldRowJ.begin();it!=oldRowJ.end();++it) {
		Index col=it->first;
		colMap_[col].erase(colMap_[col].find(j));
	}

	for (VectorIt it=oldRowI.begin();it!=oldRowI.end();++it) {
		Index col=it->first;
		rowMap_[j][col]=oldRowI[col];
		colMap_[col][j]=oldRowI[col];
	}

	for (VectorIt it=oldRowJ.begin();it!=oldRowJ.end();++it) {
		Index col=it->first;
		rowMap_[i][col]=oldRowJ[col];
		colMap_[col][i]=oldRowJ[col];
	}
}

template<class Field_>
void SparseMatrix<Field_,SparseMatrixFormat::SMM>::swapCols(Index i, Index j)
{
	VectorType oldColI=colMap_[i];
	VectorType oldColJ=colMap_[j];

	colMap_[i].clear();
	colMap_[j].clear();

	for (VectorIt it=oldColI.begin();it!=oldColI.end();++it) {
		Index row=it->first;
		rowMap_[row].erase(rowMap_[row].find(i));
	}

	for (VectorIt it=oldColJ.begin();it!=oldColJ.end();++it) {
		Index row=it->first;
		rowMap_[row].erase(rowMap_[row].find(j));
	}

	for (VectorIt it=oldColI.begin();it!=oldColI.end();++it) {
		Index row=it->first;
		colMap_[j][row]=oldColI[row];
		rowMap_[row][j]=oldColI[row];
	}

	for (VectorIt it=oldColJ.begin();it!=oldColJ.end();++it) {
		Index row=it->first;
		colMap_[i][row]=oldColJ[row];
		rowMap_[row][i]=oldColJ[row];
	}
}

template<class Field_>
typename SparseMatrix<Field_,SparseMatrixFormat::SMM>::Index SparseMatrix<Field_,SparseMatrixFormat::SMM>::nnz() const
{
	return nnz_;
}

// A -> A' = SAS^{-1}, and A' has about nnz nonzero entries.
template<class Field_>
void SparseMatrix<Field_,SparseMatrixFormat::SMM>::randomSim(Index nz, int seed)
{	typename Field::Element a;
	Index i,j;
	MersenneTwister ri;
	typename Field::RandIter r(field(),0, seed);
	//if (seed != 0) { ri.setSeed(seed); r.setSeed(seed); }
	if (seed != 0)
	{	ri.setSeed(seed);
		// ridiculous constructor only seeding!
		typename Field::RandIter s(field(), 0, seed);
		r = s;
	}
	while (nnz() < nz) {
		nonzerorandom(field(),r,a);
		i = ri.randomIntRange(0, rowdim()); j = ri.randomIntRange(0, coldim());
		if (i!=j) {
			addCol(a, j, i);
			//std::cout << nnz() << std::endl;
			field().negin(a);
			addRow(a, i, j);
		}
	}
}

// A -> A' = UAV, with U and V nonsingular, and A' has about nnz nonzero entries.
template<class Field_>
void SparseMatrix<Field_,SparseMatrixFormat::SMM>::randomEquiv(Index nz, int seed)
{	typename Field::Element a;
	Index i,j;
	MersenneTwister ri;
	typename Field::RandIter r(field(),0,seed);
	if (seed != 0)
	{	ri.setSeed(seed);
		// ridiculous seeding!
		typename Field::RandIter s(field(), 0, seed);
		r = s;
	}
	bool flip = true;
	int count=0;
	while (nnz() < nz) {
		nonzerorandom(field(),r,a);
		i = ri.randomIntRange(0, rowdim()); j = ri.randomIntRange(0, coldim());
		if (i!=j){
			if (flip) addCol(a, i, j);
			else addRow(a, i, j);
			flip = not flip;
		}
		++count;
	}
}

template<class Field_>
std::ostream& SparseMatrix<Field_,SparseMatrixFormat::SMM>::print(std::ostream& out) const
{
	for (Index i=0;i<numRows_;++i) {
		for (Index j=0;j<numCols_;++j) {
			field().write(out,getEntry(i,j));
			if (j != (numCols_-1)) {
				out << " ";
			}
		}
		out << std::endl;
	}
        return out;
}

template<class Field_>
std::ostream& SparseMatrix<Field_,SparseMatrixFormat::SMM>::write(std::ostream& out) const
{
	writeMMCoordHeader(out,*this,nnz(),"SparseMatrix<SMM>","");

        for (MapConstIt p = rowMap_.begin(); p != rowMap_.end(); ++p)
                for (VectorConstIt rp = p->second.begin(); rp != p->second.end(); ++rp)
                        field().write(out << 1+p->first << " " << 1+rp->first << " ", rp->second) << std::endl;
        out << std::endl;
        return out;
}

template<class Field_>
std::istream& SparseMatrix<Field_,SparseMatrixFormat::SMM>::read(std::istream& in) {
        Index r,c;
        Element d;
        field().init(d);
        MatrixStream<Field_> ms(field(),in);
        ms.getDimensions(r,c);
        shape(r,c);
        while (ms.nextTriple(r,c,d)) setEntry(r,c,d);
        return in;
}

template<class Field_>
typename SparseMatrix<Field_,SparseMatrixFormat::SMM>::Index SparseMatrix<Field_,SparseMatrixFormat::SMM>::rowdim() const
{
        return numRows_;
}

template<class Field_>
typename SparseMatrix<Field_,SparseMatrixFormat::SMM>::Index SparseMatrix<Field_,SparseMatrixFormat::SMM>::coldim() const
{
        return numCols_;
}

template<class Field_>
void SparseMatrix<Field_,SparseMatrixFormat::SMM>::transpose()
{
        rowMap_.swap(colMap_);
        int temp=numCols_;numCols_=numRows_;numRows_=temp;
}

template<class Field_>
template<class Matrix>
void SparseMatrix<Field_,SparseMatrixFormat::SMM>::copy(Matrix& mat) const
{
        std::stringstream ss;
        write(ss);
        mat.read(ss);
        mat.finalize();
}

template<class Field_>
template<class Matrix>
void SparseMatrix<Field_,SparseMatrixFormat::SMM>::copyFrom(Matrix& mat)
{
        std::stringstream ss;
        mat.write(ss);
        read(ss);
}

template<class Field_>
template<class Vector>
void SparseMatrix<Field_,SparseMatrixFormat::SMM>::toVector(Vector& vec) const
{
        if (numCols_ == 1) {
	        vec.resize(numRows_);
                for (Index i=0;i<numRows_;++i) {
                        vec[i]=getEntry(i,0);
                }
        } else {
	        vec.resize(numCols_);
                for (Index j=0;j<numCols_;++j) {
                        vec[j]=getEntry(0,j);
                }
        }
}

template<class Field_>
void SparseMatrix<Field_,SparseMatrixFormat::SMM>::scaleMat(const Element& k)
{
        for (size_t i=0;i<numRows_;++i) {
                timesRow(k,i);
        }
}

template<class Field_>
bool SparseMatrix<Field_,SparseMatrixFormat::SMM>::areEqual(SparseMatrix<Field_,SparseMatrixFormat::SMM>& rhs) const
{
        if (rhs.numCols_ != numCols_) {
                return false;
        }
        if (rhs.numRows_ != numRows_) {
                return false;
        }
        if (nnz_ != rhs.nnz_) {
                return false;
        }

        for (MapConstIt rowIt=rowMap_.begin();rowIt!=rowMap_.end();++rowIt) {
                for (VectorConstIt eltIt=rowIt->second.begin();eltIt!=rowIt->second.end();++eltIt) {
                        if (eltIt->second != rhs.getEntry(rowIt->first,eltIt->first)) {
                                return false;
                        }
                }
        }
        return true;
}

template<class Field_>
void SparseMatrix<Field_,SparseMatrixFormat::SMM>::generateDenseRandMat(SparseMatrix<Field_,SparseMatrixFormat::SMM>& mat,int seed)
{
        typedef typename Field::Element Element;

        size_t m=mat.rowdim(),n=mat.coldim();
	Element d;

	typename Field::RandIter ri(field(),0,seed);

        for (size_t i=0;i<m;++i) {
                for (size_t j=0;j<n;++j) {
	                ri.random(d);
                        mat.setEntry(i,j,d);
                }
        }
}

template<class Field_>
void SparseMatrix<Field_,SparseMatrixFormat::SMM>::generateRandMat(SparseMatrix<Field_,SparseMatrixFormat::SMM>& mat, int nnz, int seed)
{
        typedef typename Field::Element Element;

        
        size_t m=mat.rowdim(),n=mat.coldim();
	Element d;
	typename Field::RandIter ri(field(),0,seed);
	srand(seed);

        typedef std::pair<size_t,size_t> CoordPair;
        typedef std::set<CoordPair> PairSet;
        PairSet pairs;

	for(int i = 0; i < (int)nnz; ++i) {
                size_t row,col;
                do {
                        row = randRange(0,(int)m);
                        col = randRange(0,(int)n);
                } while (pairs.count(CoordPair(row,col))!=0);

                nonzerorandom(field(),ri,d);
                mat.setEntry(row,col,d);
                pairs.insert(CoordPair(row,col));
        }
}

template<class Field_>
void SparseMatrix<Field_,SparseMatrixFormat::SMM>::
generateScaledIdent(SparseMatrix<Field_,SparseMatrixFormat::SMM>& mat,
                    const typename Field_::Element& d)
{
        size_t m=mat.rowdim(),n=mat.coldim();
        size_t minDim=(m<n)?m:n;

        for (size_t i=0;i<minDim;++i) {
                mat.setEntry(i,i,d);
        }
}

template<class Field_>
void SparseMatrix<Field_,SparseMatrixFormat::SMM>::generateSparseNonSingular(SparseMatrix<Field_,SparseMatrixFormat::SMM>& mat, int approxNNZ, int seed)
{
        typedef typename Field::Element Element;
        int n=mat.rowdim();
        linbox_check(mat.rowdim()==mat.coldim());

        typename Field::RandIter r(mat.field(),0,seed);

	Element d;

        for (int i=0;i<n;++i) {
	        nonzerorandom(field(),r,d);
                mat.setEntry(i,i,d);
        }

        mat.randomEquiv(approxNNZ,seed);
}

template<class Field_>
void SparseMatrix<Field_,SparseMatrixFormat::SMM>::generateCompanion(SparseMatrix<Field_,SparseMatrixFormat::SMM>& mat,std::vector<typename Field::Element>& coeffs)
{
        typedef typename Field::Element Element;
        int n=mat.rowdim();
        linbox_check(mat.rowdim()==mat.coldim());
        linbox_check(n==coeffs.size());

	Element d;
	mat.field().assign(d,mat.field().one);
        for (int i=1;i<n;++i) {
                mat.setEntry(i,i-1,d);
        }
	for (int i=0;i<n;++i) {
		mat.setEntry(i,n-1,coeffs[i]);
	}
}

template<class Field_>
int SparseMatrix<Field_,SparseMatrixFormat::SMM>::randRange(int start, int end)
{
        double rval = rand();
        static const double NORMALIZING_CONSTANT = 1.0/(1.0+RAND_MAX);
        double normedRVal = rval*NORMALIZING_CONSTANT;
        double rangeSize = end-start;
        int offset = (int)(rangeSize*normedRVal);
        return start+offset;
}

template<class Field_>
typename Field_::Element&
SparseMatrix<Field_,SparseMatrixFormat::SMM>::nonzerorandom(const Field_& F,typename Field_::RandIter&r,typename Field_::Element& e)
{
	do {
		r.random(e);
	} while (F.isZero(e));
	return e;
}

}

#endif // __LINBOX_SPARSE_MAP_MAP_MATRIX_INL

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
