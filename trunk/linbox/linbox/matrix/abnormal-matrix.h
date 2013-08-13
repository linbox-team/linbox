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

/*! @file   linbox/matrix/abnormal-matrix.h
 * @ingroup matrix
 * @brief
 */

#ifndef __LINBOX_ABNORMAL_MATRIX
#define __LINBOX_ABNORMAL_MATRIX

#include <stdlib.h>
#include <fstream>

namespace LinBox
{

template <class Field>
class AbnormalHelper {
public:
        typedef typename Field::Element Element;
	typedef typename Field::Element Abnormal;

	AbnormalHelper () {}

	AbnormalHelper(const Field& F) : field_(&F) {}

	void init(const Field& field) { field_=&field; }

	inline Abnormal& mulacc(Abnormal& x, const Element& y, const Element& z) const {
		field_->mulacc(x,y,z);
		return x;
	}

	inline Element normalize(Abnormal& elt) const {
		Element d;
		field_->init(d,elt);
		return d;
	}

protected:
	Field* field_;
};




template <>
class AbnormalHelper<Modular<double> > {
public:
	typedef Modular<double> Field;
        typedef double Element;
	typedef double Abnormal;

	AbnormalHelper () {}

	AbnormalHelper(const Field& field) {init(field);}

	void init(const Field& field) {
		modulus_=field.characteristic();
		unsigned long long maxDouble = 1ULL<<52;
		bound_=(double)maxDouble;
		field_=&field;
	}

	inline Abnormal& mulacc(Abnormal& x, const Element& y, const Element& z) const
        {
		Element d;
		d=y*z;//at most (2**26-1)**2
		x=x+d;//at most (2**26-1)**2+(2**52-1)<2**53
                maybeNormalize(x);
                return x;
        }

	inline double normalize(double& elt) const
	{
		return elt=fmod(elt,modulus_);
	}

protected:
	inline double& maybeNormalize(double& elt) const
	{
		if (elt>=bound_) {
			elt=fmod(elt,modulus_);
		}
		return elt;
	}

	double modulus_;

	double bound_;

	const Field* field_;
};

template <class T, size_t ALIGNMENT>
class AlignedMat {
public:
	typedef size_t Index;

	AlignedMat() : space_(NULL), elts_(NULL),
                       rowDim_(0), colDim_(0), stride_(0) {}
		       

	AlignedMat(const AlignedMat<T,ALIGNMENT>& rhs) {
		shape(rhs.rowDim_,rhs.colDim_);
		for (Index i=0;i<rowDim_;++i) {
			for (Index j=0;j<colDim_;++j) {
				*refElt(i,j)=*rhs.refElt(i,j);
			}
		}
	}

	AlignedMat& operator=(const AlignedMat &rhs) {
		if (this == &rhs) {
			return *this;
		}
		shape(rhs.rowDim_,rhs.colDim_);
		for (Index i=0;i<rowDim_;++i) {
			for (Index j=0;j<colDim_;++j) {
				*refElt(i,j)=*rhs.refElt(i,j);
			}
		}
	}

	~AlignedMat() {
		if (space_ != NULL) {
			for (int i=0;i<rowDim_;++i) {
				for (int j=0;j<colDim_;++j) {
					refElt(i,j)->~T();
				}
			}
			delete[] space_;
		}
	}

	void shape(size_t r, size_t c) {
		stride_=roundUp(r*sizeof(T));
		rowDim_=r;
		colDim_=c;
		space_=new uint8_t[stride_*c+ALIGNMENT];
		size_t spacePtr=(size_t)space_;
		elts_=(uint8_t*)(roundUp(spacePtr));
		for (int i=0;i<rowDim_;++i) {
			for (int j=0;j<colDim_;++j) {
				new ((void*)refElt(i,j)) T();
			}
		}
	}

	inline T* refElt(size_t i, size_t j) {
		return (T*)(elts_+j*stride_+sizeof(T)*i);
	}

protected:
	inline static size_t roundUp(size_t n) {
		if (ALIGNMENT == 0) {
			return n;
		}
		return n+ALIGNMENT-(n%ALIGNMENT);
	}

	uint8_t *space_;
	uint8_t *elts_;
	Index rowDim_, colDim_, stride_;
};

template <class Field, bool ROW_MAJOR=true,size_t ALIGNMENT=0>
class AbnormalSubmatrixIterator;

template <class Field, bool ROW_MAJOR=true,size_t ALIGNMENT=0>
class AbnormalMatrix
{
public:
        class AbnormalSubmatrix;

	typedef size_t Index;
        typedef typename Field::Element Element;
	typedef typename AbnormalHelper<Field>::Abnormal Abnormal;
	typedef AbnormalMatrix<Field,ROW_MAJOR,ALIGNMENT> AbnormalMat;
        typedef AbnormalSubmatrix Submat;
        typedef AbnormalSubmatrixIterator<Field,ROW_MAJOR,ALIGNMENT> Iterator;

        AbnormalMatrix();

	AbnormalMatrix(const Field& F, Index r, Index c);

	void init(const Field& F, Index r, Index c);

	void shape(Index r, Index c);

	void setEntry(Index i, Index j, const Element& e);

	Element &getElement(Index i, Index j, Element& e);

        AbnormalSubmatrix submat(Index startRow, Index startCol,
                                 Index endRow, Index endCol);

        AbnormalSubmatrix wholeSlice();

        AbnormalSubmatrix rowSlice(Index i);

        AbnormalSubmatrix colSlice(Index j);

	class AbnormalSubmatrix {
	public:
                template <class Matrix>
                void normalize(Matrix& mat);

		template <class Matrix>
		void saxpyin(const Element& e, Matrix& mat);

                Iterator begin();

                Iterator end();

                friend Iterator;
                friend AbnormalMat;

	protected:
		AbnormalSubmatrix(AbnormalMat *mat,
				  Index startRow, Index startCol,
				  Index endRow, Index endCol) :
		mat_(mat), startRow_(startRow), startCol_(startCol),
		endRow_(endRow), endCol_(endCol) {}

                AbnormalMat *mat_;

                Index startRow_, startCol_, endRow_, endCol_;
	};

        friend Iterator;

protected:
        inline Element& refEntry(Index i, Index j);

        AlignedMat<Abnormal,ALIGNMENT> elts_;

        Index rowDim_, colDim_;

	AbnormalHelper<Field> abnormalHelper_;
};

template <class Field, bool ROW_MAJOR,size_t ALIGNMENT>
class AbnormalSubmatrixIterator
{
public:
	typedef size_t Index;
        typedef typename Field::Element Element;
	typedef typename AbnormalHelper<Field>::Abnormal Abnormal;
	typedef AbnormalMatrix<Field,ROW_MAJOR,ALIGNMENT> AbnormalMat;
        typedef typename AbnormalMat::AbnormalSubmatrix Submat;
        typedef AbnormalSubmatrixIterator<Field,ROW_MAJOR,ALIGNMENT> Iterator;

        friend AbnormalMat;
        friend Submat;
        
        inline AbnormalSubmatrixIterator& operator++() {
                if (ROW_MAJOR) {
                        ++col_;
                        if (col_>submat_->endCol_) {
                                col_=submat_->startCol_;
                                ++row_;
                        }
                } else {
                        ++row_;
                        if (row_>submat_->endRow_) {
                                row_=submat_->startRow_;
                                ++col_;
                        }
                }
                return *this;
        }

        inline Abnormal& operator*() {
                return mat_->refEntry(row_,col_);
        }

        inline Abnormal* operator->() {
                return &(mat_->refEntry(row_,col_));
        }

        inline bool operator==(Iterator rhs) {
                linbox_check(submat_==rhs.submat_);
                return ((row_==rhs.row_) && (col_==rhs.col_));
        }

        inline bool operator!=(Iterator rhs) {
                return ((row_!=rhs.row_) || (col_!=rhs.col_));
        }

protected:
        AbnormalSubmatrixIterator(Submat *submat,Index row, Index col)
                : submat_(submat), mat_(submat->mat_), row_(row), col_(col) {}

        Submat *submat_;

        AbnormalMat *mat_;

        Index row_, col_;
};
}

#include "linbox/matrix/abnormal-matrix.inl"

#endif // __LINBOX_ABNORMAL_MATRIX

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
