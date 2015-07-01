/* -*- mode: C++; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*- */

/* linbox/blackbox/quad-matrix.h
 * Copyright (C) 2006 by -bds, hui wang  
 * evolved from scalar-matrix.h
 */

#ifndef __QUAD_MATRIX_H
#define __QUAD_MATRIX_H

#include <algorithm>
#include <linbox/field/hom.h>
#include <vector>
#include <iterator>
//#include "linbox/vector/vector-traits.h"
#include "linbox/util/debug.h"
#include "linbox-config.h"
#include <linbox/blackbox/blackbox-interface.h>
#include <linbox/blackbox/scalar-matrix.h>
#include <linbox/blackbox/zo.h>
//#include <linbox/blackbox/side-by-side.h>
//#include <linbox/blackbox/over-under.h>

namespace LinBox {

template <typename Field>
class ZeroOne;

template <typename Field>
class ScalarMatrix;

template <typename Field>
class SideBySide;

template <typename Field>
class OverUnder;

// this ZOQuad class is, in effect,  a tagged union of several BB types.
// Then SideBySide does not have it's two part types as template parameters.  
// Instead the parts are always ZOQuads.

template <typename Field>
class ZOQuad
{   
	enum BBtype {zo, sm, sbs, ou};
	const void* _BBp;	
	BBtype _tag;
	size_t _r, _c;
	static const unsigned int smallThreshold = 1000;
    public:

    typedef size_t Index;
	typedef std::vector<Index> IndexVector;
    typedef std::vector<IndexVector::iterator> PointerVector;

	// constructors
	ZOQuad() : _BBp(0), _r(0), _c(0) {}
	ZOQuad(const ScalarMatrix<Field>&  A) : _BBp(&A), _tag(sm), _r(A.rowdim()), _c(A.coldim()) {} 
	ZOQuad(const SideBySide<Field>&  A) : _BBp(&A), _tag(sbs), _r(A.rowdim()), _c(A.coldim()) {} 
	ZOQuad(const OverUnder<Field>&  A) : _BBp(&A), _tag(ou), _r(A.rowdim()), _c(A.coldim()) {} 
	ZOQuad(const ZeroOne<Field>&  A) : _BBp(&A), _tag(zo), _r(A.rowdim()), _c(A.coldim()) {} 

	void init (const ZeroOne<Field>&  A)
	{
		if ( A.nnz() == 0 ) { A = NULL;	return;	}

		if ( A.nnz() <= smallThreshold ) { _BBp = &A; _tag = zo; return; }

		// if A has more rows than cols, make an OverUnder
		if ( A.coldim() <= A.rowdim() )
		{
			IndexVector firstHalf;
			PointerVector firstHalfP;
			std::back_insert_iterator< IndexVector > first(firstHalf);
			std::back_insert_iterator< PointerVector > firstP(firstHalfP);
			copy( A._index.begin(), *(A._indexP.begin() + A._indexP.size()/2), first);
			copy( A._indexP.begin(), A._indexP.begin() + A._indexP.size()/2 + 1, firstP);
			ZeroOne<Field> U(A.field(), firstHalf, firstHalfP, A._coldim, A.sorted);

			IndexVector secondHalf;
			PointerVector secondHalfP;
			std::back_insert_iterator< IndexVector > second(secondHalf);
			std::back_insert_iterator< PointerVector > secondP(secondHalfP);
			copy( *(A._indexP.begin() + A._indexP.size()/2), A._index.end(), second);
			copy( A._indexP.begin() + A._indexP.size()/2, A._indexP.end(), secondP);
			ZeroOne<Field> D(A.field(), secondHalf, secondHalfP, A._coldim, A.sorted);

			ZOQuad<Field> UU(U);
			ZOQuad<Field> DD(D);
			
			//Hui works this out!

			_BBp = new OverUnder<Field>(UU, DD); _tag = sbs; return;
		}
		// else make an sideBySide
		/*
		else
		{
			IndexVector firstHalf;
			PointerVector firstHalfP;
			std::back_insert_iterator< IndexVector > first(firstHalf);
			std::back_insert_iterator< PointerVector > firstP(firstHalfP);
			copy( A._index.begin(), *(A._indexP.begin() + A._indexP.size()/2), first);
			copy( A._indexP.begin(), A._indexP.begin() + A._indexP.size()/2 + 1, firstP);
			ZeroOne<Field> U(A.field(), firstHalf, firstHalfP, A._coldim, A.sorted);

			IndexVector secondHalf;
			PointerVector secondHalfP;
			std::back_insert_iterator< IndexVector > second(secondHalf);
			std::back_insert_iterator< PointerVector > secondP(secondHalfP);
			copy( *(A._indexP.begin() + A._indexP.size()/2), A._index.end(), second);
			copy( A._indexP.begin() + A._indexP.size()/2, A._indexP.end(), secondP);
			ZeroOne<Field> D(A.field(), secondHalf, secondHalfP, A._coldim, A.sorted);

			ZOQuad<Field> UU(U);
			ZOQuad<Field> DD(D);
			
			//Hui works this out!

		  _BBp = new SideBySide<Field>(O, U); _tag = ou; return;
		}
		*/
	} 

	~ZOQuad() 
	{	switch (_tag)
		{	case zo: 
				{	const ZeroOne<Field> *A(static_cast< const ZeroOne<Field>* >(_BBp));
					delete A ; break;
				}
		 	case sm: 
				{	const ScalarMatrix<Field> *A(static_cast< const ScalarMatrix<Field>* >(_BBp));
					delete A ; break;
				}
		 	case sbs:
				{	const SideBySide<Field> *A(static_cast< const SideBySide<Field>* >(_BBp));
					delete A ; break;
				}
		 	case ou:
				{	const OverUnder<Field> *A(static_cast< const OverUnder<Field>* >(_BBp));
					delete A ; break;
				}
		}
	}

	template <typename InVector, typename OutVector>
	OutVector& apply(OutVector& y, const InVector& x) const  
	//OutVector& apply(OutVector& y, const InVector& x)
	{	switch (_tag)
		{	case zo: 
				{
					const ZeroOne<Field> *A(static_cast< const ZeroOne<Field>* >(_BBp));
					//ZeroOne<Field> *A(static_cast< ZeroOne<Field>* >(_BBp));
					A->apply(y, x); 
					break;
				}
		 	case sm: 
				{
					const ScalarMatrix<Field> *A(static_cast< const ScalarMatrix<Field>* >(_BBp));
					//ScalarMatrix<Field> *A(static_cast< ScalarMatrix<Field>* >(_BBp));
					A->apply(y, x);
					break;
				}
		 	case sbs:
				{
					const SideBySide<Field> *A(static_cast< const SideBySide<Field>* >(_BBp));
					//SideBySide<Field> *A(static_cast< SideBySide<Field>* >(_BBp));
					A->apply(y, x);
					break;
				}
		 	case ou:
				{ 
					const OverUnder<Field> *A(static_cast< const OverUnder<Field>* >(_BBp));
					//OverUnder<Field> *A(static_cast< OverUnder<Field>* >(_BBp));
					A->apply(y, x);
					break;
				}
		}
		return y;
	}
	
	// similar applyTranspose
	template <typename InVector, typename OutVector>
	OutVector& applyTranspose(OutVector& y, const InVector& x) const
	//OutVector& applyTranspose(OutVector& y, const InVector& x)
	{	switch (_tag)
		{	case zo: 
				{
					//const ZeroOne<Field> *A(static_cast< const ZeroOne<Field>* >(_BBp));
					const ZeroOne<Field> *A(static_cast< const ZeroOne<Field>* >(_BBp));
					A->applyTranspose(y, x); 
					break;
				}
		 	case sm: 
				{
					const ScalarMatrix<Field> *A(static_cast< const ScalarMatrix<Field>* >(_BBp));
					//ScalarMatrix<Field> *A(static_cast< ScalarMatrix<Field>* >(_BBp));
					A->applyTranspose(y, x);
					break;
				}
		 	case sbs:
				{
					const SideBySide<Field> *A(static_cast< const SideBySide<Field>* >(_BBp));
					//SideBySide<Field> *A(static_cast< SideBySide<Field>* >(_BBp));
					A->applyTranspose(y, x);
					break;
				}
		 	case ou:
				{
					const OverUnder<Field> *A(static_cast< const OverUnder<Field>* >(_BBp));
					//OverUnder<Field> *A(static_cast< OverUnder<Field>* >(_BBp));
					A->applyTranspose(y, x);
					break;
				}
		}
		return y;
	}

	size_t rowdim() const {return _r;}
	size_t coldim() const {return _c;}
	const Field& field() const
	{	switch (_tag)
		{	case zo: 
				{
					const ZeroOne<Field> *A(static_cast< const ZeroOne<Field>* >(_BBp));
					return A->field();
				}
		 	case sm: 
				{
					const ScalarMatrix<Field> *A(static_cast< const ScalarMatrix<Field>* >(_BBp));
					return A->field();
				}
		 	case sbs: 
				{
					const SideBySide<Field> *A(static_cast< const SideBySide<Field>* >(_BBp));
					return A->field();
				}
		 	case ou: 
				{
					const OverUnder<Field> *A(static_cast< const OverUnder<Field>* >(_BBp));
					return A->field();
				}
		}
		// silly business to avoid a warning
		const ScalarMatrix<Field> *A(static_cast< const ScalarMatrix<Field>* >(_BBp));
		return A->field();
	}
};

template <typename Field>
class SideBySide
{
	typedef ZOQuad<Field> Quad;
	const Quad *_L, *_R;
    public:
	SideBySide(const Quad& A, const Quad& B) : _L(&A), _R(&B) {}

	template <typename InVector, typename OutVector>
	OutVector& apply(OutVector& y, const InVector& x) const
	//OutVector& apply(OutVector& y, const InVector& x)
	{
			std::vector<typename Field::Element> z(y.size());
			VectorDomain<Field> VD(field());
			std::vector<typename Field::Element> x_1( x.begin(), x.begin() + _L->coldim() );
			std::vector<typename Field::Element> x_2( x.begin() + _L->coldim(), x.end() );
			//Subvector<typename InVector::iterator, typename InVector::const_iterator> x_1(x.begin(), x.begin()+_L->coldim());
			//Subvector<typename InVector::iterator, typename InVector::const_iterator> x_2(x.begin()+_L->coldim(), x.end());
			_L->apply (y, x_1);
			_R->apply (z, x_2);
			VD.addin(y, z);

			return y;
	}

	template <typename InVector, typename OutVector>
	OutVector& applyTranspose(OutVector& y, const InVector& x) const
	//OutVector& applyTranspose(OutVector& y, const InVector& x)
	{
			std::vector<typename Field::Element> y_1( y.begin(), y.begin() + _L->coldim() );
			std::vector<typename Field::Element> y_2( y.begin() + _L->coldim(), y.end() );
			//Subvector<typename OutVector::iterator, typename OutVector::const_iterator> y_1(y.begin(), y.begin()+_L->coldim());
			//Subvector<typename OutVector::iterator, typename OutVector::const_iterator> y_2(y.begin()+_L->coldim(), y.end());
			_L->applyTranspose (y_1, x);
			_R->applyTranspose (y_2, x);
			copy(y_1.begin(), y_1.end(), y.begin());
			copy(y_2.begin(), y_2.end(), y.begin() + y_1.size());
			return y;
	}

	size_t rowdim()const{return _L->rowdim();}
	size_t coldim()const{return _L->coldim() + _R->coldim();}
	const Field& field()const {return _L->field();}	
};

template <typename Field>
class OverUnder
{
	const ZOQuad<Field> *_U, *_D;
    public:
	OverUnder(const ZOQuad<Field>& A, const ZOQuad<Field>& B) : _U(&A), _D(&B) {}

	template <typename InVector, typename OutVector>
	OutVector& apply(OutVector& y, const InVector& x) const
	//OutVector& apply(OutVector& y, const InVector& x)
	{
			std::vector<typename Field::Element> y_1( y.begin(), y.begin() + _U->rowdim() );
			std::vector<typename Field::Element> y_2( y.begin() + _U->rowdim(), y.end() );
		  	//Subvector<typename OutVector::iterator, typename OutVector::const_iterator> y_1(y.begin(), y.begin()+_U->rowdim());
			//Subvector<typename OutVector::iterator, typename OutVector::const_iterator> y_2(y.begin()+_U->rowdim(), y.end());
			//if ((_A_ptr == 0) || (_B_ptr == 0)) { throw error }
			_U->apply (y_1, x);
			_D->apply (y_2, x);
			copy(y_1.begin(), y_1.end(), y.begin());
			copy(y_2.begin(), y_2.end(), y.begin() + y_1.size());
					 
			return y;
	}

	template <typename InVector, typename OutVector>
	OutVector& applyTranspose(OutVector& y, const InVector& x) const
	//OutVector& applyTranspose(OutVector& y, const InVector& x)
	{
			std::vector<typename Field::Element> z(y.size());
			VectorDomain<Field> VD(field());
			std::vector<typename Field::Element> x_1( x.begin(), x.begin() + _U->rowdim() );
			std::vector<typename Field::Element> x_2( x.begin() + _U->rowdim(), x.end() );
			//Subvector<typename InVector::iterator, typename InVector::const_iterator> x_1(x.begin(), x.begin()+_U->rowdim());
			//Subvector<typename InVector::iterator, typename InVector::const_iterator> x_2(x.begin()+_U->rowdim(), x.end());
			_U->applyTranspose (y, x_1);
			_D->applyTranspose (z, x_2);
			VD.addin(y, z);
			return y;
	}

	size_t coldim() const{return _U->coldim();}
	size_t rowdim() const{return _U->rowdim() + _D->rowdim();}
	const Field& field() const {return _U->field();}	
};

// similar class OverUnder<Field>

}; //namespace LinBox 
#endif // __QUAD_MATRIX_H
