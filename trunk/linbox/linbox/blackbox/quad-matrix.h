/* linbox/blackbox/quad-matrix.h
 * Copyright (C) 2006 LinBox
 * Written by -bds, hui wang  
 * see COPYING for license information
 */

#ifndef __LINBOX_quad_matrix_H
#define __LINBOX_quad_matrix_H

#include <algorithm>
#include <linbox/field/hom.h>
#include <vector>
#include <iterator>
//#include "linbox/vector/vector-traits.h"
#include "linbox/util/debug.h"
#include "linbox/linbox-config.h"
#include <linbox/blackbox/blackbox-interface.h>
#include <linbox/blackbox/scalar-matrix.h>
#include <linbox/blackbox/zo.h>
//#include <linbox/blackbox/side-by-side.h>
//#include <linbox/blackbox/over-under.h>

namespace LinBox 
{/*{{{*/

	template <typename Field>
	class ZeroOne;

	template <typename Field>
	class ScalarMatrix;

	template <typename Field>
	class SideBySide;

	template <typename Field>
	class OverUnder;

/** A class of striped or block-decomposed zero-one matrices.

 Right now we implement only side-by-side blocks, which is done to improve cache performance.  This class works especially well on quite sparse matrices whose entries are rather randomly positioned.  The idea is that the vector to be multiplied by each block fits in cache.  It is being more or less randomly accessed, so cache misses are avoided for these random accesses.

 The possibility exists to set the stripe(block) width by a configure-time parameter.  Also this approach should be used on sparse matrices which are not zero-one.
 
 This ZOQuad class is like a tagged union of several BB types.
 Then SideBySide does not have it's two part types as template parameters.  
 Instead the parts are always ZOQuads.

 \ingroup blackbox

*/
	template <typename _Field>
	class ZOQuad
	{
		//friend class ZeroOne<_Field>;
		enum BBtype {zo, sbs, ou};
		const void* _BBp;	
		BBtype _tag;
		size_t _r, _c;
		static const unsigned int smallThreshold = 60000;

	protected:
		_Field _F;
    public:
		typedef _Field Field;
		
		typedef size_t Index;
		typedef std::vector<Index> IndexVector;
		typedef std::vector<IndexVector::iterator> PointerVector;
		
		// constructors
		ZOQuad() : _BBp(0), _r(0), _c(0) {}
		//ZOQuad(const ScalarMatrix<Field>&  A) : _BBp(&A), _tag(sm), _r(A.rowdim()), _c(A.coldim()) {} 
		ZOQuad(const SideBySide<Field>& A) : _BBp(&A), _tag(sbs), _r(A.rowdim()), _c(A.coldim()) {} 
		ZOQuad(const OverUnder<Field>& A) : _BBp(&A), _tag(ou), _r(A.rowdim()), _c(A.coldim()) {} 
		ZOQuad(const ZeroOne<Field>& A) //: _BBp(&A), _tag(zo), _r(A.rowdim()), _c(A.coldim())
		{
			init(A);
		}

		int init (const ZeroOne<Field>&  A)
		{
			const Field& F = A.field();
			//if ( A.nnz() == 0 )
			if ( A.coldim() <= smallThreshold )
				{
					_BBp = new ZeroOne<Field>(A);
					_tag = zo;
					_r = A.rowdim(); _c = A.coldim();
					const ZeroOne<Field> *B(static_cast< const ZeroOne<Field>* >(_BBp));
					if( !(B->sorted) )
						{
					//		std::cout << " -- init: switch sort " << A.sorted << " " << &(A.sorted) << std::endl;
							B->switch_sort();
					//		std::cout << " -- init: switch sort " << A.sorted << std::endl;
						}
					//std::cout << " -- quad init: " << A.nnz() << " ( zo: " << &A << ", " << _BBp << ", quad: " << this << ")" << std::endl;
 					return 0;
				}
			// if A has more rows than cols, make an OverUnder
			/*
			if ( A.coldim() <= A.rowdim() )
				{
					if( !A.sorted ) A.switch_sort();

					IndexVector firstHalf;
					PointerVector firstHalfP;
					std::back_insert_iterator< IndexVector > first(firstHalf);
					std::back_insert_iterator< PointerVector > firstP(firstHalfP);
					copy( A._index.begin(), *(A._indexP.begin() + A._indexP.size()/2), first);
					copy( A._indexP.begin(), A._indexP.begin() + A._indexP.size()/2 + 1, firstP);
					ZeroOne<Field> U(F, firstHalf, firstHalfP, firstHalfP.size()-1, A._coldim, A.sorted);
					ZOQuad<Field> *UU = new ZOQuad<Field>(U);
					
					IndexVector secondHalf;
					PointerVector secondHalfP;
					std::back_insert_iterator< IndexVector > second(secondHalf);
					std::back_insert_iterator< PointerVector > secondP(secondHalfP);
					copy( *(A._indexP.begin() + A._indexP.size()/2), A._index.end(), second);
					copy( A._indexP.begin() + A._indexP.size()/2, A._indexP.end(), secondP);
					ZeroOne<Field> D(F, secondHalf, secondHalfP, secondHalfP.size()-1, A._coldim, A.sorted);
					ZOQuad<Field> *DD = new ZOQuad<Field>(D);
					
					//Hui works this out!
					
					_BBp = new OverUnder<Field>(UU, DD);

					_tag = ou;
					_r = A.rowdim(); _c = A.coldim();
 					
					return 1;
				}
			*/
			// else make an sideBySide
			///*
			else
				{
					if( A.sorted ) A.switch_sort();

					IndexVector firstHalf;
					PointerVector firstHalfP;
					std::back_insert_iterator< IndexVector > first(firstHalf);
					std::back_insert_iterator< PointerVector > firstP(firstHalfP);
					copy( A._index.begin(), *(A._indexP.begin() + A._indexP.size()/2), first);
					copy( A._indexP.begin(), A._indexP.begin() + A._indexP.size()/2 + 1, firstP);
					ZeroOne<Field> L(F, firstHalf, firstHalfP, A._rowdim, firstHalfP.size()-1, A.sorted);
					ZOQuad<Field> *LL = new ZOQuad<Field>(L);
					
					IndexVector secondHalf;
					PointerVector secondHalfP;
					std::back_insert_iterator< IndexVector > second(secondHalf);
					std::back_insert_iterator< PointerVector > secondP(secondHalfP);
					copy( *(A._indexP.begin() + A._indexP.size()/2), A._index.end(), second);
					copy( A._indexP.begin() + A._indexP.size()/2, A._indexP.end(), secondP);
					ZeroOne<Field> R(F, secondHalf, secondHalfP, A._rowdim, secondHalfP.size()-1, A.sorted);
					ZOQuad<Field> *RR = new ZOQuad<Field>(R);
					
					//Hui works this out!
					
					_BBp = new SideBySide<Field>(LL, RR);

					_tag = sbs;
					_r = A.rowdim(); _c = A.coldim();
 					
					return 2;
				  }
			//*/
		} 
		std::ostream & write(std::ostream & out) const
		{
			switch (_tag)
				{
				case zo: 
					{
						const ZeroOne<Field> *A(static_cast< const ZeroOne<Field>* >(_BBp));
						out << "zo(" << A->rowdim() << ", " << A->coldim() << ", " << A->nnz() << ") ";
						break;
					}
					//case sm: 
					//{
					//const ScalarMatrix<Field> *A(static_cast< const ScalarMatrix<Field>* >(_BBp));
					//delete A ; break;
					//}
				case sbs:
					{
						out << "sidebyside(";
						const SideBySide<Field> *A(static_cast< const SideBySide<Field>* >(_BBp));
						(A->_L)->write(out) << ", ";
						(A->_R)->write(out) << ") ";
						break;
					}
				case ou:
					{
						out << "overunder(";
						const OverUnder<Field> *A(static_cast< const OverUnder<Field>* >(_BBp));
						(A->_U)->write(out) << ", ";
						(A->_D)->write(out) << ") ";
						break;
					}
				}
			return out;
		}

		~ZOQuad() 
		{/*	switch (_tag)
			{
			case zo: 
				{
					const ZeroOne<Field> *A(static_cast< const ZeroOne<Field>* >(_BBp));
					delete A ; 
					break;
				}
				//case sm: 
				//{
					//const ScalarMatrix<Field> *A(static_cast< const ScalarMatrix<Field>* >(_BBp));
					//delete A ; break;
				//}
		    case sbs:
				{
					const SideBySide<Field> *A(static_cast< const SideBySide<Field>* >(_BBp));
					delete A ; break;
				}
		 	case ou:
				{
					const OverUnder<Field> *A(static_cast< const OverUnder<Field>* >(_BBp));
					delete A ; break;
				}
			}*/
		}
		
		template <typename InVector, typename OutVector>
		OutVector& apply(OutVector& y, const InVector& x) const  
			//OutVector& apply(OutVector& y, const InVector& x)
		{
			//std::cout << " zo-quad apply size of x: " << x.size() << " " << " size of y: " << y.size() << endl;
			switch (_tag)
				{
				case zo: 
					{
						const ZeroOne<Field> *A(static_cast< const ZeroOne<Field>* >(_BBp));
						//ZeroOne<Field> *A(static_cast< ZeroOne<Field>* >(_BBp));
						//std::cout << "\n apply zo quad -- zo " << "( quad: " << this << ", zo: " << _BBp << "), (" << _r << ", " << _c << "), (" << A->rowdim() << ", " << A->coldim() << ")" << endl;
						A->apply(y, x); 
						break;
					}
					//case sm: 
					//{
					//const ScalarMatrix<Field> *A(static_cast< const ScalarMatrix<Field>* >(_BBp));
						//ScalarMatrix<Field> *A(static_cast< ScalarMatrix<Field>* >(_BBp));
						//A->apply(y, x);
						//break;
					//}
				case sbs:
					{
						const SideBySide<Field> *A(static_cast< const SideBySide<Field>* >(_BBp));
						//SideBySide<Field> *A(static_cast< SideBySide<Field>* >(_BBp));
						//std::cout << "\n apply zo quad -- sbs " << _r << " " << _c << endl;
						A->apply(y, x);
						break;
					}
				case ou:
					{ 
						const OverUnder<Field> *A(static_cast< const OverUnder<Field>* >(_BBp));
						//OverUnder<Field> *A(static_cast< OverUnder<Field>* >(_BBp));
						//std::cout << "\n apply zo quad -- ou " << _r << " " << _c << endl;
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
		{
			switch (_tag)
				{
				case zo: 
					{
						//const ZeroOne<Field> *A(static_cast< const ZeroOne<Field>* >(_BBp));
						const ZeroOne<Field> *A(static_cast< const ZeroOne<Field>* >(_BBp));
						A->applyTranspose(y, x); 
						break;
					}
					//case sm: 
					//{
					//const ScalarMatrix<Field> *A(static_cast< const ScalarMatrix<Field>* >(_BBp));
						//ScalarMatrix<Field> *A(static_cast< ScalarMatrix<Field>* >(_BBp));
						//A->applyTranspose(y, x);
						//break;
					//}
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
		
		size_t rowdim() const
		{
			switch (_tag)
				{
				case zo: 
					{
						const ZeroOne<Field> *A(static_cast< const ZeroOne<Field>* >(_BBp));
						return A->rowdim();
					}
				case sbs: 
					{
						const SideBySide<Field> *A(static_cast< const SideBySide<Field>* >(_BBp));
						return A->rowdim();
					}
				case ou: 
					{
						const OverUnder<Field> *A(static_cast< const OverUnder<Field>* >(_BBp));
						return A->rowdim();
					}
				}
			return 0;
		}
		size_t coldim() const
		{
			switch (_tag)
				{
				case zo: 
					{
						const ZeroOne<Field> *A(static_cast< const ZeroOne<Field>* >(_BBp));
						return A->coldim();
					}
				case sbs: 
					{
						const SideBySide<Field> *A(static_cast< const SideBySide<Field>* >(_BBp));
						return A->coldim();
					}
				case ou: 
					{
						const OverUnder<Field> *A(static_cast< const OverUnder<Field>* >(_BBp));
						return A->coldim();
					}
				}
			return 0;
		}
		const Field& field() const
		{
			switch (_tag)
				{
				case zo: 
					{
						const ZeroOne<Field> *A(static_cast< const ZeroOne<Field>* >(_BBp));
						return A->field();
					}
					//case sm: 
					//{
					//const ScalarMatrix<Field> *A(static_cast< const ScalarMatrix<Field>* >(_BBp));
					//return A->field();
					//}
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
			// silly business to avoid a warning ( no return statement warning )
			const ScalarMatrix<Field> *A(static_cast< const ScalarMatrix<Field>* >(_BBp));
			return A->field();
		}
	}; //ZOQuad
	
	template <typename Field>
	class SideBySide
	{
		typedef ZOQuad<Field> Quad;
	public://temp
		const Quad *_L, *_R;
    public:
		SideBySide(const Quad* A, const Quad* B) : _L(A), _R(B) {}
		//~SideBySide() {delete _L; delete _R;}
		
		template <typename InVector, typename OutVector>
		OutVector& apply(OutVector& y, const InVector& x) const
			//OutVector& apply(OutVector& y, const InVector& x)
		{
			std::vector<typename Field::Element> z(y.size());
			VectorDomain<Field> VD(field());
			//std::vector<typename Field::Element> x_1( x.begin(), x.begin() + _L->coldim() );
			//std::vector<typename Field::Element> x_2( x.begin() + _L->coldim(), x.end() );
			Subvector<typename InVector::const_iterator> x_1(x.begin(), x.begin()+_L->coldim());
			Subvector<typename InVector::const_iterator> x_2(x.begin()+_L->coldim(), x.end());
			//std::cout << " side-by-side apply size of x: " << x.size() << " " << " size of y: " << y.size() << endl;
			//std::cout << " side-by-side apply size of x_1: " << x_1.size() << " " << " size of x_2: " << x_2.size() << endl;
			_L->apply (y, x_1);
			_R->apply (z, x_2);
			VD.addin(y, z);
			
			return y;
		}
		
		template <typename InVector, typename OutVector>
		OutVector& applyTranspose(OutVector& y, const InVector& x) const
			//OutVector& applyTranspose(OutVector& y, const InVector& x)
		{
			//std::vector<typename Field::Element> y_1( y.begin(), y.begin() + _L->coldim() );
			//std::vector<typename Field::Element> y_2( y.begin() + _L->coldim(), y.end() );
			Subvector<typename OutVector::iterator, typename OutVector::const_iterator> y_1(y.begin(), y.begin()+_L->coldim());
			Subvector<typename OutVector::iterator, typename OutVector::const_iterator> y_2(y.begin()+_L->coldim(), y.end());
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
		typedef ZOQuad<Field> Quad;
	public://temp
		const Quad *_U, *_D;
    public:
		OverUnder(const Quad* A, const Quad* B) : _U(A), _D(B) {}
		//~OverUnder() {delete _U; delete _D;}

		template <typename InVector, typename OutVector>
		OutVector& apply(OutVector& y, const InVector& x) const
			//OutVector& apply(OutVector& y, const InVector& x)
		{
			//std::vector<typename Field::Element> y_1( y.begin(), y.begin() + _U->rowdim() );
			//std::vector<typename Field::Element> y_2( y.begin() + _U->rowdim(), y.end() );
		  	Subvector<typename OutVector::iterator, typename OutVector::const_iterator> y_1(y.begin(), y.begin()+_U->rowdim());
			Subvector<typename OutVector::iterator, typename OutVector::const_iterator> y_2(y.begin()+_U->rowdim(), y.end());
			//if ((_A_ptr == 0) || (_B_ptr == 0)) { throw error }
			//std::cout << " over-under apply size of x: " << x.size() << " " << " size of y: " << y.size() << endl;
			//std::cout << " over-under apply size of y_1: " << y_1.size() << " " << " size of y_2: " << y_2.size() << endl;
			_U->apply (y_1, x);
			_D->apply (y_2, x);
			//copy(y_1.begin(), y_1.end(), y.begin());
			//copy(y_2.begin(), y_2.end(), y.begin() + y_1.size());
			
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
	
}; //namespace LinBox /*}}}*/

#endif // __LINBOX_quad_matrix_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen
