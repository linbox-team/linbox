/* linbox/blackbox/block-hankel.h
 * Copyright (C) 2005 Pascal Giorgi
 *
 * Written by Pascal Giorgi pgiorgi@uwaterlo.ca
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

#ifndef __LINBOX_bb_block_hankel_H
#define __LINBOX_bb_block_hankel_H

#include <vector>
#include "linbox/matrix/dense-matrix.h"
#include "linbox/vector/vector-domain.h"
#include "linbox/matrix/matrix-domain.h"
#include "linbox/util/debug.h"

//#define BHANKEL_TIMER


namespace LinBox
{

	class BlockHankelTag {
	public:
		typedef enum{low,up,plain} shape;
	};


	// note that P is done as a mirror
	// here we compute a^d.P(1/a^d) where d is the degree of P
	template<class Field>
	void MatPolyHornerEval (const Field                             &F,
				BlasMatrix<Field>                       &R,
				const std::vector<BlasMatrix<Field> >   &P,
				const typename Field::Element           &a)
	{
		R= P[0];
		typename BlasMatrix<Field>::Iterator it_R;
		typename BlasMatrix<Field>::ConstIterator it_P;
		for (size_t k=1; k<P.size(); ++k){
			it_R = R.Begin();
			it_P = P[k].Begin();
			for (;it_R != R.End(); ++it_R, ++it_P){
				F.mulin(*it_R, a);
				F.addin(*it_R, *it_P);
			}
		}
	}


	// P is done normally
	template<class Field>
	void VectHornelEval (const Field                                    &F,
			     std::vector<typename Field::Element>           &E,
			     const std::vector<typename Field::Element>     &P,
			     size_t                                      block,
			     const typename Field::Element                  &a)
	{
		linbox_check((P.size()% block) == 0);
		E.resize(block);
		size_t numblock = P.size()/block;
		size_t idx;
		for (size_t i=0; i<block;++i)
			F.assign(E[i], P[P.size()-block+i]);

		for (size_t i=numblock-1; i--;)
			for (size_t j= (size_t)block*i; j< (size_t) block*(i+1); ++j){
				idx= j%block;
				F.mulin(E[idx], a);
				F.addin(E[idx], P[j]);
			}
	}

	template<class Field>
	void BlockHankelEvaluation(const Field                            &F,
				   std::vector<BlasMatrix<Field> >        &E,
				   const std::vector<BlasMatrix<Field> >  &P,
				   size_t                                 k)
	{
		// do the evaluation of the Block Hankel Matrix using Horner Rules
		// at k differents points (0,1,2,..,k-1)

		E.resize(k);
		E[0]=P.back();
		typename Field::Element a ;
		F.assign(a,F.one);
		for (size_t i=1;i<k;++i){
			MatPolyHornerEval(F, E[i], P, a);
			F.addin(a, F.one);
		}
	}

	template<class Field>
	void BHVectorEvaluation(const Field                                          &F,
				std::vector<std::vector<typename Field::Element> >   &E,
				const std::vector<typename Field::Element>           &P,
				size_t                                            block)
	{
		size_t k=E.size();
		E[0]=std::vector<typename Field::Element> (P.begin(), P.begin()+block);
		typename Field::Element a;
		F.assign(a,F.one);

		for (size_t i=1;i<k;++i){
			VectHornelEval (F, E[i], P, block, a);
			F.addin(a, F.one);
		}
	}


	// compute the lagrange polynomial with k points (0,1,2,..k-1)
	template <class Field>
	void BHVectorLagrangeCoeff(const Field                            &F,
				   std::vector<std::vector<typename Field::Element> >  &P,
				   size_t                                  k)
	{
		typename Field::Element  a;
		F.assign(a,F.zero);

		// compute L:= (x)(x-1)(x-2)(x-3)...(x-k+1) = a1x+a2x^2+...+a(k-1)x^(k-1)
		std::list<typename Field::Element> L(2);
		F.assign(L.front(), F.zero);
		F.assign(*(++L.begin()), F.one);
		for (size_t i=1;i<k;++i){
			F.subin(a, F.one);
			L.push_front(F.zero);
			typename std::list<typename Field::Element>::iterator it_next = L.begin();++it_next;
			typename std::list<typename Field::Element>::iterator it = L.begin();
			for (;it_next != L.end();++it, ++it_next)
				F.axpyin(*it, a, *it_next);
		}

		// compute P[i]:= L/(x-i)
		P[0]= std::vector<typename Field::Element>(++L.begin(), L.end());
		size_t deg=L.size();
		F.assign(a,F.zero);
		for (size_t i=1;i<k;++i){
			F.addin(a,F.one);
			P[i].resize(deg-1);
			typename std::list<typename Field::Element>::const_reverse_iterator rit=L.rbegin();
			F.assign(P[i][deg-2],*rit);
			++rit;
			for (int j= (int) deg-3; j>=0;--j, ++rit)
				F.axpy(P[i][j], a, P[i][j+1], *rit);
		}

		// compute P[i]= P[i] / Prod((i-j), j<>i)
		typename Field::Element prod, ui, uj, tmp;
		F.assign(ui,F.mOne);
		for (size_t i=0;i<k;++i){
			F.assign(prod,F.one);
			F.addin(ui,F.one);
			F.assign(uj,F.zero);
			for (size_t j=0;j<k;++j){
				if (j != i){
					F.sub(tmp,ui,uj);
					F.mulin(prod,tmp);
				}
				F.addin(uj,F.one);
			}
			F.invin(prod);
			//std::cout<<"coeff: ";F.write(std::cout, prod)<<"\n";
			for (size_t l=0;l<P[i].size();++l)
				F.mulin(P[i][l], prod);
		}
	}

	template<class Field>
	void BHVectorInterpolation(const Field                                                  &F,
				   std::vector<typename Field::Element>                         &x,
				   const std::vector<std::vector<typename Field::Element> >     &E,
				   const std::vector<std::vector<typename Field::Element> >     &P,
				   size_t                                                    shift)
	{
		size_t block = E[0].size();
		size_t numblock= x.size()/block;
		std::vector<typename Field::Element> acc(block);
		VectorDomain<Field> VD(F);

		for (size_t i=shift;i<numblock+shift;++i){
			VD.mul(acc, E[0], P[0][i]);//F.write(std::cout,P[0][i])<<"*",VD.write(std::cout, E[0]);
			for (size_t j=1;j<E.size();++j){
				//F.write(std::cout,P[j][i])<<"*",VD.write(std::cout, E[j]);
				VD.axpyin(acc, P[j][i], E[j]);
			}
			//VD.write(std::cout,acc)<<"\n";;
			for (size_t j=0;j<block;++j)
				F.assign(x[x.size()-((i-shift+1)*block) +j], acc[j]);
		}
	}


	template <class _Field>
	class BlockHankel {

	public:
		typedef _Field Field;
		typedef typename Field::Element Element;

		//! is this used ?
		BlockHankel() {}

		// Constructor from a stl vector of BlasMatrix reprenting
		// all different elements in the Hankel representation
		// order of element will depend on first column and/or  last row
		// (plain->[column|row];  up -> [column]; low -> [row];)
		BlockHankel (Field &F, const std::vector<BlasMatrix<Field> > &H, BlockHankelTag::shape s= BlockHankelTag::plain) :
			_field(&F), _BMD(F)
		{
			linbox_check( H.begin()->rowdim() != H.begin()->coldim());


			switch (s) {
			case BlockHankelTag::plain :
				{
					linbox_check(H.size()&0x1);
					_deg= H.size();
					_block = H.begin()->coldim();
					_rowblock= ((_deg+1)>>1);
					_colblock= _rowblock;
					_row = _rowblock*_block;
					_col = _row;
					_shape = s;
					BlockHankelEvaluation( field(), _matpoly, H, _deg+_colblock-1);
				}
				break;
			case BlockHankelTag::up :
				{
					_deg   = H.size();
					_block = H.begin()->coldim();
					_rowblock= _deg;
					_colblock= _rowblock;
					_row   = _rowblock*_block;
					_col   = _row;
					_shape = s;
					BlockHankelEvaluation( field(), _matpoly, H, _deg+_colblock-1);
				}
				break;
			case BlockHankelTag::low :
				{
					_deg   = H.size();
					_block = H.begin()->coldim();
					_rowblock= _deg;
					_colblock= _rowblock;
					_row   = _rowblock*_block;
					_col   = _row;
					_shape = s;
					BlockHankelEvaluation( field(), _matpoly, H, _deg+_colblock-1);
				}
				break;
			}
			_numpoints = _deg+_colblock-1;
			integer prime;
			F.characteristic(prime);
			if (integer(_numpoints) > prime){
				std::cout<<"LinBox ERROR: prime ("<<prime<<") is too small for number of block ("<< _numpoints <<") in block Hankel blackbox\n";
				throw LinboxError("LinBox ERROR: prime too small in block Hankel blackbox\n");
			}

			_vecpoly.resize(_numpoints, std::vector<Element>(_block));
			_veclagrange.resize(_numpoints);
			BHVectorLagrangeCoeff(field(), _veclagrange, _numpoints);


			_vander     = BlasMatrix<Field> (_field,_numpoints,_numpoints);
			_inv_vander = BlasMatrix<Field> (_field,_numpoints,_numpoints);

			std::vector<Element> points(_numpoints);
			for (size_t i=0;i<_numpoints;++i){
				F.init(points[i],i);
				_vander.setEntry(i,0, F.one);
			}


			for (size_t j=1;j<_numpoints; ++j){
				for (size_t i=0;i<_numpoints;++i){
					F.mul(_vander.refEntry(i,j), _vander.refEntry(i,j-1), points[i]);
				}
			}


			_BMD.inv(_inv_vander, _vander);

			//! @warning memory wasted
			_partial_vander= BlasMatrix<Field> (_vander, 0, 0, _numpoints, _colblock);
			size_t shift=_colblock-1;
			if ( _shape == BlockHankelTag::up)
				shift=0;

			_partial_inv_vander= BlasMatrix<Field> (_inv_vander, shift, 0, _colblock, _numpoints);

			_x = BlasMatrix<Field> (_field,_numpoints, _block);
			_y = BlasMatrix<Field> (_field,_colblock, _block);


			_Tapply.clear();
			_Teval.clear();
			_Tinterp.clear();
		}

		// Copy construtor
		BlockHankel (const BlockHankel<Field> &H) :
			_field(H._field()), _matpoly (H._matpoly), _deg(H._deg),
			_row(H._row), _col(H._col), _rowblock(H._rowblock), _colblock(H._colblock), _block(H._block), _shape(H._shape)
			// dummy defaults
			,_vander(BlasMatrix<Field>(*_field))
			,_partial_vander(BlasMatrix<Field>(*_field))
			,_inv_vander(BlasMatrix<Field>(*_field))
			,_partial_inv_vander(BlasMatrix<Field>(*_field))
			,_y(BlasMatrix<Field>(*_field))
			,_x(BlasMatrix<Field>(*_field))
			,_numpoints(0)
		{}

		// get the column dimension
		size_t coldim() const {return _col;}

		// get the row dimension
		size_t rowdim() const {return _row;}

		const Field& field() const { return *_field;}

		// get the block dimension
		size_t blockdim() const {return _block;}


		// apply the blackbox to a vector
		template<class Vector1, class Vector2>
		Vector1& apply(Vector1 &x, const Vector2 &y) const
		{
			linbox_check(this->coldim() == y.size());
			linbox_check(this->rowdim() == x.size());
			BlasMatrixDomain<Field> BMD(field());
#ifdef BHANKEL_TIMER
			_chrono.clear();
			_chrono.start();
#endif
			// evaluation of the vector seen as a vector polynomial in
			//BHVectorEvaluation(field(), _vecpoly, y, _block);

			for (size_t i=0;i<_colblock;++i)
				for (size_t j=0;j<_block;++j)
					_y.setEntry(i,j, y[i*_block+j]);

			_BMD.mul(_x, _partial_vander, _y);

			for (size_t i=0;i<_numpoints;++i){
				for (size_t j=0;j<_block; ++j)
					field().assign(_vecpoly[i][j], _x.getEntry(i,j));
			}
#ifdef BHANKEL_TIMER
			_chrono.stop();
			_Teval+=_chrono;
			_chrono.clear();
			_chrono.start();
#endif
			std::vector<std::vector<Element> > x_vecpoly(_vecpoly.size(), std::vector<Element>(_block));
			// perform the apply componentwise
			for (size_t i=0;i<_vecpoly.size();++i)
				BMD.mul(x_vecpoly[i], _matpoly[i], _vecpoly[i]);

#ifdef BHANKEL_TIMER
			_chrono.stop();
			_Tapply+=_chrono;
			_chrono.clear();
			_chrono.start();
#endif
#if 0
			// get the result according to the right part of the polynomial
			size_t shift=_colblock-1;
			if ( _shape == BlockHankelTag::up)
				shift=0;

			// interpolation to get the result vector
			BHVectorInterpolation(field(), x, x_vecpoly, _veclagrange, shift);
#endif
			for (size_t i=0;i<_numpoints;++i)
				for (size_t j=0;j<_block;++j)
					_x.setEntry(i,j, x_vecpoly[i][j]);

			_BMD.mul(_y, _partial_inv_vander, _x);

			for (size_t i=0;i<_colblock;++i)
				for (size_t j=0;j<_block;++j)
					field().assign( x[x.size() - (i+1)*_block +j], _y.getEntry(i,j));

#ifdef BHANKEL_TIMER
			_chrono.stop();
			_Tinterp +=_chrono;
#endif

			return x;
		}


		~BlockHankel() {}

		// apply the transposed of the blackbox to a vector
		template<class Vector1, class Vector2>
		Vector1& applyTranspose(Vector1 &x, const Vector2 &y) const
		{
			return apply(x,y);
		}

	private:
		const Field  *_field;
		std::vector<BlasMatrix<Field> >                _matpoly;
		mutable std::vector<std::vector<Element> >       _vecpoly;
		std::vector<std::vector<Element> >           _veclagrange;
		BlasMatrix<Field>                               _vander;
		BlasMatrix<Field>                       _partial_vander;
		BlasMatrix<Field>                           _inv_vander;
		BlasMatrix<Field>                   _partial_inv_vander;
		mutable BlasMatrix<Field>                            _y;
		mutable BlasMatrix<Field>                            _x;
		BlasMatrixDomain<Field>                              _BMD;
		size_t _deg;
		size_t _row;
		size_t _col;
		size_t _rowblock;
		size_t _colblock;
		size_t _block;
		size_t _numpoints;
		BlockHankelTag::shape _shape;
		mutable Timer _Tapply, _Teval, _Tinterp, _chrono;
	};

} // end of namespace LinBox

#endif //__LINBOX_bb_block_hankel_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
