/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
/* linbox/blackbox/blas-blackbox.h
 * Copyright (C) 2004 Pascal Giorgi 
 *
 * Written by :
 *               Pascal Giorgi  pascal.giorgi@ens-lyon.fr
 *               
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */


#ifndef __BLAS_BLACKBOX_H
#define __BLAS_BLACKBOX_H

#include <linbox/matrix/blas-matrix.h>
#include <linbox/fflas/fflas.h>
#include <linbox/matrix/matrix-domain.h>
#include <linbox/field/hom.h>


namespace LinBox {

	/// \ingroup blackbox
	template <class _Field>
	class BlasBlackbox : public BlasMatrix<typename _Field::Element> 
	{

	public:

		typedef _Field Field;
		typedef typename Field::Element Element;
                typedef BlasBlackbox<_Field> Self_t;

		BlasBlackbox (const Field& F) :  _F(F), _MD(F) { _F.init(_One,1UL), _F.init(_Zero,0UL);}

		BlasBlackbox (const Field& F, size_t m, size_t n) 
			: BlasMatrix<Element> (m,n),  _F(F), _MD(F), _row(m) , _col(n) { _F.init(_One,1UL), _F.init(_Zero,0UL);}

		BlasBlackbox (const Field& F, BlasMatrix<Element>& M) 
			: BlasMatrix<Element> (M),  _F(F), _MD(F) , _row(M.rowdim()), _col(M.coldim()) { _F.init(_One,1UL), _F.init(_Zero,0UL);}

 		template< class Blackbox >
 		BlasBlackbox (const Blackbox& M)
 			: BlasMatrix<Element> (M), _F(M.field()), _MD(M.field()), _row(M.rowdim()), _col(M.coldim()) {_F.init( _One, 1UL ); _F.init( _Zero, 0UL );}

		BlasBlackbox (const BlasBlackbox<Field>& M)
			: BlasMatrix< Element> (M), _F(M._F), _MD(M._F) , _row(M._row), _col(M._col), _One(M._One), _Zero(M._Zero) {}
		


		template <class Vector1, class Vector2> 
		Vector1&  apply (Vector1& y, const Vector2& x) const 
		{
			_MD. vectorMul (y, *this, x);
			return y;
		}

		template <class Vector1, class Vector2>
		Vector1&  applyTranspose (Vector1& y, const Vector2& x) const  
		{
			_MD.vectorMul (y, TransposeMatrix<BlasBlackbox<Field> > (*this), x);
      
			return y;
		}  

            
		template<typename _Tp1>
		struct rebind
		{
                    typedef BlasBlackbox<_Tp1> other; 
                    
                    void operator() (other *& Ap, const Self_t& A, const _Tp1& F) {
                        Ap = new other(F, A.rowdim(), A.coldim());
                        typename Self_t::ConstRawIterator A_p;
                        typename other::RawIterator Ap_p;
                        Hom<Field, _Tp1> hom(A. field(), F);
                        for (A_p = A. rawBegin(), Ap_p = Ap -> rawBegin();
                             A_p != A. rawEnd(); ++ A_p, ++ Ap_p) 
                            hom.image (*Ap_p, *A_p);
                    }
                };

		size_t rowdim() const {return _row;}

		size_t coldim() const {return _col;}


		const Field &field() const  {return _F;}
		
		/*
		  std::vector<Element>& apply(std::vector<Element>& y, const std::vector<Element>& x) const {
   
		  FFLAS::fgemv( _F, FFLAS::FflasNoTrans, 
		  this->_row, this->_col,
		  this->_One,
		  _ptr, _stride,
		  &x[0],1,
		  this->_Zero,
		  &y[0],1);  
		  return y;
		  }
		*/
	protected:
		
		const Field                 & _F;  
		const MatrixDomain<Field>   & _MD; 
		size_t                _row,_col;
		Element              _One,_Zero;


	}; // end of class BlasBlackbox

	template <class Field>
	struct MatrixTraits< BlasBlackbox<Field> >
	{
		typedef BlasBlackbox<Field> MatrixType;
		typedef typename MatrixCategories::RowColMatrixTag MatrixCategory;
	};

	template <class Field>
	struct MatrixTraits< const BlasBlackbox<Field> >
	{
		typedef const BlasBlackbox<Field> MatrixType;
		typedef typename MatrixCategories::RowColMatrixTag MatrixCategory;
	};
    
} // end of namespace LinBox

#endif
