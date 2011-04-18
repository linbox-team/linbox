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


#ifndef __LINBOX_blas_blackbox_H
#define __LINBOX_blas_blackbox_H

#include <linbox/matrix/blas-matrix.h>
#include <linbox/fflas/fflas.h>
#include <linbox/matrix/matrix-domain.h>
#include <linbox/field/hom.h>
#include <linbox/field/multimod-field.h>
#include <linbox/util/matrix-stream.h>

namespace LinBox 
{
	
	template <class Field>
	bool checkBlasApply(const Field &F, size_t n) {

		integer chara, card;
		F.characteristic(chara);
		F.cardinality(card);
		
		if ((chara != card) || chara == 0)
			return false;
		else
			if (n*chara*chara < integer("9007199254740992"))
				return true;
			else
				return false;
	}


	/** \brief dense matrix representation for BLAS based elimination. \ingroup blackbox
	 *
	 *  A BlasBlackbox can be constructed from any blackbox matrix. This costs n blackbox
	 * matrix vector products in general, but is efficiently done from a DenseMatrix 
	 * or SparseMatrix.
	 */
	template <class _Field>
	class BlasBlackbox : public BlasMatrix<typename _Field::Element> 
	{

	public:

		typedef _Field Field;
		typedef typename Field::Element Element;
		typedef BlasBlackbox<_Field> Self_t;

		//BlasBlackbox () {} // BB: ??

		BlasBlackbox (const Field& F) :  _F(F), _MD(F), _VD(F) 
		{ _F.init(_One,1UL), _F.init(_Zero,0UL);_use_fflas=false;}

		BlasBlackbox (const Field& F, const size_t m, const size_t n) 
			: BlasMatrix<Element> (m,n),  _F(F), _MD(F), _VD(F), _row(m) , _col(n) 
		{ 
			_F.init(_One,1UL), _F.init(_Zero,0UL);
			typename BlasMatrix<Element>::RawIterator it = this->rawBegin();
			for (; it != this->rawEnd(); ++it)
				_F.init(*it, 0);
			_use_fflas= checkBlasApply(_F, _col);
		}

		BlasBlackbox(MatrixStream<Field> &ms)
			: BlasMatrix<Element> (ms), _F(ms.getField()), _MD(ms.getField()), _VD(ms.getField()) 
		{
			ms.getRows(_row);
			ms.getColumns(_col);
			_use_fflas= checkBlasApply(_F, _col);
		}
		


		BlasBlackbox (const Field& F, BlasMatrix<Element>& M) 
			: BlasMatrix<Element> (M),  _F(F), _MD(F) , _VD(F),  _row(M.rowdim()), _col(M.coldim()) 
		{ _F.init(_One,1UL), _F.init(_Zero,0UL); _use_fflas= checkBlasApply(_F, _col); }


 		template< class Blackbox >
 		BlasBlackbox (const Blackbox& M)
 			: BlasMatrix<Element> (M), _F(M.field()), _MD(M.field()), _VD(M.field()), _row(M.rowdim()), _col(M.coldim()) 
		{_F.init( _One, 1UL ); _F.init( _Zero, 0UL ); _use_fflas= checkBlasApply(_F, _col);}
		

		BlasBlackbox (const BlasBlackbox<Field>& M)
			: BlasMatrix< Element> (M), _F(M._F), _MD(M._F), _VD(M._F), 
			  _row(M._row), _col(M._col), _One(M._One), _Zero(M._Zero) {_use_fflas= checkBlasApply(_F, _col);}

		BlasBlackbox (const BlasBlackbox<Field>& M, const size_t i0, const size_t j0, const size_t m, const size_t n)
			: BlasMatrix< Element> (M,i0,j0,m,n), _F(M._F), _MD(M._F), _VD(M._F), 
			  _row(m), _col(n), _One(M._One), _Zero(M._Zero) {_use_fflas= checkBlasApply(_F, _col);}
		

		BlasBlackbox (const Field &F, const BlasBlackbox<Field>& M)
			: BlasMatrix< Element> (M), _F(M._F), _MD(M._F), _VD(F),
			  _row(M._row), _col(M._col), _One(M._One), _Zero(M._Zero) {_use_fflas= checkBlasApply(_F, _col);}

		template <class Vector1, class Vector2> 
		Vector1&  apply (Vector1& y, const Vector2& x) const 
		{
			
			if (_use_fflas){
       
				FFLAS::fgemv( _F, FFLAS::FflasNoTrans, 
					      this->_row, this->_col,
					      this->_One,
					      this->_ptr, this->_stride,
					      &x[0],1,
					      this->_Zero,
					      &y[0],1);  
			}			
			else {
				_MD. vectorMul (y, *this, x);

				//typename BlasMatrix<Element>::ConstRowIterator i = this->rowBegin ();
				//typename Vector1::iterator j = y.begin ();
				
				//for (; j != y.end (); ++j, ++i)
				//	_VD.dot (*j, *i, x);
			}
			return y;
		}

		template <class Vector1, class Vector2>
		Vector1&  applyTranspose (Vector1& y, const Vector2& x) const  
		{

			if (_use_fflas)
				FFLAS::fgemv( this->_F, FFLAS::FflasTrans, 
					      this->_row, this->_col,
					      this->_One,
					      this->_ptr, this->_stride,
					      &x[0],1,
					      this->_Zero,
					      &y[0],1);  
			
			else {
				typename BlasMatrix<Element>::ConstColIterator i = this->colBegin ();
				typename Vector1::iterator j = y.begin ();
				
				for (; j != y.end (); ++j, ++i)
					_VD.dot (*j, x, *i);
			}
      
			return y;
		}  

            
            template<typename _Tp1>
            struct rebind
            {
                typedef BlasBlackbox<_Tp1> other; 
                
                void operator() (other & Ap, const Self_t& A, const _Tp1& F) {
                    typedef typename BlasMatrix<Element>::ConstRawIterator ConstRawIterator ;
                    ConstRawIterator A_p;
                    typename other::RawIterator Ap_p;
                    Hom<Field, _Tp1> hom(A. field(), F);
                    for (A_p = A. rawBegin(), Ap_p = Ap.rawBegin();
                         A_p != A. rawEnd(); ++ A_p, ++ Ap_p) 
                        hom.image (*Ap_p, *A_p);
                }
            };
            
            
            template<typename _Tp1>
            BlasBlackbox(const BlasBlackbox<_Tp1>& M, const Field& F) 
                    : BlasMatrix<Element>(M.rowdim(),M.coldim()), 
                      _F(F),_MD(F),_VD(F), 
                      _row(M.rowdim()), _col(M.coldim()),
                      _One(F.one), _Zero(F.zero) {
                _use_fflas = checkBlasApply(F, M.coldim());
                typename BlasBlackbox<_Tp1>::template rebind<Field>() (*this, M, F);
            }
            
                      

            size_t rowdim() const {return _row;}

		size_t coldim() const {return _col;}


		const Field &field() const  {return _F;}
		Field &field() {return const_cast<Field&>(_F);}
		


		/** Read the blackbox from an input stream
		 * @param file Input stream from which to read
		 */
		std::istream &read (std::istream &file){
			return BlasMatrix<Element>::read(file, _F);

		}

		/** Write the blackbox to an output stream
		 * @param os Output stream to which to write		 
		 */
		std::ostream &write (std::ostream &os) const {
			return DenseSubmatrix<Element>::write(os, _F);
		}

		

	protected:
		
		const Field                 & _F;  
		MatrixDomain<Field>    _MD; 
		VectorDomain<Field>    _VD;  
		size_t                _row,_col;
		Element              _One,_Zero;
		bool                  _use_fflas;


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

	template <class Field>
	class MatrixContainerTrait<BlasBlackbox<Field> > {
	public:
		typedef MatrixContainerCategory::BlasContainer Type;
	};

	template <class Field>
	class MatrixContainerTrait<const BlasBlackbox<Field> > {
	public:
		typedef MatrixContainerCategory::BlasContainer Type;
	};
	
    



	template<>
	class BlasBlackbox<MultiModDouble> {


	public:

		typedef MultiModDouble         Field;
		typedef std::vector<double>  Element;
                typedef BlasBlackbox<MultiModDouble> Self_t;

		//BlasBlackbox () {}
		
		BlasBlackbox (const MultiModDouble& F) :  _F(F) , _rep(F.size()), _entry(F.size())
		{}

		BlasBlackbox (const Field& F, size_t m, size_t n, bool alloc=true) 
			:  _F(F), _row(m) , _col(n) , _rep(F.size()),  _entry(F.size())
		{ 			
			for (size_t i=0;i<_rep.size();++i)
				if (alloc)_rep[i]       =  new BlasBlackbox<Modular<double> > (F.getBase(i), m, n);
			}
		
		BlasBlackbox (const BlasBlackbox<MultiModDouble> & A): _F(A._F),_row(A._row), _col(A._col), 
								       _rep(A._rep.size()), _entry(A._entry) {

			for (size_t i=0;i<_rep.size();++i)
				_rep[i]= new  BlasBlackbox<Modular<double> > (const_cast<BlasBlackbox<Modular<double> >& >( *A._rep[i]));
		}


		const BlasBlackbox<MultiModDouble>& operator=(const BlasBlackbox<MultiModDouble> & A){			
			_F   = A._F;		
			_row = A._row;
			_col = A._col;
			_rep = std::vector<BlasBlackbox<Modular<double> >* >(A._rep.size());		
			_entry = A._entry;
			for (size_t i=0;i<_rep.size();++i)
				_rep[i]= new  BlasBlackbox<Modular<double> > (const_cast<BlasBlackbox<Modular<double> >& >( *A._rep[i]));
			return *this;
		}
  

		~BlasBlackbox() {for (size_t i=0; i< _rep.size();++i) {delete _rep[i];} }

		template <class Vector1, class Vector2> 
		Vector1&  apply (Vector1& y, const Vector2& x) const 
		{	
			for (size_t i=0;i<_rep.size();++i) {				
				std::vector<double> x_tmp(x.size()), y_tmp(y.size());
				for (size_t j=0;j<x.size();++j)
					x_tmp[j]= x[j][i];		
								
				_rep[i]->apply(y_tmp, x_tmp);
			
				for (size_t j=0;j<y.size();++j){
					y[j][i]=y_tmp[j];				

				}
			}
			
			return y;
		}

		template <class Vector1, class Vector2>
		Vector1&  applyTranspose (Vector1& y, const Vector2& x) const  
		{
			for (size_t i=0;i<_rep.size();++i) {
				std::vector<double> x_tmp(x.size()), y_tmp(y.size());
				for (size_t j=0;j<x.size();++j)
					x_tmp[i]= x[j][i];

				_rep[i]->applyTranspose(y_tmp, x_tmp);
				
				for (size_t j=0;j<y.size();++j)
					y[j][i]=y_tmp[i];
			}
			
			return y;
		}  

            /*
		template<typename _Tp1>
		struct rebind
		{
                    typedef BlasBlackbox<_Tp1> other; 
                    
                    void operator() (other *& Ap, const Self_t& A, const _Tp1& F) {
                        Ap = new other(F, A.rowdim(), A.coldim());
			Hom<Field, _Tp1> hom(A. field(), F);
			
			hom.image (*Ap_p, *A_p);
                    }
                };
	    */

		size_t rowdim() const {return _row;}

		size_t coldim() const {return _col;}


		const Field &field() const  {return _F;}
		
		
		std::ostream& write(std::ostream& os) const {
			for (size_t i=0;i<_rep.size();++i)
				_rep[i]->write(os);
			return os;
		}


		void setEntry (size_t , size_t j, const Element &a_ij){
			for (size_t i=0; i< _rep.size();++i)
				_rep[i]->setEntry(i,j,a_ij[i]);
		}
		
		
		const Element& getEntry (size_t , size_t j){		
			for (size_t i=0; i< _rep.size();++i)
				_entry[i]=_rep[i]->getEntry(i,j);
			return _entry;
		}
	
		BlasBlackbox<Modular<double> >*& getMatrix(size_t i) {return _rep[i];}

	protected:
		
		MultiModDouble                 _F;  
		const std::vector<MatrixDomain<Modular<double> > >   _MD; 
		size_t                  _row,_col;
		Element                _One,_Zero;		
		std::vector<BlasBlackbox<Modular<double> >* > _rep;
		std::vector<double>       _entry;     
	};

		
		
	

	template <>
	class MatrixContainerTrait<BlasBlackbox<MultiModDouble> > {
	public:
		typedef MatrixContainerCategory::Blackbox Type;
	};



} // end of namespace LinBox

#endif //__LINBOX_blas_blackbox_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
