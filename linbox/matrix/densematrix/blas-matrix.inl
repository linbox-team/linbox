/*
 * Copyright (C) 2004 Pascal Giorgi, Clément Pernet
 *               2013, 2014 the LinBox group
 *               2018 revamped by Pascal Giorgi 
 *
 * Written by :
 *               Pascal Giorgi  <pascal.giorgi@lirmm.fr>
 *               Clément Pernet <clement.pernet@imag.fr>
 *               Brice Boyer (briceboyer) <boyer.brice@gmail.com>
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

/*!@internal
 * @file matrix/densematrix/blas-matrix.inl
 * @ingroup densematrix
 * A \c BlasMatrix<\c _Field > represents a matrix as an array of
 * <code>_Field</code>s.
 */

#ifndef __LINBOX_densematrix_blas_matrix_INL
#define __LINBOX_densematrix_blas_matrix_INL


#include "linbox/field/hom.h"
#include "linbox/vector/vector-domain.h"

#include "fflas-ffpack/fflas/fflas.h"

///////////////////
//   PROTECTED   //
///////////////////

namespace LinBox
{

	template<class _Field, class _Storage>
    template <class constIterator>
	void BlasMatrix< _Field, _Storage >::createBlasMatrix (constIterator v)
	{
		constIterator v_end = v+(_col*_row) ;
		Element_ptr iter_addr = getPointer();
		for (; v != v_end ; ++v, ++iter_addr)
            {
                field().assign(*iter_addr,*v);
            }
	}

	template<class _Field, class _Storage>
	template <class _Matrix>
	void BlasMatrix< _Field, _Storage >::createBlasMatrix (const _Matrix& A,
                                                           const size_t i0,const size_t j0,
                                                           const size_t m, const size_t n,
                                                           MatrixContainerCategory::Container)
	{
        linbox_check( areFieldEqual(A.field(),field() ) );

		typename _Matrix::ConstIterator         iter_value = A.Begin();
		typename _Matrix::ConstIndexedIterator  iter_index = A.IndexedBegin();

        // PG -> BUG if we use iter_value !=A.End() 
		for (;iter_index != A.IndexedEnd(); ++iter_value,++iter_index){
			int64_t i,j;
			i=(int64_t)iter_index.rowIndex()-(int64_t)i0;
			j=(int64_t)iter_index.colIndex()-(int64_t)j0;
			if ( (i>=0) && (j>=0) && (i< (int64_t)m) && (j < (int64_t)n))
				setEntry((size_t)i, (size_t)j, *iter_value);
		}
	}

    template<class _Field, class _Storage>
	template <class _Matrix>
	void BlasMatrix< _Field, _Storage >::createBlasMatrix (const _Matrix& A,
                                                           const size_t i0,const size_t j0,
                                                           const size_t m, const size_t n,
                                                           MatrixContainerCategory::BlasContainer)
	{
        linbox_check( areFieldEqual(A.field(),field() ) );
        FFLAS::fassign(field(), m, n, A.getPointer(), A.getStride(), getPointer(), getStride());
	}

	template<class _Field, class _Storage>
	template <class Matrix>
	void BlasMatrix< _Field, _Storage >::createBlasMatrix (const Matrix& A,
                                                           const size_t i0,const size_t j0,
                                                           const size_t m, const size_t n,
                                                           MatrixContainerCategory::Blackbox)
	{
		linbox_check( areFieldEqual(A.field(),field() ) );

        BlasVector<_Field, _Storage> e(A.field(),A.coldim(), field().zero), tmp(A.field(),A.rowdim());
		ColIterator col_p;
        
		typename BlasMatrix<_Field, _Storage >::Col::iterator elt_p;
		typename BlasVector<_Field, _Storage>::iterator e_p, tmp_p;


		for (col_p = colBegin(), e_p = e.begin();
		     e_p != e.end(); ++ col_p, ++ e_p)
            {

                field().assign(*e_p, field().one);

                A.apply (tmp, e);

                for (tmp_p = tmp.begin(), elt_p = col_p -> begin();
                     tmp_p != tmp.end(); ++ tmp_p, ++ elt_p)

                    field().assign(*elt_p, *tmp_p);

                field().assign(*e_p, field().zero);
            }
	}
            	
    ////////////
    // MEMORY //
    ////////////

    
	template < class _Field, class _Storage >
	void BlasMatrix< _Field, _Storage >::init(const size_t & m, const size_t & n)
	{
        _row = m; _col = n;
		_rep.resize(m*n);
	}

    template < class _Field, class _Storage >
	void BlasMatrix< _Field, _Storage >::resize (const size_t & m, const size_t & n, const Element& val )
	{
#ifndef NDEBUG
		if (_col > 0 && _col != n)
			std::cerr << " ***Warning*** you are resizing a matrix, possibly loosing data. " << std::endl;
#endif
		_rep.resize (m * n, val); 
		_row = m;
		_col = n;
	}


    
    //////////////////
    // CONSTRUCTORS //
    //////////////////

    
	template < class _Field, class _Storage >
	BlasMatrix< _Field, _Storage >::BlasMatrix (const _Field &F) :
		_row(0),_col(0),_rep(F){}


	template < class _Field, class _Storage >
	BlasMatrix< _Field, _Storage >::BlasMatrix ( const _Field &F, const size_t & m, const size_t & n) :
		_row(m),_col(n),_rep(F,_row*_col){}


    template < class _Field, class _Storage >
    BlasMatrix< _Field, _Storage >::BlasMatrix(MatrixStream<_Field>& ms) :
        _row(0),_col(0),_rep(ms.getField())
    {
        if( !ms.getArray(_rep) || !ms.getDimensions(_row, _col) )
            throw ms.reportError(__FUNCTION__,__LINE__);
    }


	template < class _Field, class _Storage >
	template <class Matrix>
	BlasMatrix< _Field, _Storage >::BlasMatrix (const Matrix &A) : 
		_row(A.rowdim()),_col(A.coldim()),_rep(A.field(),_row*_col)
    {
        createBlasMatrix(A,0,0,_row,_col,typename MatrixContainerTrait<Matrix>::Type());
    }

    template < class _Field, class _Storage >
    template <class Matrix>
    BlasMatrix< _Field, _Storage >::BlasMatrix (const Matrix& A,const size_t &i0, const size_t &j0,const size_t &m,  const size_t &n) :
        _row(m),_col(n),_rep(A.field(),_row*_col)
    {
                   
        createBlasMatrix(A, i0, j0, m, n,typename MatrixContainerTrait<Matrix>::Type());
    }

    template < class _Field, class _Storage >
    template<class constIterator>
    BlasMatrix< _Field, _Storage >::BlasMatrix (const _Field &F,
                                                const size_t & m, const size_t & n,
                                                const constIterator& it) :
        _row(m), _col(n),_rep(F, _row*_col)
    {
        createBlasMatrix(it);
    }

    template < class _Field, class _Storage >
    template <class StreamVector>
    BlasMatrix< _Field, _Storage >::BlasMatrix (const Field &F, VectorStream<StreamVector> &stream) :
        _row(stream.size ()), _col(stream.dim ()), _rep(F,_row*_col)
    {
        VectorDomain<Field>    VD(F);
        StreamVector tmp(F);
        typename BlasMatrix<Field,_Storage>::RowIterator p;
        VectorWrapper::ensureDim (tmp, stream.dim ());
        for (p = BlasMatrix<Field,_Storage>::rowBegin (); p != BlasMatrix<Field,_Storage>::rowEnd (); ++p) {
            stream >> tmp;

            VD.copy (*p, tmp);
        }
    }

    template < class _Field, class _Storage >
    template<class OtherMatrix>
    BlasMatrix< _Field, _Storage >::BlasMatrix (const OtherMatrix &A, const _Field &F) :
        _row(A.rowdim()), _col(A.coldim()), _rep(F, _row*_col)
    {
        //std::cout<<"GIORGI: BlasMatrix reducing mod \n";
		typename OtherMatrix::template rebind<_Field>()(*this,A);        
    }

    template < class _Field, class _Storage >
    BlasMatrix< _Field, _Storage >::BlasMatrix (const  BlasMatrix< _Field, _Storage>  &A) :
		_row(A.rowdim()), _col(A.coldim()),_rep(A.field(),_row*_col)
	{
        createBlasMatrix(A,0,0,_row,_col,MatrixContainerCategory::BlasContainer());
	}
    
    
    template < class _Field, class _Storage >
    BlasMatrix< _Field, _Storage >& BlasMatrix< _Field, _Storage >::operator= (const BlasMatrix< _Field, _Storage >& A)
    {
        if ( &A == this)
            return *this;

        _col = A.coldim();
        _row = A.rowdim();
        _rep = Storage(A.field(),_row*_col);
        createBlasMatrix(A,0,0,_row,_col,MatrixContainerCategory::BlasContainer());
        return *this;
    }
   
    template < class _Field, class _Storage >
    template < class Matrix>
    BlasMatrix< _Field, _Storage >& BlasMatrix< _Field, _Storage >::operator= (const Matrix& A)
    {
        _col = A.coldim();
        _row = A.rowdim();
        _rep = Storage(A.field(),_row*_col);
        createBlasMatrix(A,0,0,_row,_col,typename MatrixContainerTrait<Matrix>::Type());
        return *this;
    }


    //! @bug other rep
    template < class _Field, class _Storage >
    template<typename _Tp1, typename _Rep2>
    struct BlasMatrix< _Field, _Storage >::rebind {
		typedef BlasMatrix<_Tp1, _Rep2> other;

		void operator() (other & Ap, const Self_t& A) {            
			// typedef Self_t::ConstIterator ConstSelfIterator ;
			typedef typename BlasMatrix< _Field, _Storage >::ConstIterator ConstSelfIterator ;
			typedef typename other::Iterator OtherIterator ;
			OtherIterator    Ap_i = Ap.Begin();
			ConstSelfIterator A_i = A.Begin();
			Hom<Field, _Tp1> hom(A. field(), Ap. field()) ;
			for ( ; A_i != A. End(); ++ A_i, ++ Ap_i)
				hom.image (*Ap_i, *A_i);
            //Ap.write(std::cout);
		} 
    };    

    ///////////////////
    //      I/O      //
    ///////////////////

    template < class _Field, class _Storage >
    std::istream &BlasMatrix< _Field, _Storage >::read (std::istream &file)
    {
        MatrixStream<Field> ms(field(), file);
        if( !ms.getArray(_rep) || !ms.getDimensions(_row, _col) ){
            throw ms.reportError(__FUNCTION__,__LINE__);
        }
        
        return file;       
    }       

    template < class _Field, class _Storage >
    std::ostream &BlasMatrix< _Field, _Storage >::write (std::ostream &os, Tag::FileFormat f) const
    {
        //std::cout<<"BlasMatrix write: "<<&(*_rep.getPointer())<<" "<<_rep<<std::endl;
        constSubMatrixType B(*this);
        return B.write(os,f);
    }


} // end of LinBox namespace


#endif // __LINBOX_densematrix_blas_matrix_INL

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
