/*
 * Copyright (C) 2004 Pascal Giorgi, Clément Pernet
 *               2013, 2014 the LinBox group
 *               2018 revamped by Pascal Giorgi 
 * Written by :
 *               Pascal Giorgi  pascal.giorgi@lirmm.fr
 *               Clément Pernet clement.pernet@imag.fr
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

/*! @file matrix/densematrix/new-blas-submatrix.inl
 * @ingroup densematrix
 */

#include "linbox/util/write-mm.h"
#include "linbox/field/hom.h"

namespace LinBox {

    
    //////////////////
    // CONSTRUCTORS //
    //////////////////

    /*  constructors */
    
    template < class _Matrix >
    BlasSubmatrix<_Matrix>::BlasSubmatrix (typename BlasSubmatrix<_Matrix>::matrixType &M, size_t rowbeg, size_t colbeg, size_t Rowdim, size_t Coldim) :
        _ptr(M.getPointer()+rowbeg*M.getStride()+colbeg),
        _row (Rowdim),
        _col(Coldim),
        _stride(M.getStride()),
        _field(M.field())
	{ 
        linbox_check ( rowbeg  <= M.rowdim() ); // allow for NULL matrix
		linbox_check ( colbeg  <= M.coldim() );
        linbox_check ( rowbeg+Rowdim <= M.rowdim() );
        linbox_check ( colbeg+Coldim <= M.coldim() );
    }

    template < class _Matrix >
    BlasSubmatrix<_Matrix>::BlasSubmatrix (BlasSubmatrix<_Matrix>  &M, size_t rowbeg, size_t colbeg, size_t Rowdim, size_t Coldim) :
        _ptr(M.getPointer()+rowbeg*M.getStride()+colbeg),
        _row (Rowdim),
        _col(Coldim),
        _stride(M.getStride()),
        _field(M.field())
	{
        linbox_check ( rowbeg  <= M.rowdim() ); // allow for NULL matrix
		linbox_check ( colbeg  <= M.coldim() );
        linbox_check ( rowbeg+Rowdim <= M.rowdim() );
        linbox_check ( colbeg+Coldim <= M.coldim() );
    }

    template < class _Matrix >
    BlasSubmatrix<_Matrix>::BlasSubmatrix (_Matrix &M):        
        _ptr(M.getPointer()),
        _row (M.rowdim()),
        _col (M.coldim()),
		_stride(M.getStride()),
        _field(M.field())
    {
        //std::cout<<"BlasSubMatrix @# constructor: "<<getPointer()<<" ="<<M.getPointer()<<"--->"<<&M<<std::endl;
        //std::cout<<M<<std::endl;
        //std::cout<<*this<<std::endl;

    }
    template < class _Matrix >
    BlasSubmatrix<_Matrix>::BlasSubmatrix (const typename BlasSubmatrix<_Matrix>::Field& F,
                                           typename BlasSubmatrix<_Matrix>::pointer ptr,
                                           size_t Rowdim,
                                           size_t Coldim,
                                           size_t stride):        
        _ptr(ptr),
        _row (Rowdim),
        _col(Coldim),
		_stride(stride),
        _field(F)
    {
    }

     template < class _Matrix >
     BlasSubmatrix<_Matrix>::BlasSubmatrix(const constSubMatrixType &M) :
         _ptr(M.getPointer()),
         _row (M.rowdim()),
         _col(M.coldim()),
         _stride(M.getStride()),
         _field(M.field())
     {
         //std::cout<<"BlasSubMatrix Copy constructor: "<<getPointer()<<" ="<<M.getPointer()<<"--->"<<&M<<std::endl;
         // std::cout<<M<<std::endl;
         // std::cout<<*this<<std::endl;
     }
    
    template < class _Matrix >
    BlasSubmatrix<_Matrix>::BlasSubmatrix(nonconstSubMatrixType &M) :
        _ptr(M.getPointer()),
        _row (M.rowdim()),
        _col(M.coldim()),
        _stride(M.getStride()),
        _field(M.field())
    {
        //std::cout<<"BlasSubMatrix Copy constructor 2: "<<getPointer()<<" ="<<M.getPointer()<<"--->"<<&M<<std::endl;
        // std::cout<<M<<std::endl;
        // std::cout<<*this<<std::endl;
    }


    template < class _Matrix >
    BlasSubmatrix<_Matrix>::BlasSubmatrix(const nonconstSubMatrixType &M) :
        _ptr(M.getPointer()),
        _row (M.rowdim()),
        _col(M.coldim()),
        _stride(M.getStride()),
        _field(M.field())
    {
        //std::cout<<"BlasSubMatrix Copy constructor 3: "<<getPointer()<<" ="<<M.getPointer()<<"--->"<<&M<<std::endl;
        // std::cout<<M<<std::endl;
        // std::cout<<*this<<std::endl;
    }

    
    template < class _Matrix >
    template < class _AnyMatrix >
    BlasSubmatrix<_Matrix>& BlasSubmatrix<_Matrix>::copy (const _AnyMatrix & M){
        //std::cout<<"BlasSubMatrix Copy Method\n";
        if (_row == M.rowdim() && _col == M.coldim()){
            // for (size_t i=0;i<_row;i++)
            //     for(size_t j=0;j<_col;j++)
            //         setEntry(i,j,M.getEntry(i,j));
            FFLAS::fassign(_field, _row, _col, M.getPointer(), M.getStride() , _ptr, _stride); 
        }
        else
            throw LinBoxError("Calling copy from BlasSubMatrix with matrices of different dimension ... not allowed");
        return *this;
    }
    
    
    //////////////////
    //   ELEMENTS   //
    //////////////////

    template < class _Matrix >
    void BlasSubmatrix<_Matrix>::setEntry (size_t i, size_t j, const Element &a_ij) {
        field().assign(_ptr[i*_stride+j],a_ij);
    }

    template < class _Matrix >
    typename BlasSubmatrix<_Matrix>::Element& BlasSubmatrix<_Matrix>::refEntry (size_t i, size_t j) {
        return _ptr[i*_stride+j];
    }

    template < class _Matrix >
    const typename BlasSubmatrix<_Matrix>::Element& BlasSubmatrix<_Matrix>::getEntry (size_t i, size_t j) const {
        return _ptr[i*_stride+j];
    }

    template < class _Matrix >
    typename BlasSubmatrix<_Matrix>::Element& BlasSubmatrix<_Matrix>::getEntry (Element &x, size_t i, size_t j) const {
        return field().assign(x,_ptr[i*_stride+j]);
    }




    ///////////////////
    //      I/O      //
    ///////////////////

    template < class _Matrix >
    std::istream& BlasSubmatrix<_Matrix>::read (std::istream &file)
    {        
        MatrixStream<Field> ms(field(),file);
        size_t i,j;
		Element v;
        do {
            ms.nextTriple(i,j,v);
            if ( (i<0 or i> _row) and (j<0 or j> _col)){
                throw ms.reportError(__FUNCTION__,__LINE__);
            }
            setEntry(i,j,v);
        }  while( ms.getError() <= END_OF_MATRIX);
			        
     	return file;
    }

    template < class _Matrix >
    std::ostream& BlasSubmatrix<_Matrix>::write (std::ostream &os, Tag::FileFormat f) const {
        
		ConstRowIterator p;
		switch(f) {
        case (Tag::FileFormat::MatrixMarket ) : /* Matrix Market */
            {
                writeMMArray(os, *this, "BlasSubmatrix");
            }
            break;
        case (Tag::FileFormat::Plain) : /*  raw output */
            {
                integer c;
                field().cardinality (c);

                if (c >0){
                    int wid = (int) ceil (log ((double) c) / M_LN10);
                    for (p = rowBegin (); p != rowEnd ();++p) {
                    typename ConstRow::const_iterator pe;

                    os << "  [ ";

                    for (pe = p->begin (); pe != p->end (); ++pe) {
                        os.width (wid);
                        field().write (os, *pe);
                        os << ' ';
                    }
                    os << ']' << std::endl;
                    }
                }
                else {
                    for (p = rowBegin (); p != rowEnd ();++p) {
                        os<<*p<<std::endl;
                    }
                }
            }
            break;
        case (Tag::FileFormat::Pretty) : /* pretty is close to Maple */
        case (Tag::FileFormat::Maple) : /*  maple format */
            {
                const bool isMaple(f == Tag::FileFormat::Maple ? true:false);
                if (isMaple) 
                    os << "Matrix( " << rowdim() << ',' << coldim() << ",\n[" ;
                for (p = rowBegin (); p != rowEnd (); ) {
                    typename ConstRow::const_iterator pe;
                    if ( (!isMaple) || (p!=rowBegin())) os << ' ';
                    os << "[ ";
                        
                    for (pe = p->begin (); pe != p->end (); ) {
                        field().write (os, *pe);
                        ++pe ;
                        if (pe != p->end()) {
                            if (isMaple)
                                os << ',';
                            os << ' ';
                        }
                    }

                    os << "]" ;
                    ++p ;
                    if (p != rowEnd() ) {
                        if (isMaple)
                            os << ',';
                        os << '\n';
                    }

                }
                if (isMaple)
                    os << "])" ;
            }
            break;
        case (Tag::FileFormat::Guillaume) : /*  sms format */
            {

                os << rowdim() << ' ' << coldim() << " M\n" ;
                size_t i(0);
                for (p = rowBegin (); p != rowEnd (); ++i,++p) {
                    typename ConstRow::const_iterator pe;
                    size_t j(0);
                    for (pe = p->begin (); pe != p->end (); ++j,++pe) {
                        if ( !(field().isZero(*pe))) {
                            field().write (os << (i+1) << ' ' << (j+1) << ' ', *pe) << '\n';
                        }
                    }
                }
                os << "0 0 0" ;
            }
            break;
        case (Tag::FileFormat::HTML) : /*  HTML format */
            {

                os << "<table border=\"1\">" ;
                for (p = rowBegin (); p != rowEnd (); ) {
                    typename ConstRow::const_iterator pe;

                    os << "<tr>";

                    for (pe = p->begin (); pe != p->end (); ) {
                        field().write (os<< "<td>", *pe)<<"</td>";
                        ++pe ;
                    }

                    os << "</tr>" << std::endl;
                    ++p ;

                }
                os << "</table>" ;
            }
            break;
        case (Tag::FileFormat::LaTeX) : /*  LaTex format */
            {
                os << "\\begin{pmatrix} " << std::endl;
                for (p = rowBegin (); p != rowEnd (); ) {
                    typename ConstRow::const_iterator pe;

                    for (pe = p->begin (); pe != p->end (); ) {
                        field().write (os, *pe);
                        ++pe ;
                        if (pe != p->end())
                            os << "& ";
                    }

                    os << "\\\\" << std::endl;
                    ++p ;

                }
                os << "\\end{pmatrix}" ;
            }
            break;
        default : /*  this is an error */
            {
                throw LinBoxError("unknown format to write matrix in");
            }
		}
		return os;
	}

    


    ///////////////////
    //   ITERATORS   //
    ///////////////////



	/*! Raw Iterators.
	 * @ingroup iterators
	 *
	 * The raw iterator is a method for accessing all entries in the matrix
	 * in some unspecified order. This can be used, e.g. to reduce all
	 * matrix entries modulo a prime before passing the matrix into an
	 * algorithm.
	 */
	template < class _Matrix >
	class BlasSubmatrix<_Matrix>::Iterator {
	public:
		Iterator (){}

		Iterator (const typename BlasSubmatrix<_Matrix>::pointer& cur,
                  const size_t c_dim,
                  const size_t stride,
                  const size_t c_idx) :
			_cur (cur), _c_dim (c_dim), _stride(stride), _c_idx (c_idx)	{}

		/*! @internal
		 * @brief copy operator.
		 * @param r Iterator to copy.
		 */
		Iterator& operator = (const Iterator& r) {
			_cur    = r._cur;
			_c_dim  = r._c_dim;
			_stride = r._stride;
			_c_idx  = r._c_idx;
			return *this;
		}

		Iterator& operator ++() {
			if (_c_idx < _c_dim - 1){
				++_cur; ++_c_idx;
			}
			else {
				_cur = _cur + _stride - _c_dim + 1;
				_c_idx = 0;
			}

			return *this;
		}

		Iterator operator++ (int) {
            Iterator tmp(*this);
            ++(*this);
			return tmp;
		}


		/*! @internal
		 * @brief  operator !=.
		 * @param r Iterator to test inequaltity from.
		 */
		bool operator != (const Iterator& r) const {
			return (_cur != r._cur || _c_dim != r._c_dim) || (_stride != r._stride) || (_c_idx != r._c_idx);
		}

		//! @internal operator *.
		typename _Matrix::Element& operator * () {
			return *_cur;
		}

		//! @internal operator *.
		const typename _Matrix::Element& operator * () const {
			return *_cur;
		}

	protected:
        typename BlasSubmatrix<_Matrix>::pointer  _cur;
		size_t _c_dim;
		size_t _stride;
		size_t _c_idx;
	};

	/*! Raw Iterators (const version).
	 * @ingroup iterators
	 * The raw iterator is a method for accessing all entries in the matrix
	 * in some unspecified order. This can be used, e.g. to reduce all
	 * matrix entries modulo a prime before passing the matrix into an
	 * algorithm.
	 */
	template < class _Matrix >
	class BlasSubmatrix<_Matrix>::ConstIterator {
	public:
		//! @internal Null constructor
		ConstIterator (){}


		/*! @internal
		 * @brief NO DOC
		 */
		ConstIterator (const typename BlasSubmatrix<_Matrix>::pointer& cur,
                       const size_t c_dim,
                       const size_t stride,
                       const size_t c_idx) :
			_cur (cur), _c_dim (c_dim), _stride(stride), _c_idx (c_idx)
		{}

		/*! @internal
		 * @brief copy operator.
		 * @param r Iterator to copy.
		 */
		ConstIterator& operator = (const ConstIterator& r)
		{
			_cur = r._cur;
			_c_dim = r._c_dim;
			_stride = r._stride;
			_c_idx = r._c_idx;
			return *this;
		}

		/*! @internal
		 * increment.
		 * ??
		 */
		ConstIterator& operator ++()
		{
			if (_c_idx < _c_dim - 1){
				++_cur; ++_c_idx;
			}
			else {
				linbox_check(_stride >= _c_dim);
				_cur = _cur + (ptrdiff_t)(_stride - _c_dim + 1);
				_c_idx = 0;
			}

			return *this;
		}

		/*! @internal
		 * increment.
		 * ??
		 */
		ConstIterator& operator++ (int)
		{
			return this->operator++ ();
		}

		/*! @internal
		 * @brief  operator !=.
		 * @param r Iterator to test inequaltity from.
		 */
		bool operator != (const ConstIterator& r) const
		{
			return (_cur != r._cur) || (_c_dim != r._c_dim) || (_stride != r._stride) || (_c_idx != r._c_idx);
		}

		//! @internal operator *.
		const typename BlasSubmatrix<_Matrix>::Element& operator * () const
		{
			return *_cur;
		}

	protected:
        typename BlasSubmatrix<_Matrix>::pointer _cur;
		size_t _c_dim;
		size_t _stride;
		size_t _c_idx;
	};



   
	template <class _Matrix>
	template<typename _Tp1, typename _Rep2>
	struct BlasSubmatrix< _Matrix>::rebind {
		typedef BlasMatrix<_Tp1,_Rep2> other;

		void operator() (other & Ap, const Self_t& A) {
			typedef typename BlasSubmatrix<_Matrix>::ConstIterator ConstSelfIterator ;
			typedef typename other::Iterator OtherIterator ;
			OtherIterator    Ap_i;
			ConstSelfIterator A_i;
			Hom<Field, _Tp1> hom(A. field(), Ap. field());
			for (A_i = A. Begin(), Ap_i = Ap.Begin();
			     A_i != A. End(); ++ A_i, ++ Ap_i)
				hom.image (*Ap_i, *A_i);
		}
	};

} // LinBox



/* // Local Variables: */
/* // mode: C++ */
/* // tab-width: 4 */
/* // indent-tabs-mode: nil */
/* // c-basic-offset: 4 */
/* // End: */
/* // vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s */

