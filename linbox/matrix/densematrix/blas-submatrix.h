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

/*! @file matrix/densematrix/new-blas-submatrix.h
 * @ingroup densematrix
 */

#ifndef __LINBOX_densematrix_blas_submatrix_H
#define __LINBOX_densematrix_blas_submatrix_H

#include "fflas-ffpack/fflas-ffpack.h"
#include "linbox/field/rebind.h"
#include "linbox/matrix/densematrix/blas-matrix-iterator.h"

namespace LinBox
{
    // Forward Declaration
    template <class _Field, class _Storage>
    class BlasMatrix;

    /* Blas Submatrix */
  /*! Dense Submatrix representation.
   * @ingroup matrix
   * A @ref BlasSubmatrix is a matrix of \p _Field::Element, with the structure of BLAS matrices.
   * It is basically a read/write view on a vector of \p _Field::Element.
   * In the Mother model, a @ref BlasSubmatrix is not allocated.
   * <p>
   * This matrix type conforms to the same interface as @ref BlasMatrix,
   * except that you cannot resize it. It represents a submatrix of a dense
   * matrix. Upon construction, one can freely manipulate the entries in the
   * DenseSubmatrix, and the corresponding entries in the underlying
   * @ref BlasMatrix will be modified.
   */

    template <typename _Matrix>
    class MatrixEltPointer {
    public:
        typedef typename _Matrix::Field::Element_ptr     pointer;
    };
    template <typename _Matrix>
    class MatrixEltPointer <const _Matrix> {
    public:
        typedef typename _Matrix::Field::ConstElement_ptr pointer;
    };


    template <class _Matrix>
    class BlasSubmatrix {
    public :
        typedef typename _Matrix::Field                     Field;
        typedef typename Field::Element                   Element;    //!< Element type
        typedef typename _Matrix::Storage                 Storage;
        typedef typename _Matrix::RawStorage           RawStorage;
        typedef BlasSubmatrix<_Matrix>                     Self_t;    //!< Self type
        typedef typename MatrixEltPointer<_Matrix>::pointer                    pointer;    //!< pointer type to elements
        typedef typename MatrixEltPointer<const _Matrix>::pointer        const_pointer;    //!< const pointer type to elements
        typedef  _Matrix                               matrixType;    //!< matrix type
        typedef Self_t                              subMatrixType;
        typedef BlasSubmatrix<const _Matrix>   constSubMatrixType;
        typedef BlasSubmatrix<typename std::remove_const<_Matrix>::type>  nonconstSubMatrixType;

        // row and col types are unified to be a BlasSubvector
        typedef typename Storage::subVectorType           subVectorType;
        typedef typename Storage::constSubVectorType constSubVectorType;
        typedef subVectorType Row ;
        typedef subVectorType Col;
        typedef constSubVectorType ConstRow;
        typedef constSubVectorType ConstCol;


    protected:
        pointer         _ptr;  //!< pointer to the first elt of the submatrix
        size_t          _row;  //!< row dimension of Submatrix
        size_t          _col;  //!< col dimension of Submatrix
        size_t       _stride;  //!< stride of the original matrix
        const Field&  _field;

    public:

        //////////////////
        // CONSTRUCTORS //
        //////////////////

        /*  constructors */

        /** Constructor from an existing @ref BlasMatrix and dimensions.
         * \param M Pointer to @ref BlasMatrix of which to construct submatrix
         * \param rowbeg Starting row
         * \param colbeg Starting column
         * \param Rowdim Row dimension
         * \param Coldim Column dimension
         */
        BlasSubmatrix (matrixType &M,
                       size_t rowbeg,
                       size_t colbeg,
                       size_t Rowdim,
                       size_t Coldim);

        /** Constructor from an existing submatrix and dimensions
         * @param SM Constant reference to BlasSubmatrix from which to
         *           construct submatrix
         * @param rowbeg Starting row
         * @param colbeg Starting column
         * @param Rowdim Row dimension
         * @param Coldim Column dimension
         */
        BlasSubmatrix (Self_t    &SM,
                       size_t rowbeg,
                       size_t colbeg,
                       size_t Rowdim,
                       size_t Coldim);


        /** Constructor from an existing @ref BlasMatrix
         * \param M Pointer to @ref BlasMatrix of which to construct submatrix
         */
        BlasSubmatrix (matrixType &M);

        /** Constructor from a raw pointer
         * \param M Pointer to @ref BlasMatrix of which to construct submatrix
         */
        BlasSubmatrix (const Field& F,
                       pointer ptr,
                       size_t Rowdim,
                       size_t Coldim,
                       size_t stride);


        /** Copy constructor (virtual copy)
         */
        BlasSubmatrix(const constSubMatrixType &M);
        BlasSubmatrix(nonconstSubMatrixType &M);
        BlasSubmatrix(const nonconstSubMatrixType &M);

         /** Move constructor / Move operator
         */
        BlasSubmatrix(Self_t &&M)=default;
        Self_t& operator=(Self_t&& M)=default;

        //! (copying data) -> works only if dimensions are the same
        template<class _AnyMatrix>
        Self_t& copy (const _AnyMatrix& A);

            /** (swaping local data to same dimension (sub)part of A)
             *  -> works only if both dimensions of A are larger
             */
        template<class _AnyMatrix>
        Self_t& swap (_AnyMatrix& A);

        Self_t& operator=(const Self_t& M)=delete;
        // {
        //     if (&M!=this){
        //         _ptr = M.getPointer();
        //         _row = M.rowdim();
        //         _col = M.coldim();
        //         _stride= M.getStride();
        //         _field = M.field();
        //         std::cout<<"BlasSubmatrix operator=\n";
        //     }
        //     return *this;
        // }

        template <typename _TP1, typename _Rep2 = typename Rebind<RawStorage, _TP1>::other >
        struct rebind;
        /*  Members  */

        //////////////////
        //  ACCESSORS  //
        //////////////////

        /** Get the number of rows in the matrix
         * @return Number of rows in matrix
         */
        size_t rowdim () const { return _row;}

        /** Get the number of columns in the matrix
         * @return Number of columns in matrix
         */
        size_t coldim () const { return _col;}

        /** Get the stride of the matrix.
         * @return stride of submatrix (number of cols of dense base matrix)
         */
        size_t getStride() const { return _stride;}

        /** Get the field of the BlasSubMatrix
            @return Const reference to Field
        */
        const Field& field() const { return _field ;}

        /*! @internal
         * Get pointer to the matrix data (read/write access will depend on the type of the template parameter _Matrix, i.e. const or not)
         */
        pointer getPointer() {return _ptr;}
        const_pointer getPointer() const {return _ptr;}
        const_pointer getConstPointer() const {return _ptr;}

        ///////////////////
        //      I/O      //
        ///////////////////

        /** Read the matrix from an input stream.
         * @param file Input stream from which to read
         * @bug reading a submatrix should not be allowed !!
         */
        // template<class Field>
        std::istream& read (std::istream &file); // autodetect ?

        /** Write the matrix to an output stream.
         * @param os Output stream to which to write
         * @param f write in some format (@ref Tag::FileFormat::Format). Default is MM's.
         */
        //std::ostream &write (std::ostream &os,Tag::FileFormat f = Tag::FileFormat::MatrixMarket )const;
        std::ostream &write (std::ostream &os, Tag::FileFormat f = Tag::FileFormat::Plain) const;
        //////////////////
        //   ELEMENTS   //
        //////////////////


        /** Set the entry at (i, j).
         * @param i Row number, 0...rowdim () - 1
         * @param j Column number 0...coldim () - 1
         * @param a_ij Element to set
         */
        void setEntry (size_t i, size_t j, const Element &a_ij) ;

        /** Get a writeable reference to an entry in the matrix.
         * @param i Row index of entry
         * @param j Column index of entry
         * @return Reference to matrix entry
         */
        Element& refEntry (size_t i, size_t j) ;

        /** Get a read-only individual entry from the matrix.
         * @param i Row index
         * @param j Column index
         * @return Const reference to matrix entry
         */
        const Element& getEntry (size_t i, size_t j) const ;

        /** Get an entry and store it in the given value.
         * This form is more in the Linbox style and is provided for interface
         * compatibility with other parts of the library
         * @param x Element in which to store result
         * @param i Row index
         * @param j Column index
         * @return Reference to x
         */
        Element& getEntry (Element &x, size_t i, size_t j) const ;



        ///////////////////
        //   BLACK BOX   //
        ///////////////////

        //!@bug every vector we use here should have a stride/be blas vectors so it's not really templated by Vector1 Vector2 in general
        template <class Vector1, class Vector2>
        Vector1&  apply (Vector1& y, const Vector2& x) const
        {
            FFLAS::fgemv(_field, FFLAS::FflasNoTrans, _row,_col, _field.one, _ptr, _stride,
                         x.getConstPointer(), x.getInc(), _field.zero, y.getPointer(),y.getInc());
            return y;
        }
        // NEED TO BE REMOVED -> just to be done after new-blas-subvector is merged
        std::vector<Element>&  apply (std::vector<Element>& y, const std::vector<Element>& x) const
        {
            FFLAS::fgemv(_field, FFLAS::FflasNoTrans, _row,_col, _field.one, _ptr, _stride,
                         x.data(), 1, _field.zero, y.data(),1);
            return y;
        }

        //! @bug use Matrix domain
        //!@bug since removal of ApplyDomain this does not handle the case where Field if Givaro::Extension needed for charpoly computation
        template <class Vector1, class Vector2>
        Vector1&  applyTranspose (Vector1& y, const Vector2& x) const
        {
            FFLAS::fgemv(_field, FFLAS::FflasTrans, _row,_col, _field.one, _ptr, _stride,
                         x.getConstPointer(), x.getInc(), _field.zero, y.getPointer(),y.getInc());
            return y;
        }
        // NEED TO BE REMOVED -> just to be done after new-blas-subvector is merged
        std::vector<Element>&  applyTranspose (std::vector<Element>& y, const std::vector<Element>& x) const
        {
            FFLAS::fgemv(_field, FFLAS::FflasTrans, _row,_col, _field.one, _ptr, _stride,
                         x.data(), 1, _field.zero, y.data(),1);
            return y;
        }



        void zero() {
            FFLAS::fzero(field(),_row,_col,_ptr,_stride);
        }
		// init to random field elements
		void random()
		{
            typename Field::RandIter G(field());
            random(G);
        }
        template<class RandIter>
		void random(RandIter& G)
		{
            FFLAS::frand(field(),G, _row,_col,_ptr,_stride);
        }


        ///////////////////
        //   ITERATORS   //
        ///////////////////

        //! @name Forward declaration of Raw Iterators.
        //@{
        class Iterator  ;
        class ConstIterator ;
        //@}


        /** @name typedef'd Row Iterators.
         *\brief
         * The row iterator gives the rows of the
         * matrix in ascending order. Dereferencing the iterator yields
         * a row vector in dense format
         * @{
         */
        using RowIterator      = BlasMatrixIterator<Field, Storage, subVectorType>;
        using ConstRowIterator = BlasMatrixIterator<Field, Storage, constSubVectorType>;

        RowIterator      rowBegin ()       { return      RowIterator (field(), _ptr, _col, 1, _stride);}
        ConstRowIterator rowBegin () const { return ConstRowIterator (field(), _ptr, _col, 1, _stride);}
        RowIterator      rowEnd ()         { return      RowIterator (field(), _ptr+_row*_stride, _col, 1, _stride);}
        ConstRowIterator rowEnd   () const { return ConstRowIterator (field(), _ptr+_row*_stride, _col, 1, _stride);}
        //@} Row Iterators

        /** @name typedef'd Column Iterators.
         *\brief
         * The columns iterator gives the columns of the
         * matrix in ascending order. Dereferencing the iterator yields
         * a column vector in dense format
         * @{
         */
        using ColIterator      = BlasMatrixIterator<Field, Storage, subVectorType>;
        using ConstColIterator = BlasMatrixIterator<Field, Storage, constSubVectorType>;

        ColIterator      colBegin ()       { return       ColIterator (field(), _ptr, _row, _stride, 1);}
        ConstColIterator colBegin () const { return  ConstColIterator (field(), _ptr, _row, _stride, 1);}
        ColIterator      colEnd ()         { return       ColIterator (field(), _ptr+_col, _row, _stride, 1);}
        ConstColIterator colEnd ()   const { return  ConstColIterator (field(), _ptr+_col, _row, _stride, 1);}
        //@} // Column Iterators


        Iterator      Begin ()        { return      Iterator (_ptr, _col, _stride, 0); }
        ConstIterator Begin () const  { return ConstIterator (_ptr, _col, _stride, 0); }
        Iterator      End ()          { return      Iterator (_ptr+_row * _stride, _col, _stride, 0); }
        ConstIterator End ()   const  { return ConstIterator (_ptr+_row * _stride, _col, _stride, 0); }


        using IndexedIterator      = BlasMatrixIndexedIterator<Field,       pointer,       Element>;
        using ConstIndexedIterator = BlasMatrixIndexedIterator<Field, const_pointer, const Element>;

        IndexedIterator      IndexedBegin ()        { return      IndexedIterator (coldim (), _stride, 0, 0, _ptr);}
        ConstIndexedIterator IndexedBegin () const  { return ConstIndexedIterator (coldim (), _stride, 0, 0, _ptr);}
        IndexedIterator      IndexedEnd   ()        { return      IndexedIterator (coldim (), _stride, rowdim (), 0, _ptr);}
        ConstIndexedIterator IndexedEnd   () const  { return ConstIndexedIterator (coldim (), _stride, rowdim (), 0, _ptr);}

        /*!  operator[].
         * Retrieve a reference to a row
         * @param i Row index
         */
        subVectorType      operator[] (size_t i)        { return      subVectorType (field(), _ptr+i*_stride, 1, _col);}
        constSubVectorType operator[] (size_t i) const  { return constSubVectorType (field(), _ptr+i*_stride, 1, _col);}

    };

}
#include "linbox/matrix/densematrix/blas-submatrix.inl"
#endif

/* // Local Variables: */
/* // mode: C++ */
/* // tab-width: 4 */
/* // indent-tabs-mode: nil */
/* // c-basic-offset: 4 */
/* // End: */
/* // vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s */
