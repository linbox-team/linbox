/*
 * Copyright (C) 2014 the LinBox
 *
 * Written by :
 * BB <bbboyer@ncsu.edu>
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

/*! @file matrix/SparseMatrix/sparse-csc-matrix.h
 * @ingroup sparsematrix
 * @brief
 */


#ifndef __LINBOX_sparse_matrix_sparse_csc_matrix_H
#define __LINBOX_sparse_matrix_sparse_csc_matrix_H

#include <utility>
#include <iostream>
#include <algorithm>

#include "linbox/linbox-config.h"
#include "linbox/util/debug.h"
#include "linbox/field/hom.h"
#include "sparse-domain.h"

#ifndef LINBOX_CSC_TRANSPOSE
#define LINBOX_CSC_TRANSPOSE 1000
#endif

namespace LinBox
{


	/** Sparse matrix, Coordinate storage.
	 *
	 * \ingroup matrix
	 * \ingroup sparse
	 */
	template<class _Field>
	class SparseMatrix<_Field, SparseMatrixFormat::CSC > {
	public :
		typedef _Field                             Field ; //!< Field
		typedef typename _Field::Element         Element ; //!< Element
		typedef const Element               constElement ; //!< const Element
		typedef SparseMatrixFormat::CSC          Storage ; //!< Matrix Storage Format
		typedef SparseMatrix<_Field,Storage>      Self_t ; //!< Self type
		typedef typename Vector<Field>::SparseSeq    Row ; //!< @warning this is not the row type. Just used for streams.
		// typedef Vector<_Field,VectorStorage::Sparse> Rep ;

		/*! Constructors.
		 * @todo convert from other matrix.
		 *
		 */
		//@{
#if 0 /*  No empty CSTOR */
		SparseMatrix<_Field, SparseMatrixFormat::CSC> () :
			_rownb(0),_colnb(0)
			,_nbnz(0)
			,_start(1,0),_rowid(0),_data(0)
			, _field()
			, _helper()
		{
			_start[0] = 0 ;
		}
#endif

		SparseMatrix<_Field, SparseMatrixFormat::CSC> (const _Field & F) :
			_rownb(0),_colnb(0)
			,_nbnz(0)
			,_start(1,0)
			,_rowid(0)
			,_data(0)
			, _field(F)
			, _helper()
		{
			_start[0] = 0 ;
		}

		SparseMatrix<_Field, SparseMatrixFormat::CSC> (const _Field & F, size_t m, size_t n) :
			_rownb(m),_colnb(n)
			,_nbnz(0)
			,_start(n+1,0)
			,_rowid(0)
			,_data(0)
			, _field(F)
			, _helper()
		{
			_start[0] = 0 ;
		}

		SparseMatrix<_Field, SparseMatrixFormat::CSR> (const _Field & F,
							       size_t m, size_t n,
							       size_t z) :
			_rownb(m),_colnb(n)
			, _nbnz(z)
			,_start(n+1,0)
			, _rowid(z)
			,_data(z)
			, _field(F)
			, _helper()
		{
			_start[0] = 0 ;
		}

		SparseMatrix<_Field, SparseMatrixFormat::CSR> (const SparseMatrix<_Field, SparseMatrixFormat::CSR> & S) :
			_rownb(S._rownb),_colnb(S._colnb)
			,_nbnz(S._nbnz)
			,_start(S._start)
			, _rowid(S._rowid)
			,_data(S._data)
			, _field(S._field)
			, _helper()
		{
		}

#if 0
		template<class _OtherField>
		SparseMatrix<_Field, SparseMatrixFormat::COO> (const SparseMatrix<_OtherField, SparseMatrixFormat::COO> & S) :
			_rownb(S._rownb),_colnb(S._colnb),
			_nbnz(S._nbnz),
			_rowid(S._rowid),_rowid(S._rowid),_data(S._data),
			_field(S._field)
		{}
#endif

		template<typename _Tp1, typename _Rw1 = SparseMatrixFormat::CSR>
		struct rebind {
			typedef SparseMatrix<_Tp1, _Rw1> other;
		private:

			//! @todo indexed iterator
			template<class _Rw>
			void rebindMethod(SparseMatrix<_Tp1, _Rw> & Ap, const Self_t & A  /*, IndexedCategory::HasNext */)
			{
				typename _Tp1::Element e;
				Hom<typename Self_t::Field, _Tp1> hom(A.field(), Ap.field());

				size_t i, j ;
				Element f ;
				A.firstTriple();
				while ( A.nextTriple(i,j,f) ) {
					linbox_check(i < A.rowdim() && j < A.coldim()) ;
					hom. image ( e, f) ;
					if (! Ap.field().isZero(e) )
						Ap.setEntry(i,j,e);
				}
				A.firstTriple();
				Ap.finalize();
			}

			void rebindMethod(SparseMatrix<_Tp1, SparseMatrixFormat::CSC>  & Ap, const Self_t & A /*,  IndexedCategory::HasNext*/)
			{
				// we don't use nextTriple because we can do better.

				typename _Tp1::Element e;

				Hom<typename Self_t::Field, _Tp1> hom(A.field(), Ap.field());
				size_t j = 0 ;
				Ap.setStart(A.getStart());
				std::vector<size_t> offset(A.coldim()+1,0UL);
				bool changed = false ;
				for (size_t i = 0 ; i < A.coldim() ; ++i) {
					for (size_t k = A.getStart(i) ; k < A.getEnd(i) ; ++k) {
						hom. image ( e, A.getData(k) );
						if (!Ap.field().isZero(e)) {
							Ap.setRowid(j,A.getRowid(k));
							Ap.setData(j,e);
							++j;
						}
						else
						{
							offset[i+1] += 1 ;
							changed = true ;
						}
					}
				}
				if (changed) {
					for (size_t i = 0 ; i < A.coldim() ; ++i) {
						offset[i+1] += offset[i] ;
					}
					for (size_t i = 1 ; i <= A.coldim() ; ++i) {
						Ap.setStart(i,A.getStart(i)-offset[i]);
					}
				}


				if (j != Ap.size())
					Ap.resize(j);

			}

		public:

			void operator() (other & Ap, const Self_t& A)
			{
				rebindMethod(Ap, A );

			}

		};

		template<typename _Tp1, typename _Rw1>
		SparseMatrix (const SparseMatrix<_Tp1, _Rw1> &S, const Field& F) :
			_rownb(S.rowdim()),_colnb(S.coldim())
			,_nbnz(S.size())
			, _start(S.coldim()+1,0)
			, _rowid(S.size())
			,_data(S.size())
			, _field(F)
		{
			typename SparseMatrix<_Tp1,_Rw1>::template rebind<Field,Storage>()(*this, S);
			finalize();
		}


#if 0
		void merge(const SparseMatrix<_Field, SparseMatrixFormat::CSR> & S)
		{
			for (size_t i_ex = 0 ; i_ex < S.rowdim() ; ++i_ex)
				for (size_t k = S.getStart(i_ex) ; k < S.getEnd(i_ex) ; ++k) {
					// can be faster by iterating over both matrices
					setEntry(i_ex,S.getColid(k),S.getData(k));
				}
		}
#endif


		template<class VectStream>
		SparseMatrix<_Field, SparseMatrixFormat::CSC> (const _Field & F, VectStream & stream) :
			_rownb(stream.size()),_colnb(stream.dim())
			, _nbnz(0)
			, _start(_colnb+1,0)
			, _rowid(0)
			,_data(0)
			, _field(F)
		{
			SparseMatrix<_Field,SparseMatrixFormat::CSR> tmp(F,stream);
			importe(tmp);
			finalize();
		}

		SparseMatrix<_Field, SparseMatrixFormat::CSR> ( MatrixStream<Field>& ms ):
			_rownb(0),_colnb(0)
			,_nbnz(0)
			,_start(1,0)
			,_rowid(0)
			,_data(0)
			,_field(ms.field())
		{
			SparseMatrix<_Field,SparseMatrixFormat::CSR> tmp(ms);
			importe(tmp);
			finalize();
		}

		void resize(size_t zz)
		{
#ifndef NDEBUG
			if (zz < _nbnz) {
				std::cerr << "*** Warning *** you are possibly loosing data (CSR resize)" << std::endl;
				// could be a commentator()...
			}
#endif
			_rowid.resize(zz);
			_data.resize(zz);
			_nbnz = zz ;
		}

		void resize(const size_t & mm, const size_t & nn, const size_t & zz = 0)
		{
			_rownb = mm ;
			_colnb = nn ;

			_start.resize(nn+1);
			resize(zz);
		}

		/*! Default converter.
		 * @param S a sparse matrix in any storage.
		 */
		template<class _OtherStorage>
		SparseMatrix<_Field, SparseMatrixFormat::CSR> (const SparseMatrix<_Field, _OtherStorage> & S) :
			_rownb(S.rowdim()),_colnb(S.coldim()),
			_start(S.coldim()+1,0),_rowid(S.size()),_data(S.size()),
			_field(S.field())
		{
			this->importe(S); // convert Temp from anything
			finalize();
		}



		//@}
		/*! Conversions.
		 * Any sparse matrix has a converter to/from CSR.
		 * A specialisation can skip the temporary CSR matrix created.
		 */
		//@{
		/*! Import a matrix in COO format to CSR.
		 * @param S COO matrix to be converted in CSR
		 */
		void importe(const SparseMatrix<_Field,SparseMatrixFormat::COO> &S)
		{
			resize(S.rowdim(), S.coldim(), S.size());

			for (size_t i = 0 ; i < _nbnz ; ++i)
				_start[S.colid(i)+1] += 1 ;

			throw NotImplementedYet("import CSC from COO");

			finalize();

		}

		/*! Import a matrix in CSR format to CSR.
		 * @param S CSR matrix to be converted in CSR
		 */

		void importe(const SparseMatrix<_Field,SparseMatrixFormat::CSR> &S)
		{
		}
			throw NotImplementedYet("import CSC from CSR");

		/*! Export a matrix in CSR format from COO.
		 * @param S CSR matrix to be converted from COO
		 */
		SparseMatrix<_Field,SparseMatrixFormat::CSR > &
		exporte(SparseMatrix<_Field,SparseMatrixFormat::CSR> &S)
		{
			linbox_check(consistent());
			throw NotImplementedYet("export CSC to CSR");


			return S ;
		}

		SparseMatrix<_Field,SparseMatrixFormat::COO > &
		exporte(SparseMatrix<_Field,SparseMatrixFormat::COO> &S)
		{
			throw NotImplementedYet("export CSC to CSR");
			return S ;
		}


		//@}

		/*! In place transpose. Not quite...
		*/
		void transposeIn()
		{
			Self_t Temp(*this);
			Temp.transpose(*this);
		}

		/*! Transpose the matrix.
		 *  @param S [out] transpose of self.
		 *  @return a reference to \p S.
		 */
		Self_t &
		transpose(Self_t &S) const
		{


			throw NotImplementedYet("tranpose CSC ");

			return S;


		}

		/*! number of rows.
		 * @return row dimension.
		 */
		size_t rowdim() const
		{
			return _rownb ;
		}

		/*! number of columns.
		 * @return column dimension
		 */
		size_t coldim() const
		{
			return _colnb ;
		}

		/*! Number of non zero elements in the matrix.
		 * or at least the size of the _data if
		 * @return size of the _data.
		 */
		size_t size() const
		{
			// return _data.size();
			return _nbnz ;
		}

		void setSize(const size_t & z)
		{
			_nbnz = z ;
		}


		/** Get a read-only individual entry from the matrix.
		 * @param i Row _rowid
		 * @param j Column _rowid
		 * @return Const reference to matrix entry
		 */
		constElement &getEntry(const size_t &i, const size_t &j) const
		{
			throw NotImplementedYet("getEntry");
			linbox_check(i<_rownb);
			linbox_check(j<_colnb);



		}

		Element      &getEntry (Element &x, size_t i, size_t j) const
		{
			return x = getEntry (i, j);
		}

		void appendEntry(const size_t &i, const size_t &j, const Element& e)
		{
			throw NotImplementedYet("appendEntry");
			return;

		}

		/// make matrix ready to use after a sequence of setEntry calls.
		void finalize()
		{
			if (_start[coldim()] != _nbnz) { /* if it is so, then all before are 0 and we are fine... */
				for (size_t i = 2 ; i <= rowdim() ; ++i)
					_start[i] += _start[i-1];
				linbox_check(_start[rowdim()] == _nbnz);
			}
			_triples.reset();

		} // end construction after a sequence of setEntry calls.

		/** Set an individual entry.
		 * Setting the entry to 0 will not remove it from the matrix
		 * @param i Row _rowid of entry
		 * @param j Column _rowid of entry
		 * @param value Value of the new entry
		 * @todo make it faster if i is 0 or m-1 ?
		 * @warning if this is used to build a matrix and this matrix is "well formed",
		 * it can be sped up (no checking that the entry already exists).
		 */
		void setEntry(const size_t &i, const size_t &j, const Element& e
			     )
		{
			throw NotImplementedYet("setEntry");
			linbox_check(i<_rownb);
			linbox_check(j<_colnb);
			linbox_check(_start.size() == _colnb+1);
			linbox_check(_start[0] == 0);

			if (field().isZero(e)) {
				return clearEntry(i,j);
			}

				}

#if 0
		/** Get a writeable reference to an entry in the matrix.
		 * If there is no entry at the position (i, j), then a new entry
		 * with a value of zero is inserted and a reference  to it is
		 * returned.
		 * @param i Row _rowid of entry
		 * @param j Column _rowid of entry
		 * @return Reference to matrix entry
		 */
		Element &refEntry(const size_t &i, const size_t&j)
		{
			linbox_check(i<_rownb);
			linbox_check(j<_colnb);
			// Could be improved by adding an initial guess j/rowdim*size()
			typedef typename std::vector<size_t>::iterator myIterator ;

			size_t ibeg = _start[i];
			size_t iend = _start[i+1];
			if (ibeg==iend) {
				for (size_t k = i+1 ; k <= _rownb ; ++k) _start[k] +=1 ;
				_rowid.insert(_rowid.begin()+ibeg,j);
				_data.insert( _data.begin() +ibeg,field().zero);
				return _data[ibeg];
			}
			myIterator beg = _rowid.begin() ;
			myIterator low = std::lower_bound (beg+(ptrdiff_t)ibeg, beg+(ptrdiff_t)iend, j);
			if (low == beg+(ptrdiff_t)iend) {
				for (size_t k = i+1 ; k <= _rownb ; ++k) _start[k] +=1 ;
				_rowid.insert(_rowid.begin()+ibeg,j);
				_data.insert( _data.begin() +ibeg,field().zero);
				return _data[ibeg];
			}
			else {
				size_t la = low-_rowid.begin() ;
				return _data[la] ;
			}
		}
#endif

		/** Write a matrix to the given output stream using field read/write.
		 * @param os Output stream to which to write the matrix
		 * @param format Format with which to write
		 */
		std::ostream & write(std::ostream &os
				     , LINBOX_enum(Tag::FileFormat) format  = Tag::FileFormat::MatrixMarket) const
		{
			return SparseMatrixWriteHelper<Self_t>::write(*this,os,format);
		}


		/** Read a matrix from the given input stream using field read/write
		 * @param is Input stream from which to read the matrix
		 * @param format Format of input matrix
		 * @return ref to \p is.
		 */
		std::istream& read (std::istream &is
				    , LINBOX_enum(Tag::FileFormat) format = Tag::FileFormat::Detect)
		{
			return SparseMatrixReadHelper<Self_t>::read(*this,is,format);
		}

		/*! @internal
		 * @brief Deletes the entry.
		 * Deletes \c A(i,j) if it exists.
		 * @param i row _rowid
		 * @param j col _rowid
		 */
		void clearEntry(const size_t &i, const size_t &j)
		{
			linbox_check(i<_rownb);
			linbox_check(j<_colnb);
			throw NotImplementedYet("clearEntry");

		}

		/*! @internal
		 * @brief cleans 0 entries.
		 */
		void clean()
		{
			throw NotImplementedYet("clean");
		}

		// y= Ax
		// y[i] = sum(A(i,j) x(j)
		// start(i)<k < start(i+1) : _delta[k] = A(i,colid(k))
		template<class inVector, class outVector>
		outVector& apply(outVector &y, const inVector& x, const Element & a ) const
		{
			// linbox_check(consistent());
			prepare(field(),y,a);


			// std::cout << "apply" << std::endl;
			for (size_t i = 0 ; i < _colnb ; ++i) {
				for (size_t k = _start[i] ; k < _start[i+1] ; ++k)
					// field().axpyin( y[i], _data[k], x[_rowid[k]] ); //! @todo delay !!!
					field().axpyin(y[colid[k]],data[k],x[k]);
			}

			return y;
		}


		// y= A^t x
		// y[i] = sum(A(j,i) x(j)
		template<class inVector, class outVector>
		outVector& applyTranspose(outVector &y, const inVector& x, const Element & a) const
		{
			for (size_t i = 0 ; i < _colnb ; ++i) {
				for (size_t k = _start[i] ; k < _start[i+1] ; ++k)
					// field().axpyin( y[i], _data[k], x[_rowid[k]] ); //! @todo delay !!!
					field().axpyin(y[i],data[k],x[colid[k]]);
			}
			return y;
		}



		template<class inVector, class outVector>
		outVector& apply(outVector &y, const inVector& x ) const
		{
			return apply(y,x,field().zero);
		}

		template<class inVector, class outVector>
		outVector& applyTranspose(outVector &y, const inVector& x ) const
		{
			return applyTranspose(y,x,field().zero);
		}

		const Field & field()  const
		{
			return _field ;
		}

		//! @todo
		bool consistent() const
		{
			return true ;
		}

	private :

		class Helper {
			bool _useable ;
			bool _optimized ;
			bool blackbox_usage ;
			Self_t *_AT ;
		public:

			Helper() :
				_useable(false)
				,_optimized(false)
				, blackbox_usage(true)
				, _AT(NULL)
			{}

			~Helper()
			{
				if ( _AT ) {
					delete _AT ;
				}
			}

			bool optimized(const Self_t & A)
			{
				if (!_useable) {
					getHelp(A);
					_useable = true;
				}
				return	_optimized;
			}

			void getHelp(const Self_t & A)
			{
				if ( A.size() > LINBOX_CSR_TRANSPOSE ) { // and/or A.rowDensity(), A.coldim(),...
					// std::cout << "optimizing..." ;
					_optimized = true ;
					_AT = new Self_t(A.field(),A.coldim(),A.rowdim());
					A.transpose(*_AT); // or convert to CSC
					// std::cout << "done!" << std::endl;
				}
			}

			const Self_t & matrix() const
			{
				return *_AT ;
			}

		};

	public:
		// pseudo iterators
		size_t getStart(const size_t & i) const
		{
			return _start[i];
		}

		size_t getEnd(const size_t & i) const
		{
			return _start[i+1];
		}


		void setStart(const size_t &i, const size_t & j)
		{
			if (i > _rownb) this->resize(i,_colnb,_nbnz);
			_start[i] = j ;
		}

		void setStart(const std::vector<size_t> &  new_start)
		{
			// linbox_check(_start.size() == new_start.size());
			_start = new_start ;
		}

		std::vector<size_t>  getStart( ) const
		{
			return _start ;
		}

		size_t rowLength(const size_t & i)
		{
			return getEnd(i) - getStart(i) ;
		}

		// remove empty last lines
		void compress()
		{
			size_t i = rowdim();
			for ( ; i != 0 ; --i) {
				if ( _start[i] != _start[i-1] )
					break;
			}
			_rownb = i ;

		}

		size_t getRowid(const size_t & i) const
		{
			return _rowid[i];
		}

		void setRowid(const size_t & i, const size_t & j)
		{
			if (i>=_nbnz) this->resize(i+1);
			linbox_check(i <= _rowid.size())
			_rowid[i]=j;
		}

		void setRowid(std::vector<size_t> new_rowid)
		{
			_rowid = new_rowid ;
		}

		std::vector<size_t>  getRowid( ) const
		{
			return _rowid ;
		}

		const Element & getData(const size_t & i) const
		{
			return _data[i];
		}

		void setData(const size_t & i, const Element & e)
		{
			if (i>=_nbnz) this->resize(i+1);
			linbox_check(i <= _data.size());
			field().assign(_data[i],e);
		}

		void setData(const std::vector<Element> & new_data)
		{
			_data = new_data ;
		}

		std::vector<Element>  getData( ) const
		{
			return _data ;
		}

		void firstTriple() const
		{
			_triples.reset();
		}

		bool nextTriple(size_t & i, size_t &j, Element &e) const
		{
			thow NotImplementedYet("nextTriple");

			ptrdiff_t idx =_triples.next( _start );
			i = _triples._row ;
			if (idx >= _nbnz || i >= _rownb ) {
				_triples.reset() ;
				return false;
			}


			j = _rowid[idx];
			e = _data[idx];

			return true;
		}
#if 0 /*  not implemented yet */
		template<class element_iterator, class Element>
		class _Iterator {
		private :
			element_iterator _data_it ;
			const element_iterator _data_beg ;
			const element_iterator _data_end ;
		public:
			typedef Element value_type ;
			_Iterator(element_iterator e_beg, element_iterator e_end) :
				_data_it(e_beg)
				, _data_beg(e_beg)
				, _data_end(e_end)
			{}

			_Iterator (const _Iterator &iter) :
				_data_it(iter._data_it)
				, _data_beg(iter._data_beg)
				, _data_end(iter._data_end)
			{}

			_Iterator &operator = (const _Iterator &iter)
			{
				_data_it  = iter._data_it  ;
				_data_beg  = iter._data_beg  ;
				_data_end  = iter._data_end  ;

				return *this;
			}

			bool operator == (const _Iterator &i) const
			{
				return  (_data_it == i._data_it) && (_data_beg == i._data_beg) && (_data_end == i._data_end);
			}

			bool operator != (const _Iterator &i) const
			{
				return  (_data_it != i._data_it) || (_data_beg != i._data_beg) || (_data_end != i._data_end);
			}

			_Iterator &operator ++ ()
			{
				++_data_it ;

				if (_data_it == _data_end)
					return *this ;
			}

			_Iterator operator ++ (int)
			{
				_Iterator tmp = *this;
				++(*this);
				return tmp;
			}

			_Iterator &operator -- ()
			{
				--_data_it ;

				if (_data_it == _data_beg)
					return *this ;
			}

			_Iterator operator -- (int)
			{
				_Iterator tmp = *this;
				--(*this);
				return tmp;
			}

			value_type &operator * ()
			{
				return *_data_it;
			}

			value_type *operator -> ()
			{
				return _data_it ;
			}

			const value_type &operator*() const
			{
				return *_data_it;
			}

			const value_type *operator -> () const
			{
				return _data_it ;
			}

			const value_type &value() const
			{
				return *_data_it;
			}

		};

		template<class index_iterator, class element_iterator, class Element>
		class _IndexedIterator {
		private :
			typedef  index_iterator    index_it ;
			typedef  element_iterator  data_it ;
			index_it _start_it ;
			index_it _start_beg ;
			index_it _rowid_it ;
			data_it _data_it ;
			const data_it _data_beg ;
			const data_it _data_end ;
		public:
			typedef Element value_type ;
			_IndexedIterator( const index_it &i
					  , const index_it &j
					  , const data_it &e
					  , const data_it &e_e) :
				_start_it(i)
				, _start_beg(i)
				, _rowid_it(j)
				, _data_it(e)
				, _data_beg(e)
				, _data_end(e_e)
			{}

			_IndexedIterator (const _IndexedIterator &iter) :
				_start_it(iter._start_it)
				, _start_beg(iter._start_beg)
				, _rowid_it(iter._rowid_it)
				, _data_it(iter._data_it)
				, _data_beg(iter._data_beg)
				, _data_end(iter._data_end)
			{}

			_IndexedIterator &operator = (const _IndexedIterator &iter)
			{
				_start_it = iter._start_it ;
				_start_beg = iter._start_beg ;
				_rowid_it = iter._rowid_it ;
				_data_it  = iter._data_it  ;
				_data_beg = iter._data_beg ;
				_data_end  = iter._data_end  ;

				return *this;
			}

			bool operator == (const _IndexedIterator &i) const
			{
				// we assume consistency
				return  (_data_it == i._data_it);
			}

			bool operator != (const _IndexedIterator &i) const
			{
				// we assume consistency
				return  (_data_it != i._data_it) ;
			}

			_IndexedIterator &operator ++ ()
			{

				++_data_it  ;
				if (_data_it == _data_end) {
					return *this ;
				}
				++_rowid_it ;

				while (std::distance(_data_beg, _data_it) >= *(_start_it+1)) {
					++_start_it ;
				}
				return *this;
			}

			_IndexedIterator operator ++ (int)
			{
				_IndexedIterator tmp = *this;
				++(*this);
				return tmp;
			}

			_IndexedIterator &operator -- ()
			{
				throw NotImplementedYet("not sure");
				--_data_it  ;
				if (_data_it == _data_beg) {
					return *this ;
				}
				--_rowid_it ;
				while (std::distance(_data_beg, _data_it) > *(_start_it)) {
					--_start_it ;
				}

				return *this;
			}

			_IndexedIterator operator -- (int)
			{
				_IndexedIterator tmp = *this;
				--(*this);
				return tmp;
			}

			value_type &operator * ()
			{
				return *_data_it;
			}

			value_type *operator -> ()
			{
				return _data_it ;
			}

			const value_type &operator*() const
			{
				return *_data_it;
			}

			const value_type *operator -> () const
			{
				return _data_it ;
			}

			size_t rowIndex () const
			{
				return std::distance(_start_beg, _start_it);
			}

			size_t colIndex () const
			{
				return *_rowid_it;
			}

			const value_type &value() const
			{
				return *_data_it;
			}


		};

		typedef _Iterator<typename std::vector<Element>::iterator, Element> Iterator;
		typedef _Iterator<typename std::vector<Element>::const_iterator, constElement> ConstIterator;

		typedef _IndexedIterator<std::vector<size_t>::iterator, typename std::vector<Element>::iterator, Element> IndexedIterator;
		typedef _IndexedIterator<std::vector<size_t>::const_iterator, typename std::vector<Element>::const_iterator, constElement> ConstIndexedIterator;


		Iterator      Begin ()
		{
			return Iterator(_data.begin(),_data.end()) ;
		}

		Iterator      End   ()
		{
			return Iterator(_data.end(), _data.end()) ;
		}

		ConstIterator      Begin () const
		{
			return ConstIterator(_data.begin(),_data.end()) ;
		}

		ConstIterator      End   () const
		{
			return ConstIterator(_data.end(), _data.end()) ;
		}

		IndexedIterator      IndexedBegin ()
		{
			return IndexedIterator(_start.begin(), _rowid.begin(), _data.begin(),_data.end()) ;
		}

		IndexedIterator      IndexedEnd   ()
		{
			return IndexedIterator(_start.end(), _rowid.end(), _data.end(),_data.end()) ;
		}

		ConstIndexedIterator      IndexedBegin () const
		{
			return ConstIndexedIterator(_start.begin(), _rowid.begin(), _data.begin(),_data.end()) ;
		}

		ConstIndexedIterator      IndexedEnd   () const
		{
			return ConstIndexedIterator(_start.end(), _rowid.end(), _data.end(),_data.end()) ;
		}
#endif


	protected :
		friend class SparseMatrixWriteHelper<Self_t >;
		friend class SparseMatrixReadHelper<Self_t >;


		size_t              _rownb ;
		size_t              _colnb ;
		size_t               _nbnz ;

		std::vector<size_t> _start ;
		std::vector<size_t> _rowid ;
		std::vector<Element> _data ;

		const _Field & _field;

		mutable Helper _helper ;

		mutable struct _triples {
			// XXX not done
			ptrdiff_t _row ;
			ptrdiff_t _nnz ;
			_triples() :
				_row(-1)
				, _nnz(-1)
			{}

			ptrdiff_t next( const std::vector<size_t> & start)
			{
				_nnz +=1 ;
				while (_row+1 < start.size() && _nnz >= start[_row+1]) {
					_row += 1;
				}
				return _nnz ;
			}

			void reset()
			{
				_row = -1 ;
				_nnz = -1 ;
			}
		}_triples;
	};





} // namespace LinBox

#endif // __LINBOX_sparse_matrix_sparse_csc_matrix_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
