/* linbox/matrix/sparse-csr-matrix.h
 * Copyright (C) 2013 the LinBox
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

/*! @file matrix/SparseMatrix/sparse-ell-matrix.h
 * @ingroup sparsematrix
 * @brief
 */


#ifndef __LINBOX_matrix_sparsematrix_sparse_ell_matrix_H
#define __LINBOX_matrix_sparsematrix_sparse_ell_matrix_H

#include <utility>
#include <iostream>
#include <algorithm>

#include "linbox/linbox-config.h"
#include "linbox/util/debug.h"
#include "linbox/field/hom.h"
#include "sparse-domain.h"

#ifndef LINBOX_ELL_TRANSPOSE
#define LINBOX_ELL_TRANSPOSE 1000
#endif

namespace LinBox
{


	/** Sparse matrix, Coordinate storage.
	 *
	 * \ingroup matrix
	 * \ingroup sparse
	 */
	template<class _Field>
	class SparseMatrix<_Field, SparseMatrixFormat::ELL > {
	public :
		typedef _Field                             Field ; //!< Field
		typedef typename _Field::Element         Element ; //!< Element
		typedef const Element               constElement ; //!< const Element
		typedef SparseMatrixFormat::ELL         Storage ; //!< Matrix Storage Format
		typedef SparseMatrix<_Field,Storage>     Self_t ; //!< Self type
		typedef typename Vector<Field>::SparseSeq    Row ; //!< @warning this is not the row type. Just used for streams.
		// typedef Vector<_Field,VectorStorage::Sparse> Rep ;

		/*! Constructors.
		 * @todo convert from other matrix.
		 *
		 */
		//@{
#if 0 /*  No empty CSTOR */
		SparseMatrix<_Field, SparseMatrixFormat::ELL> () :
			_rownb(0),_colnb(0),_maxc(0)
			,_nbnz(0)
			,_colid(0),_data(0)
			, _field()
			, _helper()
		{
		}
#endif

		SparseMatrix<_Field, SparseMatrixFormat::ELL> (const _Field & F) :
			_rownb(0),_colnb(0)
			,_maxc(0)
			,_nbnz(0)
			,_colid(0)
			,_data(0)
			, _field(F)
			, _helper()
		{
		}

		SparseMatrix<_Field, SparseMatrixFormat::ELL> (const _Field & F, size_t m, size_t n) :
			_rownb(m),_colnb(n)
			,_maxc(0)
			,_nbnz(0)
			,_colid(0)
			,_data(0)
			, _field(F)
			, _helper()
		{
		}

		SparseMatrix<_Field, SparseMatrixFormat::ELL> (const _Field & F,
							       size_t m, size_t n,
							       size_t z) :
			_rownb(m),_colnb(n)
			, _maxc(0)
			, _nbnz(z)
			, _colid(z)
			,_data(z)
			, _field(F)
			, _helper()
		{
		}

		SparseMatrix<_Field, SparseMatrixFormat::ELL> (const SparseMatrix<_Field, SparseMatrixFormat::CSR> & S) :
			_rownb(S._rownb),_colnb(S._colnb)
			,_maxc(S._maxc)
			,_nbnz(S._nbnz)
			, _colid(S._colid)
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
			_rowid(S._rowid),_colid(S._colid),_data(S._data),
			_field(S._field)
		{}
#endif

		template<typename _Tp1, typename _Rw1 = SparseMatrixFormat::ELL>
		struct rebind {
			typedef SparseMatrix<_Tp1, _Rw1> other;
		private:

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
						Ap.appendEntry(i,j,e);
				}
				A.firstTriple();
				Ap.finalize();
			}

			void rebindMethod(SparseMatrix<_Tp1, Storage>  & Ap, const Self_t & A /*,  IndexedCategory::HasNext*/)
			{
				// we don't use nextTriple because we can do better.
				linbox_check(A.consistent());
				// Ap = new other(F, A.rowdim(), A.coldim());
				Ap.resize(A.rowdim(),A.coldim(),A.size(),A.ld());

				typename _Tp1::Element e;

				Hom<typename Self_t::Field, _Tp1> hom(A.field(), Ap.field());
				size_t newz = 0 ;
				for (size_t i = 0 ; i < A.rowdim() ; ++i) {
					size_t j = 0 ;
					for (size_t k = 0 ; k < A.ld() ; ++k) {
						if (A.field().isZero(A.getData(i,k))) {
							break;
						}
						hom. image ( e, A.getData(i,k) );
						if (!Ap.field().isZero(e)) {
							Ap.setColid(i,j,A.getColid(i,k));
							Ap.setData(i,j,e);
							++j;
						}
						else {
							++newz ;
						}
					}
				}
				Ap.setSize(Ap.size() - newz) ;
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
			,_maxc(0)
			,_nbnz(S.size())
			, _colid(0)
			,_data(0)
			, _field(F)
		{
			typename SparseMatrix<_Tp1,_Rw1>::template rebind<Field,Storage>()(*this, S);
			finalize();
		}




		template<class VectStream>
		SparseMatrix<_Field, SparseMatrixFormat::ELL> (const _Field & F, VectStream & stream) :
			_rownb(stream.size()),_colnb(stream.dim())
			, _maxc(0)
			, _nbnz(0)
			, _colid(0)
			,_data(0)
			, _field(F)
		{
			//! @todo
			SparseMatrix<_Field,SparseMatrixFormat::CSR> Tmp(F,stream);
			importe(Tmp);
			finalize();
		}

		SparseMatrix<_Field, SparseMatrixFormat::ELL> ( MatrixStream<Field>& ms ):
			_rownb(0),_colnb(0)
			,_maxc(0)
			,_nbnz(0)
			,_colid(0)
			,_data(0)
			,_field(ms.field())
		{
			firstTriple();

			Element val;
			size_t i, j;
			while( ms.nextTriple(i,j,val) ) {
				if (! field().isZero(val)) {
					if( i >= _rownb ) {
						// _rownb = i + 1;
						resize(i+1,_colnb,_nbnz,(size_t)_maxc);
					}
					if( j >= _colnb ) {
						// _colnb = j + 1;
						resize(_rownb,j+1,_nbnz,(size_t)_maxc);
					}
					appendEntry(i,j,val);
				}
			}
			if( ms.getError() > END_OF_MATRIX )
				throw ms.reportError(__func__,__LINE__);
			if( !ms.getDimensions( i, j ) )
				throw ms.reportError(__func__,__LINE__);
#ifndef NDEBUG
			if( i != _rownb  || j != _colnb) {
				std::cout << " ***Warning*** the sizes got changed" << __func__ << ',' << __LINE__ << std::endl;
			}
#endif

			firstTriple();
			finalize();
			linbox_check(consistent());
		}


		void resize(const size_t   mm, const size_t   nn, const size_t   zz = 0, const size_t   ll = 0)
		{
			// linbox_check(_rownb*_maxc == _colid.size());
			// attention RowMajor/ColMajor
			if (!_maxc || mm == _maxc) {
				_colid.resize(mm*ll,0);
				_data .resize(mm*ll,field().zero);
			}
			else if ( ll > _maxc) {
				reshape(ll);
			}
			else {
#ifndef NDEBUG
				if (_maxc && ll != _maxc)
					std::cout << " ***Warning*** possibly loosing data in ELL resize " << std::endl;
#endif
				_colid.resize(mm*ll,0);
				_data .resize(mm*ll,field().zero);
			}
			_rownb = mm ;
			_colnb = nn ;
			_nbnz  = zz;
			_maxc  = ll;

			linbox_check(_rownb*_maxc == _colid.size());
		}

		/*! Default converter.
		 * @param S a sparse matrix in any storage.
		 */
		template<class _OtherStorage>
		SparseMatrix<_Field, SparseMatrixFormat::ELL> (const SparseMatrix<_Field, _OtherStorage> & S) :
			_rownb(S.rowdim()),_colnb(S.coldim())
			,_nbnz(S.size())
			,_maxc(0)
			,_colid(0)
			,_data(0)
			,_field(S.field())
		{
			// _nbnz is set there:
			this->importe(S); // convert Temp from anything
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
		void importe(const SparseMatrix<_Field,SparseMatrixFormat::CSR> &S)
		{
			// can be sped up on multicores.
			for (size_t i = 0 ; i < S.rowdim() ; ++i)
				_maxc = std::max(_maxc, (size_t)(S.getEnd(i)-S.getStart(i)));

			resize(S.rowdim(), S.coldim(), S.size(),(size_t)_maxc);

			for (size_t i = 0 ; i < S.rowdim() ; ++i) {
				size_t k = 0 ;
				for (size_t j = S.getStart(i) ; j < S.getEnd(i) ; ++j, ++k) {
					setColid(i,k,S.getColid(j));
					setData(i,k,S.getData(j));
				}
			}

		}

		/*! Import a matrix in CSR format to CSR.
		 * @param S CSR matrix to be converted in CSR
		 */

		void importe(const SparseMatrix<_Field,SparseMatrixFormat::ELL> &S)
		{
			resize( S.rowdim(), S.coldim(), S.size() , S.ld());

			setColid(S.getColid());
			setData(S.getData());

		}

		template<class _OtherStorage>
		void importe(const SparseMatrix<_Field,_OtherStorage> &S)
		{
			SparseMatrix<_Field,SparseMatrixFormat::CSR> Tmp(S);
			this->importe(S);

		}

		/*! Export a matrix in CSR format from COO.
		 * @param S CSR matrix to be converted from COO
		 */
		SparseMatrix<_Field,SparseMatrixFormat::CSR > &
		exporte(SparseMatrix<_Field,SparseMatrixFormat::CSR> &S) const
		{
			linbox_check(consistent());

			S.resize(_rownb, _colnb, _nbnz);
			S.setStart(0,0);
			size_t k = 0 ;
			for (size_t i = 0 ; i < S.rowdim() ; ++i) {
				for (size_t j = 0  ; j < _maxc; ++j ) {
					if (field().isZero(getData(i,j)))
						break;
					S.setColid(k,getColid(i,j));
					S.setData(k,getData(i,j));
					++k;
				}
				S.setStart(i+1,k);
			}
			linbox_check(k == _nbnz);

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
			// outStart
			size_t maxc = 0;
			{
				std::vector<size_t> start(coldim(),0);

				for (size_t i = 0 ; i < rowdim() ; ++i) {
					for (size_t j = 0 ; j < _maxc ; ++j) {
						if (!field().isZero(getData(i,j)))
							start[getColid(i,j)] += 1 ;
						else
							break;
					}
				}

				maxc = start[0] ;
				for (size_t i = 1 ; i < coldim() ; ++i)
					maxc= std::max(maxc, start[i]);
			}



			S.resize((size_t)_colnb, (size_t)_rownb, (size_t)_nbnz, (size_t)maxc ); // necessary copy to temp, no const ref

			{
				std::vector<size_t> done_col(coldim(),0);


				for (size_t i = 0 ; i < rowdim() ; ++i) {
					for (size_t j = 0 ; j < _maxc ; ++j) {
						if (!field().isZero(getData(i,j))) {
							size_t loc_col = getColid(i,j) ;
							S.setData (loc_col,done_col[loc_col],getData(i,j));
							S.setColid(loc_col,done_col[loc_col],i);

							linbox_check(done_col[loc_col] < maxc);
							done_col[loc_col] += 1 ;

						}
						else
							break;
					}
				}
			}

			S.finalize();


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

		void setSize(const size_t  z)
		{
			_nbnz = z ;
		}


		/** Get a read-only individual entry from the matrix.
		 * @param i Row _colid
		 * @param j Column _colid
		 * @return Const reference to matrix entry
		 */
		constElement &getEntry(const size_t &i, const size_t &j) const
		{
			ptrdiff_t off = _triples.next(_maxc);
			size_t ii = i ;
			while ( field().isZero(_data[ii*_maxc+off]) )  {
				_triples._off = off = 0;
				_triples._row = ii  = ii+1 ;
			}
			if (  i == ii &&  _colid[ii*_maxc+off]  == j  ) { /* sort of nextTriple */
				linbox_check(!field().isZero(_data[ii*_maxc+off]));
				return _data[ii*_maxc+off];
			}
			else { /* searching */

				linbox_check(consistent());

				linbox_check(i<_rownb);
				linbox_check(j<_colnb);

				const size_t * beg = &_colid[i*_maxc];
				const Element * dat = &_data[i*_maxc];
				for (size_t k = 0 ; k < _maxc ; ++k) {
					if ( field().isZero(dat[k])) {
						return field().zero;
					}
					// replace
					if (beg[k] == j) {
						_triples._off = k;
						_triples._row = i;
						return dat[k];
					}
					if (beg[k] > j) {
						return field().zero;
					}
				}
				return field().zero;
			}

		}

		Element      &getEntry (Element &x, size_t i, size_t j) const
		{
			return x = getEntry (i, j);
		}

		void appendEntry(const size_t &i, const size_t &j, const Element& e)
		{
			linbox_check(i < rowdim());
			linbox_check(j < coldim());
			linbox_check(_colid.size() == _data.size());
			linbox_check(_rownb*_maxc == _colid.size());

			if (field().isZero(e)) {
				return ;
			}
			ptrdiff_t row = _triples._row ;
			ptrdiff_t off = _triples._off ;
			if (row != (ptrdiff_t)i) { /* new row */
				linbox_check((ptrdiff_t)i>row);
				_triples._row = i ;
				_triples._off = 0 ;
				if (_maxc == 0) {
					insert(i,_maxc,j,e);
				}
				else {
					linbox_check(i*_maxc < _colid.size());
					_colid[i*_maxc] = j ;
					_data [i*_maxc] = e ;
					++_nbnz;
				}
			}
			else { /* same row */
				linbox_check(_triples._row == i);
				_triples._off = off = off + 1 ;
				if (off == (ptrdiff_t)_maxc) {
					insert(i,_maxc,j,e);
				}
				else {
					linbox_check(i*_maxc+off < _colid.size());
					_colid[i*_maxc+off] = j ;
					_data [i*_maxc+off] = e ;
					++_nbnz;
				}

			}
			linbox_check(_rownb*_maxc == _colid.size());
		}

		/// make matrix ready to use after a sequence of setEntry calls.
		void finalize(){
			// could check that maxc is not too large and shrink ? Is is optimize job ?
			_triples.reset();
		} // end construction after a sequence of setEntry calls.

		/** Set an individual entry.
		 * Setting the entry to 0 will not remove it from the matrix
		 * @param i Row _colid of entry
		 * @param j Column _colid of entry
		 * @param value Value of the new entry
		 * @todo make it faster if i is 0 or m-1 ?
		 * @warning if this is used to build a matrix and this matrix is "well formed",
		 * it can be sped up (no checking that the entry already exists).
		 */
		void setEntry(const size_t &i, const size_t &j, const Element& e
			     )
		{
			linbox_check(i<_rownb);
			linbox_check(j<_colnb);

			linbox_check(consistent());

			if (field().isZero(e)) {
				return clearEntry(i,j);
			}

			size_t * beg = &_colid[i*_maxc];
			Element * dat = &_data[i*_maxc];
			bool found = false;
			for (size_t k = 0 ; k < _maxc ; ++k) {
				if ( field().isZero(dat[k])) {
					field().assign(dat[k], e) ;
					beg[k] = j;
					found = true;
					++_nbnz ;
					break;
				}
				if (beg[k] == j) {
					if (field().isZero(dat[k])) {
						++_nbnz ;
					}
					field().assign(dat[k], e) ;
					beg[k] = j;
					found = true;
					break;
				}
				if (beg[k] > j) {
					found = true;
					insert(i,k,j,e);
					break;
				}
			}
			if (!found) {
				insert(i,_maxc,j,e);
			}
		}

#if 0
		/** Get a writeable reference to an entry in the matrix.
		 * If there is no entry at the position (i, j), then a new entry
		 * with a value of zero is inserted and a reference  to it is
		 * returned.
		 * @param i Row _colid of entry
		 * @param j Column _colid of entry
		 * @return Reference to matrix entry
		 */
		Element &refEntry(const size_t &i, const size_t&j)
		{
			linbox_check(i<_rownb);
			linbox_check(j<_colnb);
			// Could be improved by adding an initial guess j/rowdim*size()

			size_t ibeg = _start[i];
			size_t iend = _start[i+1];
			if (ibeg==iend) {
				for (size_t k = i+1 ; k <= _rownb ; ++k) _start[k] +=1 ;
				_colid.insert(_colid.begin()+ibeg,j);
				_data.insert( _data.begin() +ibeg,field().zero);
				return _data[ibeg];
			}
			typedef typename std::vector<size_t>::iterator myIterator ;
			myIterator beg = _colid.begin() ;
			myIterator low = std::lower_bound (beg+(ptrdiff_t)ibeg, beg+(ptrdiff_t)iend, j);
			if (low == beg+(ptrdiff_t)iend) {
				for (size_t k = i+1 ; k <= _rownb ; ++k) _start[k] +=1 ;
				_colid.insert(_colid.begin()+ibeg,j);
				_data.insert( _data.begin() +ibeg,field().zero);
				return _data[ibeg];
			}
			else {
				size_t la = low-_colid.begin() ;
				return _data[la] ;
			}
		}
#endif

		/** Write a matrix to the given output stream using field read/write.
		 * @param os Output stream to which to write the matrix
		 * @param format Format with which to write
		 */
		std::ostream & write(std::ostream &os
				     , LINBOX_enum(Tag::FileFormat) format = Tag::FileFormat::MatrixMarket) const
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
		 * @param i row _colid
		 * @param j col _colid
		 * @bug rows are increasing order => faster search.
		 */
		void clearEntry(const size_t &i, const size_t &j)
		{
			linbox_check(i<_rownb);
			linbox_check(j<_colnb);

			size_t k = 0 ;
			for ( ; k < _maxc ; ++k) {
				if (_colid[i*_maxc+k] == j)
					break;
				if (_colid[i*_maxc+k] > j)
					return;
			}
			if (k == _maxc)
				return; // not found

			if (field().isZero(_data[i*_maxc+k]) && _colid[i*_maxc+k] == 0)
				return;

			for (size_t l = k ; l < _maxc-1 ; ++l) {
				field().assign(_data[i*_maxc+l],_data[i*_maxc+l+1]);
				_colid[i*_maxc+l] = _colid[i*_maxc+l+1];
			}
			field().assign(_data[i*_maxc+_maxc-1],field().zero) ;
			_colid[i*_maxc+_maxc-1] = 0;
			--_nbnz;
		}

		/*! @internal
		 * @brief cleans 0 entries.
		 */
		void clean()
		{
#if 0
			size_t i = 0 ;
			while(i < _data.size()) {
				if ( field().isZero(_data[i]) ) {
					for (size_t k = i+1 ; k <= _rownb ; ++k) _start[k] -= 1 ;
					_colid.erase(_colid.begin()+i);
					_data. erase(_data. begin()+i);
				}
				else
					++i ;
			}
			return ;
#endif
		}

		// y= Ax
		// y[i] = sum(A(i,j) x(j)
		template<class Vector>
		Vector& apply(Vector &y, const Vector& x, const Element & a ) const
		{
			// linbox_check(consistent());
			prepare(field(),y,a);


			FieldAXPY<Field> accu(field());
			for (size_t i = 0 ; i < _rownb ; ++i) {
				accu.reset();
				for (size_t k = 0   ; k < _maxc ; ++k)
					if (!field().isZero(getData(i,k)))
						// field().axpyin( y[i], getData(i,k), x[getColid(i,k)] ); //! @todo delay !!!
						accu.mulacc( getData(i,k), x[getColid(i,k)] );
					else {
						break;
					}
				accu.get(y[i]);
			}

			return y;
		}


		// y= A^t x
		// y[i] = sum(A(j,i) x(j)
		template<class Vector>
		Vector& applyTranspose(Vector &y, const Vector& x, const Element & a ) const
		{
			linbox_check(consistent());
			if (_helper.optimized(*this)) {
				return _helper.matrix().apply(y,x,a) ; // NEVER use applyTranspose on that thing.
			}

			prepare(field(),y,a);


#if 0
			for (size_t i = 0 ; i < _rownb ; ++i)
				for (size_t k = 0   ; k < _maxc ; ++k)
					if (!field().isZero(getData(i,k)))
						field().axpyin(y[getColid(i,k)], getData(i,k), x[i] ); //! @todo delay !!!
					else
						break;
#else
			const FieldAXPY<Field> accu0(field());
			std::vector<FieldAXPY<Field> > Y(_colnb, accu0);


			for (size_t i = 0 ; i < _rownb ; ++i)
				for (size_t k = 0   ; k < _maxc ; ++k)
					if (!field().isZero(getData(i,k)))
						Y[getColid(i,k)].mulacc( getData(i,k), x[i] );
					else
						break;

			for (size_t i = 0 ; i < _colnb ; ++i)
				Y[i].get(y[i]) ;
#endif

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

		/** @todo Non element marker.
		 * We could end up a line with a marker.
		 * A field F would contain an element that does not belong to
		 * it. eg a nan for a Modular<double>. It could act as a
		 * marker.
		 */
		bool consistent() const
		{
			size_t nbnz = 0 ;
			if (_rownb * _maxc != _data.size()) {
				std::cout << "bad data size" << std::endl;
				return false;
			}
			if (_data.size() != _colid.size()) {
				std::cout << "bad colid size" << std::endl;
				return false;
			}
			// bool ok = true ;
			for (size_t i = 0 ; i < _rownb ; ++i) {
				// bool row_ok = true ;
				bool zero = false ;
				for (size_t j = 0 ; j < _maxc ; ++j) {
					if (field().isZero(getData(i,j))) {
						zero = true ;
					}
					else {
						if (zero) {
							std::cout << "non zero after a 0" << std::endl;
							return false ; // elements after a 0...
						}
						++nbnz ;

					}
				}
			}
			if (nbnz != _nbnz)
				std::cout << "bad number of non zero" << std::endl;
			return (nbnz == _nbnz);

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
				if ( A.size() > LINBOX_ELL_TRANSPOSE ) { // and/or A.rowDensity(), A.coldim(),...
					// std::cout << "optimizing..." ;
					_optimized = true ;
					_AT = new Self_t(A.field(),A.coldim(),A.rowdim());
					A.transpose(*_AT);
					// std::cout << "done!" << std::endl;
				}
			}

			const Self_t & matrix() const
			{
				return *_AT ;
			}

		};



		void reshape(const size_t ll)
		{
			linbox_check(_rownb*_maxc == _colid.size());
			linbox_check(ll > _maxc);

			_data .resize(ll*_rownb);
			_colid.resize(ll*_rownb);
			for (size_t i = _rownb ; i>0 ; --i ) {
				// std::copy(_data.begin()+(i-1)*_maxc, _data.begin()+(i)*_maxc, _data.begin()+(i-1)*ll);
				std::copy_backward(_data.begin()+(i-1)*_maxc, _data.begin()+(i)*_maxc, _data.begin()+(i-1)*ll+_maxc);
				for (size_t j = (i-1)*ll+_maxc ; j < (i)*ll ; ++j )
					_data[j] = field().zero ;
			}
			for (size_t i = _rownb ; i>0 ; --i ) {
				// std::copy(_colid.begin()+(i-1)*_maxc, _colid.begin()+(i)*_maxc, _colid.begin()+(i-1)*ll);
				std::copy_backward(_colid.begin()+(i-1)*_maxc, _colid.begin()+(i)*_maxc, _colid.begin()+(i-1)*ll+_maxc);
				for (size_t j = (i-1)*ll+_maxc ; j < (i)*ll ; ++j )
					_colid[j] = 0 ;
			}

			_maxc = ll ;
			linbox_check(_rownb*_maxc == _colid.size());
		}

	public:
		// pseudo iterators

		size_t getColid(const size_t & i, const size_t j) const
		{
			linbox_check(i < _rownb);
			linbox_check(_maxc*_rownb == _colid.size());
			return _colid[i*_maxc+j];
		}

		void setColid(const size_t & i, const size_t & j, const size_t & k)
		{
			linbox_check(i < _rownb);
			linbox_check(i*_maxc+j < _colid.size());
			linbox_check(_maxc*_rownb == _colid.size());

			_colid[i*_maxc+j]=k;
		}

		void setColid(std::vector<size_t> new_colid)
		{
			_colid = new_colid ;
		}

		std::vector<size_t>  getColid( ) const
		{
			return _colid ;
		}

		const Element & getData(const size_t & i, const size_t & j) const
		{
			linbox_check(i < _rownb);
			linbox_check(j < _colnb);
			linbox_check(_maxc*_rownb == _colid.size());
			return _data[ i*_maxc+j ];
		}

		void setData(const size_t & i, const size_t & j, const Element & e)
		{
			linbox_check(i*_maxc+j < _data.size());
			linbox_check(_maxc*_rownb == _colid.size());
			field().assign(_data[i*_maxc+j],e);
		}

		void setData(const std::vector<Element> & new_data)
		{
			_data = new_data ;
		}

		std::vector<Element>  getData( ) const
		{
			return _data ;
		}

		size_t ld() const
		{
			return _maxc;
		}

		void firstTriple() const
		{
			_triples.reset();
		}

		// a zero means end of line
		bool nextTriple(size_t & i, size_t &j, Element &e) const
		{
			ptrdiff_t off =_triples.next(_maxc);

			i = _triples._row ;

			if (  i >= rowdim() ) {
				_triples.reset() ;
				return false ;
			}

			while ( field().isZero(_data[i*_maxc+off]) )  {
				_triples._off = off = 0;
				_triples._row = i   = i+1 ;
				if (i >= rowdim()) {
					_triples.reset();
					return false;
				}
			}

			j = _colid[i*_maxc+off];
			e = _data[i*_maxc+off];

			return true;
		}

		template<class element_iterator, class Field>
		class _Iterator {
		private :
			element_iterator _data_it ;
			const element_iterator _data_beg ;
			const element_iterator _data_end ;
			const Field & _field ;
			typedef typename Field::Element Element;
		public:
			typedef Element value_type ;
			_Iterator(const Field & F , const element_iterator & e_beg, const element_iterator & e_end) :
				_data_it(e_beg)
				, _data_beg(e_beg)
				, _data_end(e_end)
				,_field(F)
			{}

			_Iterator (const _Iterator &iter) :
				_data_it(iter._data_it)
				, _data_beg(iter._data_beg)
				, _data_end(iter._data_end)
				,_field(iter._field)

			{}

			_Iterator &operator = (const _Iterator &iter)
			{
				_data_it  = iter._data_it  ;
				_data_beg  = iter._data_beg  ;
				_data_end  = iter._data_end  ;
				_field  = iter._field  ;

				return *this;
			}

			bool operator == (const _Iterator &i) const
			{
				return  (_data_it == i._data_it) ;
			}

			bool operator != (const _Iterator &i) const
			{
				return  (_data_it != i._data_it) ;
			}

			_Iterator &operator ++ ()
			{
				do {
					++_data_it ;
				} while (_data_it  != _data_end && _field.isZero(*_data_it));
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
				throw NotImplementedYet("not sure");
				do {
					--_data_it ;
				} while (_data_it  != _data_beg && _field.isZero(*_data_it));
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

		template<class index_iterator, class element_iterator, class Field>
		class _IndexedIterator {
		private :
			typedef  index_iterator    index_it ;
			typedef  element_iterator  data_it ;
			index_it _colid_beg ;
			index_it _colid_it ;
			data_it _data_it ;
			const data_it _data_beg ;
			const data_it _data_end ;
			const Field & _field ;
			const size_t & _ld ;
			size_t  _row ;
			typedef typename Field::Element Element;
		public:
			typedef Element value_type ;
			_IndexedIterator( const Field & F
					  , const size_t & ld
					  , const index_it &j
					  , const data_it &e
					  , const data_it &e_e) :
				_colid_beg(j)
				, _colid_it(j)
				, _data_it(e)
				, _data_beg(e)
				, _data_end(e_e)
				, _field(F)
				, _ld(ld)
				, _row(0)
			{}

			_IndexedIterator (const _IndexedIterator &iter) :
				_colid_beg(iter._colid_beg)
				, _colid_it(iter._colid_it)
				, _data_it(iter._data_it)
				, _data_beg(iter._data_beg)
				, _data_end(iter._data_end)
				, _field(iter._field)
				, _ld(iter._ld)
				, _row(iter._row)
			{}

			_IndexedIterator &operator = (const _IndexedIterator &iter)
			{
				_colid_beg = iter._colid_beg ;
				_colid_it  = iter._colid_it ;
				_data_it   = iter._data_it  ;
				_data_beg  = iter._data_beg ;
				_data_end  = iter._data_end  ;
				_field     = iter._field ;
				_ld        = iter._ld ;
				_row       = iter._row ;

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
				if (std::distance(_data_beg,_data_it) % _ld == 0)
					++_row ;

				++_colid_it ;
				while (_field.isZero(*_data_it)) {
					++_row ;
					_data_it = _data_beg + _row * _ld ;
					_colid_it = _colid_beg + _row * _ld ;
					if (_data_it == _data_end) {
						return *this;
					}
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
				if (std::distance(_data_beg,_data_it) % _ld == 0)
					--_row ;

				--_colid_it ;
				while (_field.isZero(*_data_it)) {
					--_row ;
					_data_it = _data_beg + _row * _ld ;
					_colid_it = _colid_beg + _row * _ld ;
					if (_data_it == _data_beg) {
						return *this;
					}
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
				return _row;
			}

			size_t colIndex () const
			{
				return *_colid_it;
			}

			const value_type &value() const
			{
				return *_data_it;
			}


		};

		typedef _Iterator<typename std::vector<Element>::iterator, Element> Iterator;
		typedef _Iterator<typename std::vector<Element>::const_iterator, constElement> ConstIterator;

		typedef _IndexedIterator<std::vector<size_t>::iterator, typename std::vector<Element>::iterator, Field> IndexedIterator;
		typedef _IndexedIterator<std::vector<size_t>::const_iterator, typename std::vector<Element>::const_iterator, const Field> ConstIndexedIterator;


		Iterator      Begin ()
		{
			return Iterator(field(),_data.begin(),_data.end()) ;
		}

		Iterator      End   ()
		{
			return Iterator(field(),_data.end(), _data.end()) ;
		}

		ConstIterator      Begin () const
		{
			return ConstIterator(field(),_data.begin(),_data.end()) ;
		}

		ConstIterator      End   () const
		{
			return ConstIterator(field(),_data.end(), _data.end()) ;
		}

		IndexedIterator      IndexedBegin ()
		{
			return IndexedIterator(field(), _maxc , _colid.begin(), _data.begin(),_data.end()) ;
		}

		IndexedIterator      IndexedEnd   ()
		{
			return IndexedIterator(field(), _maxc, _colid.end(), _data.end(),_data.end()) ;
		}

		ConstIndexedIterator      IndexedBegin () const
		{
			return ConstIndexedIterator(field(), _maxc, _colid.begin(), _data.begin(),_data.end()) ;
		}

		ConstIndexedIterator      IndexedEnd   () const
		{
			return ConstIndexedIterator(field(), _maxc, _colid.end(), _data.end(),_data.end()) ;
		}



	private:

		void insert (const size_t  i, const size_t  k, const size_t  j, const Element  e)
		{
			linbox_check(_rownb*_maxc == _colid.size());
			if (k == _maxc) {
				resize(_rownb,_colnb,_nbnz,_maxc+1);
				linbox_check(_rownb*_maxc == _colid.size());
				_colid[_maxc*i+k] = j;
				field().assign(_data[_maxc*i+k],e);
				++_nbnz;
				return;
			}
			size_t l = k;
			for ( ; l < _maxc ; ++l)
				if (field().isZero(_data[i*_maxc+l]))
					break;
			if (l == _maxc)
				resize(_rownb,_colnb,_nbnz,_maxc+1);

			for (size_t u = l ; u > k ; --u) {
				_colid[_maxc*i+u] = _colid[_maxc*i+u-1] ;
				field().assign(_data[_maxc*i+u],_data[_maxc*i+u-1]);
			}
			_colid[_maxc*i+k] = j;
			field().assign(_data[_maxc*i+k],e);
			++_nbnz;
			linbox_check(_rownb*_maxc == _colid.size());
			return;
		}

		void write_raw() const
		{
			std::cout << "colids" << std::endl;
			for (size_t i = 0 ; i < _rownb ; ++i){
				for (size_t j = 0 ; j < _maxc ; ++j)
					std::cout << _colid[i*_maxc+j] << ' ' ;
				std::cout << std::endl;
			}

			std::cout << "data" << std::endl;
			for (size_t i = 0 ; i < _rownb ; ++i){
				for (size_t j = 0 ; j < _maxc ; ++j)
					std::cout << _data[i*_maxc+j] << ' ' ;
				std::cout << std::endl;
			}
		}

	protected :
		friend class SparseMatrixWriteHelper<Self_t >;
		friend class SparseMatrixReadHelper<Self_t >;


		size_t              _rownb ;
		size_t              _colnb ;
		size_t               _maxc ; //!< longest row
		size_t               _nbnz ;

		std::vector<size_t> _colid ; //!< \p _colid is \p _rownb x \p _maxc in RowMajor
		std::vector<Element> _data ; //!< \p _data  is \p _rownb x \p _maxc in RowMajor

		const _Field            & _field;

		mutable Helper _helper ;

		mutable struct _triples {
			ptrdiff_t _row ;
			ptrdiff_t _off ;
			_triples() :
				_row(-1)
				, _off(-1)
			{}

			ptrdiff_t next(size_t maxc)
			{
				++ _off ;
				if (_row == -1) {
					++_row ;
					return _off ;
				}
				if (_off >= (ptrdiff_t)maxc) {
					_row += 1 ;
					_off = 0 ;
				}
				return _off;
			}

			void reset()
			{
				_row = -1 ;
				_off = -1 ;
			}
		}_triples;
	};





} // namespace LinBox

#endif // __LINBOX_matrix_sparsematrix_sparse_ell_matrix_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
