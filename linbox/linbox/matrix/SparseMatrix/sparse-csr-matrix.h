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

/*! @file matrix/sparse-matrix.h
 * @ingroup matrix
 * @ingroup sparse
 * A <code>SparseMatrix2<_Field ></code> ....
 */


#ifndef __LINBOX_sparse_matrix_sparse_csr_matrix_H
#define __LINBOX_sparse_matrix_sparse_csr_matrix_H

#include <utility>
#include <iostream>
#include <algorithm>

#include "linbox-config.h"
#include "linbox/util/debug.h"
#include "linbox/matrix/sparse.h"
// #include "linbox/blackbox/factory.h"
// #include "linbox/vector/vector-traits.h"
// #include "linbox/matrix/matrix-domain.h"
// #include "linbox/util/matrix-stream.h"


namespace LinBox
{


	/** Sparse matrix, Coordinate storage.
	 *
	 * \ingroup matrix
	 * \ingroup sparse
	 */
	template<class _Field>
	class SparseMatrix2<_Field, SparseMatrix2Format::CSR > {
	public :
		typedef _Field                             Field ; //!< Field
		typedef typename _Field::Element         Element ; //!< Element
		typedef const Element               constElement ; //!< const Element
		typedef SparseMatrix2Format::CSR         Storage ; //!< Matrix Storage Format
		typedef SparseMatrix2<_Field,Storage>     Self_t ; //!< Self type
		// typedef Vector<_Field,VectorStorage::Sparse> Rep ;

		/*! Constructors.
		 * @todo convert from other matrix.
		 *
		 */
		//@{
		SparseMatrix2<_Field, SparseMatrix2Format::CSR> () :
			_rownb(0),_colnb(0),
			_nbnz(0),
			_start(0),_colid(0),_data(0)
			, _field()
		{}

		SparseMatrix2<_Field, SparseMatrix2Format::CSR> (const _Field & F) :
			_rownb(0),_colnb(0),
			_nbnz(0),
			_start(0),_colid(0),_data(0)
			, _field(F)
		{}

		SparseMatrix2<_Field, SparseMatrix2Format::CSR> (const _Field & F, size_t m, size_t n) :
			_rownb(m),_colnb(n),
			_nbnz(0),
			_start(0),_colid(0),_data(0)
			, _field(F)
		{}

		SparseMatrix2<_Field, SparseMatrix2Format::CSR> (const _Field & F,
							       size_t m, size_t n,
							       size_t z) :
			_rownb(m),_colnb(n),
			_nbnz(z),
			_start(m+1),_colid(z),_data(z),
			_field(F)
		{}

		SparseMatrix2<_Field, SparseMatrix2Format::CSR> (const SparseMatrix2<_Field, SparseMatrix2Format::CSR> & S) :
			_rownb(S._rownb),_colnb(S._colnb),
			_nbnz(S._nbnz),
			_start(S._start),_colid(S._colid),_data(S._data),
			_field(S._field)
		{}

#if 0
		template<class _OtherField>
		SparseMatrix2<_Field, SparseMatrix2Format::COO> (const SparseMatrix2<_OtherField, SparseMatrix2Format::COO> & S) :
			_rownb(S._rownb),_colnb(S._colnb),
			_nbnz(S._nbnz),
			_rowid(S._rowid),_colid(S._colid),_data(S._data),
			_field(S._field)
		{}
#endif

		// XXX only for CSR
		template<typename _Tp1, typename _Rw1 = SparseMatrix2Format::CSR>
		struct rebind ;

		template<typename _Tp1>
		struct rebind<_Tp1/*  ,SparseMatrix2Format::COO */ > {
			typedef SparseMatrix2<_Tp1, SparseMatrix2Format::CSR> other;

			void operator() (other & Ap, const Self_t& A)
			{
				// Ap = new other(F, A.rowdim(), A.coldim());

				typename _Tp1::Element e;

				Hom<typename Self_t::Field, _Tp1> hom(A.field(), Ap.field());
				size_t j = 0 ;
				for (size_t i = 0 ; i < A.rowdim()+1 ; ++i) {
					Ap.setStart(i,A.getStart(i));
				}
				std::vector<size_t> offset(A.rowdim()+1,0UL);
				bool changed = false ;
				for (size_t i = 0 ; i < A.rowdim() ; ++i) {
					for (size_t k = A.getStart(i) ; k < A.getStart(i+1) ; ++k) {
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
					for (size_t i = 0 ; i < A.rowdim() ; ++i) {
						offset[i+1] += offset[i] ;
					}
					for (size_t i = 0 ; i < A.rowdim()+1 ; ++i) {
						Ap.setStart(i,offset[i]);
					}
				}


				if (j != Ap.size())
					Ap.resize(j);
			}
		};

		template<typename _Tp1, typename _Rw1>
		SparseMatrix2 (const SparseMatrix2<_Tp1, _Rw1> &S, const Field& F) :
			_rownb(S.rowdim()),_colnb(S.coldim()),
			_nbnz(S.size()),
			_start(S.rowdim()+1),_colid(S.size()),_data(S.size()),
			_field(F)
		{
			typename SparseMatrix2<_Tp1,_Rw1>::template rebind<Field,SparseMatrix2Format::CSR>()(*this, S);
		}




		template<class VectStream>
		SparseMatrix2<_Field, SparseMatrix2Format::CSR> (const _Field & F, VectStream & stream) :
			_rownb(stream.size()),_colnb(stream.size()),
			_start(0),_colid(0),_data(0)
			, _nbnz(0)
			, _field(F)
		{
			// VectStream == RandomSparseStream<Field, typename Vector<Field>::SparseSeq>
			{

			}
			for (size_t i = 0 ; i< _rownb ; ++i) {
				typename Vector<Field>::SparseSeq lig_i ;
				stream >> lig_i ;
				_start[i+1] = lig_i.size();
				for (size_t j = 0 ; j < lig_i.size() ; ++j) {
					_nbnz++;
					resize(_nbnz);
					_colid[j] = lig_i[j].first ;
					_data[j] = lig_i[j].second ;
				}
			}
			for (size_t i = 0 ; i < rowdim() ; ++i)
				_start[i+1] += _start[i];
		}

		void resize(size_t nn)
		{
#ifndef NDEBUG
			if (nn < _nbnz) {
				std::cerr << "*** Warning *** you are possibly loosing data (CSR resize)" << std::endl;
				// could be a commentator()...
			}
#endif
			// _rowid.resize(nn);
			_colid.resize(nn);
			_data.resize(nn);
			_nbnz = nn ;
		}

		/*! Default converter.
		 * @param S a sparse matrix in any storage.
		 */
		template<class _OtherStorage>
		SparseMatrix2<_Field, SparseMatrix2Format::CSR> (const SparseMatrix2<_Field, _OtherStorage> & S) :
			_rownb(S._rownb),_colnb(S._colnb),
			_start(S.size()),_colid(S.size()),_data(S.size()),
			_field(S._field)
		{
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
		void importe(const SparseMatrix2<_Field,SparseMatrix2Format::COO> &S)
		{
			_rownb = S.rowdim() ;
			_colnb = S.coldim() ;

			_data = S.getData();
			_colid.resize(S.size());
			for (size_t i = 0 ; i < _rownb ; ++i)
				for (size_t j = S.getStart(i) ; j < S.getStart(i+1); ++j)
					_colid[j] = i ;

		}

		/*! Import a matrix in CSR format to CSR.
		 * @param S CSR matrix to be converted in CSR
		 */

		void importe(const SparseMatrix2<_Field,SparseMatrix2Format::CSR> &S)
		{
			resize(S.rowdim(), S.coldim(), S.size());

			_data = S.getData();
			_colid = S.getColid();
			_start = S.getStart();
		}

		/*! Export a matrix in CSR format from COO.
		 * @param S CSR matrix to be converted from COO
		 */
		SparseMatrix2<_Field,SparseMatrix2Format::CSR > &
		exporte(SparseMatrix2<_Field,SparseMatrix2Format::CSR> &S)
		{
			// S = *this ;

			S.resize(_rownb, _colnb, _nbnz);
			S.setData(  _data ) ;
			S.setStart( _start );
			S.setColdid (_colid );

			return S ;
		}
		//@}

		/*! In place transpose
		*/
		void transposeIn()
		{
			SparseMatrix2<_Field,SparseMatrix2Format::CSR> Temp(*this);
			resize(_colnb, _rownb, _nbnz );

			for (size_t i = 0 ; i <= _colnb ; ++i)
				_start[i] = 0 ;

			for (size_t i = 0 ; i < size() ; ++i)
				_start[_colid[i]+1] += 1 ;

			for (size_t i = 0 ; i < _colnb ; ++i)
				_start[i+1] += _start[i] ;
			{
				size_t i = 0 ;
				std::vector<size_t> done_col(_rownb,0);
				for (size_t nextlig = 1 ; nextlig <= _colnb ; ++nextlig) {
					// treating line before nextlig
					while (i < Temp._start[nextlig]){
						size_t cur_place ;
						cur_place = _start[Temp.getColid(i)] + done_col[Temp.getColid(i)] ;
						_data[ cur_place ]  = Temp.getData(i) ;
						_colid[ cur_place ] = nextlig-1 ;
						done_col[Temp.getColid(i)] += 1 ;
						++i;
					}
				}
			}
		}

		/*! Transpose the matrix.
		 *  @param S [out] transpose of self.
		 *  @return a reference to \p S.
		 */
		SparseMatrix2<_Field,SparseMatrix2Format::CSR> &
		transpose(SparseMatrix2<_Field,SparseMatrix2Format::CSR> &S)
		{
			// assert(S.rowdim() == _colnb);
			// assert(S.coldim() == _rownb);
			S.importe(*this);
			S.transposeIn();
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

		/** Get a read-only individual entry from the matrix.
		 * @param i Row _colid
		 * @param j Column _colid
		 * @return Const reference to matrix entry
		 */
		constElement &getEntry(const size_t &i, const size_t &j) const
		{
			assert(i<_rownb);
			assert(j<_colnb);
			// Could be improved by adding an initial guess j/rodim*size()
			// typedef typename std::vector<size_t>::iterator myIterator ;
			typedef typename std::vector<size_t>::const_iterator myConstIterator ;

			std::pair<myConstIterator,myConstIterator> bounds = std::equal_range (_start.begin(), _start.end(), i);
			size_t ibeg = (size_t)(bounds.first-_start.begin());
			size_t iend = (size_t)(bounds.second-_start.begin())-ibeg;
			if (!iend)
				return _field.zero;

			myConstIterator beg = _colid.begin()+(ptrdiff_t)ibeg ;
			myConstIterator low = std::lower_bound (beg, beg+(ptrdiff_t)iend, j);
			if (low == beg+(ptrdiff_t)iend)
				return _field.zero;
			else {
				// not sure
				size_t la = (size_t)(low-_colid.begin()) ;
				return _data[la] ;
			}
		}

		Element      &getEntry (Element &x, size_t i, size_t j) const
		{
			return x = getEntry (i, j);
		}

		/** Set an individual entry.
		 * Setting the entry to 0 will not remove it from the matrix
		 * @param i Row _colid of entry
		 * @param j Column _colid of entry
		 * @param value Value of the new entry
		 * @todo make it faster if i is 0 or m-1 ?
		 */
		void setEntry(const size_t &i, const size_t &j, const Element& e)
		{
			if (_field.isZero(e)) {
				return clearEntry(i,j);
			}
			typedef typename std::vector<size_t>::iterator myIterator ;
			std::pair<myIterator,myIterator> bounds = std::equal_range (_start.begin(), _start.end(), i);
			size_t ibeg = bounds.first-_start.begin();
			size_t iend = (bounds.second-_start.begin())-ibeg;
			if (!iend) {
				_start.insert(_start.begin()+ibeg,i);
				_colid.insert(_colid.begin()+ibeg,j);
				_data.insert( _data.begin() +ibeg,e);
				return ;
			}
			myIterator beg = _colid.begin()+(ptrdiff_t)ibeg ;
			myIterator low = std::lower_bound (beg, beg+(ptrdiff_t)iend, j);
			ibeg = low-_colid.begin();
			if (low == beg+(ptrdiff_t)iend) {
				_start.insert(_start.begin()+(ptrdiff_t)ibeg,i);
				_colid.insert(_colid.begin()+(ptrdiff_t)ibeg,j);
				_data.insert( _data.begin() +(ptrdiff_t)ibeg,e);
				return ;
			}
			else {
				_data[ibeg] = e ;
				return ;
			}
		}


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
			assert(i<_rownb);
			assert(j<_colnb);
			// Could be improved by adding an initial guess j/rodim*size()
			typedef typename std::vector<size_t>::iterator myIterator ;

			std::pair<myIterator,myIterator> bounds = std::equal_range (_start.begin(), _start.end(), i);
			size_t ibeg = bounds.first-_start.begin();
			size_t iend = (bounds.second-_start.begin())-ibeg;
			if (!iend) {
				_start.insert(_start.begin()+ibeg,i);
				_colid.insert(_colid.begin()+ibeg,j);
				_data.insert( _data.begin() +ibeg,_field.zero);
				return _data[ibeg];
			}
			myIterator beg = _colid.begin()+ibeg ;
			myIterator low = std::lower_bound (beg, beg+(ptrdiff_t)iend, j);
			if (low == beg+(ptrdiff_t)iend) {
				_start.insert(_start.begin()+ibeg,i);
				_colid.insert(_colid.begin()+ibeg,j);
				_data.insert( _data.begin() +ibeg,_field.zero);
				return _data[ibeg];
			}
			else {
				size_t la = low-_colid.begin() ;
				return _data[la] ;
			}
		}

		/** Write a matrix to the given output stream using field read/write.
		 * @param os Output stream to which to write the matrix
		 * @param format Format with which to write
		 */
		template<class Format>
		std::ostream & write(std::ostream &os,
				     Format = SparseFileFormat::CSR())
		{
			return this->writeSpecialized(os,Format());
		}

		/** Read a matrix from the given input stream using field read/write
		 * @param file Input stream from which to read the matrix
		 * @param format Format of input matrix
		 * @return ref to \p file.
		 */
		template<class Format>
		std::istream& read (std::istream &file,
				    Format = SparseFileFormat::CSR())
		{
			return readSpecialized(file,Format());
		}

		/*! @internal
		 * @brief Deletes the entry.
		 * Deletes \c A(i,j) if it exists.
		 * @param i row _colid
		 * @param j col _colid
		 */
		void clearEntry(const size_t &i, const size_t &j)
		{
			assert(i<_rownb);
			assert(j<_colnb);
			// Could be improved by adding an initial guess j/rodim*size()
			typedef typename std::vector<size_t>::iterator myIterator ;

			std::pair<myIterator,myIterator> bounds = std::equal_range (_start.begin(), _start.end(), i);
			size_t ibeg = bounds.first-_start.begin();
			size_t iend = (bounds.second-_start.begin())-ibeg;
			if (!iend)
				return ;

			myIterator beg = _colid.begin()+ibeg ;
			myIterator low = std::lower_bound (beg, beg+(ptrdiff_t)iend, j);
			if (low == beg+(ptrdiff_t)iend)
				return ;
			else {
				// not sure
				size_t la = low-_colid.begin() ;
				_start.erase(_start.begin()+la);
				_colid.erase(_colid.begin()+la);
				_data. erase(_data. begin()+la);
				return  ;
			}
		}

		/*! @internal
		 * @brief cleans 0 entries.
		 */
		void clean()
		{
			for (size_t i = 0 ; i < _data.size() ; ++i) {
				if ( _field.isZero(_data[i]) ) {
					_start.erase(_start.begin()+i);
					_colid.erase(_colid.begin()+i);
					_data. erase(_data. begin()+i);
				}
			}
			return ;
		}

		template<class Vector>
		Vector& apply(Vector &y, const Vector& x) const
		{
			//! @bug why always zero-assign ?
			for (size_t i = 0 ; i < y.size() ; ++i)
				_field.assign(y[i],_field.zero);

			for (size_t i = 0 ; i < _nbnz ; ++i)
				_field.axpyin( y[_start[i]], _data[i], x[_colid[i]] );

			return y;
		}

		template<class Vector>
		Vector& applyTranspose(Vector &y, const Vector& x) const
		{
			//! @bug why always zero-assign ?
			for (size_t i = 0 ; i < y.size() ; ++i)
				_field.assign(y[i],_field.zero);

			for (size_t i = 0 ; i < _nbnz ; ++i)
				_field.axpyin( y[_colid[i]], _data[i], x[_start[i]] );

			return y;
		}

		const Field & field()  const
		{
			return _field ;
		}


	private :

		/*! @internal
		 * write for CSR format.
		 */
		std::ostream & writeSpecialized(std::ostream &os,
						SparseFileFormat::CSR) const
		{
			os << _rownb << ' ' << _colnb  << ' ' << size() << std::endl;
			size_t lig = 0 ;
			size_t i = 0 ;
			while(i < size()) {
				while(lig < _start[i]) {
					os << "-1" << std::endl;
					++lig ;
				}
				while (lig == _start[i]) {
					_field.write(_data[i], os << _colid[i] << ' ') << std::endl;
					++i;
				}
				++lig ;
			}
			return os << "0 0 0" << std::endl;
		}

		/*! @internal
		 * write for COO format.
		 */
		std::ostream & writeSpecialized(std::ostream &os,
						SparseFileFormat::COO) const
		{
			os << _rownb << ' ' << _colnb  << ' ' << size() << std::endl;
			size_t i = 0 ;
			while(i < size()) {
				_field.write(_data[i], os << _start[i] << ' ' << _colid[i] << ' ') << std::endl;
				++i;
			}
			return os << "0 0 0" << std::endl;
		}

		std::ostream & writeSpecialized(std::ostream &os,
						FileFormatTag format) const
		{
			switch (format) {
			case FORMAT_MAPLE:
				{

					linbox_check(_colnb > 0);
					os << "[";
					bool firstrow=true;
					size_t idx = 0 ;

					for (size_t i = 0 ; i < _rownb ; ++i ) {
						if (firstrow) {
							os << "[";
							firstrow =false;
						}
						else
							os << ", [";


						for (size_t j = 0; j < _colnb ; ++j) {
							if (_colid[idx] == j && _start[idx] ==i) {
								_field.write (os, _data[idx]);
								++idx;
							}
							else {
								_field.write (os, _field.zero);
							}

							if (j < _colnb - 1)
								os << ", ";
						}

						os << " ]";
					}

					os << "]";
					linbox_check(idx == _nbnz);

					break;
				}
			default :
				os << "I don't know" << std::endl;

			}
			return os ;

		}


		/*! @internal
		 * Read for CSR format.
		 */
		std::istream & readSpecialized(std::istream &is,
					       SparseFileFormat::CSR)
		{
			size_t nnz = 0 ;
			bool sms = true ;
			std::string firstLine ;
			std::string x ;
			getline(is, firstLine);
			std::istringstream line(firstLine);
			line >> _rownb >> _colnb >> x ;
			size_t mem = 10 ;
			if (!_rownb || _colnb)
				throw LinBoxError("bad input");
			if (x.empty() || x.compare("M")) {  /* SMS */
				// mem = m ;
				_start.reserve(mem);
				_colid.reserve(mem);
				_data.reserve(mem);
			}
			else { /* SMF */
				sms = false ;
				std::istringstream (x) >> nnz ;
				if (!nnz)
					throw LinBoxError("bad input");
				mem = nnz ;
				_start.reserve(nnz);
				_colid.reserve(nnz);
				_data.reserve(nnz);
			}
			Element z ;
			if (sms) { /*  SMS */
				size_t lig = 0 ;
				nnz = 0 ;
				int n ;
				while (is>>n) {
					if (n == 0)
						break;
					while (n == -1) {
						++lig ;
						is >> n ;
					}
					_field.read(is,z)  ;
					if (n<0 || lig >=_rownb || n >> _colnb)
						throw LinBoxError("bad input");
					if (!_field.isZero(z)){
						if (mem == nnz) {
							mem+=20 ;
							_start.resize(mem);
							_colid.resize(mem);
							_data.resize (mem);
						}

						_start[nnz]= lig ;
						_colid[nnz]= n ;
						_data[nnz] = z ;
						++nnz ;
					}
				}
				_start.resize(nnz);
				_colid.resize(nnz);
				_data.resize (nnz);

			}
			else { /*  SMF */
				size_t lig = 0 ;
				int n ;
				size_t loc = 0;
				while (is>>n) {
					if (n == 0)
						break;
					while (n == -1) {
						++lig ;
						is >> n ;
					}
					_field.read(is,z)  ;
					if (n<0 || lig >=_rownb || n >> _colnb)
						throw LinBoxError("bad input");
					if (!_field.isZero(z)){
						_start[loc]= lig ;
						_colid[loc]= n ;
						_data[loc] = z ;
						++loc ;
					}
				}
				if (loc > nnz)
					throw LinBoxError("bad input");
				_start.resize(loc);
				_colid.resize(loc);
				_data.resize (loc);
			}
			return is ;
		}

		/*! @internal
		 * Read for COO format.
		 */
		std::istream & readSpecialized(std::istream &is,
					       SparseFileFormat::COO)
		{
			size_t nnz = 0;
			bool sms = true ;
			std::string firstLine ;
			getline(is, firstLine);
			std::istringstream line(firstLine);
			std::string x ;
			line >> _rownb >> _colnb >> x ;
			size_t mem  = 10 ;
			if (!_rownb || _colnb)
				throw LinBoxError("bad input");
			if (x.empty() || x.compare("M")) {  /* SMS */
				// mem = m ;
				_start.reserve(mem);
				_colid.reserve(mem);
				_data.reserve(mem);
			}
			else { /* SMF */
				sms = false ;
				std::istringstream (x) >> nnz ;
				if (!nnz)
					throw LinBoxError("bad input");
				mem = nnz ;
				_start.reserve(nnz);
				_colid.reserve(nnz);
				_data.reserve(nnz);
			}
			Element z ;
			if (sms) { /*  SMS */
				// size_t lig = 0 ;
				nnz = 0 ;
				int m,n ;
				while (is>>m >> n) {
					if (m == 0 && n == 0)
						break;
					if (n<0 || m<0 ||  m >=_rownb || n >> _colnb)
						throw LinBoxError("bad input");
					_field.read(is,z)  ;
					if (!_field.isZero(z)){
						if (mem == nnz) {
							mem+=20 ;
							_start.resize(mem);
							_colid.resize(mem);
							_data.resize (mem);
						}
						_start[nnz]= m ;
						_colid[nnz]= n ;
						_data[nnz] = z ;
						++nnz ;
					}
				}
				_start.resize(nnz);
				_colid.resize(nnz);
				_data.resize (nnz);

			}
			else { /*  SMF */
				size_t loc = 0 ;
				int m,n ;
				while (is>>m >> n) {
					if (m == 0 && n == 0)
						break;
					if (n<0 || m<0 ||  m >=_rownb || n >> _colnb)
						throw LinBoxError("bad input");
					_field.read(is,z)  ;
					if (!_field.isZero(z)){
						_start[loc]= m ;
						_colid[loc]= n ;
						_data[loc] = z ;
						++loc ;
					}
				}

				if (loc > nnz)
					throw LinBoxError("bad input");
				_start.resize(loc);
				_colid.resize(loc);
				_data.resize (loc);

			}
			return is ;
		}


	public:
		// pseudo iterators
		size_t getRowid(const size_t i) const
		{
			return _start[i];
		}

		void setRowid(const size_t i, const size_t j)
		{
			if (i>=_nbnz) this->resize(i);
			_start[i]=j;
		}

		size_t getColid(const size_t i) const
		{
			return _colid[i];
		}

		void setColid(const size_t i, const size_t j)
		{
			if (i>=_nbnz) this->resize(i);
			_start[i]=j;
		}

		const Element & getData(const size_t i) const
		{
			return _data[i];
		}

		void setData(const size_t i, const Element & e)
		{
			if (i>=_nbnz) this->resize(i);
			field().assign(_data[i],e);
		}

	protected :

		size_t              _rownb ;
		size_t              _colnb ;
		size_t               _nbnz ;

		std::vector<size_t> _start ;
		std::vector<size_t> _colid ;
		std::vector<Element> _data ;

		const _Field            & _field ;
	};





} // namespace LinBox

#endif // __LINBOX_sparse_matrix_sparse_csr_matrix_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
