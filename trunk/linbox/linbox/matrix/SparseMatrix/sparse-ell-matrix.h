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


#ifndef __LINBOX_sparse_matrix_sparse_ell_matrix_H
#define __LINBOX_sparse_matrix_sparse_ell_matrix_H

#include <utility>
#include <iostream>
#include <algorithm>

#include "linbox-config.h"
#include "linbox/util/debug.h"
#include "sparse-domain.h"


namespace LinBox
{


	/** Sparse matrix, Coordinate storage.
	 *
	 * \ingroup matrix
	 * \ingroup sparse
	 */
	template<class _Field>
	class SparseMatrix2<_Field, SparseMatrixFormat::ELL > {
	public :
		typedef _Field                             Field ; //!< Field
		typedef typename _Field::Element         Element ; //!< Element
		typedef const Element               constElement ; //!< const Element
		typedef SparseMatrixFormat::ELL         Storage ; //!< Matrix Storage Format
		typedef SparseMatrix2<_Field,Storage>     Self_t ; //!< Self type
		// typedef Vector<_Field,VectorStorage::Sparse> Rep ;

		/*! Constructors.
		 * @todo convert from other matrix.
		 *
		 */
		//@{
		SparseMatrix2<_Field, SparseMatrixFormat::ELL> () :
			_rownb(0),_colnb(0),_maxc(0)
			,_nbnz(0)
			,_colid(0),_data(0)
			, _field()
		{
		}

		SparseMatrix2<_Field, SparseMatrixFormat::ELL> (const _Field & F) :
			_rownb(0),_colnb(0)
			,_maxc(0)
			,_nbnz(0),
			_colid(0),_data(0)
			, _field(F)
		{
		}

		SparseMatrix2<_Field, SparseMatrixFormat::ELL> (const _Field & F, size_t m, size_t n) :
			_rownb(m),_colnb(n)
			,_maxc(0)
			,_nbnz(0)
			,_colid(0),_data(0)
			, _field(F)
		{
		}

		SparseMatrix2<_Field, SparseMatrixFormat::ELL> (const _Field & F,
							       size_t m, size_t n,
							       size_t z) :
			_rownb(m),_colnb(n),
			_maxc(0),
			_nbnz(z),
			_colid(z),_data(z),
			_field(F)
		{
		}

		SparseMatrix2<_Field, SparseMatrixFormat::ELL> (const SparseMatrix2<_Field, SparseMatrixFormat::CSR> & S) :
			_rownb(S._rownb),_colnb(S._colnb),
			_maxc(S._maxc),
			_nbnz(S._nbnz),
			_colid(S._colid),_data(S._data),
			_field(S._field)
		{}

#if 0
		template<class _OtherField>
		SparseMatrix2<_Field, SparseMatrixFormat::COO> (const SparseMatrix2<_OtherField, SparseMatrixFormat::COO> & S) :
			_rownb(S._rownb),_colnb(S._colnb),
			_nbnz(S._nbnz),
			_rowid(S._rowid),_colid(S._colid),_data(S._data),
			_field(S._field)
		{}
#endif

		// XXX only for CSR
		template<typename _Tp1, typename _Rw1 = SparseMatrixFormat::ELL>
		struct rebind ;

		template<typename _Tp1>
		struct rebind<_Tp1/*  ,SparseMatrixFormat::COO */ > {
			typedef SparseMatrix2<_Tp1, SparseMatrixFormat::ELL> other;

			void operator() (other & Ap, const Self_t& A)
			{
				linbox_check(A.consistent());
				// Ap = new other(F, A.rowdim(), A.coldim());
				Ap.resize(A.rowdim(),A.coldim(),A.size(),A.ld());

				typename _Tp1::Element e;

				Hom<typename Self_t::Field, _Tp1> hom(A.field(), Ap.field());
				size_t newz = 0 ;
				for (size_t i = 0 ; i < A.rowdim() ; ++i) {
					size_t j = 0 ;
					for (size_t k = 0 ; k < A.ld() ; ++k) {
						if (A.field().isZero(A.getData(i,k)))
							break;
						hom. image ( e, A.getData(i,k) );
						if (!Ap.field().isZero(e)) {
							Ap.setColid(i,j,A.getColid(i,k));
							Ap.setData(i,j,e);
							++j;
						}
						else
							++newz ;
					}
				}
				Ap.setSize(Ap.size() - newz) ;
			}
		};

		template<typename _Tp1, typename _Rw1>
		SparseMatrix2 (const SparseMatrix2<_Tp1, _Rw1> &S, const Field& F) :
			_rownb(S.rowdim()),_colnb(S.coldim()),
			_maxc(0),
			_nbnz(S.size()),
			_colid(0),_data(0),
			_field(F)
		{
			typename SparseMatrix2<_Tp1,_Rw1>::template rebind<Field,SparseMatrixFormat::ELL>()(*this, S);
		}




		template<class VectStream>
		SparseMatrix2<_Field, SparseMatrixFormat::ELL> (const _Field & F, VectStream & stream) :
			_rownb(stream.size()),_colnb(stream.dim()),
			_colid(0),_data(0)
			,_maxc(0)
			, _nbnz(0)
			, _field(F)
		{
			SparseMatrix2<_Field,SparseMatrixFormat::CSR> Tmp(F,stream);
			importe(Tmp);
		}



		void resize(const size_t & mm, const size_t & nn, const size_t & zz, const size_t & ll)
		{
#ifndef NDEBUG
			if (_maxc && ll != _maxc)
				std::cout << " ***Warning*** loosing data in ELL resize " << std::endl;
#endif
			_rownb = mm ;
			_colnb = nn ;
			_nbnz  = zz;
			_maxc  = ll;

			_colid.resize(mm*ll,0);
			_data .resize(mm*ll,0);
		}

		/*! Default converter.
		 * @param S a sparse matrix in any storage.
		 */
		template<class _OtherStorage>
		SparseMatrix2<_Field, SparseMatrixFormat::ELL> (const SparseMatrix2<_Field, _OtherStorage> & S) :
			_rownb(S._rownb),_colnb(S._colnb),_nbnz(S.size()),
			_maxc(0),_colid(0),_data(0),
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
		void importe(const SparseMatrix2<_Field,SparseMatrixFormat::CSR> &S)
		{
			// can be sped up on multicores.
			for (size_t i = 0 ; i < S.rowdim() ; ++i)
				_maxc = std::max(_maxc, S.getEnd(i)-S.getStart(i));

			resize(S.rowdim(), S.coldim(), S.size(),_maxc);

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

		void importe(const SparseMatrix2<_Field,SparseMatrixFormat::ELL> &S)
		{
			resize( S.rowdim(), S.coldim(), S.size() , S.ld());

			setColid(S.getColid());
			setData(S.getData());

		}

		template<class _OtherStorage>
		void importe(const SparseMatrix2<_Field,_OtherStorage> &S)
		{
			SparseMatrix2<_Field,SparseMatrixFormat::CSR> Tmp(S);
			this->importe(S);

		}

		/*! Export a matrix in CSR format from COO.
		 * @param S CSR matrix to be converted from COO
		 */
		SparseMatrix2<_Field,SparseMatrixFormat::CSR > &
		exporte(SparseMatrix2<_Field,SparseMatrixFormat::CSR> &S) const
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
			SparseMatrix2<_Field,SparseMatrixFormat::ELL> Temp(*this);
			Temp.transposeIn();
			importe(Temp);
		}

		/*! Transpose the matrix.
		 *  @param S [out] transpose of self.
		 *  @return a reference to \p S.
		 */
		SparseMatrix2<_Field,SparseMatrixFormat::ELL> &
		transpose(SparseMatrix2<_Field,SparseMatrixFormat::ELL> &S)
		{
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

		void setSize(const size_t & z)
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
#if 0
			// std::cout << "get entry : " << i << ',' << j << std::endl;
			linbox_check(i<_rownb);
			linbox_check(j<_colnb);
			typedef typename std::vector<size_t>::const_iterator myConstIterator ;

			size_t ibeg = _start[i] ;
			size_t iend = _start[i+1] ;
			if (ibeg == iend) {
				// std::cout << "get entry : " << 0 << std::endl;
				return field().zero;
			}

			myConstIterator beg = _colid.begin() ;
			myConstIterator low = std::lower_bound (beg+(ptrdiff_t)ibeg, beg+(ptrdiff_t)iend, j);
			if (low == beg+(ptrdiff_t)iend) {
				// std::cout << "get entry : " << 0 << std::endl;
				return field().zero;
		}
			else {
				// not sure
				size_t la = (size_t)(low-beg) ;
				// std::cout << "get entry : " << _data[la] << std::endl;
				return _data[la] ;
			}
#endif
		}

		Element      &getEntry (Element &x, size_t i, size_t j) const
		{
			return x = getEntry (i, j);
		}

		void finalize(){} // end construction after a sequence of setEntry calls.
		/** Set an individual entry.
		 * Setting the entry to 0 will not remove it from the matrix
		 * @param i Row _colid of entry
		 * @param j Column _colid of entry
		 * @param value Value of the new entry
		 * @todo make it faster if i is 0 or m-1 ?
		 */
		void setEntry(const size_t &i, const size_t &j, const Element& e)
		{
#if 0
			linbox_check(i<_rownb);
			linbox_check(j<_colnb);

			if (field().isZero(e)) {
				return clearEntry(i,j);
			}

			size_t ibeg = _start[i];
			size_t iend = _start[i+1];
			// element does not exist, insert
			if (ibeg == iend) {
				for (size_t k = i+1 ; k <= _rownb ; ++k) _start[k] += 1 ;
				_colid.insert(_colid.begin()+ibeg,j);
				_data.insert( _data.begin() +ibeg,e);
				return ;
			}
			// element may exist
			typedef typename std::vector<size_t>::iterator myIterator ;
			myIterator beg = _colid.begin() ;
			myIterator low = std::lower_bound (beg+(ptrdiff_t)ibeg, beg+(ptrdiff_t)iend, j);
			ibeg = low-beg;
			// insert
			if (low == beg+(ptrdiff_t)iend) {
				for (size_t k = i ; k <= _rownb ; ++k) _start[k] += 1 ;
				_colid.insert(_colid.begin() +(ptrdiff_t)ibeg,j);
				_data.insert( _data. begin() +(ptrdiff_t)ibeg,e);
				return ;
			}
			// replace
			else {
				_data[ibeg] = e ;
				return ;
			}
#endif
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
#if 0
			linbox_check(i<_rownb);
			linbox_check(j<_colnb);
			// Could be improved by adding an initial guess j/rodim*size()

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
#endif
		}

		/** Write a matrix to the given output stream using field read/write.
		 * @param os Output stream to which to write the matrix
		 * @param format Format with which to write
		 */
		template<class Format>
		std::ostream & write(std::ostream &os,
				     Format = SparseFileFormat::CSR()) const
		{
			return this->writeSpecialized(os,Format());
		}

		std::ostream & write(std::ostream &os,
				     enum LINBOX_enum(Tag::FileFormat) ff  = Tag::FileFormat::Maple) const
		{
			return this->writeSpecialized(os,ff);
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
#if 0
			linbox_check(i<_rownb);
			linbox_check(j<_colnb);

			size_t ibeg = _start[i];
			size_t iend = _start[i+1];
			if (ibeg == iend)
				return ;

			typedef typename std::vector<size_t>::iterator myIterator ;
			myIterator beg = _colid.begin() ;
			myIterator low = std::lower_bound (beg+(ptrdiff_t)ibeg, beg+(ptrdiff_t)iend, j);
			if (low == beg+(ptrdiff_t)iend)
				return ;
			else {
				// not sure
				size_t la = low-beg ;
				for (size_t k = i+1 ; k <= _rownb ; ++k) _start[k] -= 1 ;
				_colid.erase(_colid.begin()+la);
				_data. erase(_data. begin()+la);
				return  ;
			}
#endif
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
			linbox_check(consistent());
			prepare(field(),y,a);

			for (size_t i = 0 ; i < _rownb ; ++i) {
				for (size_t k = 0   ; k < _maxc ; ++k)
					if (!field().isZero(getData(i,k)))
						field().axpyin( y[i], getData(i,k), x[getColid(i,k)] ); //! @todo delay !!!
					else
						break;
			}

			return y;
		}

		class Helper ; // transpose

		// y= A^t x
		// y[i] = sum(A(j,i) x(j)
		template<class Vector>
		Vector& applyTranspose(Vector &y, const Vector& x, const Element & a ) const
		{
			linbox_check(consistent());
			//! @bug if too big, create transpose.
			prepare(field(),y,a);

			for (size_t i = 0 ; i < _rownb ; ++i)
				for (size_t k = 0   ; k < _maxc ; ++k)
					if (!field().isZero(getData(i,k)))
						field().axpyin( y[getColid(i,k)], getData(i,k), x[i] ); //! @todo delay !!!
					else
						break;

			return y;
		}

		template<class Vector>
		Vector& apply(Vector &y, const Vector& x ) const
		{
			return apply(y,x,field().zero);
		}
		template<class Vector>
		Vector& applyTranspose(Vector &y, const Vector& x ) const
		{
			return apply(y,x,field().zero);
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
			bool ok = true ;
			for (size_t i = 0 ; i < _rownb ; ++i) {
				bool row_ok = true ;
				bool zero = false ;
				for (size_t j = 0 ; j < _maxc ; ++j) {
					if (field().isZero(getData(i,j))) {
						zero = true ;
					}
					else {
						if (zero)
							return false ; // elements after a 0...
						++nbnz ;

					}
				}
			}
			return (nbnz == _nbnz);

		}


	private :

		std::ostream & writeSpecialized(std::ostream &os,
						LINBOX_enum(Tag::FileFormat) format) const
		{
			SparseMatrix2<Field,SparseMatrixFormat::CSR> Temp(field());
			this->exporte(Temp);
			Temp.write(os,format);
			return os ;

		}

#if 0 /*  not updated and to be cleaned */
		/*! @internal
		 * write for CSR format.
		 * @bug wrong.
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
					field().write(_data[i], os << _colid[i] << ' ') << std::endl;
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
			for (size_t i = 0 ; i < rowdim() ; ++i)
				for (size_t j = _start[i] ; j < _start[j+1] ; ++j)
					field().write(_data[j], os << i << ' ' << _colid[j] << ' ') << std::endl;

			return os << "0 0 0" << std::endl;
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
					field().read(is,z)  ;
					if (n<0 || lig >=_rownb || n >> _colnb)
						throw LinBoxError("bad input");
					if (!field().isZero(z)){
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
					field().read(is,z)  ;
					if (n<0 || lig >=_rownb || n >> _colnb)
						throw LinBoxError("bad input");
					if (!field().isZero(z)){
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
					field().read(is,z)  ;
					if (!field().isZero(z)){
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
					field().read(is,z)  ;
					if (!field().isZero(z)){
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
#endif

		void reshape(const size_t &ll)
		{
			linbox_check(ll > _maxc);

			_data .resize(ll*_rownb);
			_colid.resize(ll*_rownb);
			for (size_t i = _rownb ; i>1 ; --i ) {
				std::copy(_data.begin()+(i-1)*_maxc, _data.begin()+(i)*_maxc+1, _data.begin()+(i-1)*ll);
				for (size_t j = (i-1)*_maxc ; j < (i-1)*ll ; ++j ) _data[j] = field().zero ;
			}
			for (size_t i = _rownb ; i>1 ; --i ) {
				std::copy(_colid.begin()+(i-1)*_maxc, _colid.begin()+(i)*_maxc+1, _colid.begin()+(i-1)*ll);
				for (size_t j = (i-1)*_maxc ; j < (i-1)*ll ; ++j ) _colid[j] = 0 ;
			}

			_maxc = ll ;
		}

	public:
		// pseudo iterators

		size_t getColid(const size_t & i, const size_t j) const
		{
			return _colid[i*_maxc+j];
		}

		void setColid(const size_t & i, const size_t & j, const size_t & k)
		{
			linbox_check(i <= _rownb);
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
			return _data[ i*_maxc+j ];
		}

		void setData(const size_t & i, const size_t & j, const Element & e)
		{
			field().assign(_data[i*_maxc+j],e);
		}

		void setData(std::vector<size_t> & new_data)
		{
			_data = new_data ;
		}

		std::vector<size_t>  getData( ) const
		{
			return _data ;
		}

		size_t ld() const
		{
			return _maxc;
		}
	protected :

		size_t              _rownb ;
		size_t              _colnb ;
		size_t               _maxc ;
		size_t               _nbnz ;

		std::vector<size_t> _colid ;
		std::vector<Element> _data ;

		const _Field            & _field ;
	};





} // namespace LinBox

#endif // __LINBOX_sparse_matrix_sparse_ell_matrix_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
