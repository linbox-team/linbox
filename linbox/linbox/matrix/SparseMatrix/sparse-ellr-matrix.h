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

/*! @file matrix/SparseMatrix/sparse-ellr-matrix.h
 * @ingroup sparsematrix
 */


#ifndef __LINBOX_matrix_sparsematrix_sparse_ellr_matrix_H
#define __LINBOX_matrix_sparsematrix_sparse_ellr_matrix_H

#include <utility>
#include <iostream>
#include <algorithm>

#include "linbox/linbox-config.h"
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
	class SparseMatrix<_Field, SparseMatrixFormat::ELL_R > {
	public :
		typedef _Field                             Field ; //!< Field
		typedef typename _Field::Element         Element ; //!< Element
		typedef const Element               constElement ; //!< const Element
		typedef SparseMatrixFormat::ELL_R         Storage ; //!< Matrix Storage Format
		typedef SparseMatrix<_Field,Storage>     Self_t ; //!< Self type
		// typedef Vector<_Field,VectorStorage::Sparse> Rep ;

		/*! Constructors.
		 * @todo convert from other matrix.
		 *
		 */
		//@{
		SparseMatrix<_Field, SparseMatrixFormat::ELL_R> () :
			_rownb(0),_colnb(0),_maxc(0)
			,_nbnz(0)
			,_colid(0),_data(0)
			,_rowid(0)
			, _field()
		{
		}

		SparseMatrix<_Field, SparseMatrixFormat::ELL_R> (const _Field & F) :
			_rownb(0),_colnb(0)
			,_maxc(0),_nbnz(0)
			,_rowid(0),_colid(0)
			,_data(0), _field(F)
		{
		}

		SparseMatrix<_Field, SparseMatrixFormat::ELL_R> (const _Field & F, size_t m, size_t n) :
			_rownb(m),_colnb(n)
			,_maxc(0),_nbnz(0)
			,_rowid(m,0),_colid(0)
			,_data(0), _field(F)
		{
		}

		SparseMatrix<_Field, SparseMatrixFormat::ELL_R> (const _Field & F,
							       size_t m, size_t n,
							       size_t z) :
			_rownb(m),_colnb(n)
			,_maxc(0), _nbnz(z)
			,_rowid(m,0), _colid(z)
			,_data(z), _field(F)
		{
		}

		SparseMatrix<_Field, SparseMatrixFormat::ELL_R> (const SparseMatrix<_Field, SparseMatrixFormat::CSR> & S) :
			_rownb(S._rownb),_colnb(S._colnb)
			,_maxc(S._maxc), _nbnz(S._nbnz)
			,_rowid(S._rownb,0), _colid(S._colid)
			,_data(S._data), _field(S._field)
		{}

#if 0
		template<class _OtherField>
		SparseMatrix<_Field, SparseMatrixFormat::COO> (const SparseMatrix<_OtherField, SparseMatrixFormat::COO> & S) :
			_rownb(S._rownb),_colnb(S._colnb),
			_nbnz(S._nbnz),
			_rowid(S._rowid),_colid(S._colid),_data(S._data),
			_field(S._field)
		{}
#endif

		// XXX only for CSR
		template<typename _Tp1, typename _Rw1 = SparseMatrixFormat::ELL_R>
		struct rebind ;

		template<typename _Tp1>
		struct rebind<_Tp1/*  ,SparseMatrixFormat::COO */ > {
			typedef SparseMatrix<_Tp1, SparseMatrixFormat::ELL_R> other;

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
						if (A.field().isZero(A.getData(i,k))) {
							Ap.setRowid(i,  k) ;
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
					if (Ap.getRowid(i) == 0 && j>0)
						Ap.setRowid(i,  A.ld()) ;
				}
				Ap.setSize(Ap.size() - newz) ;
			}
		};

		template<typename _Tp1, typename _Rw1>
		SparseMatrix (const SparseMatrix<_Tp1, _Rw1> &S, const Field& F) :
			_rownb(S.rowdim()),_colnb(S.coldim()),
			_maxc(0),
			_nbnz(S.size())
			, _rowid(S.rowdim(),0)
			, _colid(0),_data(0)
			, _field(F)
		{
			typename SparseMatrix<_Tp1,_Rw1>::template rebind<Field,SparseMatrixFormat::ELL_R>()(*this, S);
		}




		template<class VectStream>
		SparseMatrix<_Field, SparseMatrixFormat::ELL_R> (const _Field & F, VectStream & stream) :
			_rownb(stream.size()),_colnb(stream.dim())
			, _maxc(0)
			, _nbnz(0)
			, _rowid(_rownb,0)
			, _colid(0),_data(0)
			, _field(F)
		{
			SparseMatrix<_Field,SparseMatrixFormat::CSR> Tmp(F,stream);
			importe(Tmp);
		}



		void resize(const size_t & mm, const size_t & nn, const size_t & zz = 0, const size_t & ll = 0)
		{
			// attention RowMajor/ColMajor
			if (!_maxc || mm == _maxc) {
				_colid.resize(mm*ll,0);
				_data .resize(mm*ll,field().zero);
				_rowid.resize(mm   ,0);
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
				_rowid.resize(mm   ,0);
			}
			_rownb = mm ;
			_colnb = nn ;
			_nbnz  = zz;
			_maxc  = ll;
		}

		/*! Default converter.
		 * @param S a sparse matrix in any storage.
		 */
		template<class _OtherStorage>
		SparseMatrix<_Field, SparseMatrixFormat::ELL_R> (const SparseMatrix<_Field, _OtherStorage> & S) :
			_rownb(S._rownb),_colnb(S._colnb),_nbnz(S.size()),
			_maxc(0),_colid(0),_data(0),
			_field(S._field)
			,_rowid(S.rowdim(),0)
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
		void importe(const SparseMatrix<_Field,SparseMatrixFormat::CSR> &S)
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
				_rowid[i] = S.getEnd(i)-S.getStart(i);
			}

		}

		/*! Import a matrix in CSR format to CSR.
		 * @param S CSR matrix to be converted in CSR
		 */

		void importe(const SparseMatrix<_Field,SparseMatrixFormat::ELL_R> &S)
		{
			resize( S.rowdim(), S.coldim(), S.size() , S.ld());

			setRowid(S.getRowid());
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
				for (size_t j = 0  ; j < _rowid[i]; ++j ) {
					S.setColid(k,getColid(i,j));
					S.setData(k,getData(i,j));
					++k;
				}
				S.setStart(i+1,S.getStart(i)+_rowid[i]);
			}
			linbox_check(k == _nbnz);

			return S ;

		}


		//@}

		/*! In place transpose. Not quite...
		*/
		void transposeIn()
		{
			SparseMatrix<_Field,SparseMatrixFormat::ELL_R> Temp(*this);
			Temp.transposeIn();
			importe(Temp);
		}

		/*! Transpose the matrix.
		 *  @param S [out] transpose of self.
		 *  @return a reference to \p S.
		 */
		SparseMatrix<_Field,SparseMatrixFormat::ELL_R> &
		transpose(SparseMatrix<_Field,SparseMatrixFormat::ELL_R> &S)
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

		void appendEntry(size_t i, size_t j, const Element & value) { setEntry(i,j,value) ;}
		// end construction after a sequence of setEntry calls.
		void finalize(){
			// could check that maxc is not too large and shrink ? Is is optimize job ?
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
			linbox_check(i<_rownb);
			linbox_check(j<_colnb);

			if (field().isZero(e)) {
				return clearEntry(i,j);
			}

			size_t * beg = &_colid[i*_maxc];
			Element * dat = &_data[i*_maxc];
			bool found = false;
			for (size_t k = 0 ; k < _maxc ; ++k) {
				// std::cout << "beg = " << beg[k] << std::endl;
				// std::cout << "dat = " << dat[k] << std::endl;
				// std::cout << "nz  = " << _nbnz << std::endl;
				if (dat[k]==0 && field().isZero(dat[k])) {
					// std::cout << "===2===" << std::endl;
					field().assign(dat[k], e) ;
					beg[k] = j;
					found = true;
					_rowid[i] +=1;
					++_nbnz ;
					// write_raw();
					break;
				}
				if (beg[k] == j) {
					// std::cout << "===1===" << std::endl;
					if (field().isZero(dat[k])) {
						_rowid[i] +=1;
						++_nbnz ;
					}
					field().assign(dat[k], e) ;
					beg[k] = j;
					found = true;
					// write_raw();
					break;
				}
				if (beg[k] > j) {
					// std::cout << "===3===" << std::endl;
					found = true;
					insert(i,k,j,e);
					// write_raw();
					break;
				}
			}
			if (!found) {
				// std::cout << "===4===" << std::endl;
				insert(i,_maxc,j,e);
				// write_raw();
			}
			return;
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
#endif
		}

		/** Write a matrix to the given output stream using field read/write.
		 * @param os Output stream to which to write the matrix
		 * @param format Format with which to write
		 */
		std::ostream & write(std::ostream &os,
				     LINBOX_enum(Tag::FileFormat) format = Tag::FileFormat::MatrixMarket) const
		{
			return SparseMatrixWriteHelper<Self_t>::write(*this,os,format);
		}


		/** Read a matrix from the given input stream using field read/write
		 * @param is Input stream from which to read the matrix
		 * @param format Format of input matrix
		 * @return ref to \p is.
		 */
		std::istream& read (std::istream &is,
				    LINBOX_enum(Tag::FileFormat) format = Tag::FileFormat::Detect)
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
			size_t k = 0 ;
			for ( ; k < _rowid[i] ; ++k) {
				if (_colid[i*_maxc+k] == j)
					break;
				if (_colid[i*_maxc+k] > j)
					return;
			}
			if (k == _rowid[i])
				return; // not found

			linbox_check(!field().isZero(_data[i*_maxc+k]));
				return;

			for (size_t l = k ; l <_rowid[i]-1  ; ++l) {
				field().assign(_data[i*_maxc+l],_data[i*_maxc+l+1]);
				_colid[i*_maxc+l] = _colid[i*_maxc+l+1];
			}
			field().assign(_data[i*_maxc+_rowid[i]-1],field().zero) ;
			_colid[i*_maxc+_rowid[i]] = 0;
			_rowid[i] -= 1 ;
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
			linbox_check(consistent());
			prepare(field(),y,a);

			FieldAXPY<Field> accu(field());
			for (size_t i = 0 ; i < _rownb ; ++i) {
				accu.reset();
				for (size_t k = 0   ; k < _rowid[i] ; ++k)
						// field().axpyin( y[i], getData(i,k), x[getColid(i,k)] ); //! @todo delay !!!
						 accu.mulacc( getData(i,k), x[getColid(i,k)] );
				accu.get(y[i]);
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
				for (size_t k = 0   ; k < _rowid[i] ; ++k)
						field().axpyin( y[getColid(i,k)], getData(i,k), x[i] ); //! @todo delay !!!

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
			return applyTranspose(y,x,field().zero);
		}


		const Field & field()  const
		{
			return _field ;
		}



		bool consistent() const
		{
			size_t nbnz = 0 ;
			// bool ok = true ;
			for (size_t i = 0 ; i < _rownb ; ++i) {
				// bool row_ok = true ;
				bool zero = false ;
				for (size_t j = 0 ; j < _maxc ; ++j) {
					if (field().isZero(getData(i,j))) {
						if (_rowid[i] != j && !zero) {
							std::cout << "@" << i << " : " << _rowid[i] <<"!=" << j << std::endl;
							return false ;
						}
						zero = true ;
					}
					else {
						if (zero) {
							std::cout <<  "@" << i << " : " << _rowid[i] <<"!= 0 (should be)"  << std::endl;
							return false ; // elements after a 0...
						}
						++nbnz ;
					}
				}
				if(!zero)
					if (_rowid[i] != _maxc) {
							std::cout <<  "@" << i << " : " << _rowid[i] <<"!= " << _maxc << " (should be)"  << std::endl;
						return false;
					}
			}
			return (nbnz == _nbnz);

		}



	private :

		std::ostream & writeSpecialized(std::ostream &os,
						LINBOX_enum(Tag::FileFormat) format) const
		{
			SparseMatrix<Field,SparseMatrixFormat::CSR> Temp(field());
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
		size_t getRowid(const size_t & i) const
		{
			return _rowid[i];
		}

		void setRowid(const size_t & i, const size_t & k)
		{
			linbox_check(i <= _rownb);
			_rowid[i]=k;
		}

		void setRowid(std::vector<size_t> new_rowid)
		{
			_rowid = new_rowid ;
		}

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

		bool nextTriple(size_t & i, size_t &j, Element &e) const
		{
			ptrdiff_t off =_triples.next(_rowid);

			i = _triples._row ;

			if (  i >= rowdim() ) {
				_triples.reset() ;
				return false ;
			}

			j = _colid[i*_maxc+off];
			e = _data[i*_maxc+off];

			return true;
		}

	private:

		void insert (const size_t &i, const size_t &k, const size_t &j, const Element& e)
		{
			// std::cout << "before inserting at " << k << std::endl; write_raw();
			if (k == _maxc) {
				resize(_rownb,_colnb,_nbnz,_maxc+1);
				_colid[_maxc*i+k-1] = j;
				field().assign(_data[_maxc*i+k-1],e);
				++_nbnz;
				_rowid[i] += 1;
				// std::cout << "after (1)" << std::endl; write_raw();
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
			_rowid[i] += 1;
			// std::cout << "after (2)" << std::endl; write_raw();
			return;
		}

		void write_raw() {
			std::cout << "rowids" << std::endl;
			std::cout << _rowid << std::endl;
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

		size_t              _rownb ;
		size_t              _colnb ;
		size_t               _maxc ; //!< longest column
		size_t               _nbnz ;

		std::vector<size_t> _rowid ;
		std::vector<size_t> _colid ; //!< \p _colid is \p _rownb x \p _maxc in RowMajor
		std::vector<Element> _data ; //!< \p _data  is \p _rownb x \p _maxc in RowMajor

		const _Field            & _field ;

		mutable struct _triples {
			ptrdiff_t _row ;
			ptrdiff_t _off ;
			_triples() :
				_row(-1)
				, _off(-1)
			{}
			ptrdiff_t next(const std::vector<size_t> & rowid) {
				++ _off ;
				if (_row == -1) {
					++_row ;
					return _off ;
				}
				while (_off >= rowid[_row]) {
					_row += 1 ;
					_off = 0 ;
					if (_row >= rowid.size())
						break;
				}
				return _off;
			}
			void reset() {
				_row = -1 ;
				_off = -1 ;
			}
		}_triples;
	};





} // namespace LinBox

#endif // __LINBOX_matrix_sparsematrix_sparse_ellr_matrix_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
