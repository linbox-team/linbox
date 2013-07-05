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


#ifndef __LINBOX_sparse_matrix_sparse_hyb_matrix_H
#define __LINBOX_sparse_matrix_sparse_hyb_matrix_H

#include <utility>
#include <iostream>
#include <algorithm>

#include "linbox-config.h"
#include "linbox/util/debug.h"
#include "linbox/matrix/sparse.h"

/*! @todo benchmark me */
#define HYB_ELL_THRESHOLD 0.9
/*! @todo benchmark me */
#define HYB_ELL_COO_THRESHOLD 0.1

namespace LinBox
{

	template<class Field>
	class Stats {
		const SparseMatrix2<Field,SparseMatrix2Format::CSR> & _Mat;
		const Field & _field ;
		size_t one ;
		size_t mone ;
		size_t avg ;
		// size_t avg_one ;
		// size_t avg_mone ;
		std::vector<size_t> row ;
		std::vector<size_t> row_one ;
		std::vector<size_t> row_mone ;
		std::vector<size_t> null_row ;
		// size_t mean ;
		size_t ell ;
		size_t ell_one ;
		size_t ell_mone ;
	public :
		Stats( const SparseMatrix2<Field,SparseMatrix2Format::CSR> & Mat) :
			_Mat(Mat),_field(_Mat.field())
			,one(0),mone(0)
			,row(Mat.rowdim(),0)
			,row_one(Mat.rowdim(),0)
			,row_mone(Mat.rowdim(),0)
			,null_row(0)
			,avg(0)
			,ell(0)
			// ,_off_ell(Mat.rowdim(),0);
		{
			if (Mat.size() > 100) { /*! @todo what is small ?  */
				for (size_t i = 0 ; i < Mat.rowdim() ; ++i) {
					if (Mat.getStart(i) == Mat.getEnd(i)) {
						null_row.push_back(i);
					}
					for (size_t k = Mat.getStart(i) ; k < Mat.getEnd(i) ; ++k) {
						if (_field.isOne(Mat.getEntry(k))) {
							++one ;
							row_one[i] += 1 ;
						}
						else if (_field.isMone(Mat.getEntry(k))) {
							++mone ;
							row_mone[i] += 1 ;
						}
						row[i] += 1;
					}
				}
			}
			avg = std::ceil((Mat.rowdim()-null_row.size())/Mat.size());
		}

		void getOptimsedFormat(SparseMatrix2<Field,SparseMatrix2Format::HYB> & hyb)
		{
			std::vector<SparseMatrix2Format::ANY>  try1 ;
			// std::sort(row.start(), row.end());
			// std::vector<size_t>::iterator up;
			// up= std::upper_bound (row.begin(), row.end(), 0); //
			// linbox_check(null_row.size() == (size_t)(up-row.begin()));
			size_t t1 = 0 , t2 = 0 , t3 = 0 ;
			size_t e1 = avg-1 , e2 = avg , e3 = avg+1 ;
			for (size_t i = 0 ; i < row.size() ; ++i) {
				t1 += std::min(e1,row[i]);
				t2 += std::min(e2,row[i]);
				t3 += std::min(e3,row[i]);
			}
			double r1 =  (double)t1/double(hyb.size()) ;
			double r2 =  (double)t2/double(hyb.size()) ;
			double r3 =  (double)t3/double(hyb.size()) ;
			ell = 0 ;
			if ( r1 > HYB_ELL_THRESHOLD) {
				ell = e1 ;
			}
			if (r2 > std::min(r1, HYB_ELL_THRESHOLD)) { /*! @todo benchmark me */
				ell = e2 ;
			}
			if (r3 > std::min( std::min(r2,r1), HYB_ELL_THRESHOLD) ) { /*! @todo benchmark me */
				ell = e3 ;
			}

			size_t choose_ell = 0 ; // 0 is nothing, 1 is ell, 2 is ellr
			size_t choose_coo = 0 ; // 0 is nothing, 1 is coo, 2 is csr

			if (ell != 0) {
				if (ell*hyb.rowdim() == hyb.size())
					choose_ell = 1 ;
				else
					choose_ell = 2 ;
			}

			size_t rem = (hyb.size()-ell*hyb.rowdim()) ;
			if (rem > 0) {
				double r4 = rem/hyb.rowdim();
				if (r4 < HYB_ELL_COO_THRESHOLD)
					choose_coo = 1 ;
				else
					choose_coo = 2 ;
			}

			switch(choose_ell) {
				case (1) :
					// hyb.newELL();
				case (2) :
					hyb.newELL_R();
				default :
					{}
			}

			switch(choose_coo) {
				case (1) :
					hyb.newCOO();
				case (2) :
					hyb.newCSR();
				default :
					{}
			}


			//! @todo ±1 !

		}
	};

	/** Sparse matrix, Coordinate storage.
	 *
	 * \ingroup matrix
	 * \ingroup sparse
	 */
	template<class _Field>
	class SparseMatrix2<_Field, SparseMatrix2Format::HYB > {
	public :
		typedef _Field                             Field ; //!< Field
		typedef typename _Field::Element         Element ; //!< Element
		typedef const Element               constElement ; //!< const Element
		typedef SparseMatrix2Format::HYB         Storage ; //!< Matrix Storage Format
		typedef SparseMatrix2<_Field,Storage>     Self_t ; //!< Self type
		// typedef Vector<_Field,VectorStorage::Sparse> Rep ;

		/*! Constructors.
		 * @todo convert from other matrix.
		 *
		 */
		//@{
		SparseMatrix2<_Field, SparseMatrix2Format::HYB> () :
			_rownb(0),_colnb(0) ,_nbnz(0)
			,_coo(NULL),_csr(NULL),_ell_r(NULL)
			, _field()
		{
		}

		SparseMatrix2<_Field, SparseMatrix2Format::HYB> (const _Field & F) :
			_rownb(0),_colnb(0) ,_nbnz(0)
			,_coo(NULL),_csr(NULL),_ell_r(NULL)
			, _field(F)
		{
		}

		SparseMatrix2<_Field, SparseMatrix2Format::HYB> (const _Field & F, size_t m, size_t n) :
			_rownb(m),_colnb(n) ,_nbnz(0)
			,_coo(NULL),_csr(NULL),_ell_r(NULL)
			, _field(F)
		{
		}

		SparseMatrix2<_Field, SparseMatrix2Format::HYB> (const SparseMatrix2<_Field, SparseMatrix2Format::CSR> & S) :
			_rownb(S._rownb),_colnb(S._colnb), _nbnz(S._nbnz)
			,_coo(NULL),_csr(NULL),_ell_r(NULL)
			,_field(S._field)
		{
			importe(S);

		}


		// XXX only for CSR
		template<typename _Tp1, typename _Rw1 = SparseMatrix2Format::HYB>
		struct rebind ;

		template<typename _Tp1>
		struct rebind<_Tp1/*  ,SparseMatrix2Format::COO */ > {
			typedef SparseMatrix2<_Tp1, SparseMatrix2Format::HYB> other;

			void operator() (other & Ap, const Self_t& A)
			{
				Ap.resize(A.rowdim(),A.coldim(),A.size());
				// size might decrease...

				if (A.have_reader()) {
					Ap.newReader();
					typename SparseMatrix2<_Tp1,SparseMatrix2Format::CSR>::template rebind<Field,SparseMatrix2Format::COO>()(Ap.reader(),A.reader());
				}

				if (A.have_coo()) {
					Ap.newCOO();
					typename SparseMatrix2<_Tp1,SparseMatrix2Format::COO>::template rebind<Field,SparseMatrix2Format::COO>()(Ap.coo(),A.coo());
				}

				if (A.have_csr()) {
					Ap.newCSR();
					typename SparseMatrix2<_Tp1,SparseMatrix2Format::CSR>::template rebind<Field,SparseMatrix2Format::CSR>()(Ap.csr(),A.csr());
				}

				if (A.have_ell_r()) {
					Ap.newELL_R();
					typename SparseMatrix2<_Tp1,SparseMatrix2Format::ELL_R>::template rebind<Field,SparseMatrix2Format::ELL_R>()(Ap.ell_r(),A.ell_r());
				}
			}
		};

		template<typename _Tp1, typename _Rw1>
		SparseMatrix2 (const SparseMatrix2<_Tp1, _Rw1> &S, const Field& F) :
			_rownb(S.rowdim()),_colnb(S.coldim()), _nbnz(S.size())
			,_coo(NULL)
			,_csr(NULL)
			,_ell_r(NULL)
		{
			typename SparseMatrix2<_Tp1,_Rw1>::template rebind<Field,SparseMatrix2Format::HYB>()(*this, S);
		}




		template<class VectStream>
		SparseMatrix2<_Field, SparseMatrix2Format::HYB> (const _Field & F, VectStream & stream) :
			_rownb(stream.size()),_colnb(stream.dim()), _nbnz(0)
			,_coo(NULL),_csr(NULL),_ell_r(NULL)
			, _field(F)
		{
			_csr = (new SparseMatrix2<Field,SparseMatrix2Format::CSR>(F,stream));
		}


		void resize(const size_t & mm, const size_t & nn, const size_t & zz)
		{
			_rownb = mm ;
			_colnb = nn ;
			_nbnz  = zz;
		}

		/*! Default converter.
		 * @param S a sparse matrix in any storage.
		 */
		template<class _OtherStorage>
		SparseMatrix2<_Field, SparseMatrix2Format::HYB> (const SparseMatrix2<_Field, _OtherStorage> & S) :
			_rownb(S._rownb),_colnb(S._colnb),_nbnz(S._nbnz)
			,_coo(NULL),_csr(NULL),_ell_r(NULL)
			,_field(S._field)
		{
			this->importe(S); // convert Temp from anything
		}

		~SparseMatrix2<_Field, SparseMatrix2Format::HYB>()
		{
			if(have_coo()) {
				delete _coo ;
				_coo = NULL ;
			}
			if(have_csr()) {
				delete _csr;
				_csr = NULL ;
			}
			if(have_reader()) {
				delete _reader;
				_reader = NULL ;
			}

			if(have_ell_r()) {
				delete _ell_r;
				_ell_r = NULL ;
			}
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
		void importe(const SparseMatrix2<_Field,SparseMatrix2Format::CSR> &S)
		{
			optimise();
		}

		void optimise()
		{
			Stats<Field> stats(*this);
		}

		/*! Import a matrix in CSR format to CSR.
		 * @param S CSR matrix to be converted in CSR
		 */

		void importe(const SparseMatrix2<_Field,SparseMatrix2Format::HYB> &A)
		{
			resize( A.rowdim(), A.coldim(), A.size() );

			if (A.have_coo()) {
				newCOO();
				coo().importe(A.coo());
			}
			if (A.have_csr()) {
				newCSR();
				csr().importe(A.csr());
			}
			if (A.have_ell_r()) {
				newELL_R();
				ell_r().importe(A.ell_r());
			}


		}

		/*! Export a matrix in CSR format from COO.
		 * @param S CSR matrix to be converted from COO
		 */
		SparseMatrix2<_Field,SparseMatrix2Format::HYB > &
		exporte(SparseMatrix2<_Field,SparseMatrix2Format::HYB> &S) const
		{

			// could save memory by doing that iteratively and resize to 0,0,0
			SparseMatrix2<_Field,SparseMatrix2Format::CSR> S1(field());
			SparseMatrix2<_Field,SparseMatrix2Format::CSR> S2(field());
			SparseMatrix2<_Field,SparseMatrix2Format::CSR> S3(field());
			if (have_coo()) {
				S1.importe(coo());
			}
			if (have_csr()) {
				S1.importe(csr());
			}
			if (have_ell_r()) {
				S1.importe(ell_r());
			}
			S.merge(S1);
			S.merge(S2);
			S.merge(S3);

			return S ;

		}


		//@}

		/*! In place transpose. Not quite...
		*/
		void transposeIn()
		{
			if (have_coo()) {
				coo().transposeIn();
			}
			if (have_csr()) {
				csr().transposeIn();
			}
			if (have_ell_r()) {
				ell_r().transposeIn();
			}
			if (have_reader()) {
				reader().transposeIn();
			}

		}

		/*! Transpose the matrix.
		 *  @param S [out] transpose of self.
		 *  @return a reference to \p S.
		 */
		SparseMatrix2<_Field,SparseMatrix2Format::HYB> &
		transpose(SparseMatrix2<_Field,SparseMatrix2Format::HYB> &S)
		{
			if (have_coo()) {
				S.newCoo();
				coo().transpose(S.coo());
			}
			if (have_csr()) {
				S.newCSR();
				csr().transpose(S.csr());
			}
			if (have_ell_r()) {
				S.newELLR();
				ell_r().transpose(S.ell_r());
			}
			if (have_reader()) {
				S.newCSR();
				reader().transpose(S.reader());
			}

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
			// get Entry in various opt, then reader
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
			// only setentry on _reader, destroy other ?
			// or add in coo/ wherever possible ? (seems doable)
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
			// dangerous if ref is 0, but doable.
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
				    Format fmt = SparseFileFormat::CSR())
		{
			if ( ! have_reader() )
				newReader();
			return reader().read(file,fmt);
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

			if (have_ell_r())
				ell_r().apply(y,x,field().one);
			if (have_coo())
				coo()  .apply(y,x,field().one);
			if (have_csr())
				csr()  .apply(y,x,field().one);



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

			if (have_ell_r())
				ell_r().applyTranspose(y,x,field().one);
			if (have_coo())
				coo()  .applyTranspose(y,x,field().one);
			if (have_csr())
				csr()  .applyTranspose(y,x,field().one);


			return y;
		}

		template<class Vector>
	Vector& apply(Vector &y, const Vector& x ) const
		{
			return apply(y,x,field().zero());
		}
		template<class Vector>
		Vector& applyTranspose(Vector &y, const Vector& x ) const
		{
			return apply(y,x,field().zero());
		}

		const Field & field()  const
		{
			return _field ;
		}


	protected:
		/** @todo Non element marker.
		 * We could end up a line with a marker.
		 * A field F would contain an element that does not belong to
		 * it. eg a nan for a Modular<double>. It could act as a
		 * marker.
		 */
		bool consistent() const
		{
			size_t nbnz = 0 ;
			size_t col = 0 ;
			size_t row = 0 ;
			if (have_coo()) {
				if (!coo().consistent()) return false ;
				nbnz += coo().size();
				col = std::max(col,coo().coldim());
				row = std::max(row,coo().rowdim());
			}
			if (have_csr()) {
				if (!csr().consistent()) return false ;
				nbnz += csr().size();
				col = std::max(col,csr().coldim());
				row = std::max(row,csr().rowdim());
			}
			if (have_ell_r()) {
				if (!ell_r().consistent()) return false ;
				nbnz += ell_r().size();
				col = std::max(col,ell_r().coldim());
				row = std::max(row,ell_r().rowdim());
			}
			if (nbnz != size())
				return false ;
			if (col > coldim())
				return false;
			if (row > rowdim())
				return false;

			return true ;
		}


	private :

		std::ostream & writeSpecialized(std::ostream &os,
						LINBOX_enum(Tag::FileFormat) format) const
		{
			if (have_reader())
				return reader().write(os,format);

			// SparseMatrix2<Field,SparseMatrix2Format::CSR> Temp(field());
			// this->exporte(Temp);
			// Temp.write(os,format);
			return os << "no export yet in HYB" << __FILE__ << ','  __LINE__ ;

		}


	public:



		void newCOO ()
		{
			linbox_check (!have_coo())
			_coo= new SparseMatrix2<Field,SparseMatrix2Format::COO>(_field);

		}

		void newCSR ()
		{
			linbox_check (!have_csr()) ;
			_csr= new SparseMatrix2<Field,SparseMatrix2Format::CSR>(_field);

		}

		void newELL_R ()
		{
			linbox_check (!have_ell_r()) ;
			_ell_r= new SparseMatrix2<Field,SparseMatrix2Format::ELL_R>(_field);

		}

		void newReader ()
		{
			linbox_check(!have_reader()) ;
			_reader= new SparseMatrix2<Field,SparseMatrix2Format::CSR>(_field);
		}


		SparseMatrix2<Field,SparseMatrix2Format::COO> & coo(void)
		{
			linbox_check(have_coo());
			return _coo[0] ;
		}

		bool have_coo() const
		{
			return _coo != NULL ;
		}
		bool have_csr() const
		{
			return _csr != NULL ;
		}
		bool have_ell_r() const
		{
			return _ell_r != NULL ;
		}
		bool have_reader() const
		{
			return _reader != NULL ;
		}


		SparseMatrix2<Field,SparseMatrix2Format::CSR> & csr(void)
		{
			linbox_check(have_csr());
			return _csr[0] ;
		}

		SparseMatrix2<Field,SparseMatrix2Format::ELL_R> & ell_r(void)
		{
			linbox_check(have_ell_r());
			return _ell_r[0] ;
		}

		SparseMatrix2<Field,SparseMatrix2Format::ELL_R> & reader(void)
		{
			linbox_check(have_reader());
			return _reader[0] ;
		}

		const	SparseMatrix2<Field,SparseMatrix2Format::ELL_R> & reader(void) const
		{
			linbox_check(have_reader());
			return _reader[0] ;
		}



	protected :


		size_t              _rownb ;
		size_t              _colnb ;
		size_t               _nbnz ;

		SparseMatrix2<Field,SparseMatrix2Format::COO>   * _coo  ;
		SparseMatrix2<Field,SparseMatrix2Format::CSR>   * _csr  ;
		SparseMatrix2<Field,SparseMatrix2Format::ELL_R> * _ell_r ;

		SparseMatrix2<Field,SparseMatrix2Format::CSR>   * _reader  ;

		const _Field            & _field ;

	};





} // namespace LinBox

#endif // __LINBOX_sparse_matrix_sparse_hyb_matrix_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
