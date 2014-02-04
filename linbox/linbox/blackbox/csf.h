/* linbox/blackbox/csr.h
 * Compressed Sparse Format BB
 * Author: Bryan Youse
 * overhaul of ./zo.h, a CSR formatted BB for {0,1}-Matrices
 * ------------------------------------
 *
 * Copyright (c) LinBox
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 *.
 */

#ifndef __LINBOX_CSF_H
#define __LINBOX_CSF_H

#include "linbox/integer.h"
//#include "linbox/vector/vector-traits.h"
#include "linbox/util/debug.h"
#include "linbox/util/matrix-stream.h"
#include "linbox/field/modular.h"
#include "linbox/blackbox/blackbox-interface.h"

// For STL pair in IndexIterator
#include <utility>
#include <iterator>
#include <vector> // For vectors in _col2row and _row2col
#include <cstdlib> // For randomness in randomized quicksort
#include <ctime>

namespace LinBox
{
	/** \brief Space efficient representation of sparse matrices.
	 *
	 * Compressed Sparse Row/Column (CSR/C) involves two arrays of length NNZ (number of non-zeros)
	 * One of column indices and another of matrix values.
	 * A third array contains pointers indicating the division of the two NNZ arrays into rows.
	 * (or vice versa w/r/t rows & columns)
	 *
	 \ingroup blackbox
	 */

	/*
	bool revLexLess(const std::pair<size_t,size_t>& a, const std::pair<size_t,size_t> b)
	{ return a.second < b.second || (b.second == a.second && a.first < b.first); }
	*/

	template<class _Field>
	class CSF : public BlackboxInterface {
	public:
		//  if these are ints, SuperLU can use the data directly
		//  otherwise they need converted from size_t (8byte v. 4byte)
		typedef size_t Index;//int64_t Index;
		typedef CSF<_Field> Self_t;
		typedef _Field Field;
		typedef typename _Field::Element Element;
		typedef std::vector<Index> IndexVector;
		typedef std::vector<Element> ElementVector;
		typedef IndexVector PtrVector;
		typedef std::pair<Index, Index> IndexPair;
		typedef std::pair<IndexPair, Element> Triple;
		typedef std::vector<Triple> Data;

		//to denote by which way we sort our matrix
		enum csformat { csr, csc }; //( true / false )

		// DEFAULT CONSTRUCTOR, do nothing. Matrix will be uninitialized.
		CSF(){};

		// Destructor, once again do nothing
		~CSF(){};

		// basic constructor, can be used with subsequent read.
		CSF(const Field& F) :
			_field(&F), sorted(true)
			,_rowdim(0),_coldim(0),isCSR(false)
		{}

		/* The real constructor /TODO give docs here
		  assuming entries are sorted in lexicographic order by (row,col) pair.
		  This depends on an accurate nnz being passed in.
		  Probably should change to inputting vectors sometime...
		  */
		CSF
		(Field& F, Index* rowP, Index* colP, Element* valP, Index rows, Index cols, Index nnz) :
			_field(&F), _rowdim(rows), _coldim(cols), sorted(true), isCSR(true)
		{
			Data data;
			for (Index i = 0; i < nnz; ++i, ++rowP, ++colP, ++valP)
				data.push_back(Triple(IndexPair(static_cast<Index>(*rowP), static_cast<Index>(*colP)), static_cast<Element>(*valP)));
			init(data);
		}

		/* constructor from a MatrixStream */
		CSF( MatrixStream<Field>& ms ) :
			_field(&(ms.getField()))
			,isCSR(false),sorted(false)
		{
			read(ms);
		}

		//  copy constructor, easy enough
		CSF(const CSF<Field>& A) :
			_field(A._field), _inds(A._inds), _vals(A._vals), _ptrs(A._ptrs), _rowdim(A._rowdim), _coldim(A._coldim), sorted(A.sorted)
		{ }

		// TODO
#if 0
		//switching the way in which the matrix is sorted
		void switch_sort() const
		{
			//std::cout << " -- switch_sort: " << std::endl;
			Index dim;

			if( sorted ) dim = _coldim;
			else dim = _rowdim;

			//std::cout << " -- in switch sort, before allocating temp -- " << std::endl;
			std::vector< IndexVector > temp(dim);
			//std::cout << temp.size() << std::endl;
			//IndexVector temp[dim]; //maybe this needs toooooo much memory space when dim is very large
			//std::cout << " -- in switch sort, after allocating temp -- " << std::endl;

			for( PtrVector::iterator i = _ptrs.begin(); i < _ptrs.end() - 1; ++i )
				for( IndexVector::iterator j = *i; j != *(i+1); ++j )
					temp[*j].push_back( (Index)(i - _ptrs.begin()) );

			_inds.clear(); _ptrs.clear();
			std::back_insert_iterator < std::vector<Index> > colend( _inds) ;
			for( size_t k = 0; k < dim; ++k )
			{
				_ptrs.push_back( _inds.end() );
				copy( temp[k].begin(), temp[k].end(), colend );
			}

			_ptrs.push_back( _inds.end() );
			sorted = !sorted;

			return;
		}
#endif

/*
		template<class OutVector, class InVector>
		OutVector& applyTranspose(OutVector& y, const InVector& x) const; // y = ATx
		//OutVector& applyTranspose(OutVector& y, const InVector& x); // y = ATx

		template<class OutVector, class InVector>
		OutVector& apply(OutVector& y, const InVector& x) const; // y = Ax;
		//OutVector& apply(OutVector& y, const InVector& x); // y = Ax;
*/
		template<class OutVector, class InVector>
		OutVector & apply(OutVector & y, const InVector & x) const {
			linbox_check((y.size()==rowdim())&&(x.size()==coldim()));

			FieldAXPY<Field> accum (field());

			typename OutVector::iterator yp;
			typename InVector::const_iterator xp;
			// PtrVector::const_iterator ip;

			/*
			   if( !sorted )
			   switch_sort();
			 */

			xp=x.begin();
			yp=y.begin();
			accum.reset();

			for(Index i = _ptrs[0]; (size_t)i < _ptrs.size()-1; ++i, ++yp) {
				for(Index j = _ptrs[i]; j < _ptrs[i+1]; ++j) {
					accum.mulacc(_vals[j], *(xp + _inds[j]) ); // y = a*x
				}
				accum.get(*yp);
				accum.reset();
			}

			return y;
		}

		template<class OutVector, class InVector>
		OutVector & applyTranspose(OutVector & y, const InVector & x) const {
			linbox_check((y.size()==coldim())&&(x.size()==rowdim()));

			for(size_t i = 0; i < y.size(); ++i) y[i] = field().zero;

			for(Index i = _ptrs[0]; (size_t)i < _ptrs.size()-1; ++i) {
				for(Index j = _ptrs[i]; j < _ptrs[i+1]; ++j) {
				// process row i:  yj += xi Aij , yindsj += xi valsj
					field().axpyin(y[i], x[_inds[j]], _vals[j]);
				}
			}

			return y;
		}

		Element & getEntry(Element& x, Index i, Index j) {
			size_t k;
			for (k = _ptrs[i]; k < _ptrs[i+1]; ++k)
				if (_inds[k] == j) break;
			if (k == _ptrs[i+1]) return field().assign(x, field().zero);
			else return field().copy(x, _vals[k]);
		}

		Element & setEntry(Index i, Index j, Element& x) {
			// data must exist.
			data.push_back(i, j, x);
		}

		void finalize() { // from data to csf
			init(_data);
		}

		double &d00norm(double &norm){
			norm = 0;
			double  t;
			//  maximal row inf (OO) norm
			for(Index i = _ptrs[0]; (size_t)i < _ptrs.size()-1; ++i) {
				double old;
				old = norm;
				for(Index j = _ptrs[i]; j < _ptrs[i+1]; ++j) {
					field().convert(t, _vals[j]);
					norm += abs(t);
				}
				if (norm < old) norm = old;
			}

			return norm;
		}

		integer &hadamardBound(integer &res){
			res = 1L;
			integer tmp;

			//  product of ||A_i||_2 norms for rows A_i
			for(Index i = _ptrs[0]; (size_t)i < _ptrs.size()-1; ++i) {
				tmp = 0;
				for(Index j = _ptrs[i]; j < _ptrs[i+1]; ++j) {
					tmp += static_cast<integer>(_vals[j]) * _vals[j];
				}
				res *= tmp;
			}

			res = sqrt (res);
			return res;
		}

		/** Read the matrix from a matrix stream
		 *  @param ms Stream from which to read the matrix
		 */
		void read(MatrixStream<Field> &ms)
		{
			Data d;
			size_t r, c;
			Element v;

			int64_t count = 0;

			ms.getDimensions( _rowdim, _coldim );

			while (ms.nextTriple(r, c, v) )
			{
				d.push_back(Triple(IndexPair(static_cast<Index>(r), static_cast<Index>(c)), static_cast<Element>(v)));
				++count;
			}

			init(d);
		}

		/*  Accessor methods */
		size_t rowdim() const
		{ return _rowdim; }

		size_t coldim() const
		{ return _coldim; }

		// TODO generecize for csc
		IndexVector &getRows()
		{
			return _ptrs;
		}

		IndexVector &getCols()
		{
			return _inds;
		}

		ElementVector &getVals()
		{
			return _vals;
		}

		std::ostream& write_summary(std::ostream& out =  std::cout) const
		{
			out << "CSF Matrix: _inds.size() " << _inds.size();
			out << ", _ptrs.size() " << _ptrs.size();
			out << ", _rowdim " << _rowdim;
			out << ", _coldim " << _coldim;
			return out;
		}

		//  helper to write out triples, essentially sms w/o header/footer
		std::ostream& write_sms(std::ostream& out =  std::cout) const {
			Index row = 0;
			integer val;
			for(Index i = 0; i < _ptrs.size() - 1; ++i, ++row)
				for(Index j = _ptrs[i]; j < _ptrs[i+1]; ++j){
					field().convert(val, _vals[j]);
					out << row << " " << _inds[j] << " " << val << std::endl;
				}

			return out;
		}

		// TODO generecize for CSC format
		std::ostream& write(std::ostream& out =  std::cout) const {
			integer val;
			//  for each row
			for(Index i=0, k=0; (size_t)i < _ptrs.size() - 1; ++i){
				k = 0;
				out << "  [";
				//  j will be the index in _inds and _vals of data
				for(Index j = _ptrs[i]; j < _ptrs[i+1]; ++j){
					//  print zeros up to data
					for(; k<_inds[j]; ++k) out << " 0";
					field().convert(val, _vals[j]);
					//  print data
					out << " " <<  val;
					++k;  //  adjust zero printing counter by one
				}
				for(; (size_t)k<_coldim; ++k) out << " 0";  // print zeros to end
				out << " ]";
			}
			out << std::endl;
			return out;
		}

		const Field& field() const
		{
			return *_field;
		}

		/* Non blackbox function.  Tells the number of nonzero entries */
		size_t nnz() const
		{
			return _inds.size();
		};

		typedef MatrixCategories::BlackboxTag MatrixCategory;

// TODO come back to these things later.
#if 0
		template<typename _Tp1>
		struct rebind
		{
			typedef CSF<_Tp1> other;
			void operator() (other *& Ap,
					 const Self_t& A,
					 const _Tp1& F) {
				Ap = new other(F, A._inds, A._ptrs, A._rowdim, A._coldim, A.sorted);
			}
		};
#endif

	protected:
		const Field *_field; // The field used by this class
		IndexVector _inds; // The nnz indices sorted by row or by col
		ElementVector _vals; // The values corresonding to nnz indices
		PtrVector _ptrs; // the pointers to beginning of each (row/col)
		size_t _rowdim, _coldim;

		// if data is empty, matrix is ready to use.
		Data _data;
		bool sorted, isCSR;

		class sort_data_by_col{
			private:
				Data _d;
			public:
				sort_data_by_col(Data d) : _d(d) {}
				//  col sorted
				bool operator()(const Triple &a, const Triple &b){
					return (a.first.second < b.first.second ||
						(a.first.second == b.first.second && a.first.first < b.first.first));
				}
		};

		//  INITIAL
		void init(Data& d)
		{
			Index nnz = d.size();
			sort(d.begin(), d.end());
			//sort(d.begin(), d.end(), sort_data_by_col(d));

			/*
			typename Data::iterator di = d.begin();
			for(;di != d.end(); ++di){
				std::cout << (*di).first.first << " " << (*di).first.second << " " << (*di).second << std::endl;
			}
			*/

			// set up _inds and _vals
			// if CSR, matrix cols represented by _inds and rows rep. by _ptrs
			for (Index i = 0; i < nnz; ++i){
				if(isCSR) _inds.push_back(d[i].first.second);
				else _inds.push_back(d[i].first.first);
				_vals.push_back(d[i].second);
			}

			// p represents position in the index/value arrays
			Index p = 0;

			// set up _ptrs
			_ptrs.push_back(p);

			typename Data::iterator q = d.begin();
			//  for all the data we have
			for (Index i = 0; q != d.end(); ++q, ++p){
				//  if we see a new row
				if (i != q->first.first) {
					//  loop to encapsulate all possible "zero rows"
					for (Index j = i; j < q->first.first; j++)
						_ptrs.push_back(p); //add the index for the new row
					i = q->first.first;  // now we're searching against the new row
				}
			}
			// lastly, a pointer to AFTER last elt.
			_ptrs.push_back(nnz);

		}  // init()

#if 0
			// TODO decide between 2 copies and sort switching.
			 keep another copy is not needed if we can switch sort between row and col
			// sort by cols first, then do the same as above
			sort(ip.begin(), ip.end(), revLexLess);

			// set up _row
			for (Index i = 0; i < NNz; ++i)
			_row.push_back(ip[i].first);

			// set up _colP
			p = _row.begin();
			_colP.push_back(p);

			q =ip.begin();
			i = q->second;
			p++;q++;

			for ( ; q != ip.end(); ++q, ++p)
			if (i != q->second)
			{
			for (Index j = i; j < (q+1)->second; j++)
			_colP.push_back(p);
			i = q->second;
			}
			_colP.push_back(_row.end());
		}
#endif

	}; //CSF

}//End of LinBox

//#include "csf.inl"

#endif // __LINBOX_CSF


// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
