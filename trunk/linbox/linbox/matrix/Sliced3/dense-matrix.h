/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* linbox/matrix/dense-matrix.h
 * Copyright (C) 2011 B. David Saunders,
 * See COPYING for license information
 *
 * evolved from dense-submatrix and blas-matrix by B. David Saunders <saunders@cis.udel.edu>, BSY
 */

/*! @file matrix/Sliced3/dense-submatrix.h
 * @ingroup matrix
 * @brief Representation of a submatrix of a dense matrix, not resizeable.
 * This matrix type conforms to the \c LinBox::DenseMatrixBase interface.
 * \c LinBox::BlasMatrix is an example of DenseSubmatrix.
 */

#ifndef __LINBOX_dense_matrix_H
#define __LINBOX_dense_matrix_H

#include <utility>
#include "linbox/linbox-config.h"

//#include "linbox/util/debug.h"
//#include "linbox/matrix/matrix-domain.h"
#include "linbox/vector/subvector.h"
#include "submat-iterator.h"

namespace LinBox
{

	/** @brief to be used in standard matrix domain
	
	 * Matrix variable declaration, sizing, entry initialization may involve one to 3 steps.
	 * Matrix ops are container ops. (sizing, copying)  
	 * Mathematically meaningful operations are to be found only in an associated matrix domain
	 *
	 * A matrix may be allocated or not.  A matrix initialized by a submatrix() call is not allocated.
	 * When an allocated matrix goes out of scope or is reinitialized with init(), the memory is released 
	 * and all submatrices of it become invalid.
	 *
	 * Allocating:
	 * DenseMatrix A(2, 3); // allocation of mem for 6 entries at construction
	 * DenseMatrix B; B.init(10, 10); // default constr and subsequent allocation for 100 entries.
	 *
	 * Allocation of memory plus entry initialization:
	 * // a meaningful value of DenseMatrix::Entry x is set by a field or matrix domain.
	 * DenseMatrix A(10, 10, x); 
	 * DenseMatrix B: B.init(10, 10, x); 
	 * DenseMatrix C(A); // allocation at copy construction.  A could be a submatrix of another.
	 * A.read(istream)
	 *
	 * Nonallocation sizing:
	 * // assume D declared, A initialized, n = A.coldim().
	 * D.submatrix(A, 0, 1, 2, n-1); // D is 2 by n-1 in upper right of A.
	 *
	 * Entry initialization (and overwriting) in already sized matrices:
	 * A.setEntry(i, j, x);
	 * A = B; // A and B must have the same shape.

	 * Entry read access. OK on const matrices
	 * getEntry, write

	 * Under consideration:  Require A.clear() on an allocated matrix before any action that would abandon
	 * the allocated mem (init or submatrix).
	 \ingroup matrix
	 */
	template<class _Element>
	class DenseMatrix {
		typedef _Element Entry;
	public: // protected:  //TODO scope
		Entry *_rep; // matrix entries on the heap.
		bool _alloc; // iff _alloc, I am responsible for _rep allocation/deallocation.
		size_t _rows; 
		size_t _cols;
		size_t _stride; // stride from row to row.  Entries in a row are contiguous.

	public:
		size_t rowdim() const { return _rows; } 
		size_t coldim() const { return _cols; }

		Entry& getEntry(Entry& x, size_t i, size_t j) const { 
			linbox_check(i < _rows && j < _cols);
			return x = *(_rep + i*_stride + j); 
		}

		Entry& setEntry(size_t i, size_t j, const Entry& x ) { 
			//linbox_check((i < _rows && j < _cols));
			return *(_rep + i*_stride + j) = x; 
		}
		// no refEntry - not well supportable in some matrix domains.

		DenseMatrix() : _rep(NULL), _alloc(false), _rows(0), _cols(0), _stride(0) {}

		void init(size_t m = 0, size_t n = 0) {
			//std::cerr << m << " " << n << " <<<<<<" << std::endl;
			if (_alloc) delete _rep; // abandon any prior def
			if (m*n != 0) { 
				_rep = new Entry[m*n]; _alloc = true;
			} else { 
				_rep = NULL; _alloc = false; 
			}
			size(m, n);
		}

		void init(size_t m, size_t n, const Entry& filler) {
			if (_alloc) delete _rep; // abandon any prior def
			if (m*n != 0) { 
				_rep = new Entry[m*n]; _alloc = true;
				for (size_t i = 0; i < m*n; ++i) _rep[i] = filler;
			} else { 
				_rep = NULL; _alloc = false; 
			}
			size(m, n);
		}

		/* Copy construction makes this a completely distinct copy of A.  
		 * Entries of this are stored contiguously even if A is not contiguous (is a submatrix of another).
		 * 
		 * If memcpy is not valid for your entry type, specialize DenseMatrix for it.
		 */
		DenseMatrix(const DenseMatrix& A)
		: _rep(new Entry[A._rows*A._cols]), _alloc(true), _rows(A._rows), _cols(A._cols), _stride(A._cols) 
		{	
			//std::cout << "copy construction " << _rep << std::endl;
			//std::cout << "copy cons " << _rows << " " << _cols << std::endl;
			*this = A; // copy A's data into _rep
		}

		~DenseMatrix() {
			if (_alloc) delete _rep;
		}

		// For assignment, the matrices must have same size and not overlap in memory.
		// This restriction could be loosened...
		DenseMatrix& operator=(const DenseMatrix& B) {
			linbox_check(_rows == B._rows && _cols == B._cols);
			if (_cols == _stride && B._cols == B._stride) // both are contiguous
				memcpy(_rep, B._rep, sizeof(Entry)*_rows*_cols);
			else
				for (size_t i = 0; i < _rows; ++i) // copy row by row
					memcpy(_rep + i*_stride, B._rep + i*B._stride, sizeof(Entry)*_cols);
			return *this;	
		}
		
		/** Set this to be an m by n submatrix of A with upper left corner at i,j position of A.
		 * Requires i+m <= A.rowdim(), j+n <= A.coldim(). 
		 *
		 * For instance, B.submatrix(A, i, 0, 1, A.coldim()) makes B the i-th row of A.
		 */
		void submatrix(const DenseMatrix & A, size_t i, size_t j, size_t m, size_t n) {
			linbox_check(i+m <= A._rows && j+n <= A._cols);
			if (_alloc) delete _rep; // abandon any prior def
			_rep = A._rep + (i*A._stride + j);
			_alloc = false;
			_rows = m; _cols = n; _stride = A._stride;
		}

		/*  SLICED helper functions... maybe shouldn't be public */

		//  size a matrix w/o allocating memory
		void size(size_t m, size_t n){
			_rows = m; _cols = _stride = n;
		}

		template<class vp>
		int did_swapback(vp &swaps, size_t hot, size_t l, size_t r){
			for(size_t i=l; i<r; ++i){
				if(swaps[i].second == hot)
					return swaps[i].first;
			}
			return -1;
		}

		template<class vp>
		void flocs(vp & tos, vp & swaps){
			std::vector<size_t> pos;
			for(size_t i =0; i < coldim(); ++i)  pos.push_back(i);

			for(int i = coldim()-1; i>= 0; --i){
				//std::cerr << swaps[i] << std::endl;
				swap(pos[swaps[i].first], pos[swaps[i].second]);
			}

			for(size_t i = 0; i < coldim(); ++i){
				std::pair<size_t, size_t> p(i, pos[i]);
				tos.push_back(p);
			}
		}

		template<class vp>
		size_t final_loc(vp &swaps, size_t j){
			//  first check if we were swapped back prior to being reached
			int to = did_swapback(swaps, j, 0, j+1);
			if(to >= 0) return to;
			//  next check if we get swapped back from our immediate neightbor
			//size_t nbor = swaps[j].second;
			//to = did_swapback(swaps, nbor, j+1, nbor);
			//if(to >= 0) return to;
			//  bad luck, we are chained forward, so 
			size_t i = j;
			while(to < 0){
				size_t lc = swaps[i].first;
				size_t rc = swaps[i].second;
				to = did_swapback(swaps, rc, lc+1, rc+1);
				//std::cerr << j << ": " << to << "= ds(swps, " << rc << "," << lc+1 << "," << rc+1 << ")" << std::endl;

				i = rc;
				if(to < 0 && rc == swaps.size() - 1){
					to = swaps.size() - 1;
				}
			}
			return to;
		}

		void randomColPermutation() {
			typedef std::pair<size_t, size_t> pair;
			typedef std::vector<pair> vp;
			vp swaps;
			vp tos;
			vp tos2;
			//  create pairs
			for (size_t j = 0; j < coldim(); ++j){
				// Each iteration swap col j with a random col in range [j..n-1].
				int k = j + rand()%(coldim()-j);
				//std::cerr << j << "->" << k << std::endl;
				pair p(j,k);
				swaps.push_back(p);
				//for (size_t i = 0; i < rowdim(); ++i)
					//swap( _rep[i*_stride + j], _rep[i*_stride + k]);
			}
			/*
			//flocs(tos2, swaps);
			//  find each final location
			for(size_t j = 0; j < coldim(); ++j){
				size_t to = final_loc(swaps, j);
				pair p(j, to);
				if(p != tos2[j]){
					std::cerr << "HOUSTON: PROBLEM" << std::endl;
				}
				else{
					std::cout << "YAY ";
				}
				tos.push_back(p);
			}
				*/
			flocs(tos, swaps);
			//  tos is the correct mapping now permute
			Entry *perm_row = new Entry[coldim()];
			for(size_t i = 0; i < rowdim(); ++i){
				for(vp::iterator ti = tos.begin(); ti!= tos.end(); ++ti)
					perm_row[(*ti).second] = _rep[i*_stride + (*ti).first];
				memcpy(&(_rep[i*_stride]), perm_row, coldim()*sizeof(Entry));	
			}
			delete[] perm_row;
			
			//  CHECKING CODE
			/*
			std::cerr << "swaps: ";
			for(vp::iterator i = swaps.begin(); i!= swaps.end(); ++i){
				std::cerr << *i << " ";
			}
			std::cerr << std::endl;
			std::cerr << "tos: ";
			for(vp::iterator i = tos.begin(); i!= tos.end(); ++i){
				std::cerr << *i << " ";
			}
			std::cerr << std::endl;
			std::vector<int> orig;
			for(size_t i = 0; i < coldim(); ++i) orig.push_back((int)i);
			std::vector<int> out(orig);
			std::vector<int> out2(orig);
			for(vp::iterator i = swaps.begin(); i!= swaps.end(); ++i){
				std::swap(out[(*i).first], out[(*i).second]);
			}
			for(vp::iterator i = tos.begin(); i!= tos.end(); ++i){
				out2[(*i).second] = orig[(*i).first];
			}
			for(std::vector<int>::iterator i=out.begin(), i2=out2.begin(); i!=out.end(); ++i, ++i2){
				if(*i != *i2){
					std::cerr << "LOSER: ";
					std::cerr << *i << " v " << *i2 << std::endl;
				}
			}
			*/
		}

		void randomLowerTriangularColTransform() {
			for (size_t j = 0; j < coldim(); ++j){
				// Each iteration swap col j with a random col in range [j..n-1].
				//int l = 1, k = 0;
				//do { l *= RAND_MAX; k = k*RAND_MAX + rand(); } while (l < j);
				//k = k%j;
				int k = rand()%j;
				for (size_t i = 0; i < rowdim(); ++i)
					 _rep[i*_stride + j] += _rep[i*_stride + k];
			}
		}
		/*  Iterators */

		typedef SubMatIterator<_Element> RawIterator;
		typedef ConstSubMatIterator<_Element> ConstRawIterator;

		RawIterator rawBegin() { 
			return RawIterator(_rep, _cols, _stride); };
		RawIterator rawEnd() { 
			return RawIterator(_rep + _rows*_stride); }
		RawIterator rowBegin(size_t i) { 
			//return rawBegin() + i*_cols;
			return RawIterator(_rep+i*_stride, _cols, _stride);
		}
		RawIterator rowEnd(size_t i) { 
			return rowBegin(i + 1); 
		}
		/*
		RawIterator rowBegin(size_t i) { 
			return rawBegin() + i*_cols; }
		RawIterator rowEnd(size_t i) { 
			return rowBegin(i + 1); }
			*/

		ConstRawIterator rawBegin() const { 
			return ConstRawIterator(_rep, _cols, _stride); }
		ConstRawIterator rawEnd() const { 
			return ConstRawIterator(_rep + _rows*_stride); }
		ConstRawIterator rowBegin(size_t i) const { 
			return rawBegin() + i*_cols; }
		ConstRawIterator rowEnd(size_t i) const { 
			return rowBegin(i + 1); }

		typedef DenseMatrix<_Element>   Self_t;       //!< Self type

#if 0  // TODO factor out
		/** @name typedef'd Row Iterators.
		 *\brief
		 * The row iterator gives the rows of the
		 * matrix in ascending order. Dereferencing the iterator yields
		 * a row vector in dense format
		 * @{
		 */
		typedef Entry* RowIterator;
		typedef const Entry* ConstRowIterator;
		typedef Subvector<Entry*> Row;
		typedef Subvector<const Entry*> ConstRow;
		 //@} Row Iterators

		/** @name typedef'd Column Iterators.
		 *\brief
		 * The columns iterator gives the columns of the
		 * matrix in ascending order. Dereferencing the iterator yields
		 * a column vector in dense format
		 * @{
		 */
		typedef Subiterator<Entry> ColIterator;
		typedef Subiterator<const Entry> ConstColIterator;
		typedef Subvector<ColIterator>                   Col;
		typedef Subvector<ConstColIterator>                   ConstCol;
		//@} // Column Iterators
#endif

		template<typename _Tp1>
		struct rebind {
			typedef DenseMatrix<typename _Tp1::Element> other;
		};

		/** Read the matrix from an input stream.
		 * @param file Input stream from which to read
		 * @param field
		 */
		template<class Field>
		std::istream& read (std::istream &file, const Field& field);

		/** Write the matrix to an output stream.
		 * @param os Output stream to which to write
		 * @param field
		 * @param mapleFormat write in Maple(r) format ?
		 */
		template<class Field>
		std::ostream& write (std::ostream &os, const Field& field,
				     bool mapleFormat = false) const;

		/** Write the matrix to an output stream.
		 * This a raw version of \c write(os,F) (no field is given).
		 * @param os Output stream to which to write
		 * @param mapleFormat write in Maple(r) format ?
		 */
		std::ostream& write (std::ostream &os,
				     bool mapleFormat = false) const;

		/*  TODO factor out row/col iterator
		RowIterator rowBegin ();
		RowIterator rowEnd ();
		ConstRowIterator rowBegin () const;
		ConstRowIterator rowEnd () const;

		ColIterator colBegin ();
		ColIterator colEnd ();
		ConstColIterator colBegin () const;
		ConstColIterator colEnd () const;
		*/

};
	/*! Write a matrix to a stream.
	 * The C++ way using <code>operator<<</code>
	 * @param o output stream
	 * @param Mat matrix to write.
	 */
	template<class T>
	std::ostream& operator<< (std::ostream & o, const DenseMatrix<T> & Mat)
	{
		return Mat.write(o);
	}

	/*! @internal
	 * @brief MatrixTraits
	 */
	/*  necessary?
	template <class Element>
	struct MatrixTraits< DenseMatrix<Element> > {
		typedef DenseMatrix<Element> MatrixType;
		typedef typename MatrixCategories::RowColMatrixTag MatrixCategory;
	};
	*/

	template <class _Element>
	template <class Field>
	std::istream& DenseMatrix<_Element>::read (std::istream &file, const Field& field)
	{
		RawIterator p;

		for (p = rawBegin (); p != rawEnd (); ++p) {
			// each entry is seperated by one space.
			file.ignore (1);
			field.read (file, *p);
		}

		return file;
	}

	template <class _Element>
	template <class Field>
	std::ostream &DenseMatrix<_Element>::write (std::ostream &os, const Field& field,
						       bool mapleFormat) const
	{
		return os;
	}
/*  TODO factor out RowIterator
		ConstRowIterator p;

		// integer c;
		//int wid;

		// field.cardinality (c);
		//wid = (int) ceil (log ((double) c) / M_LN10); //BB : not used !

		typename ConstRow::const_iterator pe;

		if (mapleFormat) os << "[";

		for (p = rowBegin (); p != rowEnd (); ++p) {
			if (mapleFormat && (p != rowBegin()))
				os << ',';
			if (mapleFormat) os << "[";

			for (pe = p->begin (); pe != p->end (); ++pe) {
				if (mapleFormat && (pe != p->begin())) os << ',';
				// matrix base does not provide this field(), maybe should?
				//_M.field ().write (os, *pe);
				//os << *pe;
				//fixed by using extra field

				field.write (os, *pe);
				os << " ";
			}

			if (!mapleFormat)
				os << std::endl;
			else os << ']';
		}

		if (mapleFormat) os << ']';
		return os;
	}
	*/

	template <class _Element>
	std::ostream &DenseMatrix<_Element>::write (std::ostream &os, bool mapleFormat) const
	{
		return os;
	}
	/*
		ConstRowIterator p;



		typename ConstRow::const_iterator pe;

		if (mapleFormat) os << "[";

		for (p = rowBegin (); p != rowEnd (); ++p) {
			if (mapleFormat && (p != rowBegin()))
				os << ',';
			if (mapleFormat) os << "[";

			for (pe = p->begin (); pe != p->end (); ++pe) {
				if (mapleFormat && (pe != p->begin())) os << ',';

				os << *pe;
				os << " ";
			}

			if (!mapleFormat)
				os << std::endl;
			else os << ']';
		}

		if (mapleFormat) os << ']';
		return os;
	}
	*/

} // namespace LinBox

#endif // __LINBOX_dense_matrix_H

