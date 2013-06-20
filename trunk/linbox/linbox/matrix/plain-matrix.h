/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* linbox/matrix/plain-matrix.h
 * -bds 2013 
 * See COPYING for license information
 *
 * evolved from Dense-submatrix and blas-matrix 
 */

/*! @file matrix/plain-matrix.h
 * @ingroup matrix
 * @brief Reference representation of a PlainMatrix (dense, memory allocating) class
 * and PlainSubmatrix (dense, non-allocating) class.
 * \c other dense submatrix classes such as LinBox::BlasSubmatrix exihibit this functionality.
 */

#ifndef __LINBOX_plain_matrix_h
#define __LINBOX_plain_matrix_h

#include "linbox/util/error.h"

#include "linbox/util/debug.h"

#include "linbox/matrix/plain-domain.h"

namespace LinBox
{

	/** @brief to be used in reference matrix domain (PlainDomain).
	
	 * Matrix variable declaration, sizing, entry initialization may involve one to 3 steps.
	 * Matrix ops are container ops. (sizing, copying)  
	 *
	 * Mathematical operations are to be found only in an associated matrix domain ).
	 * (exceptions are some use of domain scalars in, eg., zero(), random(), setEntry(), getEntry().
	 *
	 * A Submatrix does not allocate heap memory.  It shares (subset of) the memory of a (memory allocating) DenseMatrix.
	 * When a DenseMatrix goes out of scope or is reinitialized with init(), the memory is released 
	 * and all Submatrices of it become invalid.
	 *
	 * Allocating:
	 * Given a matrix domain, MatDom MD,

	 * MatDom::Matrix A(MD, 2, 3); // allocation of mem for 6 entries at construction
	 * MatDom::Matrix B; B.init(MD, 10, 10); // default constr and subsequent allocation.
	 *
	 * Allocation of memory plus entry initialization:
	 * // a meaningful value of DenseMatrix::Entry x is set by a field.
	 * MatDom::Matrix B(A); // allocation at copy construction.  A could be a submatrix of another.
	 * MatDom::Matrix A; A.read(stream); // allocation at read time.
	 * MatDom::Submatrix A(MD, n, m); A.read(stream); // no allocation at read time. Shape must match.
	 *
	 * Nonallocation sizing:
	 * MatDom::Submatrix S,T;
	 * S.submatrix(A, 1, 0, 2, A.coldim()); // S is second 2 rows of A
	 * T.submatrix(S, 0, S.coldim()-2, 2, 2); // T is 2by2 at right end of S, shares mem with S and A.
	 *
	 * Entry initialization (and overwriting) in already sized matrices:
	 * S.setEntry(i, j, x);
	 * S.copy(B); S = B; // A and B must have the same shape.
	 * S.read(stream); // A and matrix in stream must have the same shape.

	 * Entry read access. OK on const matrices
	 * S.getEntry(x,i,j), S.write(stream)

	 \ingroup matrix
	 */
	 /* defined here:
	 class PlainSubmatrix<Dom>;
	 class PlainMatrix<Dom>: PlainSubmatrix<Dom>;
	 */

	template<class MatDom>
	class PlainSubmatrix /*: public DenseMatrixInterface<MatDom>*/ {
	public:
		typedef PlainSubmatrix<MatDom> Self_t;
		typedef MatDom MatrixDomain;
		typedef size_t Index;
		typedef typename MatrixDomain::Scalar Entry;
	protected:
		Entry *rep_; // matrix entries on the heap.
		const MatrixDomain* domain_; // scalar, vector, matrix arithmetic context
		Index rows_; 
		Index cols_;
		Index row_stride_; // stride from row to row.  Entries in a row are contiguous.
	public:
		Index rowdim() const 
		{ return rows_; } 
		Index coldim() const 
		{ return cols_; }
		inline const MatrixDomain& domain() const 
		{ return *domain_; }
		const MatrixDomain& field() const // transitional
		{ return *domain_; }

		Entry& getEntry(Entry& x, Index i, Index j) const 
		{ return x = rep_[i*row_stride_ + j]; }
		void setEntry(Index i, Index j, const Entry& x ) 
		{ rep_[i*row_stride_ + j] = x; }

		void submatrix(const Self_t & A, Index i, Index j, Index m, Index n) 
		{	rep_ = A.rep_ + i*row_stride_ + j;
			rows_ = m; cols_= n, row_stride_ = A.row_stride_;
			domain_ = A.domain(); 
		}

		Self_t& zero() // set to zeroes, no shape change
		{	for (Index i = 0; i < rowdim(); ++ i)
				for (Index j = 0; j < coldim(); ++ j)
					setEntry(i, j, field().zero);
		}
		Self_t& identity() // set to I, must be square, no shape change
		{	this->zero();
			for (Index i = 0; i < rowdim(); ++i)
				setEntry(i, i, field().one);
		}
		Self_t& random() // set to random entries, no shape change
		{	Entry x; field().init(x,0);
			typename MatrixDomain::RandIter r(field());
			for (Index i = 0; i < rowdim(); ++ i)
				for (Index j = 0; j < coldim(); ++ j)
					setEntry(i, j, r.random(x));
		}

		std::istream& read (std::istream &is) // The matrix read must have the same shape.
		{ throw(LinboxError("no PlainSubmatrix read yet")); return is; }
		std::ostream& write (std::ostream &os) const
		{ throw(LinboxError("no PlainSubmatrix write yet")); return os; }

		PlainSubmatrix() :rep_(0), rows_(0), cols_(0), row_stride_(1)
		{}
		PlainSubmatrix(const Self_t& A) 
		{ submatrix(A, 0, 0, A.rowdim(), A.coldim()); }

		Self_t& copy(const Self_t& B) // deep copy
		{	linbox_check(rowdim() == B.rowdim() and coldim() == B.coldim());
			linbox_check(&(domain()) == &(B.domain()));
			Entry x; domain().init(x, 0);
			for (Index i = 0; i < rowdim(); ++ i)
				for (Index j = 0; j < coldim(); ++ j)
					setEntry(i, j, B.getEntry(x, i,j));
			return *this;
		}
		Self_t& operator=(const Self_t& B) 
		{	return copy(B); }

		// can have trace, rank, det, etc.

	}; // PlainSubmatrix

	template<class Domain_>
	class PlainMatrix : public PlainSubmatrix<Domain_> {
	public:
		typedef PlainMatrix<Domain_> Self_t;
		typedef PlainSubmatrix<Domain_> Father_t;
		typedef typename Father_t::MatrixDomain MatrixDomain;
		typedef typename Father_t::Index Index;
		typedef typename Father_t::Entry Entry;
	protected:
		using Father_t::domain_;
		using Father_t::rep_;
		using Father_t::rows_;
		using Father_t::cols_;
		using Father_t::row_stride_;
	public:
		PlainMatrix() 
		: Father_t() {}
		//~PlainMatrix() 
		//{} // default dstor works
		PlainMatrix(const MatrixDomain& D, Index m, Index n)
		: Father_t() { init(D, m, n); }
		PlainMatrix(const PlainMatrix& A) // deep copy
		: Father_t() { init(A.domain(), A.rowdim(), A.coldim()); this->copy(A); }

		void init(const MatrixDomain& D, Index m, Index n)
		{	domain_ = &D;
			if (rows_*cols_ != m*n) // must replace current mem.
			{	if (rep_) delete(rep_); 
				rep_ = new(Entry[m*n]);  
			}
			rows_ = m; cols_ = n; row_stride_ = 1; 
		}
		PlainMatrix& operator=(const PlainMatrix& A) // deep copy, reallocate if necessary
		{ init(A.domain(), A.rowdim(), A.coldim()); copy(A); }

		std::istream& read (std::istream &is) // Will be reshaped to match the matrix read.
		{ throw(LinboxError("no PlainMatrix read yet")); return is; }

		using Father_t::rowdim;
		using Father_t::coldim;
		using Father_t::field;
		using Father_t::domain;

		using Father_t::getEntry;
		using Father_t::setEntry;
		using Father_t::submatrix;
		using Father_t::zero;
		using Father_t::random;
		using Father_t::write;
		using Father_t::copy;
	}; //PlainMatrix
} //LinBox 
#endif //__LINBOX_plain_matrix_h
