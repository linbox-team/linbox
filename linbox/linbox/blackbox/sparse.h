/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/blackbox/sparse.h
 * Copyright (C) 1999-2001 William J Turner,
 *               2001-2002 Bradford Hovinen
 *
 * Written by William J Turner <wjturner@math.ncsu.edu>,
 *            Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * ------------------------------------
 * 2002-08-06  Bradford Hovinen  <hovinen@cis.udel.edu>
 *
 * Renamed to sparse.h from sparse0.h
 * ------------------------------------
 * Modified by Bradford Hovinen <hovinen@cis.udel.edu>
 * Modified by Jean-Guillaume Dumas <Jean-Guillaume.Dumas@imag.fr>
 * 	       28.08.2002 : added back : field()
 *
 * Refactoring:
 *   - Eliminated SparseMatrixAux and moved that functionality into Sparse0
 *   - Made SparseMatrixBase parameterized only on the element type
 *   - New read/write implementations for SparseMatrixBase, supporting multiple
 *     formats
 *   - Eliminated Gaussian elimination code
 *   - Added iterators, including ColOfRowsIterator, RawIterator, and
 *     RawIndexIterator
 *   - Eliminated operator []; added getEntry; changed put_value to setEntry
 * ------------------------------------
 * 
 * See COPYING for license information.
 */

#ifndef __BLACKBOX_SPARSE_H
#define __BLACKBOX_SPARSE_H

#include "linbox/blackbox/archetype.h"
#include "linbox/matrix/sparse.h"
#include "linbox/field/vector-domain.h"
#include "linbox/vector/vector-traits.h"
#include "linbox/vector/stream.h"
#include "linbox/util/field-axpy.h"

// Namespace in which all LinBox library code resides
namespace LinBox
{

/** Sparse matrix
 * This is a generic black box for a sparse matrix. It inherits
 * \ref{SparseMatrixBase}, which implements all of the underlying
 * accessors and iterators.
 */
template <class Field,
	  class _Vector = typename LinBox::Vector<Field>::Dense,
	  class _Row    = typename LinBox::Vector<Field>::Sparse,
	  class Trait   = typename VectorTraits<_Vector>::VectorCategory>
class SparseMatrix : public SparseMatrixBase<typename Field::Element, _Row>, public BlackboxArchetype<_Vector>
{
    public:

	typedef _Vector Vector;
	typedef typename Field::Element Element;
	typedef typename SparseMatrixBase<typename Field::Element, _Row>::Row Row;
	typedef typename SparseMatrixBase<typename Field::Element, _Row>::Format Format;
	typedef typename SparseMatrixBase<typename Field::Element, _Row>::RawIterator RawIterator;
	typedef typename SparseMatrixBase<typename Field::Element, _Row>::RawIndexedIterator RawIndexedIterator;

	/** Constructor.
	 * Builds a zero m x n matrix
	 * Note: the copy constructor and operator= will work as intended
	 *       because of STL's container design
	 * @param  F  Field over which entries exist
	 * @param  m  Row dimension
	 * @param  n  Column dimension
	 */
	SparseMatrix (const Field &F, size_t m, size_t n);

	/** Constructor from a vector stream
	 * @param  F  Field over which entries exist
	 * @param  stream  Stream with which to generate row vectors
	 */
	SparseMatrix (const Field &F, VectorStream<Row> &stream); 

	/** Copy constructor
	 */
	SparseMatrix (const SparseMatrix<Field, Row, Vector> &B); 

	/** Destructor. */
	~SparseMatrix ();

	/** Create a clone of the matrix
	 */
	BlackboxArchetype<Vector> *clone () const
		{ return new SparseMatrix (*this); }

	/** Matrix-vector product
	 * y = A x.
	 * @return reference to output vector y
	 * @param  x input vector
	 */
	Vector &apply (Vector &y, const Vector &x) const;

	/** Transpose matrix-vector product
	 * y = A^T x.
	 * @return reference to output vector y
	 * @param  x input vector
	 */
	Vector &applyTranspose (Vector &y, const Vector &x) const
		{ return applyTransposeSpecialized (y, x, VectorTraits<Row>::VectorCategory ()); }

	/** Retreive row dimensions of Sparsemat matrix.
	 * @return integer number of rows of SparseMatrix0Base matrix.
	 */
	size_t rowdim () const { return _m; }

	/** Retreive column dimensions of Sparsemat matrix.
	 * @return integer number of columns of SparseMatrix0Base matrix.
	 */
	size_t coldim () const { return _n; }

	/** Read the matrix from a stream in the given format
	 * @param is Input stream from which to read the matrix
	 * @param format Format of input matrix
	 * @return Reference to input stream
	 */
	std::istream &read (std::istream &is, Format format = FORMAT_DETECT);

	/** Write the matrix to a stream in the given format
	 * @param os Output stream to which to write the matrix
	 * @param format Format of output
	 * @return Reference to output stream
	 */
	std::ostream &write (std::ostream &os, Format format = FORMAT_PRETTY);

	// JGD 28.08.2002
	/** Access to the base field
	 */
	const Field& field () const { return _F;}


    private:

	const Field                            &_F;      // Field used for all arithmetic
	VectorDomain<Field>                     _VD;     // Vector domain for vector operations
	mutable std::vector<FieldAXPY<Field> >  _faxpy;  // FieldAXPY objects used for applyTranspose

	inline Vector &applyTransposeSpecialized (Vector &y, const Vector &x) const;

	template <class RowTrait>
	inline Vector &applyTransposeSpecialized (Vector &y, const Vector &x,
						  VectorCategories::SparseSequenceVectorTag<RowTrait> tag) const
		{ return applyTransposeSpecialized (y, x); }
	template <class RowTrait>
	inline Vector &applyTransposeSpecialized (Vector &y, const Vector &x,
						  VectorCategories::SparseAssociativeVectorTag<RowTrait> tag) const
		{ return applyTransposeSpecialized (y, x); }
	template <class RowTrait>
	inline Vector &applyTransposeSpecialized (Vector &y, const Vector &x,
						  VectorCategories::SparseParallelVectorTag<RowTrait> tag) const;
};

template <class Field, class _Vector, class _Row, class VectorTrait>
class SparseMatrix<Field, _Vector, _Row, VectorCategories::DenseVectorTag<VectorTrait> >
	: public SparseMatrixBase<typename Field::Element, _Row>, public BlackboxArchetype<_Vector>
{
    public:

	typedef _Vector Vector;
	typedef typename Field::Element Element;
	typedef typename SparseMatrixBase<typename Field::Element, _Row>::Row Row;
	typedef typename SparseMatrixBase<typename Field::Element, _Row>::Format Format;
	typedef typename SparseMatrixBase<typename Field::Element, _Row>::RawIterator RawIterator;
	typedef typename SparseMatrixBase<typename Field::Element, _Row>::RawIndexedIterator RawIndexedIterator;

	SparseMatrix (const Field &F, size_t m, size_t n) 
		: SparseMatrixBase<Element, Row> (m, n), _F (F), _VD (F) {}
	SparseMatrix (const Field &F, VectorStream<Row> &stream)
		: SparseMatrixBase<typename Field::Element, Row> (stream.m (), stream.n ()), _F (F), _VD (F)
	{
		linbox_check (stream.m () > 0);

		typename Rep::iterator i = _A.begin ();

		while (stream) {
			stream.next (*i);
			i++;
		}
	}

	SparseMatrix (const SparseMatrix<Field, Row, Vector> &B)
		: SparseMatrixBase<Element, Row> (B), _F (B._F), _VD (B._F), _faxpy (B._faxpy) {}
	~SparseMatrix () {}
	BlackboxArchetype<Vector> *clone () const
		{ return new SparseMatrix (*this); }
	Vector &apply (Vector &y, const Vector &x) const;
	Vector &applyTranspose (Vector &y, const Vector &x) const
		{ return applyTransposeSpecialized (y, x, VectorTraits<Row>::VectorCategory ()); }

	size_t rowdim () const { return _m; }
	size_t coldim () const { return _n; }

	std::istream &read (std::istream &is, Format format = FORMAT_DETECT)
		{ return SparseMatrixBase<Element, Row>::read (is, _F, format); }
	std::ostream &write (std::ostream &os, Format format = FORMAT_PRETTY)
		{ return SparseMatrixBase<Element, Row>::write (os, _F, format); }

	// JGD 28.08.2002
	/** Access to the base field
	 */
	const Field& field () const { return _F;}


    private:

	const Field &_F;
	VectorDomain<Field> _VD;
	mutable std::vector<FieldAXPY<Field> > _faxpy;

	inline Vector &applyTransposeSpecialized (Vector &y, const Vector &x) const;

	template <class RowTrait>
	inline Vector &applyTransposeSpecialized (Vector &y, const Vector &x,
						  VectorCategories::SparseSequenceVectorTag<RowTrait> tag) const
		{ return applyTransposeSpecialized (y, x); }
	template <class RowTrait>
	inline Vector &applyTransposeSpecialized (Vector &y, const Vector &x,
						  VectorCategories::SparseAssociativeVectorTag<RowTrait> tag) const
		{ return applyTransposeSpecialized (y, x); }
	template <class RowTrait>
	inline Vector &applyTransposeSpecialized (Vector &y, const Vector &x,
						  VectorCategories::SparseParallelVectorTag<RowTrait> tag) const;
};
	  
template <class Field, class _Vector, class _Row, class VectorTrait>
class SparseMatrix<Field, _Vector, _Row, VectorCategories::SparseSequenceVectorTag<VectorTrait> >
	: public SparseMatrixBase<typename Field::Element, _Row>, public BlackboxArchetype<_Vector>
{
    public:

	typedef _Vector Vector;
	typedef typename Field::Element Element;
	typedef typename SparseMatrixBase<typename Field::Element, _Row>::Row Row;
	typedef typename SparseMatrixBase<typename Field::Element, _Row>::Format Format;
	typedef typename SparseMatrixBase<typename Field::Element, _Row>::RawIterator RawIterator;
	typedef typename SparseMatrixBase<typename Field::Element, _Row>::RawIndexedIterator RawIndexedIterator;

	SparseMatrix (const Field &F, size_t m, size_t n) 
		: SparseMatrixBase<Element, Row> (m, n), _F (F), _VD (F) {}
	SparseMatrix (const Field &F, VectorStream<Row> &stream)
		: SparseMatrixBase<typename Field::Element, Row> (stream.m (), stream.n ()), _F (F), _VD (F)
	{
		linbox_check (stream.m () > 0);

		typename Rep::iterator i = _A.begin ();

		while (stream) {
			stream.next (*i);
			i++;
		}
	}

	SparseMatrix (const SparseMatrix<Field, Row, Vector> &B)
		: SparseMatrixBase<Element, Row> (B), _F (B._F), _VD (B._F), _faxpy (B._faxpy) {}
	~SparseMatrix () {}
	BlackboxArchetype<Vector> *clone () const
		{ return new SparseMatrix (*this); }
	Vector &apply (Vector &y, const Vector &x) const;
	Vector &applyTranspose (Vector &y, const Vector &x) const
		{ return applyTransposeSpecialized (y, x, VectorTraits<Row>::VectorCategory ()); }

	size_t rowdim () const { return _m; }
	size_t coldim () const { return _n; }

	std::istream &read (std::istream &is, Format format = FORMAT_DETECT)
		{ return SparseMatrixBase<Element, Row>::read (is, _F, format); }
	std::ostream &write (std::ostream &os, Format format = FORMAT_PRETTY)
		{ return SparseMatrixBase<Element, Row>::write (os, _F, format); }

	// JGD 28.08.2002
	/** Access to the base field
	 */
	const Field& field () const { return _F;}


    private:

	const Field &_F;
	VectorDomain<Field> _VD;
	mutable std::vector<FieldAXPY<Field> > _faxpy;

	inline Vector &applyTransposeSpecialized (Vector &y, const Vector &x) const;

	template <class RowTrait>
	inline Vector &applyTransposeSpecialized (Vector &y, const Vector &x,
						  VectorCategories::SparseSequenceVectorTag<RowTrait> tag) const
		{ return applyTransposeSpecialized (y, x); }
	template <class RowTrait>
	inline Vector &applyTransposeSpecialized (Vector &y, const Vector &x,
						  VectorCategories::SparseAssociativeVectorTag<RowTrait> tag) const
		{ return applyTransposeSpecialized (y, x); }
	template <class RowTrait>
	inline Vector &applyTransposeSpecialized (Vector &y, const Vector &x,
						  VectorCategories::SparseParallelVectorTag<RowTrait> tag) const;
};

template <class Field, class _Vector, class _Row, class VectorTrait>
class SparseMatrix<Field, _Vector, _Row, VectorCategories::SparseAssociativeVectorTag<VectorTrait> >
	: public SparseMatrixBase<typename Field::Element, _Row>, public BlackboxArchetype<_Vector>
{
    public:

	typedef _Vector Vector;
	typedef typename Field::Element Element;
	typedef typename SparseMatrixBase<typename Field::Element, _Row>::Row Row;
	typedef typename SparseMatrixBase<typename Field::Element, _Row>::Format Format;
	typedef typename SparseMatrixBase<typename Field::Element, _Row>::RawIterator RawIterator;
	typedef typename SparseMatrixBase<typename Field::Element, _Row>::RawIndexedIterator RawIndexedIterator;

	SparseMatrix (const Field &F, size_t m, size_t n) 
		: SparseMatrixBase<Element, Row> (m, n), _F (F), _VD (F) {}
	SparseMatrix (const Field &F, VectorStream<Row> &stream)
		: SparseMatrixBase<typename Field::Element, Row> (stream.m (), stream.n ()), _F (F), _VD (F)
	{
		linbox_check (stream.m () > 0);

		typename Rep::iterator i = _A.begin ();

		while (stream) {
			stream.next (*i);
			i++;
		}
	}

	SparseMatrix (const SparseMatrix<Field, Row, Vector> &B)
		: SparseMatrixBase<Element, Row> (B), _F (B._F), _VD (B._F), _faxpy (B._faxpy) {}
	~SparseMatrix () {}
	BlackboxArchetype<Vector> *clone () const
		{ return new SparseMatrix (*this); }
	Vector &apply (Vector &y, const Vector &x) const;
	Vector &applyTranspose (Vector &y, const Vector &x) const
		{ return applyTransposeSpecialized (y, x, VectorTraits<Row>::VectorCategory ()); }

	size_t rowdim () const { return _m; }
	size_t coldim () const { return _n; }

	std::istream &read (std::istream &is, Format format = FORMAT_DETECT)
		{ return SparseMatrixBase<Element, Row>::read (is, _F, format); }
	std::ostream &write (std::ostream &os, Format format = FORMAT_PRETTY)
		{ return SparseMatrixBase<Element, Row>::write (os, _F, format); }

	// JGD 28.08.2002
	/** Access to the base field
	 */
	const Field& field () const { return _F;}


    private:

	const Field &_F;
	VectorDomain<Field> _VD;
	mutable std::vector<FieldAXPY<Field> > _faxpy;

	inline Vector &applyTransposeSpecialized (Vector &y, const Vector &x) const;

	template <class RowTrait>
	inline Vector &applyTransposeSpecialized (Vector &y, const Vector &x,
						  VectorCategories::SparseSequenceVectorTag<RowTrait> tag) const
		{ return applyTransposeSpecialized (y, x); }
	template <class RowTrait>
	inline Vector &applyTransposeSpecialized (Vector &y, const Vector &x,
						  VectorCategories::SparseAssociativeVectorTag<RowTrait> tag) const
		{ return applyTransposeSpecialized (y, x); }
	template <class RowTrait>
	inline Vector &applyTransposeSpecialized (Vector &y, const Vector &x,
						  VectorCategories::SparseParallelVectorTag<RowTrait> tag) const;
};

template <class Field, class _Vector, class _Row, class VectorTrait>
class SparseMatrix<Field, _Vector, _Row, VectorCategories::SparseParallelVectorTag<VectorTrait> >
	: public SparseMatrixBase<typename Field::Element, _Row>, public BlackboxArchetype<_Vector>
{
    public:

	typedef _Vector Vector;
	typedef typename Field::Element Element;
	typedef typename SparseMatrixBase<typename Field::Element, _Row>::Row Row;
	typedef typename SparseMatrixBase<typename Field::Element, _Row>::Format Format;
	typedef typename SparseMatrixBase<typename Field::Element, _Row>::RawIterator RawIterator;
	typedef typename SparseMatrixBase<typename Field::Element, _Row>::RawIndexedIterator RawIndexedIterator;

	SparseMatrix (const Field &F, size_t m, size_t n) 
		: SparseMatrixBase<Element, Row> (m, n), _F (F), _VD (F) {}
	SparseMatrix (const Field &F, VectorStream<Row> &stream)
		: SparseMatrixBase<typename Field::Element, Row> (stream.m (), stream.n ()), _F (F), _VD (F)
	{
		linbox_check (stream.m () > 0);

		typename Rep::iterator i = _A.begin ();

		while (stream) {
			stream.next (*i);
			i++;
		}
	}

	SparseMatrix (const SparseMatrix<Field, Row, Vector> &B)
		: SparseMatrixBase<Element, Row> (B), _F (B._F), _VD (B._F), _faxpy (B._faxpy) {}
	~SparseMatrix () {}
	BlackboxArchetype<Vector> *clone () const
		{ return new SparseMatrix (*this); }
	Vector &apply (Vector &y, const Vector &x) const;
	Vector &applyTranspose (Vector &y, const Vector &x) const
		{ return applyTransposeSpecialized (y, x, VectorTraits<Row>::VectorCategory ()); }

	size_t rowdim () const { return _m; }
	size_t coldim () const { return _n; }

	std::istream &read (std::istream &is, Format format = FORMAT_DETECT)
		{ return SparseMatrixBase<Element, Row>::read (is, _F, format); }
	std::ostream &write (std::ostream &os, Format format = FORMAT_PRETTY)
		{ return SparseMatrixBase<Element, Row>::write (os, _F, format); }

	// JGD 28.08.2002
	/** Access to the base field
	 */
	const Field& field () const { return _F;}


    private:

	const Field &_F;
	VectorDomain<Field> _VD;
	mutable std::vector<FieldAXPY<Field> > _faxpy;

	inline Vector &applyTransposeSpecialized (Vector &y, const Vector &x) const;

	template <class RowTrait>
	inline Vector &applyTransposeSpecialized (Vector &y, const Vector &x,
						  VectorCategories::SparseSequenceVectorTag<RowTrait> tag) const
		{ return applyTransposeSpecialized (y, x); }
	template <class RowTrait>
	inline Vector &applyTransposeSpecialized (Vector &y, const Vector &x,
						  VectorCategories::SparseAssociativeVectorTag<RowTrait> tag) const
		{ return applyTransposeSpecialized (y, x); }
	template <class RowTrait>
	inline Vector &applyTransposeSpecialized (Vector &y, const Vector &x,
						  VectorCategories::SparseParallelVectorTag<RowTrait> tag) const;
};

/** Sparse matrix factory
 * This class inherits \ref{BlackboxFactory} and provides a method for using a
 * \ref{SparseMatrixBase} object with integer or rational data type as input to
 * the high-level integer and rational solutions functions.
 */

template <class Field,
	  class BElement = typename Field::Element,
	  class _Vector  = typename LinBox::Vector<Field>::Dense,
	  class Row      = typename LinBox::Vector<Field>::Sparse,
	  class BRow     = typename LinBox::RawVector<BElement>::Sparse>
class SparseMatrixFactory : public BlackboxFactory<Field, _Vector> 
{
	const SparseMatrixBase<BElement, BRow> &_A;

    public:

	typedef _Vector Vector;

	SparseMatrixFactory (const SparseMatrixBase<BElement, BRow> &A)
		: _A (A) 
	{}

	BlackboxArchetype<Vector> *makeBlackbox (const Field &F);

	// FIXME: This function assumes basically that the matrix is over the integers
	integer &maxNorm (integer &res)
	{
		typename SparseMatrixBase<BElement, BRow>::ConstRawIterator i;

		res = 0L;

		integer tmp;

		for (i = _A.rawBegin (); i != _A.rawEnd (); ++i) {
			tmp = abs (*i);

			if (res < tmp)
				res = tmp;
		}

		return res;
	}

	size_t rowdim ()
		{ return _A.rowdim (); }
	size_t coldim ()
		{ return _A.coldim (); }
};

template <class Field, class _Vector, class _Row, class Trait>
struct MatrixTraits< SparseMatrix<Field, _Vector, _Row, Trait> >
{ 
	typedef SparseMatrix<Field, _Vector, _Row, Trait> MatrixType;
	typedef typename MatrixCategories::RowMatrixTag<MatrixTraits<MatrixType> > MatrixCategory; 
};

} // namespace LinBox

#include "linbox/blackbox/sparse.inl"

#endif // __BLACKBOX_SPARSE_H
