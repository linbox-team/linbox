/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/blackbox/sparse0.h
 * Copyright (C) 1999-2001 William J Turner,
 *               2001-2002 Bradford Hovinen
 *
 * Written by William J Turner <wjturner@math.ncsu.edu>,
 *            Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * ------------------------------------
 * Modified by Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * Refactoring:
 *   - Eliminated SparseMatrix0Aux and moved that functionality into Sparse0
 *   - Made SparseMatrix0Base parameterized only on the element type
 *   - New read/write implementations for SparseMatrix0Base, supporting multiple
 *     formats
 *   - Eliminated Gaussian elimination code
 *   - Added iterators, including ColOfRowsIterator, RawIterator, and
 *     RawIndexIterator
 *   - Eliminated operator []; added getEntry; changed put_value to setEntry
 * ------------------------------------
 * 
 * See COPYING for license information.
 */

#ifndef __SPARSE0_H
#define __SPARSE0_H

#include "linbox/blackbox/archetype.h"
#include "linbox/blackbox/sparse0-base.h"
#include "linbox/field/vector-domain.h"
#include "linbox/vector/vector-traits.h"
#include "linbox/util/vector-factory.h"
#include "linbox/util/field-axpy.h"

// Namespace in which all LinBox library code resides
namespace LinBox
{

/** Sparse matrix
 * This is a generic black box for a sparse matrix. It inherits
 * \ref{SparseMatrix0Base}, which implements all of the underlying
 * accessors and iterators.
 */
template <class Field, class Vector, class Row = std::vector<std::pair<size_t, typename Field::Element> >, class Trait = typename VectorTraits<Vector>::VectorCategory>
class SparseMatrix0 : public SparseMatrix0Base<typename Field::Element, Row>, public BlackboxArchetype<Vector>
{
    public:

	typedef typename Field::Element Element;
	typedef typename SparseMatrix0Base<typename Field::Element, Row>::Format Format;
	typedef typename SparseMatrix0Base<typename Field::Element, Row>::RawIterator RawIterator;
	typedef typename SparseMatrix0Base<typename Field::Element, Row>::RawIndexIterator RawIndexIterator;

	/** Constructor.
	 * Builds a zero m x n matrix
	 * Note: the copy constructor and operator= will work as intended
	 *       because of STL's container design
	 * @param  F  Field over which entries exist
	 * @param  m  Row dimension
	 * @param  n  Column dimension
	 */
	SparseMatrix0 (const Field &F, size_t m, size_t n);

	/** Constructor from a vector factory
	 * @param  F  Field over which entries exist
	 * @param  factory  Factory with which to generate row vectors
	 */
	SparseMatrix0 (const Field &F, VectorFactory<Row> &factory); 

	/** Constructor from a SparseMatrix0Base
	 * This constructor initializes all elements of the matrix from those
	 * present in the SparseMatrix0Base container. It can therefore be used
	 * to reduce an integer matrix modulo a prime before performing
	 * computations over a finite field. It requires that the element type
	 * have an implicit conversion to the type \ref{integer}.
	 * @param  F  Field over which entries exist
	 * @param  B  Container
	 */
	template <class BElement, class BRow>
	SparseMatrix0 (const Field &F, const SparseMatrix0Base<BElement, BRow> &B); 

	/** Copy constructor
	 */
	SparseMatrix0 (const SparseMatrix0<Field, Row, Vector> &B); 

	/** Destructor. */
	~SparseMatrix0 ();

	/** Create a clone of the matrix
	 */
	BlackboxArchetype<Vector> *clone () const
		{ return new SparseMatrix0 (*this); }

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
	Vector &applyTranspose (Vector &y, const Vector &x) const;

	/** Retreive row dimensions of Sparsemat matrix.
	 * @return integer number of rows of SparseMatrix00Base matrix.
	 */
	size_t rowdim (void) const { return _m; }

	/** Retreive column dimensions of Sparsemat matrix.
	 * @return integer number of columns of SparseMatrix00Base matrix.
	 */
	size_t coldim (void) const { return _n; }

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
	std::ostream &write (std::ostream &os, Format format = FORMAT_GUILLAUME);

    private:

	const Field                            &_F;      // Field used for all arithmetic
	VectorDomain<Field>                     _VD;     // Vector domain for vector operations
	mutable std::vector<FieldAXPY<Field> >  _faxpy;  // FieldAXPY objects used for applyTranspose
};

template <class Field, class Row, class Vector, class VectorTrait>
class SparseMatrix0<Field, Vector, Row, VectorCategories::DenseVectorTag<VectorTrait> >
	: public SparseMatrix0Base<typename Field::Element, Row>, public BlackboxArchetype<Vector>
{
    public:

	typedef typename Field::Element Element;
	typedef typename SparseMatrix0Base<typename Field::Element, Row>::Format Format;
	typedef typename SparseMatrix0Base<typename Field::Element, Row>::RawIterator RawIterator;
	typedef typename SparseMatrix0Base<typename Field::Element, Row>::RawIndexIterator RawIndexIterator;

	SparseMatrix0 (const Field &F, size_t m, size_t n) 
		: SparseMatrix0Base<Element, Row> (m, n), _F (F), _VD (F) {}
	SparseMatrix0 (const Field &F, VectorFactory<Row> &factory)
		: SparseMatrix0Base<typename Field::Element, Row> (factory.m (), factory.n ()), _F (F), _VD (F)
	{
		linbox_check (factory.m () > 0);

		typename Rep::iterator i = _A.begin ();

		while (factory) {
			factory.next (*i);
			i++;
		}
	}

	template <class BElement, class BRow>
	SparseMatrix0 (const Field &F, SparseMatrix0Base<BElement, BRow> &B)
		: SparseMatrix0Base<typename Field::Element, Row> (B.rowdim (), B.coldim ()), _F (F), _VD (F)
	{
		typename SparseMatrix0Base<BElement, BRow>::RawIterator i;
		typename SparseMatrix0Base<BElement, BRow>::RawIndexIterator j;

		for (i = B.rawBegin (), j = B.indexBegin (); i != B.rawEnd (); i++, j++)
			_F.init (getEntry (j->first, j->second), *i);
	}

	SparseMatrix0 (const SparseMatrix0<Field, Row, Vector> &B)
		: SparseMatrix0Base<Element, Row> (B), _F (B._F), _VD (B._F), _faxpy (B._faxpy) {}
	~SparseMatrix0 () {}
	BlackboxArchetype<Vector> *clone () const
		{ return new SparseMatrix0 (*this); }
	Vector &apply (Vector &y, const Vector &x) const;
	Vector &applyTranspose (Vector &y, const Vector &x) const;

	size_t rowdim (void) const { return _m; }
	size_t coldim (void) const { return _n; }

	std::istream &read (std::istream &is, Format format = FORMAT_DETECT)
		{ return SparseMatrix0Base<Element, Row>::read (is, _F, format); }
	std::ostream &write (std::ostream &os, Format format = FORMAT_GUILLAUME)
		{ return SparseMatrix0Base<Element, Row>::write (os, _F, format); }

    private:

	const Field &_F;
	VectorDomain<Field> _VD;
	mutable std::vector<FieldAXPY<Field> > _faxpy;
};
	  
template <class Field, class Row, class Vector, class VectorTrait>
class SparseMatrix0<Field, Vector, Row, VectorCategories::SparseSequenceVectorTag<VectorTrait> >
	: public SparseMatrix0Base<typename Field::Element, Row>, public BlackboxArchetype<Vector>
{
    public:

	typedef typename Field::Element Element;
	typedef typename SparseMatrix0Base<typename Field::Element, Row>::Format Format;
	typedef typename SparseMatrix0Base<typename Field::Element, Row>::RawIterator RawIterator;
	typedef typename SparseMatrix0Base<typename Field::Element, Row>::RawIndexIterator RawIndexIterator;

	SparseMatrix0 (const Field &F, size_t m, size_t n) 
		: SparseMatrix0Base<Element, Row> (m, n), _F (F), _VD (F) {}
	SparseMatrix0 (const Field &F, VectorFactory<Row> &factory)
		: SparseMatrix0Base<typename Field::Element, Row> (factory.m (), factory.n ()), _F (F), _VD (F)
	{
		linbox_check (factory.m () > 0);

		typename Rep::iterator i = _A.begin ();

		while (factory) {
			factory.next (*i);
			i++;
		}
	}

	template <class BElement, class BRow>
	SparseMatrix0 (const Field &F, SparseMatrix0Base<BElement, BRow> &B)
		: SparseMatrix0Base<typename Field::Element, Row> (B.rowdim (), B.coldim ()), _F (F), _VD (F)
	{
		typename SparseMatrix0Base<BElement, BRow>::RawIterator i;
		typename SparseMatrix0Base<BElement, BRow>::RawIndexIterator j;

		for (i = B.rawBegin (), j = B.indexBegin (); i != B.rawEnd (); i++, j++)
			_F.init (getEntry (j->first, j->second), *i);
	}

	SparseMatrix0 (const SparseMatrix0<Field, Row, Vector> &B)
		: SparseMatrix0Base<Element, Row> (B), _F (B._F), _VD (B._F), _faxpy (B._faxpy) {}
	~SparseMatrix0 () {}
	BlackboxArchetype<Vector> *clone () const
		{ return new SparseMatrix0 (*this); }
	Vector &apply (Vector &y, const Vector &x) const;
	Vector &applyTranspose (Vector &y, const Vector &x) const;

	size_t rowdim (void) const { return _m; }
	size_t coldim (void) const { return _n; }

	std::istream &read (std::istream &is, Format format = FORMAT_DETECT)
		{ return SparseMatrix0Base<Element, Row>::read (is, _F, format); }
	std::ostream &write (std::ostream &os, Format format = FORMAT_GUILLAUME)
		{ return SparseMatrix0Base<Element, Row>::write (os, _F, format); }

    private:

	const Field &_F;
	VectorDomain<Field> _VD;
	mutable std::vector<FieldAXPY<Field> > _faxpy;
};

template <class Field, class Row, class Vector, class VectorTrait>
class SparseMatrix0<Field, Vector, Row, VectorCategories::SparseAssociativeVectorTag<VectorTrait> >
	: public SparseMatrix0Base<typename Field::Element, Row>, public BlackboxArchetype<Vector>
{
    public:

	typedef typename Field::Element Element;
	typedef typename SparseMatrix0Base<typename Field::Element, Row>::Format Format;
	typedef typename SparseMatrix0Base<typename Field::Element, Row>::RawIterator RawIterator;
	typedef typename SparseMatrix0Base<typename Field::Element, Row>::RawIndexIterator RawIndexIterator;

	SparseMatrix0 (const Field &F, size_t m, size_t n) 
		: SparseMatrix0Base<Element, Row> (m, n), _F (F), _VD (F) {}
	SparseMatrix0 (const Field &F, VectorFactory<Row> &factory)
		: SparseMatrix0Base<typename Field::Element, Row> (factory.m (), factory.n ()), _F (F), _VD (F)
	{
		linbox_check (factory.m () > 0);

		typename Rep::iterator i = _A.begin ();

		while (factory) {
			factory.next (*i);
			i++;
		}
	}

	template <class BElement, class BRow>
	SparseMatrix0 (const Field &F, SparseMatrix0Base<BElement, BRow> &B)
		: SparseMatrix0Base<typename Field::Element, Row> (B.rowdim (), B.coldim ()), _F (F), _VD (F)
	{
		typename SparseMatrix0Base<BElement, BRow>::RawIterator i;
		typename SparseMatrix0Base<BElement, BRow>::RawIndexIterator j;

		for (i = B.rawBegin (), j = B.indexBegin (); i != B.rawEnd (); i++, j++)
			_F.init (getEntry (j->first, j->second), *i);
	}

	SparseMatrix0 (const SparseMatrix0<Field, Row, Vector> &B)
		: SparseMatrix0Base<Element, Row> (B), _F (B._F), _VD (B._F), _faxpy (B._faxpy) {}
	~SparseMatrix0 () {}
	BlackboxArchetype<Vector> *clone () const
		{ return new SparseMatrix0 (*this); }
	Vector &apply (Vector &y, const Vector &x) const;
	Vector &applyTranspose (Vector &y, const Vector &x) const;

	size_t rowdim (void) const { return _m; }
	size_t coldim (void) const { return _n; }

	std::istream &read (std::istream &is, Format format = FORMAT_DETECT)
		{ return SparseMatrix0Base<Element, Row>::read (is, _F, format); }
	std::ostream &write (std::ostream &os, Format format = FORMAT_GUILLAUME)
		{ return SparseMatrix0Base<Element, Row>::write (os, _F, format); }

    private:

	const Field &_F;
	VectorDomain<Field> _VD;
	mutable std::vector<FieldAXPY<Field> > _faxpy;
};

} // namespace LinBox

#include "linbox/blackbox/sparse0.inl"

#endif // __SPARSE0_H
