/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/blackbox/submatrix.h
 * Copyright (C) 2001 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * ------------------------------------
 * Modified by Dmitriy Morozov <linbox@foxcub.org>. May 27, 2002.
 *
 * Added parametrization of VectorCategory tags by VectorTraits. See
 * vector-traits.h for more details.
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __SUBMATRIX_H
#define __SUBMATRIX_H

#include "linbox/vector/vector-traits.h"
#include "linbox/util/debug.h"
#include "linbox/util/error.h"

// Namespace in which all LinBox library code resides
namespace LinBox
{

	/** @memo Blackbox for a leading principal minor of existing matrix without copying.
	 * @doc
	 * leading principal minor of an existing matrix in a black box fashion.
	 *
	 * The matrix itself is not stored in memory.  Rather, its apply
	 * methods use a vector of {@link Fields field} elements, which are 
	 * used to "multiply" the matrix to a vector.
	 * 
	 * This class has three template parameters.  The first is the field in 
	 * which the arithmetic is to be done.  The second is the type of 
	 * \Ref{LinBox} vector to which to apply the matrix.  The 
	 * third is chosen be default to be the \Ref{LinBox} vector trait
	 * of the vector.  This class is then specialized for dense and sparse 
	 * vectors.
	 *
	 * @param Field \Ref{LinBox} field
	 * @param Vector \Ref{LinBox} dense or sparse vector of field elements
	 * @param Trait  Marker whether to use dense or sparse LinBox vector 
	 *               implementation.  This is chosen by a default parameter 
	 *               and partial template specialization.  */
	template <class Blackbox, class Trait = typename VectorTraits<typename LinBox::Vector<typename Blackbox::Field>::Dense >::VectorCategory>
	class Submatrix;

	/* Specialization for dense vectors */

	template <class Blackbox, class VectorTrait>
	class Submatrix<Blackbox, VectorCategories::DenseVectorTag<VectorTrait> >
	{
	    public:
		typedef typename Blackbox::Field Field; 

		typedef typename Field::Element Element;
		/** Constructor from field and dense vector of field elements.
		 * @param BB   Black box from which to extract the submatrix
		 * @param row  First row of the submatrix to extract (1.._BB->rowdim ())
		 * @param col  First column of the submatrix to extract (1.._BB->coldim ())
		 * @param rowdim Row dimension
		 * @param coldim Column dimension
		 */
		Submatrix (const Blackbox *BB,
			   size_t          row,
			   size_t          col,
			   size_t          rowdim,
			   size_t          coldim)
			: _BB (BB),
			_row (row), _col (col), _rowdim (rowdim), _coldim (coldim),
			  _z (_BB->coldim ()), _y (_BB->rowdim ())
		{
			linbox_check (row + rowdim <= _BB->rowdim ());
			linbox_check (col + coldim <= _BB->coldim ());

			BB->field().init (_zero, 0);
		}

		/** Destructor
		 */
		virtual ~Submatrix () {}

		/** Application of BlackBox matrix.
		 * y= A*x.
		 * Requires one vector conforming to the \Ref{LinBox}
		 * vector {@link Archetypes archetype}.
		 * Required by abstract base class.
		 * @return reference to vector y containing output.
		 * @param  x constant reference to vector to contain input
		 */
		template<class OutVector, class InVector>
	        OutVector& apply (OutVector &y, const InVector& x) const
	        {
			std::fill (_z.begin (), _z.begin () + _col, _zero);
			std::fill (_z.begin () + _col + _coldim, _z.end (), _zero);

			copy (x.begin (), x.end (), _z.begin () + _col);  // Copying. Yuck.
			_BB->apply (_y, _z);
			copy (_y.begin () + _row, _y.begin () + _row + _rowdim, y.begin ());
			return y;
	        }

		/** Application of BlackBox matrix transpose.
		 * y= transpose(A)*x.
		 * Requires one vector conforming to the \Ref{LinBox}
		 * vector {@link Archetypes archetype}.
		 * Required by abstract base class.
		 * @return reference to vector y containing output.
		 * @param  x constant reference to vector to contain input
		 */
		template<class OutVector, class InVector>
		OutVector& applyTranspose (OutVector &y, const InVector& x) const
		{
			std::fill (_y.begin (), _y.begin () + _row, _zero);
			std::fill (_y.begin () + _row + _rowdim, _y.end (), _zero);

			copy (x.begin (), x.end (), _y.begin () + _row);  // Copying. Yuck.
			_BB->applyTranspose (_z, _y);
			copy (_z.begin () + _col, _z.begin () + _col + _coldim, y.begin ());
			return y;
		}

		/** Retreive _row dimensions of BlackBox matrix.
		 * This may be needed for applying preconditioners.
		 * Required by abstract base class.
		 * @return integer number of _rows of black box matrix.
		 */
		size_t rowdim (void) const
			{ return _rowdim; }
    
		/** Retreive _column dimensions of BlackBox matrix.
		 * Required by abstract base class.
		 * @return integer number of _columns of black box matrix.
		 */
		size_t coldim (void) const
			{ return _coldim; }

		const Field& field() const {return _BB->field();}

	    private:

		const Blackbox *_BB;
		size_t    _row;
		size_t    _col;
		size_t    _rowdim;
		size_t    _coldim;

	        // Temporaries for reducing the amount of memory allocation we do
	        mutable std::vector<Element> _z;
	        mutable std::vector<Element> _y;

		typename Field::Element _zero;

	}; // template <Vector> class Submatrix

} // namespace LinBox

#endif // __SUBMATRIX_H
