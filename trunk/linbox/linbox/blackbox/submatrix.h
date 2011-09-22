/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
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
 * Modified by Zhendong Wan
 *
 * Added specialization for Protected::DenseMatrix.
 *
 * -------
 *
 * See COPYING for license information.
 */

#ifndef __LINBOX_bb_submatrix_H
#define __LINBOX_bb_submatrix_H

#include "linbox/vector/vector-traits.h"
#include "linbox/util/debug.h"
#include "linbox/util/error.h"
#include <linbox/matrix/dense-submatrix.h>
#include <linbox/vector/vector-domain.h>
#include <linbox/blackbox/blackbox-interface.h>


// Namespace in which all LinBox library code resides
namespace LinBox
{

	/** \brief leading principal minor of existing matrix without copying.
	 * \ingroup blackbox
	 * leading principal minor of an existing matrix in a black box fashion.
	 *
	 * The matrix itself is not stored in memory.  Rather, its apply
	 * methods use a vector of @link Fields field@endlink elements, which are
	 * used to "multiply" the matrix to a vector.
	 *
	 * This class has three template parameters.  The first is the field in
	 * which the arithmetic is to be done.  The second is the type of
	 * \ref LinBox vector to which to apply the matrix.  The
	 * third is chosen be default to be the \ref LinBox vector trait
	 * of the vector.  This class is then specialized for dense and sparse
	 * vectors.
	 *
	 * @param Field \ref LinBox field
	 * @param Vector \ref LinBox dense or sparse vector of field elements
	 * @param Trait  Marker whether to use dense or sparse LinBox vector
	 *               implementation.  This is chosen by a default parameter
	 *               and partial template specialization.
	 */
	//@{
	// Basic declaration.
	template <class Blackbox, class Trait = typename VectorTraits<typename LinBox::Vector<typename Blackbox::Field>::Dense >::VectorCategory>
	class Submatrix : public BlackboxInterface {

	private:

		Submatrix() {}

	};

	template <class Blackbox, class Trait = typename VectorTraits<typename LinBox::Vector<typename Blackbox::Field>::Dense >::VectorCategory>
	class SubmatrixOwner;

}


// Namespace in which all LinBox library code resides
namespace LinBox
{

	/** Specialization for dense vectors */
	template <class Blackbox>
	class Submatrix<Blackbox, VectorCategories::DenseVectorTag > : public BlackboxInterface {
	public:
		typedef typename Blackbox::Field Field;
		typedef typename Field::Element Element;
		typedef Blackbox Blackbox_t;
		typedef Submatrix<Blackbox, VectorCategories::DenseVectorTag > Self_t;

		/** Constructor from field and dense vector of field elements.
		 * @param BB   Black box from which to extract the submatrix
		 * @param row  First row of the submatrix to extract (1.._BB->rowdim ())
		 * @param col  First column of the submatrix to extract (1.._BB->coldim ())
		 * @param Rowdim Row dimension
		 * @param Coldim Column dimension
		 */
		Submatrix (const Blackbox *BB,
			   size_t          row,
			   size_t          col,
			   size_t          Rowdim,
			   size_t          Coldim) :
			_BB (BB),
			_row (row), _col (col), _rowdim (Rowdim), _coldim (Coldim),
			_z (_BB->coldim ()), _y (_BB->rowdim ())
		{
			linbox_check (row + Rowdim <= _BB->rowdim ());
			linbox_check (col + Coldim <= _BB->coldim ());

		}

		/** Destructor
		*/
		virtual ~Submatrix () {}

		/** Application of BlackBox matrix.
		 * <code>y= A*x</code>.
		 * Requires one vector conforming to the \ref LinBox
		 * vector @link Archetypes archetype@endlink.
		 * Required by abstract base class.
		 * @return reference to vector y containing output.
		 * @param  x constant reference to vector to contain input
		 * @param y
		 */
		template<class OutVector, class InVector>
		OutVector& apply (OutVector &y, const InVector& x) const
		{
			std::fill (_z.begin (), _z.begin () + _col, _BB->field().zero);
			std::fill (_z.begin () + _col + _coldim, _z.end (), _BB->field().zero);

			copy (x.begin (), x.end (), _z.begin () + _col);  // Copying. Yuck.
			_BB->apply (_y, _z);
			copy (_y.begin () + _row, _y.begin () + _row + _rowdim, y.begin ());
			return y;
		}

		/** Application of BlackBox matrix transpose.
		 * <code>y= transpose(A)*x</code>.
		 * Requires one vector conforming to the \ref LinBox
		 * vector @link Archetypes archetype@endlink.
		 * Required by abstract base class.
		 * @return reference to vector y containing output.
		 * @param  x constant reference to vector to contain input
		 * @param y
		 */
		template<class OutVector, class InVector>
		OutVector& applyTranspose (OutVector &y, const InVector& x) const
		{
			std::fill (_y.begin (), _y.begin () + _row, _BB->field().zero);
			std::fill (_y.begin () + _row + _rowdim, _y.end (), _BB->field().zero);

			copy (x.begin (), x.end (), _y.begin () + _row);  // Copying. Yuck.
			_BB->applyTranspose (_z, _y);
			copy (_z.begin () + _col, _z.begin () + _col + _coldim, y.begin ());
			return y;
		}


		template<typename _Tp1>
		struct rebind
		{
			typedef SubmatrixOwner<typename Blackbox_t::template rebind<_Tp1>::other, VectorCategories::DenseVectorTag> other;
			void operator() (other & Ap, const Self_t& A, const _Tp1& F) {
				typename Blackbox_t::template rebind<_Tp1> Rebinder;
				Rebinder( Ap.getData(), *(A.getPtr()), F);
			}
		};


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

		size_t rowfirst() const { return _row; }
		size_t colfirst() const { return _col; }



		const Field& field() const {return _BB->field();}

		const Blackbox * getPtr() const { return _BB; }

	private:

		const Blackbox *_BB;
		size_t    _row;
		size_t    _col;
		size_t    _rowdim;
		size_t    _coldim;

		// Temporaries for reducing the amount of memory allocation we do
		mutable std::vector<Element> _z;
		mutable std::vector<Element> _y;

	}; // template <Vector> class Submatrix


	/** Specialization for dense ZeroOne vectors */
	template <class Blackbox>
	class Submatrix<Blackbox, VectorCategories::DenseZeroOneVectorTag > : public Submatrix<Blackbox, VectorCategories::DenseVectorTag> {
	public:
		typedef typename Blackbox::Field Field;
		typedef typename Field::Element Element;
		typedef Blackbox Blackbox_t;
		typedef Submatrix<Blackbox, VectorCategories::DenseZeroOneVectorTag > Self_t;
		typedef Submatrix<Blackbox, VectorCategories::DenseVectorTag > Father_t;

		/** Constructor from field and dense vector of field elements.
		 * @param BB   Black box from which to extract the submatrix
		 * @param row  First row of the submatrix to extract (1.._BB->rowdim ())
		 * @param col  First column of the submatrix to extract (1.._BB->coldim ())
		 * @param Rowdim Row dimension
		 * @param Coldim Column dimension
		 */
		Submatrix (const Blackbox *BB,
			   size_t          row,
			   size_t          col,
			   size_t          Rowdim,
			   size_t          Coldim) :
			Father_t(BB,row,col,Rowdim,Coldim)
		{}
	};

	namespace Protected {
		template <class Field>
		class DenseMatrix;
	}

	/** special case for the submatrix of a dense matrix
	*/
	template<class _Field>
	class Submatrix<Protected::DenseMatrix<_Field>, VectorCategories::DenseVectorTag> : public Protected::DenseSubmatrix<typename _Field::Element> {
	public:

		typedef _Field Field;
		typedef Submatrix<Protected::DenseMatrix<_Field>, VectorCategories::DenseVectorTag> Self_t;
		typedef Protected::DenseSubmatrix<typename _Field::Element> Father_t;

	private:

		Field f;

		VectorDomain<Field> vd;

	public:

		typedef typename Field::Element Element;

		/** Constructor from an existing \ref Protected::DenseMatrix  and dimensions
		 * @param M Pointer to \ref Protected::DenseMatrix  of which to construct submatrix
		 * @param row Starting row
		 * @param col Starting column
		 * @param Rowdim Row dimension
		 * @param Coldim Column dimension
		 */

		Submatrix (const Protected::DenseMatrix<Field> *Mat,
			   size_t row,
			   size_t col,
			   size_t Rowdim,
			   size_t Coldim) :
			Protected::DenseSubmatrix<Element>(const_cast<Protected::DenseMatrix<Field>& >(*Mat), row, col, Rowdim, Coldim),
			f(Mat -> field()), vd(Mat -> field())
		{
		}

		/** Constructor from an existing \ref Protected::DenseMatrix  and dimensions
		 * @param M reference to \ref Protected::DenseMatrix  of which to construct submatrix
		 * @param row Starting row
		 * @param col Starting column
		 * @param Rowdim Row dimension
		 * @param Coldim Column dimension
		 */
		Submatrix (const Protected::DenseMatrix<Field> &Mat,
			   size_t row,
			   size_t col,
			   size_t Rowdim,
			   size_t Coldim) :
			Protected::DenseSubmatrix<Element>(const_cast<Protected::DenseMatrix<Field>& >(Mat), row, col, Rowdim, Coldim),
			f(Mat.field()), vd(Mat.field()) {
			}

		/** Constructor from an existing submatrix and dimensions
		 * @param SM pointer to Submatrix from which to
		 *           construct submatrix
		 * @param row Starting row
		 * @param col Starting column
		 * @param Rowdim Row dimension
		 * @param Coldim Column dimension
		 */
		Submatrix (const Submatrix<Protected::DenseMatrix<Field> > *SM,
			   size_t row,
			   size_t col,
			   size_t Rowdim,
			   size_t Coldim ) :
			Protected::DenseSubmatrix<Element> (const_cast<Submatrix<Protected::DenseMatrix<Field> >&>(*SM), row, col, Rowdim, Coldim),
			f (SM ->  field()), vd(SM -> field()){
			}

		/** Constructor from an existing submatrix and dimensions
		 * @param SM reference to Submatrix from which to
		 *           construct submatrix
		 * @param row Starting row
		 * @param col Starting column
		 * @param Rowdim Row dimension
		 * @param Coldim Column dimension
		 */
		Submatrix (const Submatrix<Protected::DenseMatrix<Field> >& SM,
			   size_t row,
			   size_t col,
			   size_t Rowdim,
			   size_t Coldim ) :
			Protected::DenseSubmatrix<Element> (const_cast<Submatrix<Protected::DenseMatrix<Field> >&>(SM), row, col, Rowdim, Coldim),
			f (SM. field()), vd(SM. field()){
			}

		const Field& field() const {

			return f;
		}

		std::istream& read (std::istream& is) {

			Protected::DenseSubmatrix<Element>::read (is, f);

			return is;
		}

		std::ostream& write (std::ostream& os) const {

			Protected::DenseSubmatrix<Element>::write (os, f);

			return os;
		}

		/** Generic matrix-vector apply
		 * <code>y = A * x</code>.
		 * This version of apply allows use of arbitrary input and output vector         * types.
		 * @param y Output vector
		 * @param x Input vector
		 * @return Reference to output vector
		 */
		template<class Vect1, class Vect2>
		Vect1 &apply (Vect1 &y, const Vect2 &x) const {

			typename Protected::DenseSubmatrix<Element>::ConstRowIterator p;

			typename Vect1::iterator p_y = y.begin ();

			for (p = this->rowBegin (); p != this->rowEnd (); ++p, ++p_y)
				vd.dot (*p_y, *p, x);

			return y;
		}

		/** Generic matrix-vector transpose apply
		 * <code>y = A^T * x</code>
		 * This version of applyTranspose allows use of arbitrary input and
		 * output vector types
		 * @param y Output vector
		 * @param x Input vector
		 * @return Reference to output vector
		 */
		template<class Vect1, class Vect2>
		Vect1 &applyTranspose (Vect1 &y, const Vect2 &x) const {

			typename Protected::DenseSubmatrix<Element>::ConstColIterator colp;

			typename Vect1::iterator p_y = y.begin ();

			for (colp = this->colBegin (); colp != this->colEnd (); ++colp, ++p_y)
				vd. dot (*p_y, *colp, x);

			return y;
		}

		template<typename _Tp1>
		struct rebind
		{
			typedef SubmatrixOwner<Protected::DenseMatrix<_Tp1>, VectorCategories::DenseVectorTag> other;

			void operator() (other & Ap, const Self_t& A, const _Tp1& F) {

				typename other::Father_t A1;
				typename Father_t::template rebind<_Tp1> () ( A1, static_cast<Father_t>(A), F);
				Ap = other(A1, A._row, A._col, A._rowdim, A._coldim);
			}

		};
	};

	//@}
} // namespace LinBox




// Namespace in which all LinBox library code resides
namespace LinBox
{

	/** Specialization for dense vectors */
	template <class Blackbox>
	class SubmatrixOwner<Blackbox, VectorCategories::DenseVectorTag > : public BlackboxInterface {
	public:
		typedef typename Blackbox::Field Field;
		typedef typename Field::Element Element;
		typedef Blackbox Blackbox_t;
		typedef SubmatrixOwner<Blackbox_t, VectorCategories::DenseVectorTag > Self_t;

		/** Constructor from field and dense vector of field elements.
		 * @param BB   Black box from which to extract the submatrix
		 * @param row  First row of the submatrix to extract (1.._BB_data.rowdim ())
		 * @param col  First column of the submatrix to extract (1.._BB_data.coldim ())
		 * @param Rowdim Row dimension
		 * @param Coldim Column dimension
		 */
		SubmatrixOwner (const Blackbox *BB,
				size_t          row,
				size_t          col,
				size_t          Rowdim,
				size_t          Coldim) :
			_BB_data (*BB),
			_row (row), _col (col), _rowdim (Rowdim), _coldim (Coldim),
			_z (_BB_data.coldim ()), _y (_BB_data.rowdim ())
		{
			linbox_check (row + Rowdim <= _BB_data.rowdim ());
			linbox_check (col + Coldim <= _BB_data.coldim ());

		}

		/** Destructor
		*/
		virtual ~SubmatrixOwner () {}

		/** Application of BlackBox matrix.
		 * <code>y= A*x</code>.
		 * Requires one vector conforming to the \ref LinBox
		 * vector @link Archetypes archetype@endlink.
		 * Required by abstract base class.
		 * @return reference to vector y containing output.
		 * @param  x constant reference to vector to contain input
		 * @param y
		 */
		template<class OutVector, class InVector>
		OutVector& apply (OutVector &y, const InVector& x) const
		{
			std::fill (_z.begin (), _z.begin () + _col, _BB_data.field().zero);
			std::fill (_z.begin () + _col + _coldim, _z.end (), _BB_data.field().zero);

			copy (x.begin (), x.end (), _z.begin () + _col);  // Copying. Yuck.
			_BB_data.apply (_y, _z);
			copy (_y.begin () + _row, _y.begin () + _row + _rowdim, y.begin ());
			return y;
		}

		/** Application of BlackBox matrix transpose.
		 * <code>y= transpose(A)*x</code>.
		 * Requires one vector conforming to the \ref LinBox
		 * vector @link Archetypes archetype@endlink.
		 * Required by abstract base class.
		 * @return reference to vector y containing output.
		 * @param  x constant reference to vector to contain input
		 * @param y
		 */
		template<class OutVector, class InVector>
		OutVector& applyTranspose (OutVector &y, const InVector& x) const
		{
			std::fill (_y.begin (), _y.begin () + _row, _BB_data.field().zero);
			std::fill (_y.begin () + _row + _rowdim, _y.end (), _BB_data.field().zero);

			copy (x.begin (), x.end (), _y.begin () + _row);  // Copying. Yuck.
			_BB_data.applyTranspose (_z, _y);
			copy (_z.begin () + _col, _z.begin () + _col + _coldim, y.begin ());
			return y;
		}


		template<typename _Tp1>
		struct rebind {
			typedef SubmatrixOwner< typename Blackbox::template rebind<_Tp1>::other, VectorCategories::DenseVectorTag> other;

			void operator() (other & Ap, const Self_t& A, const _Tp1& F) {
				typename Blackbox_t::template rebind<_Tp1> () ( Ap.getData(), A.getData(), F);

			}

		};

		template<typename _BB, typename _Vc, class Field>
		SubmatrixOwner (const Submatrix<_BB, _Vc>& T, const Field& F) :
			_BB_data(*(T.getPtr()), F),
			_row(T.rowfirst()), _col(T.colfirst()),
			_rowdim(T.rowdim()), _coldim(T.coldim()),
			_z (_BB_data.coldim ()), _y (_BB_data.rowdim ())
		{
			typename Submatrix<_BB,_Vc>::template rebind<Field>()(*this,T, F);
		}
		template<typename _BB, typename _Vc, class Field>
		SubmatrixOwner (const SubmatrixOwner<_BB,_Vc>& T, const Field& F) :
			_BB_data(T.getData(), F),
			_row(T.rowfirst()), _col(T.colfirst()),
			_rowdim(T.rowdim()), _coldim(T.coldim()),
			_z (_BB_data.coldim ()), _y (_BB_data.rowdim ())
		{
			typename SubmatrixOwner<_BB,_Vc>::template rebind<Field>()(*this,T, F);
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

		size_t rowfirst() const { return _row; }
		size_t colfirst() const { return _col; }

		const Field& field() const {return _BB_data.field();}

		const Blackbox& getData() const {return  _BB_data;}
		Blackbox& getData() {return  _BB_data;}

	private:

		Blackbox _BB_data;
		size_t    _row;
		size_t    _col;
		size_t    _rowdim;
		size_t    _coldim;

		// Temporaries for reducing the amount of memory allocation we do
		mutable std::vector<Element> _z;
		mutable std::vector<Element> _y;

	}; // template <Vector> class SubmatrixOwner


	//@}
} // namespace LinBox



#endif // __LINBOX_bb_submatrix_H

