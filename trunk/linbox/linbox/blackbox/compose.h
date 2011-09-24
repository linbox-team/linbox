/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* linbox/blackbox/compose.h
 * Copyright (C) 1999-2001 William J Turner,
 *               2001 Bradford Hovinen
 *
 * Written by William J Turner <wjturner@math.ncsu.edu>,
 *            Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin Street - Fifth Floor,
 * Boston, MA 02110-1301, USA.
 */

#ifndef __LINBOX_compose_H
#define __LINBOX_compose_H


#include "linbox/util/debug.h"
#include "linbox/linbox-config.h"
#include <linbox/blackbox/blackbox-interface.h>

namespace LinBox
{
	template <class _Blackbox1, class _Blackbox2 = _Blackbox1>
	class Compose;

	template <class _Blackbox1, class _Blackbox2 = _Blackbox1>
	class ComposeOwner;
}


namespace LinBox
{

	/**
	 * Blackbox of a product: \f$C = AB\f$, i.e \f$Cx \gets A(Bx)\f$.

	 * This is a class that multiplies two matrices by implementing an
	 * apply method that calls the apply methods of both of the consituent
	 * matrices, one after the other.
	 *
	 * This class, like the Black Box archetype from which it is derived,
	 * is templatized by the vector type to which the matrix is applied.
	 * Both constituent matrices must also use this same vector type.
	 * For specification of the blackbox members see \ref BlackboxArchetype.
	 *
	 * <b> Template parameter:</b> must meet the \ref Vector requirement.
	 \ingroup blackbox
	 */
	//@{
	/// General case
	template <class _Blackbox1, class _Blackbox2>
	class Compose : public BlackboxInterface {
		typedef Compose<_Blackbox1, _Blackbox2> Self_t;
	public:

		typedef _Blackbox1 Blackbox1;
		typedef _Blackbox2 Blackbox2;

		typedef typename Blackbox2::Field Field;
		typedef typename Field::Element Element;

		/** Constructor of C := A*B from blackbox matrices A and B.
		 * Build the product A*B of any two black box matrices of compatible dimensions.
		 * @pre <code>A.coldim() == B.rowdim()</code>.
		 * @param A blackbox
		 * @param B blackbox
		 */
		Compose (const Blackbox1 &A, const Blackbox2 &B) :
			_A_ptr(&A), _B_ptr(&B)
		{
			// Rich Seagraves - "It seems VectorWrapper somehow
			// became depricated.  Makes the assumption that
			// this vector type supports resize"
			// VectorWrapper::ensureDim (_z, _A_ptr->coldim ());
			_z.resize(_A_ptr->coldim());
		}

		/** Constructor of C := (*A_ptr)*(*B_ptr).
		 * This constructor creates a matrix that is a product of two black box
		 * matrices: A*B from pointers to them.
		 * @param A_ptr blackbox
		 * @param B_ptr blackbox
		 */
		Compose (const Blackbox1 *A_ptr, const Blackbox2 *B_ptr) :
			_A_ptr(A_ptr), _B_ptr(B_ptr)
		{
			linbox_check (A_ptr != (Blackbox1 *) 0);
			linbox_check (B_ptr != (Blackbox2 *) 0);
			linbox_check (A_ptr->coldim () == B_ptr->rowdim ());

			//			VectorWrapper::ensureDim (_z, _A_ptr->coldim ());
			_z.resize(_A_ptr->coldim());
		}

		/** Copy constructor.
		 * Copies the composed matrix (a small handle).  The underlying two matrices
		 * are not copied.
		 * @param[in] M blackbox to copy.
		 */
		Compose (const Compose<Blackbox1, Blackbox2>& Mat) :
			_A_ptr ( Mat._A_ptr), _B_ptr ( Mat._B_ptr)
			//{ VectorWrapper::ensureDim (_z, _A_ptr->coldim ()); }
		{
			_z.resize(_A_ptr->coldim());
		}

		/// Destructor
		~Compose () {}

		/*- Virtual constructor.
		 * Required because constructors cannot be virtual.
		 * Make a copy of the BlackboxArchetype object.
		 * Required by abstract base class.
		 * @return pointer to new blackbox object
		*/
#if 0
		 		BlackboxArchetype<_Vector> *clone () const
		 			{ return new Compose (*this); }
#endif

		/** Matrix * column vector product.
		 * \f$ y \gets (A\cdot B)\cdot x\f$
		 * Applies B, then A.
		 * @return reference to vector y containing output.
		 * @param  x constant reference to vector to contain input
		 * @param[out] y the result.
		 */
		template <class OutVector, class InVector>
		inline OutVector& apply (OutVector& y, const InVector& x) const
		{
			if ((_A_ptr != 0) && (_B_ptr != 0)) {
				_B_ptr->apply (_z, x);
				_A_ptr->apply (y, _z);
			}

			return y;
		}

		/** row vector * matrix product.
		 * \f$ y \gets (A\cdot B)^t  \cdot x\f$.
		 * Applies A^t then B^t.
		 * @return reference to vector y containing output.
		 * @param  x constant reference to vector to contain input
		 * @param[out] y the result.
		 */
		template <class OutVector, class InVector>
		inline OutVector& applyTranspose (OutVector& y, const InVector& x) const
		{
			if ((_A_ptr != 0) && (_B_ptr != 0)) {
				_A_ptr->applyTranspose (_z, x);
				_B_ptr->applyTranspose (y, _z);
			}

			return y;
		}

		template<typename _Tp1, typename _Tp2 = _Tp1>
		struct rebind {
			typedef ComposeOwner<
			typename Blackbox1::template rebind<_Tp1>::other,
				 typename Blackbox2::template rebind<_Tp2>::other
				 > other;

			void operator() (other & Ap, const Self_t& A, const _Tp1& F)
			{
				typename Blackbox1::template rebind<_Tp1> () ( Ap.getLeftData(), *(A.getLeftPtr()), F);
				typename Blackbox2::template rebind<_Tp2> () ( Ap.getRightData(), *(A.getRightPtr()), F);
			}

		};



		/*- Retreive row dimensions of BlackBox matrix.
		 * This may be needed for applying preconditioners.
		 * Required by abstract base class.
		 * @return integer number of rows of black box matrix.
		 */
		/// The number of rows
		size_t rowdim (void) const
		{
			if (_A_ptr != 0)
				return _A_ptr->rowdim ();
			else
				return 0;
		}

		/*- Retreive column dimensions of BlackBox matrix.
		 * Required by abstract base class.
		 * @return integer number of columns of black box matrix.
		 */
		/// The number of columns
		size_t coldim(void) const
		{
			if (_B_ptr != 0)
				return _B_ptr->coldim ();
			else
				return 0;
		}
		/// The field.
		const Field& field() const
		{return _B_ptr->field();}

		/// accessor to the blackboxes
		const Blackbox1* getLeftPtr() const
		{return  _A_ptr;}

		/// accessor to the blackboxes
		const Blackbox2* getRightPtr() const
		{return  _B_ptr;}

	protected:

		// Pointers to A and B matrices
		const Blackbox1 *_A_ptr;
		const Blackbox2 *_B_ptr;

		// local intermediate vector
		mutable std::vector<Element> _z;
	};

	/// specialization for _Blackbox1 = _Blackbox2
	template <class _Blackbox>
	class Compose <_Blackbox, _Blackbox> : public BlackboxInterface {
		typedef Compose<_Blackbox, _Blackbox> Self_t;
	public:
		typedef _Blackbox Blackbox;

		typedef typename _Blackbox::Field Field;
		typedef typename _Blackbox::Element Element;


		Compose (const Blackbox& A, const Blackbox& B) {
			_BlackboxL.push_back(&A);
			_BlackboxL.push_back(&B);

			_zl.resize(1);

			_zl.front().resize (A.coldim());
		}

		Compose (const Blackbox* Ap, const Blackbox* Bp) {
			_BlackboxL.push_back(Ap);
			_BlackboxL.push_back(Bp);

			_zl.resize(1);

			_zl.front().resize (Ap ->coldim());
		}

		/** Constructor of C := A*B from blackbox matrices A and B.
		 * Build the product A*B of any two black box matrices of compatible dimensions.
		 * Requires A.coldim() equals B.rowdim().
		 */
		template<class BPVector>
		Compose (const BPVector& v) :
			_BlackboxL(v.begin(), v.end())
		{

			linbox_check(v.size() > 0);
			_zl.resize(v.size() - 1);
			typename std::vector<const Blackbox*>::iterator b_p;
			typename std::vector<std::vector<Element> >::iterator z_p;
			// it would be good to use just 2 vectors and flip/flop.
			for ( b_p = _BlackboxL.begin(), z_p = _zl.begin();
			      z_p != _zl.end(); ++ b_p, ++ z_p)
				z_p -> resize((*b_p) -> coldim());
		}

		~Compose () {}

		template <class OutVector, class InVector>
		inline OutVector& apply (OutVector& y, const InVector& x) const
		{

			typename std::vector<const Blackbox*>::const_reverse_iterator b_p;
			typename std::vector<std::vector<Element> >::reverse_iterator z_p, pz_p;
			b_p = _BlackboxL.rbegin();
			pz_p = z_p = _zl.rbegin();

			(*b_p) -> apply(*pz_p, x);
			++ b_p;  ++ z_p;

			for (; z_p != _zl.rend(); ++ b_p, ++ z_p, ++ pz_p)
				(*b_p) -> apply (*z_p,*pz_p);

			(*b_p) -> apply(y, *pz_p);

			return y;
		}

		/*! Application of BlackBox matrix transpose.
		 * <code>y= transpose(A*B)*x</code>.
		 * Requires one vector conforming to the \ref LinBox
		 * vector @link Archetypes archetype@endlink.
		 * Required by abstract base class.
		 * @return reference to vector y containing output.
		 * @param  x constant reference to vector to contain input
		 * \param y result
		 */
		template <class OutVector, class InVector>
		inline OutVector& applyTranspose (OutVector& y, const InVector& x) const
		{
			typename std::vector<const Blackbox*>::reverse_iterator b_p;
			typename std::vector<std::vector<Element> >::reverse_iterator z_p, nz_p;

			b_p = _BlackboxL.rbegin();
			z_p = nz_p = _zl.rbegin();

			(*b_p) -> applyTranspose (*z_p, x);

			++ b_p; ++ nz_p;

			for (; nz_p != _zl.rend(); ++ z_p, ++ nz_p, ++ b_p)
				(*b_p) -> applyTranspose (*nz_p, *z_p);

			(*b_p) -> applyTranspose (y, *z_p);

			return y;
		}

		template<typename _Tp1>
		struct rebind
		{
			typedef Compose<typename Blackbox::template rebind<_Tp1>::other, typename Blackbox::template rebind<_Tp1>::other> other;

			void operator() (other *& Ap, const Self_t& A, const _Tp1& F) {
				std::vector<typename other::Blackbox *> newPtrV;
				typename std::vector<typename other::Blackbox *>::iterator np;
				typename std::vector<const Blackbox* >::const_iterator bp;
				for( bp = A._BlackboxL.begin(), np = newPtrV.begin();
				     bp != A._BlackboxL.end(); ++bp, ++np) {
					typename Blackbox::template rebind<_Tp1> () (*np, *(*bp), F);
				}
				Ap = new other(newPtrV);
			}
		};

		/*- Retreive row dimensions of BlackBox matrix.
		 * This may be needed for applying preconditioners.
		 * Required by abstract base class.
		 * @return integer number of rows of black box matrix.
		 */
		size_t rowdim (void) const
		{
			return _BlackboxL.front() -> rowdim();
		}

		/*- Retreive column dimensions of BlackBox matrix.
		 * Required by abstract base class.
		 * @return integer number of columns of black box matrix.
		 */
		size_t coldim(void) const
		{
			return _BlackboxL[_BlackboxL.size() - 1] -> coldim();
		}

		const Field& field() const
		{return _BlackboxL.front() -> field();}

		// accessors to the blackboxes

		const Blackbox* getLeftPtr() const
		{return  _BlackboxL.front();}

		const Blackbox* getRightPtr() const
		{return  _BlackboxL.back();}


	protected:

		// Pointers to A and B matrices
		std::vector<const Blackbox*> _BlackboxL;

		// local intermediate vector
		mutable std::vector<std::vector<Element> > _zl;
	};

	//@}

} // namespace LinBox

// was compose-traits.h (by Zhendong Wan)
#include "linbox/blackbox/dense.h"
#include "linbox/blackbox/blas-blackbox.h" // ???
namespace LinBox
{

	/// used in ..., for example
	template<class IMatrix>
	class ComposeTraits {
	public:
		typedef Compose<IMatrix, IMatrix> value_type;
	};

#if 0
	/// used in smith-binary, for example
	template<class Field>
	class ComposeTraits<Protected::DenseMatrix<Field> > {
	public:

		// define the return value type
		typedef Protected::DenseMatrix<Field> value_type;
	};
#endif
	// template<class _Field>
	// class BlasBlackbox ;

	/// used in smith-binary, for example
	template<class Field>
	class ComposeTraits<  BlasBlackbox<Field> > {
	public:

		// define the return value type
		typedef BlasBlackbox<Field> value_type;
	};

}


namespace LinBox
{

	/**
	 * Blackbox of a product: \f$C = AB\f$, i.e \f$Cx \gets A(Bx)\f$.

	 * This is a class that multiplies two matrices by implementing an
	 * apply method that calls the apply methods of both of the consituent
	 * matrices, one after the other.
	 *
	 * This class, like the Black Box archetype from which it is derived,
	 * is templatized by the vector type to which the matrix is applied.
	 * Both constituent matrices must also use this same vector type.
	 * For specification of the blackbox members see \ref BlackboxArchetype.
	 *
	 * <b> Template parameter:</b> must meet the \ref Vector requirement.
	 \ingroup blackbox
	 */
	//@{
	/// General case
	template <class _Blackbox1, class _Blackbox2>
	class ComposeOwner : public BlackboxInterface {
		typedef ComposeOwner<_Blackbox1, _Blackbox2> Self_t;
	public:

		typedef _Blackbox1 Blackbox1;
		typedef _Blackbox2 Blackbox2;

		typedef typename Blackbox2::Field Field;
		typedef typename Field::Element Element;

		/** Constructor of C := A*B from blackbox matrices A and B.
		 * Build the product A*B of any two black box matrices of compatible dimensions.
		 * Requires A.coldim() equals B.rowdim().
		 */
		ComposeOwner (const Blackbox1 &A, const Blackbox2 &B) :
			_A_data(A), _B_data(B)
		{
			// Rich Seagraves - "It seems VectorWrapper somehow
			// became depricated.  Makes the assumption that
			// this vector type supports resize"
			// VectorWrapper::ensureDim (_z, _A_data.coldim ());
			_z.resize(_A_data.coldim());
		}

		/** Constructor of C := (*A_data)*(*B_data).
		 * This constructor creates a matrix that is a product of two black box
		 * matrices: A*B from pointers to them.
		 */
		ComposeOwner (const Blackbox1 *A_data, const Blackbox2 *B_data) :
			_A_data(*A_data), _B_data(*B_data)
		{
			linbox_check (A_data != (Blackbox1 *) 0);
			linbox_check (B_data != (Blackbox2 *) 0);
			linbox_check (A_data.coldim () == B_data.rowdim ());

			// VectorWrapper::ensureDim (_z, _A_data.coldim ());
			_z.resize(_A_data.coldim());
		}

		/** Copy constructor.
		 * Copies the composed matrix (a small handle).  The underlying two matrices
		 * are not copied.
		 * \param M matrix to be copied.
		 */
		ComposeOwner (const ComposeOwner<Blackbox1, Blackbox2>& Mat) :
			_A_data ( Mat.getLeftData()), _B_data ( Mat.getRightData())
		{
			_z.resize(_A_data.coldim());
		}


		/// Destructor
		~ComposeOwner () {}

#if 0
		/*- Virtual constructor.
		 * Required because constructors cannot be virtual.
		 * Make a copy of the BlackboxArchetype object.
		 * Required by abstract base class.
		 * @return pointer to new blackbox object
		*/
		BlackboxArchetype<_Vector> *clone () const
		{
			return new ComposeOwner (*this);
		}
#endif

		/** Matrix * column vector product.
		 * \f$ y= (A\cdot B) \cdot x.\f$
		 * Applies B, then A.
		 * @return reference to vector y containing output.
		 * @param  x constant reference to vector to contain input
		 * \param y result
		 */
		template <class OutVector, class InVector>
		inline OutVector& apply (OutVector& y, const InVector& x) const
		{
			return _A_data.apply (y, _B_data.apply (_z, x));
		}

		/** row vector * matrix product \f$y= (A \times B)^T \cdot x\f$.
		 * Applies A^t then B^t.
		 * @return reference to vector y containing output.
		 * @param  x constant reference to vector to contain input
		 * @param y
		 */
		template <class OutVector, class InVector>
		inline OutVector& applyTranspose (OutVector& y, const InVector& x) const
		{
			return _B_data.applyTranspose (y, _A_data.applyTranspose (_z, x));
		}

		template<typename _Tp1, typename _Tp2 = _Tp1>
		struct rebind {
			typedef ComposeOwner<
			typename Blackbox1::template rebind<_Tp1>::other,
				 typename Blackbox2::template rebind<_Tp2>::other
				 > other;

			void operator() (other & Ap, const Self_t& A, const _Tp1& F) {
				typename Blackbox1::template rebind<_Tp1> () ( Ap.getLeftData(), A.getLeftData(), F);
				typename Blackbox2::template rebind<_Tp2> () ( Ap.getRightData(), A.getRightData(), F);
			}

		};


		template<typename _BBt1, typename _BBt2, typename Field>
		ComposeOwner (const Compose<_BBt1, _BBt2> &Mat, const Field& F) :
			_A_data(*(Mat.getLeftPtr()), F),
			_B_data(*(Mat.getRightPtr()), F),
			_z(_A_data.coldim())
		{
			typename Compose<_BBt1, _BBt2>::template rebind<Field>()(*this,Mat,F);
		}

		template<typename _BBt1, typename _BBt2, typename Field>
		ComposeOwner (const ComposeOwner<_BBt1, _BBt2> &Mat, const Field& F) :
			_A_data(Mat.getLeftData(), F),
			_B_data(Mat.getRightData(), F) ,
			_z(_A_data.coldim())
		{
			typename ComposeOwner<_BBt1, _BBt2>::template rebind<Field>()(*this,Mat,F);
		}


		/*- Retreive row dimensions of BlackBox matrix.
		 * This may be needed for applying preconditioners.
		 * Required by abstract base class.
		 * @return integer number of rows of black box matrix.
		 */
		/// The number of rows
		size_t rowdim (void) const
		{
			return _A_data.rowdim ();
		}

		/*- Retreive column dimensions of BlackBox matrix.
		 * Required by abstract base class.
		 * @return integer number of columns of black box matrix.
		 */
		/// The number of columns
		size_t coldim(void) const
		{
			return _B_data.coldim ();
		}
		/// The field.
		const Field& field() const
		{return _B_data.field();}


		// accessors to the blackboxes without ownership
		const Blackbox1& getLeftData() const
		{return  _A_data;}
		Blackbox1& getLeftData() {return  _A_data;}

		const Blackbox2& getRightData() const
		{return  _B_data;}
		Blackbox2& getRightData() {return  _B_data;}

	protected:

		// A and B matrices
		Blackbox1 _A_data;
		Blackbox2 _B_data;

		// local intermediate vector
		mutable std::vector<Element> _z;
	};

}


#endif // __LINBOX_compose_H

