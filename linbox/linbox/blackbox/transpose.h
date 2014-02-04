
/* linbox/blackbox/transpose.h
 * Copyright (C) 2001 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
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

#ifndef __LINBOX_transpose_H
#define __LINBOX_transpose_H

#include "linbox/blackbox/blackbox-interface.h"

namespace LinBox
{
	template <class Blackbox>
	class Transpose;


	template <class Blackbox>
	class TransposeOwner;
}




// Namespace in which all LinBox library code resides
namespace LinBox
{

	/** \brief transpose matrix without copying.

	  \ingroup blackbox

	 * @param Vector \ref LinBox dense or sparse vector of field elements
	 * @bug no write here. test-blackbox.h requires it
	 */
	template <class Blackbox>
	class Transpose : public BlackboxInterface {

	public:
		typedef Blackbox Blackbox_t;
		typedef Transpose<Blackbox> Self_t;

		typedef typename Blackbox::Field Field;
		typedef typename Blackbox::Element Element;

		/** Constructor from a black box.
		 * This constructor creates a matrix that is the transpose of a black box
		 * matrix A
		 * @param A pointer to black box matrix.
		 */
		Transpose (Blackbox& A) :
			_A_ptr(&A)
		{}

		Transpose (const Blackbox& A) :
			_A_ptr(&A)
		{}

		Transpose (const Blackbox *A_ptr) :
			_A_ptr(A_ptr)
		{
		}

		/** Copy constructor.
		 * @param Mat constant reference to compose black box matrix
		 */
		Transpose (const Transpose<Blackbox> &Mat) :
			_A_ptr(Mat._A_ptr)
		{
			// create new copies of matrices in dynamic memory
			//linbox_check (M._A_ptr != NULL);
			//_A_ptr = M._A_ptr->clone ();
		}

		/// Destructor
		~Transpose (void)
		{
		}

		template<typename _Tp1>
		struct rebind {
			typedef TransposeOwner<typename Blackbox_t::template rebind<_Tp1>::other> other;
			void operator() (other & Ap, const Self_t& A)
			{
				typename Blackbox_t::template rebind<_Tp1> () ( Ap.getData(), *(A.getPtr()));
			}
		};

		/** Application of BlackBox matrix.
		 * <code>y= (A*B)*x</code>.
		 * Requires one vector conforming to the \ref LinBox
		 * vector @link Archetypes archetype@endlink.
		 * Required by abstract base class.
		 * @return reference to vector y containing output.
		 * @param  x constant reference to vector to contain input
		 * @param y
		 */
		template <class Vector1, class Vector2>
		inline Vector1 &apply (Vector1 &y, const Vector2 &x) const
		{
			if (_A_ptr != 0) _A_ptr->applyTranspose (y, x);
			return y;
		}


		/** Application of BlackBox matrix transpose.
		 * <code>y= transpose(A*B)*x</code>.
		 * Requires one vector conforming to the \ref LinBox
		 * vector @link Archetypes archetype@endlink.
		 * Required by abstract base class.
		 * @return reference to vector y containing output.
		 * @param  x constant reference to vector to contain input
		 * @param y
		 */
		template <class Vector1, class Vector2>
		inline Vector1 &applyTranspose (Vector1 &y, const Vector2 &x) const
		{
			if (_A_ptr != 0) _A_ptr->apply (y, x);
			return y;
		}

		/** Retreive row dimensions of BlackBox matrix.
		 * This may be needed for applying preconditioners.
		 * Required by abstract base class.
		 * @return integer number of rows of black box matrix.
		 */
		size_t rowdim (void) const
		{
			if (_A_ptr != 0)
				return _A_ptr->coldim ();
			else
				return 0;
		}

		/** Retreive column dimensions of BlackBox matrix.
		 * Required by abstract base class.
		 * @return integer number of columns of black box matrix.
		 */
		size_t coldim (void) const
		{
			if (_A_ptr != 0)
				return _A_ptr->rowdim ();
			else
				return 0;
		}


		const Field& field() const {return _A_ptr->field();}

		// accessors to the blackboxes
		const Blackbox* getPtr() const {return  _A_ptr;}

		Element& getEntry(Element& x, size_t i, size_t j) const
		{
			return _A_ptr->getEntry(x, j, i);
		}

		void setEntry(size_t i, size_t j, const Element& x)
		{
			const_cast<Blackbox_t*>(_A_ptr)->setEntry( j, i, x);
		}

		std::ostream &write(std::ostream & os) const {
			return _A_ptr->write(os << "transpose of:");
		}
		std::istream &read(std::istream & is) {
			throw LinBoxError("you don't want to call read here");
			return is;
		}

	protected:

		// Pointers to the matrix
		const Blackbox *_A_ptr;


	}; // template <Vector> class Transpose

} // namespace LinBox


namespace LinBox
{

	/** \brief transpose matrix without copying.

	  \ingroup blackbox

	 * @param Vector \ref LinBox dense or sparse vector of field elements
	 */
	template <class Blackbox>
	class TransposeOwner : public BlackboxInterface {

	public:
		typedef Blackbox Blackbox_t;
		typedef TransposeOwner<Blackbox> Self_t;

		typedef typename Blackbox::Field Field;
		typedef typename Blackbox::Element Element;

		/** Constructor from a black box.
		 * This constructor creates a matrix that is the transpose of a black box
		 * matrix A
		 */
		TransposeOwner (const Blackbox& A) :
			_A_data(A)
		{}

		TransposeOwner (const Blackbox *A_data) :
			_A_data(*A_data)
		{ }

		/** Copy constructor.
		 * @param Mat constant reference to compose black box matrix
		 */
		TransposeOwner (const TransposeOwner<Blackbox> &Mat) :
			_A_data(Mat.getData())
		{
#if 0
			create new copies of matrices in dynamic memory
			linbox_check (M.getData() != NULL);
			_A_data = M.getData().clone ();
#endif
		}

		/// Destructor
		~TransposeOwner (void)
		{
		}

		template<typename _Tp1>
		struct rebind {
			typedef TransposeOwner<typename Blackbox::template rebind<_Tp1>::other> other;
			void operator() (other & Ap, const Self_t& A)
			{
				typename Blackbox_t::template rebind<_Tp1> () ( Ap.getData(), A.getData());
			}
		};

		template<typename _BB, class Field>
		TransposeOwner (const Transpose<_BB>& T, const Field& F) :
			_A_data(*(T.getPtr()), F)
		{
			typename Transpose<_BB>::template rebind<Field>()(*this,T);
		}
		template<typename _BB, class Field>
		TransposeOwner (const TransposeOwner<_BB>& T, const Field& F) :
			_A_data(T.getData(), F)
		{
			typename TransposeOwner<_BB>::template rebind<Field>()(*this,T);
		}


		/** Application of BlackBox matrix.
		 * <code>y= (A*B)*x</code>.
		 * Requires one vector conforming to the \ref LinBox
		 * vector @link Archetypes archetype@endlink.
		 * Required by abstract base class.
		 * @return reference to vector y containing output.
		 * @param  x constant reference to vector to contain input
		 * @param y
		 */
		template <class Vector1, class Vector2>
		inline Vector1 &apply (Vector1 &y, const Vector2 &x) const
		{
			return _A_data.applyTranspose (y, x);
		}


		/** Application of BlackBox matrix transpose.
		 * <code>y= transpose(A*B)*x</code>.
		 * Requires one vector conforming to the \ref LinBox
		 * vector @link Archetypes archetype@endlink.
		 * Required by abstract base class.
		 * @return reference to vector y containing output.
		 * @param  x constant reference to vector to contain input
		 * @param y
		 */
		template <class Vector1, class Vector2>
		inline Vector1 &applyTranspose (Vector1 &y, const Vector2 &x) const
		{
			return _A_data.apply (y, x);
		}

		/** Retreive row dimensions of BlackBox matrix.
		 * This may be needed for applying preconditioners.
		 * Required by abstract base class.
		 * @return integer number of rows of black box matrix.
		 */
		size_t rowdim (void) const
		{
			return _A_data.coldim ();
		}

		/** Retreive column dimensions of BlackBox matrix.
		 * Required by abstract base class.
		 * @return integer number of columns of black box matrix.
		 */
		size_t coldim (void) const
		{
			return _A_data.rowdim ();
		}


		const Field& field() const {return _A_data.field();}

		// accessors to the blackboxes without ownership
		const Blackbox& getData() const {return  _A_data;}
		Blackbox& getData() {return  _A_data;}
	private:
		// Takes ownership of the data
		Blackbox _A_data;


	};

} // namespace LinBox



#endif // __LINBOX_transpose_H

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

