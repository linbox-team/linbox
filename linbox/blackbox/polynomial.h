/* linbox/blackbox/polynomial.h
 * Copyright (C) 2005 Cl'ement Pernet
 *
 * Written by Cl'ement Pernet <Clement.Pernet@imag.fr>
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

#ifndef __LINBOX_bb_polynomial_H
#define __LINBOX_bb_polynomial_H

#include "linbox/blackbox/blackbox-interface.h"
#include "linbox/vector/vector-domain.h"
// Namespace in which all LinBox library code resides
namespace LinBox
{
	template <class Blackbox, class Poly>
	class PolynomialBB ;
	template <class Blackbox, class Poly>
	class PolynomialBBOwner ;
}


namespace LinBox
{

	/** \brief represent the matrix P(A) where A is a blackbox and P a polynomial

	  \ingroup blackbox

*/
	template <class Blackbox, class Poly>
	class PolynomialBB : public BlackboxInterface {
	public:

		typedef typename Blackbox::Field Field;
		typedef typename Blackbox::Element Element;
		typedef Poly Polynomial;
		typedef PolynomialBB<Blackbox,Polynomial> Self_t;

		/** Constructor from a black box and a polynomial.
		*/
		PolynomialBB (const Blackbox& A, const Polynomial& P) :
			_A_ptr(&A), _P_ptr(&P), _VD(A.field())
		{}

		PolynomialBB (const Blackbox *A_ptr, const Polynomial * P_ptr)
		: _A_ptr(A_ptr), _P_ptr(P_ptr), _VD(A_ptr->field())
		{
		}

		/** Copy constructor.
		 * Creates new black box objects in dynamic memory.
		 * @param Mat constant reference to compose black box matrix
		 */
		PolynomialBB (const PolynomialBB<Blackbox, Polynomial> &Mat) :
			_A_ptr(Mat._A_ptr), _P_ptr(Mat._P_ptr), _VD(Mat._VD)
		{
		}

		/// Destructor
		~PolynomialBB (void)
		{
		}


		/** Application of BlackBox matrix.
		 * <code>y = P(A)x</code>
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
			Vector2 u (x);
			Vector2 v(field(),u.size());
			_VD.mul( y, x, _P_ptr->operator[](0) );
			for (size_t i=1; i<_P_ptr->size(); ++i){
				_A_ptr->apply( v, u );
				_VD.axpyin( y, _P_ptr->operator[](i), v);
				u=v;
			}
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
			Vector2 u( x );
			Vector2 v(field(),u.size());
			_VD.mul( y, x, _P_ptr->operator[](0));
			for (size_t i=1; i<_P_ptr->size(); ++i){
				_A_ptr->applyTranspose( v, u );
				_VD.axpyin( y, _P_ptr->operator[](i), v);
				u=v;
			}
			return y;
		}


		template<typename _Tp1, class Poly1 = typename Polynomial::template rebind<_Tp1>::other>
		struct rebind {
			typedef PolynomialBBOwner<typename Blackbox::template rebind<_Tp1>::other, Poly1> other;

			void operator() (other & Ap, const Self_t& A, const _Tp1& F) {
				typename Polynomial::template rebind<_Tp1>() (Ap.getDataPolynomial(), *A.getPolynomial(), F);
				typename Blackbox::template rebind<_Tp1>() (Ap.getDataBlackbox(), *A.getBlackbox(),F);

			}
		};



		/** Retreive row dimensions of BlackBox matrix.
		 * This may be needed for applying preconditioners.
		 * Required by abstract base class.
		 * @return integer number of rows of black box matrix.
		 */
		size_t rowdim (void) const
		{
			if (_A_ptr != 0)
				return _A_ptr->rowdim ();
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
				return _A_ptr->coldim ();
			else
				return 0;
		}


		const Polynomial* getPolynomial () const  { return _P_ptr; }
		const Blackbox* getBlackbox () const { return _A_ptr; }
		const Field& field () const {return _A_ptr->field();}
	private:

		// Pointers to A and P
		const Blackbox *_A_ptr;
		const Polynomial *_P_ptr;
		const VectorDomain<Field> _VD;

	};

} // namespace LinBox

// Namespace in which all LinBox library code resides
namespace LinBox
{

	/** \brief represent the matrix P(A) where A is a blackbox and P a polynomial

	  \ingroup blackbox

*/
	template <class Blackbox, class Poly>
	class PolynomialBBOwner : public BlackboxInterface {
	public:

		typedef typename Blackbox::Field Field;
		typedef typename Blackbox::Element Element;
		typedef Poly Polynomial;
		typedef PolynomialBBOwner<Blackbox,Polynomial> Self_t;

		/** Constructor from a black box and a polynomial.
		*/
		PolynomialBBOwner (const Blackbox& A, const Polynomial& P) :
			_A_data(A), _P_data(P), _VD(A.field())
		{}

		PolynomialBBOwner (const Blackbox *A_data, const Polynomial * P_data) :
			_A_data(*A_data), _P_data(*P_data), _VD(A_data->field())
		{
		}

		/** Copy constructor.
		 * Creates new black box objects in dynamic memory.
		 * @param Mat constant reference to compose black box matrix
		 */
		PolynomialBBOwner (const PolynomialBBOwner<Blackbox, Polynomial> &Mat) :
			_A_data(Mat._A_data), _P_data(Mat._P_data), _VD(Mat._VD)
		{
		}

		/// Destructor
		~PolynomialBBOwner (void)
		{
		}


		/** Application of BlackBox matrix.
		 * <code>y = P(A)x</code>
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
			Vector2 u (x);
			Vector2 v(field(),u.size());
			_VD.mul( y, x, _P_data[0] );
			for (size_t i=1; i<_P_data.size(); ++i){
				_A_data.apply( v, u );
				_VD.axpyin( y, _P_data[i], v);
				u=v;
			}
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
			Vector2 u( x );
			Vector2 v(field(), u.size());
			_VD.mul( y, x, _P_data[0]);
			for (size_t i=1; i<_P_data.size(); ++i){
				_A_data.applyTranspose( v, u );
				_VD.axpyin( y, _P_data[i], v);
				u=v;
			}
			return y;
		}


		template<typename _Tp1, class Poly1 = typename Polynomial::template rebind<_Tp1>::other>
		struct rebind {
			typedef PolynomialBBOwner<typename Blackbox::template rebind<_Tp1>::other, Poly1> other;

			void operator() (other & Ap, const Self_t& A, const _Tp1& F) {
				typename Polynomial::template rebind<_Tp1>() (Ap.getDataPolynomial(), A.getDataPolynomial(), F);
				typename Blackbox::template rebind<_Tp1>() (Ap.getDataBlackbox(), A.getDataPolynomial(),F);

			}
		};

		template<typename _BBt, typename _Polt, typename Field>
		PolynomialBBOwner (const PolynomialBB<_BBt, _Polt> &Mat, const Field& F) :
                                _VD(F),
                                _A_data(*(Mat.getBlackbox()), F),
                                _P_data(*(Mat.getPolynomial()), F)
		{
            // Not needed anymore, rebind is done by constructor with other F
// 			typename _BBt::template rebind<Field>()(_A_data, *(Mat.getBlackbox()));
// 			typename _Polt::template rebind<Field>()(_P_data, *(Mat.getPolynomial()));
		}

		template<typename _BBt, typename _Polt, typename Field>
		PolynomialBBOwner (const PolynomialBBOwner<_BBt, _Polt> &Mat, const Field& F) :
			_A_data(Mat.getDataBlackbox(), F),
			_P_data(Mat.getDataPolynomial(), F)
		{
            // Not needed anymore, rebind is done by constructor with other F
// 			typename _BBt::template rebind<Field>()(_A_data, Mat.getDataBlackbox(), F);
// 			typename _Polt::template rebind<Field>()(_P_data, Mat.getDataPolynomial(), F);
		}



		/** Retreive row dimensions of BlackBox matrix.
		 * This may be needed for applying preconditioners.
		 * Required by abstract base class.
		 * @return integer number of rows of black box matrix.
		 */
		size_t rowdim (void) const
		{
			return _A_data.rowdim ();
		}

		/** Retreive column dimensions of BlackBox matrix.
		 * Required by abstract base class.
		 * @return integer number of columns of black box matrix.
		 */
		size_t coldim (void) const
		{
			return _A_data.coldim ();
		}


		const Polynomial& getDataPolynomial () const  { return _P_data; }
		const Blackbox& getDataBlackbox () const { return _A_data; }
		const Field& field () const {return _A_data.field();}
	private:

		const VectorDomain<Field> _VD;
		// Matrix A and polynomial P
		Blackbox _A_data;
		Polynomial _P_data;

	};

} // namespace LinBox

#endif // __LINBOX_bb_polynomial_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
