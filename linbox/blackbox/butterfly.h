/* linbox/blackbox/butterfly.h
 * Copyright (C) 1999-2001 William J Turner,
 *               2001 Bradford Hovinen
 *
 * Written by William J Turner <wjturner@math.ncsu.edu>,
 *            Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * -----------------------------------------------------------
 * 2002-09-26  Bradford Hovinen  <bghovinen@math.uwaterloo.ca>
 *
 * Refactoring: The switch object now only contains the information necessary
 * for a single 2x2 block. The butterfly black box maintains a vector of switch
 * objects that it keeps in parallel with its vector of indices. There is a new
 * lightweight class, called a SwitchFactory, that constructs switches on the
 * fly. It is defined individually for each switch type, and a instance thereof
 * is passed to the butterfly, which then uses it to construct its vector.
 *
 * This eliminates two problems: first, because switch objects are constructed
 * by the butterfly itself, there is no need to know a priori the length of the
 * vector of indices. Second, the switch object itself becomes simpler, as it
 * need only be responsible for a single 2x2 block.
 *
 * -----------------------------------------------------------
 *
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 *
 */

#ifndef __LINBOX_butterfly_H
#define __LINBOX_butterfly_H

#include "linbox/blackbox/blackbox-interface.h"
#include "linbox/vector/vector-domain.h"

/*! @file blackbox/butterfly.h
*/

namespace LinBox
{


	/** The default butterfly switch object.
	 *
	 * This is a predicate object that is applied
	 * to two elements to switch them as needed
	 * by the \ref Butterfly\ Switching\ Network\ BlackBox\ Matrix\ Object
	 * following the exchange matrix introduced in "Efficient Matrix
	 * Preconditioners for Black Box Linear Algebra" by Chen, Eberly,
	 * Kaltofen, Saunders, Turner, and Villard.
	 */
	template <class Field> class CekstvSwitch ;

	/// Alternate butterfly switch object for testing.
	class BooleanSwitch;

	/** @name Butterfly
	 * @brief Butterfly preconditioner and supporting function
	 */
	//@{
	//
	/** \brief Switching Network based BlackBox Matrix.  A good preconditioner.

	 * Implements butterfly switching network on a LinBox vector
	 * as a black box matrix through the use of a switch object.
	 *
	 * This is a blackbox matrix object, and it implements all
	 * purely virtual methods of the abstract base class.
	 * See \ref BlackboxArchetype for the specification of these methods.
	 *
	 * This matrix requires a dense vector to be used.  Sparse vectors must
	 * somehow be converted to dense vectors before this matrix may
	 * be applied to them.
	 *
	 * @param Vector LinBox dense vector type
	 * @param Switch switch object type
	 \ingroup blackbox
	 */
	template <class _Field, class Switch = CekstvSwitch<_Field>>
	class Butterfly : public BlackboxInterface {
	public:
		typedef _Field Field;
		typedef Butterfly<_Field, Switch> Self_t;
		typedef typename Field::Element Element;

		/** No-Op Constructor
		*/
		Butterfly (const Field &F, size_t n) :
			_field (&F), _VD (F), _n (n)
		{}



		/** Constructor from an integer and a switch object.
		 * The switch object is an object that is applied
		 * to two references to elements to switch them.  It must have both
		 * an apply and an applyTranspose method.
		 * It must contain all information needed by the switch other
		 * than the elements themselves.  This includes any random
		 * numbers or sequences of values.  It must also be able to
		 * be applied as many times as needed.  In particular, it must be able
		 * to create new random elements or repeat a stored sequence
		 * of values.
		 * This is not required by the abstract base class.
		 * @param n integer size of vectors to be applied to
		 * @param F
		 * @param factory switch predicate object object
		 */
		Butterfly (const Field &F, size_t n, typename Switch::Factory &factory);

		/* Destructor. */
		~Butterfly () {}


		/** Application of BlackBox matrix.
		 * <code>y = A*x</code>.
		 * Requires one vector conforming to the \ref LinBox
		 * vector @link Archetypes archetype@endlink.
		 * Required by abstract base class.
		 * For this matrix, this involves applying each switch in order to the
		 * input vector.
		 * @return reference to vector y containing output (after switching).
		 * @param  x constant reference to vector to contain input
		 * 			(before switching)
		 * @param y
		 */

		template<class OutVector, class InVector>
		OutVector& apply (OutVector& y, const InVector& x) const;

		/** Application of BlackBox matrix transpose.
		 * <code>y = transpose (A)*x</code>.
		 * Requires one vector conforming to the \ref LinBox
		 * vector @link Archetypes archetype@endlink.
		 * Required by abstract base class.
		 * For this matrix, this involves applying the transpose of each switch
		 * to the input vector in the reverse order of the apply function.
		 * @return reference to vector y containing output (after switching).
		 * @param  x constant reference to vector to contain input
		 * 			(before switching)
		 * @param y
		 */
		template<class OutVector, class InVector>
		OutVector& applyTranspose (OutVector& y, const InVector& x) const;

		template<typename _Tp1, typename _Sw1 = typename Switch::template rebind<_Tp1>::other>
		struct rebind {
			typedef Butterfly<_Tp1, _Sw1> other;

			void operator() (other & Ap, const Self_t& A) {
				//             other LAp(F,A._n);
				Ap.n_vec() = A.n_vec();
				Ap.l_vec() = A.l_vec();
				Ap.indices() = A.indices();

				typename std::vector<Switch>::const_iterator sit = A.switchesBegin();

				for( ; sit != A.switchesEnd(); ++sit) {
					_Sw1 newsw;
					typename Switch::template rebind<_Tp1>() (newsw, *sit, Ap.field(), A.field());
					Ap.switches().push_back( newsw );
				}
				//             Ap = new other(LAp);
			}
		};

		template<typename _Tp1, typename _Sw1>
		Butterfly (const Butterfly<_Tp1,_Sw1>& B, const Field &F) :
			_field (&F), _VD (F), _n (B.rowdim())
		{
			typename Butterfly<_Tp1,_Sw1>::template rebind<Field>() (*this, B);
		}



		/*- Retreive row dimensions of BlackBox matrix.
		 * This may be needed for applying preconditioners.
		 * Required by abstract base class.
		 * @return integer number of rows of black box matrix.
		 */
		size_t rowdim () const
		{ return _n; }

		/*- Retreive column dimensions of BlackBox matrix.
		 * Required by abstract base class.
		 * @return integer number of columns of black box matrix.
		 */
		size_t coldim () const
		{ return _n; }

		const Field& field() const
		{return *_field;}


		// Required for rebind
		// Don't know how to tell that rebind should be friend ...
		std::vector<size_t> n_vec() const
		{ return this->_n_vec; }
		std::vector<size_t> l_vec() const
		{ return this->_l_vec; }
		std::vector< std::pair< size_t, size_t > > indices() const
		{ return this->_indices; }
		std::vector<size_t>& n_vec() { return this->_n_vec; }
		std::vector<size_t>& l_vec() { return this->_l_vec; }
		std::vector< std::pair< size_t, size_t > >& indices() { return this->_indices; }
		typename std::vector<Switch>::const_iterator switchesBegin() const
		{ return this->_switches.begin();}
		typename std::vector<Switch>::const_iterator switchesEnd() const
		{ return this->_switches.end(); }
		std::vector<Switch>& switches() { return _switches; }


	private:


		// Field over which we are working
		const Field *_field;
		VectorDomain<Field> _VD;

		// Number of rows and columns of square matrix.
		size_t _n;

		// Vectors of sizes of sub-groups and number of levels in each
		// These may not need to be stored in general.
		// They may only be used in the constructor
		std::vector<size_t> _n_vec, _l_vec;

		// Vector of index pairs.  These are the indices to be switched with
		// a given switch.
		std::vector< std::pair< size_t, size_t > > _indices;

		// Vector of switches
		std::vector<Switch> _switches;

		// Build the vector of indices
		void buildIndices ();

	}; // template <class Field, class Vector> class Butterfly

	/** A function used with Butterfly Blackbox Matrices.
	 * This function takes an STL vector x of booleans, and returns
	 * a vector y of booleans such that setting the switches marked
	 * by true flags in y to be on (or to swap elements) the true
	 * elements x will be switched to a given contiguous block
	 * through the use of a Butterfly switching network.
	 * The integer parameter j marks where this block is to begin.
	 * If x has r true elements, the Butterfly switching network will place
	 * these elements in a contiguous block starting at j and ending at
	 * j + r - 1.
	 * Wrap around shall be considered to preserve contiguity.
	 * The value of j is defaulted to be zero, and it is only allowed to
	 * be non-zero is the size of x is a power of 2.
	 * @return vector of booleans for setting switches
	 * @param x vector of booleans marking elements to switch into
	 *	      contiguous block
	 * @param j offset of contiguous block
	 */
	inline std::vector<bool> setButterfly (const std::vector<bool>& x,
					       size_t j = 0);

} // namespace LinBox

#include "butterfly.inl"

#endif // __LINBOX_butterfly_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
