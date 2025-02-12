/* linbox/blackbox/butterfly.inl
 * Copyright (C) 1999-2001 William J Turner,
 *               2001 Bradford Hovinen
                 2015 reorg bds
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

#ifndef __LINBOX_butterfly_INL
#define __LINBOX_butterfly_INL

#include <vector>
#include "linbox/blackbox/blackbox-interface.h"
#include "linbox/vector/vector-domain.h"
#include "linbox/field/hom.h"

/** @file blackbox/butterfly.inl
 *
 * First linbox block: butterfly method implementations
 * Second LinBox block: butterfly switch and switch factory impls.
 *
 */
namespace LinBox
{
	// Implementation of Butterfly methods

	template <class Field, class Switch>
	inline Butterfly<Field, Switch>::Butterfly (const Field &F, size_t n, typename Switch::Factory &factory) :
		_field (&F), _VD (F), _n (n)
	{
		buildIndices ();

		for (unsigned int i = 0; i < _indices.size (); ++i)
			_switches.push_back (factory.makeSwitch ());
	}

	template <class Field, class Switch>
	template<class OutVector, class InVector>
	inline OutVector& Butterfly<Field, Switch>::apply (OutVector& y, const InVector& x) const
	{
		std::vector< std::pair<size_t, size_t> >::const_iterator idx_iter = _indices.begin ();
		typename std::vector<Switch>::const_iterator switch_iter = _switches.begin ();

		_VD.copy (y, x);

		for (; idx_iter != _indices.end (); ++idx_iter, ++switch_iter)
			switch_iter->apply (field(), y[idx_iter->first], y[idx_iter->second]);

		return y;
	}

	template <class Field, class Switch>
	template <class OutVector, class InVector>
	inline OutVector& Butterfly<Field, Switch>::applyTranspose (OutVector& y, const InVector& x) const
	{
		std::vector< std::pair<size_t, size_t> >::const_reverse_iterator idx_iter = _indices.rbegin ();
		typename std::vector<Switch>::const_reverse_iterator switch_iter = _switches.rbegin ();

		_VD.copy (y, x);

		for (; idx_iter != _indices.rend (); ++idx_iter, ++switch_iter)
			switch_iter->applyTranspose (field(), y[idx_iter->first], y[idx_iter->second]);

		return y;
	}

	template <class Field, class Switch>
	void Butterfly<Field, Switch>::buildIndices ()
	{
		for (size_t value (_n), l_p (0), n_p (1);
		     n_p != 0;
		     value >>= 1, ++l_p, n_p <<= 1)
		{
			if (value & 1) {
				_l_vec.push_back (l_p);
				_n_vec.push_back (n_p);
			}
		}

		// Create vector of indices to switch
		size_t n_p ;   	// size of group and number of levels in group
		size_t level (0), difference (1);	// track levels done for powers of 2

		// Vector containing indices for last level of last power of 2.
		std::vector< std::pair< size_t, size_t > > p_ind;

		// Vector and iterator used for computing p_ind.
		std::vector< std::pair< size_t, size_t > > temp_ind;
		std::vector< std::pair< size_t, size_t > >::iterator iter;

		// Loop over sub-groups of powers of two
		for (size_t p (0), start_index (0);
		     p < _n_vec.size ();
		     ++p, start_index += n_p)
		{
			// update size
			n_p = _n_vec[p];
			size_t l_p = _l_vec[p];

			// loop over levels of sub-group network
			for ( ; level < l_p; ++level, difference <<= 1) {
				// Create
				temp_ind = p_ind;

				// the second sub group is a shift of the first
				for (iter = temp_ind.begin (); iter != temp_ind.end (); ++iter) {
					iter->first += difference;
					iter->second += difference;
				}

				// add the second group to the first
				p_ind.insert (p_ind.end (), temp_ind.begin (), temp_ind.end ());

				// add switches to mix the two sub groups
				temp_ind = std::vector< std::pair<size_t, size_t> >
				(difference, std::pair<size_t, size_t> (0, 0));

				size_t i = 0;
				for (iter = temp_ind.begin (); iter != temp_ind.end (); ++i, ++iter) {
					iter->first += i;
					iter->second += i + difference;
				}

				// add the combining group to the first and second
				p_ind.insert (p_ind.end (), temp_ind.begin (), temp_ind.end ());
			}

			// Add this level to total list of indices and correct starting point
			temp_ind = p_ind;

			for (iter = temp_ind.begin (); iter != temp_ind.end (); ++iter) {
				iter->first += start_index;
				iter->second += start_index;
			}

			_indices.insert (_indices.end (), temp_ind.begin (), temp_ind.end ());

			// Combine everything so far
			temp_ind = std::vector< std::pair<size_t, size_t> > (start_index, std::pair<size_t, size_t> (0, 0));

			iter = temp_ind.begin ();
			for (size_t index = 0; index < start_index; ++index, ++iter) {
				iter->first = index;
				iter->second += index + n_p;
			}

			_indices.insert (_indices.end (), temp_ind.begin (), temp_ind.end ());
		}
	}

	inline std::vector<bool> setButterfly (const std::vector<bool>& x,
					       size_t j)
	{
		size_t n = x.size ();

		commentator().start ("Setting butterfly switches", "setButterfly");

		std::ostream &report = commentator().report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);

		report << "Called set switches with vector of size " << n
		<< " and offset " << j << std::endl;

		// return empty vector if zero or one elements in x because
		// no switching will be done.
		if (x.size () <= 1) {
			commentator().indent (report);
			report << "No switches needed. Returning with empty vector." << std::endl;

			commentator().stop ("done");
			return std::vector<bool> ();
		}

		commentator().indent (report);
		report << "Counting the number of switches that exist." << std::endl;

		// break inputs into groups of size powers of 2.
		// calculate size of groups, and powers of 2 that give sizes
		// store these values in vectors n and l, respectively
		std::vector<size_t> l_vec, n_vec;

		for (size_t value (n), l_p (0), n_p (1);
		     n_p != 0;
		     value >>= 1, ++l_p, n_p <<= 1)
		{
			commentator().indent (report);
			report << "  looping at value = " << value
			<< ", l_p = " << l_p
			<< ", n_p = " << n_p << std::endl;

			if (value & 1) {
				l_vec.push_back (l_p);
				n_vec.push_back (n_p);

				commentator().indent (report);
				report << "    inserted value = " << value
				<< ", l_p = " << l_p
				<< ", n_p = " << n_p << std::endl;
			}
		}

		// Calculate total number of switches required
		size_t s (0);

		for (size_t ii = 0; ii < n_vec.size (); ++ii)
			s += n_vec[ii] * l_vec[ii] / 2;

		for (size_t ii = 0; ii < n_vec.size () - 1; ++ii)
			for (size_t jj = 0; jj <= ii; ++jj)
				s += n_vec[jj];

		commentator().indent (report);
		report << "There are a total of " << s << " switches" << std::endl;

		// Set largest power of 2 in decomposition of n = x.size ()
		size_t n_p (*n_vec.rbegin ());

		commentator().indent (report);
		report << "Found largest power of 2 in decomposition of " << n
		<< " as n_p = " << n_p << std::endl;

		if ( (n != n_p) && (j != 0) ) {
			commentator().indent (report);
			report << "Non-zero offset " << j
			<< " used with non-power size."
			<< "Offset reset to zero." << std::endl;

			j = 0;
		}
		else
			j %= n;

		if (n == n_p) {
			n_p /= 2;	  // >> is not portable!

			commentator().indent (report);
			report << "n = " << n << " is a power of two.  "
			<< "Resetting n_p to be half of n: n_p = " << n_p << std::endl;
		}

		// count true elements not in largest power of 2 block
		size_t r_1(0);

		for (std::vector<bool>::const_iterator iter = x.begin ();
		     iter != x.begin () + (ptrdiff_t)(n - n_p);
		     ++iter)
			if (*iter) ++r_1;

		// count total number of true elements in x.
		size_t r (r_1);

		for (std::vector<bool>::const_iterator iter = x.begin () + (ptrdiff_t)(n - n_p);
		     iter != x.end ();
		     ++iter)
			if (*iter) ++r;

		commentator().indent (report);
		report << "The vector x will be broken into two sub-vectors,"
		<< "x_1 = x[0,...," << n - n_p - 1 << "] and x_2 = x["
		<< n - n_p << ",...," << n - 1 << "]."
		<< "There are a total of " << r << " true Elements in x, "
		<< r_1 << " of which occured in the first sub-vector."
		<< "The output vector will have " << s << " entries and will"
		<< "switch the true Elements of x into a contiguous block"
		<< "[" << j << "," << j + r
		<< ") = [" << j << "," << j + r - 1<< "]." << std::endl;

		if (r == 0) {
			commentator().indent (report);
			report << "There are no true Elements in x, so the recursion is"
			<< "being broken and a vector of false flags returned." << std::endl;

			commentator().stop ("done");
			return std::vector<bool> (s, false);
		}
		else if (r == n) {
			commentator().indent (report);
			report << "There are no false Elements in x, so the recursion is"
			<< "being broken and a vector of false flags returned." << std::endl;

			commentator().stop ("done");
			return std::vector<bool> (s, false);
		}

		// Calculate where the true elements are supposed to end up
		// Here, they will be in a contiguous block starting after the
		// offset.  s_1 are the true elements after the offset and in the first
		// sub-group, s_2 are the ones in the second sub group, and s_3 are the
		// elements that wrap around to the beginning.  s_1 and s_3 cannot both
		// be non-zero unless s_2 == n_p.  (I.e., the second group is full.)
		// Also, because for n != 2 n_p the offset is zero, in that case
		// s_3 must be zero.  Any of them may be zero if the corrsponding block
		// is empty.
		// s_2 is only used for tracing the program, so it is not always computed.

		size_t s_1;

		if (j < n - n_p) {
			if (j + r < n - n_p)
				s_1 = r;
			else
				s_1 = n - n_p - j;
		}
		else
			s_1 = 0;

		size_t s_2 = 0;

		if (commentator().isPrinted (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)) {
			if (j + r < n - n_p)
				s_2 = 0;
			else {
				if (j + r < n)
					s_2 = j + r;
				else
					s_2 = n;

				if (j < n - n_p)
					s_2 -= (n - n_p);
				else
					s_2 -= j;
			}
		}

		size_t s_3 = ((j + r) > n) ? j + r - n : 0;

		commentator().indent (report);
		report << "The number of Elements in each of the three blocks of "
		<< "true Elements in the end result are"
		<< "s_1 = " << s_1
		<< ", s_2 = " << s_2
		<< ", and s_3 = " << s_3 << "." << std::endl;

		// Create empty vector for output. y_temp is used to retrieve output
		// from recursion before inserting into output.
		std::vector<bool> y_1, y_2, y_3 = std::vector<bool> (n - n_p, false);

		if ((s_1 + s_3) == r_1) {
			commentator().indent (report);
			report << "Case I: s_1 + s_3 == r_1 and s_2 == r - r_1."
			<< "No Elements are moved between the two sub-vectors." << std::endl;

			if (j < (n - n_p)) {
				commentator().indent (report);
				report << "  A: j < (n - n_p).  j_1 = j = " << j << ", j_2 = 0";

				y_1 = setButterfly (std::vector<bool>(x.begin (), x.begin () + (ptrdiff_t)(n - n_p)), j);
				y_2 = setButterfly (std::vector<bool>(x.begin () + (ptrdiff_t)(n - n_p), x.end ()), 0);

			}
			else {
				commentator().indent (report);
				report << "  A: j >= (n - n_p).  j_1 = 0, j_2 = j - (n - n_p) = "
				<< j - (n - n_p) << std::endl;

				// This case cannot occur for n != 2*n_p because j != 0

				y_1 = setButterfly (std::vector<bool>(x.begin (), x.begin () + (ptrdiff_t)(n - n_p)), 0);
				y_2 = setButterfly (std::vector<bool>(x.begin () + (ptrdiff_t)(n - n_p), x.end ()), j - (n - n_p));
			}
		}
		else if ((s_1 + s_3) > r_1) {
			commentator().indent (report);
			report << "Case II: s_1 + s_3 > r_1 and s_2 < r - r_1."
			<< "Elements are moved from the right sub-vector to the left." << std::endl;

			// This means that s_2 < n_p, so either s_1 = 0 or s_3 = 0 (or both).

			if (j < (n - n_p)) {
				commentator().indent (report);
				report << "  A: j < (n - n_p).  j_1 = j, j_2 = 2*n_p + j + r_1 - n = "
				<< 2*n_p + j + r_1 - n << std::endl;

				// In this case, s_1 > 0, so s_3 = 0, and wrap-around cannot occur.

				y_1 = setButterfly (std::vector<bool>(x.begin (), x.begin () + (ptrdiff_t)(n - n_p)), j);
				y_2 = setButterfly (std::vector<bool>(x.begin () + (ptrdiff_t)(n - n_p), x.end ()), 2*n_p + j + r_1 - n);

				for (std::vector<bool>::iterator iter = (y_3.begin () + (ptrdiff_t)(j + r_1));
				     iter != (y_3.begin () + (ptrdiff_t)(n - n_p));
				     ++iter)
					*iter = true;
			}
			else {
				commentator().indent (report);
				report << "  A: j >= (n - n_p).  j_1 = j + r - n - r_1 = "
				<< j + r - n - r_1 << ", j_2 = j - (n - n_p) = "
				<< j - (n - n_p) << std::endl;

				// In this case, s_1 = 0, so s_3 >= 0, and wrap-around may occur.
				// This case cannot occur for n != 2*n_p because j != 0.

				y_1 = setButterfly (std::vector<bool>(x.begin (), x.begin () + (ptrdiff_t)(n - n_p)), j + r - n - r_1);
				y_2 = setButterfly (std::vector<bool>(x.begin () + (ptrdiff_t)(n - n_p), x.end ()), j - (n - n_p));

				for (std::vector<bool>::iterator iter = y_3.begin ();
				     iter != (y_3.begin () + (ptrdiff_t)(j + r - n - r_1));
				     ++iter)
					*iter = true;
			}
		}
		else if ((s_1 + s_3) < r_1) {
			commentator().indent (report);
			report << "Case III: s_1 + s_3 < r_1 and s_2 > r - r_1."
			<< "Elements are moved from the left sub-vector to the right." << std::endl;

			// This case also means that s_1 + s_3 < n - n_p, or the contiguous
			// block cannot encompass the entire first sub-vector.  For this
			// reason, this case is not considered when n != 2*n_p (when j = 0).

			if (j < (n - n_p)) {
				commentator().indent (report);
				report << "  A: j < (n - n_p).  j_1 = j = " << j
				<< ", j_2 = j + r_1 - n + n_p = " << j + r_1 - n + n_p << std::endl;
				// In this case, s_1 > 0, so s_3 = 0, and wrap-around cannot occur.

				y_1 = setButterfly (std::vector<bool>(x.begin (), x.begin () +(ptrdiff_t) (n - n_p)), j);
				y_2 = setButterfly (std::vector<bool>(x.begin () +(ptrdiff_t) (n - n_p), x.end ()), j + r_1 - n + n_p);

				for (std::vector<bool>::iterator iter = (y_3.begin () +(ptrdiff_t) s_3);
				     iter != (y_3.begin () + (ptrdiff_t)(j + r_1 - n + n_p));
				     ++iter)
					*iter = true;
			}
			else {
				commentator().indent (report);
				report << "  A: j >= (n - n_p).  j_1 = j + r - n_p - r_1 = "
				<< j + r - n_p - r_1 << ", j_2 = j - (n - n_p) = "
				<< j - (n - n_p) << std::endl;

				// In this case, s_1 = 0, so s_3 >= 0, and wrap-around may occur.
				// This case cannot occur for n != 2*n_p because j != 0.

				y_1 = setButterfly (std::vector<bool>(x.begin (), x.begin () + (ptrdiff_t)(n - n_p)), j + r - n_p - r_1);
				y_2 = setButterfly (std::vector<bool>(x.begin () + (ptrdiff_t)(n - n_p), x.end ()), j - (n - n_p));

				for (std::vector<bool>::iterator iter (y_3.begin () + (ptrdiff_t)(j + r - n_p - r_1));
				     iter != (y_3.begin () + (ptrdiff_t)(n - n_p));
				     ++iter)
					*iter = true;
			}
		}

		// Create output vector.
		std::vector<bool> y (y_1);
		y.insert (y.end (), y_2.begin (), y_2.end ());
		y.insert (y.end (), y_3.begin (), y_3.end ());

		commentator().indent (report);
		report << "The output vector for n = " << n << " has " << y.size ()
		<< " entries."
		<< "  " << y_1.size () << " from the first sub-vector"
		<< "  " << y_2.size () << " from the second sub-vector"
		<< "  " << y_3.size () << " from recombining the two"
		<< "And the output vector y is:"
		<< "-------------------------- " << std::endl;

		for (size_t i = 0; i < y.size (); ++i) {
			commentator().indent (report);
			report << "  " << i << ": " << y[i] << std::endl;
		}

		commentator().indent (report);
		report << "-------------------------- " << std::endl;

		commentator().stop ("done");

		return y;

	} // std::vector<bool> setButterfly (const std::vector<bool>& x, size_t j)

	//@}

    // Begin cekstv switch

	template <class Field> class CekstvSwitchFactory;

	template <class Field>
	class CekstvSwitch {
	public:

		/// Typedef
		typedef typename Field::Element Element;
		typedef CekstvSwitch<Field> Self_t;
		typedef CekstvSwitchFactory<Field> Factory;

		CekstvSwitch () {}

		/** Constructor from a field and a field element.
		 * @param a vector of switches
		 */
		CekstvSwitch (const typename Field::Element &a) :
		       	_a (a)
		{}

		/** Destructor.
		*/
		~CekstvSwitch () {}

		/** Apply switch function.
		 * Switches the elements in references according to the
		 * exchange matrix introduced in "Efficient Matrix
		 * Preconditioners for Black Box Linear Algebra" by Chen, Eberly,
		 * Kaltofen, Saunders, Turner, and Villard and the current field element
		 * specified in the switch object.
		 * @return bool true if swapped, false otherwise
		 * @param F
		 * @param x reference to first element to be switched
		 * @param y reference to second element to be switched
		 */
		bool apply (const Field &F, Element &x, Element &y) const;

		/** Apply switch transpose function.
		 * Switches the elements in references according to the
		 * transpose of the exchange matrix introduced in "Efficient Matrix
		 * Preconditioners for Black Box Linear Algebra" by Chen, Eberly,
		 * Kaltofen, Saunders, Turner, and Villard and the current field element
		 * specified in the switch object.
		 * @return bool true if swapped, false otherwise
		 * @param F
		 * @param x reference to first element to be switched
		 * @param y reference to second element to be switched
		 */
		bool applyTranspose (const Field &F, Element &x, Element &y) const;


		template<typename _Tp1>
		struct rebind
		{
			typedef CekstvSwitch<_Tp1> other;

			// special rebind operator() with two fields,
			// indeed local field is not stored in the switch
			void operator() (other & Ap, const Self_t& A, const _Tp1& T, const Field& F) {
				Hom<Field, _Tp1>(F,T).image(Ap.getData(), A.getData());
			}
		};

		typename Field::Element& getData() { return _a; }
		const typename Field::Element& getData() const { return _a; }


	private:

		// Parameter of this 2x2 block
		typename Field::Element _a;
	};

	template <class Field>
	class CekstvSwitchFactory {
	public:
		/** Constructor from an STL vector of bools
		*/
		CekstvSwitchFactory (typename Field::RandIter r) :
		       	_r (r)
		{}

		/** Construct and return a boolean switch object
		*/
		CekstvSwitch<Field> makeSwitch ()
		{ typename Field::Element a; return CekstvSwitch<Field> (_r.random (a)); }

	private:

		typename Field::RandIter _r;
	};

	template <class Field>
	inline bool CekstvSwitch<Field>::apply (const Field             &F,
						typename Field::Element &x,
						typename Field::Element &y) const
	{
		F.axpyin (x, _a, y);
		F.addin (y, x);

		return true;
	}

	template <class Field>
	inline bool CekstvSwitch<Field>::applyTranspose (const Field             &F,
							 typename Field::Element &x,
							 typename Field::Element &y) const
	{
		F.addin (x, y);
		F.axpyin (y, _a, x);

		return true;
	}

// End cekstv switch
#if 0
// Begin specialization of cekstv switch object
	template <>
	class CekstvSwitch<GF2>
	{
	public:
		typedef GF2 Field;
		/// Typedef
		typedef Field::Element Element;
		typedef CekstvSwitch<Field> Self_t;
		typedef CekstvSwitchFactory<Field> Factory;

		/** Constructor from a field and a field element.
		 * @param F field in which arithmetic is done
		 * @param switches vector of switches
		 */
		CekstvSwitch (const Field::Element &a) :
			_a (a)
		{}

		~CekstvSwitch () {}

		bool apply (const Field &F, Element &x, Element &y) const
		{
			F.axpyin (x, _a, y);
			F.addin (y, x);
			return true;
		}

		bool applyTranspose (const Field &F, Element &x, Element &y) const
		{
			F.addin (x, y);
			F.axpyin (y, _a, x);
			return true;
		}

		bool apply (const Field &F, stdBitReference x, stdBitReference y) const
		{
			F.axpyin (x, _a, y);
			F.addin (y, x);
			return true;
		}

		bool applyTranspose (const Field &F, stdBitReference x, stdBitReference y) const
		{
			F.addin (x, y);
			F.axpyin (y, _a, x);
			return true;
		}

		template<typename _Tp1>
		struct rebind
		{
			typedef CekstvSwitch<_Tp1> other;

			// special rebind operator() with two fields,
			// indeed local field is not stored in the switch
			void operator() (other *& Ap, const Self_t& A, const _Tp1& T, const Field& F) {
				typename _Tp1::Element u;
				Hom<Field, _Tp1>(F,T).image(u, A._a);
				Ap = new other(u);
			}
		};


	private:

		// Parameter of this 2x2 block
		Field::Element _a;
	};


// End, specialization of cekstv switch object
#endif
// Begin boolean switch

	class BooleanSwitchFactory;

	/** Boolean switch object.
	 * This is a switch predicate object that is applied
	 * to two references to elements to switch them as needed
	 * by the \ref Butterfly\ Switching\ Network\ BlackBox\ Matrix\ Object.
	 */
	class BooleanSwitch {
	public:

		typedef BooleanSwitch Self_t;
		typedef BooleanSwitchFactory Factory;

		/** Constructor from an STL vector of booleans.
		 * The switch is applied using the vector of booleans.
		 * A true value means to swap the two elements, and a false
		 * value means not to.
		 * The apply function starts at the beginning of the vector moving
		 * forward through it, and applyTranspose function starts at the end
		 * moving backwards.  Both repeat the vector after they pass through it.
		 * @param s vector of switches
		 */
		BooleanSwitch (const bool s) :
		       	_s (s)
		{}

		/** Destructor.
		*/
		~BooleanSwitch () {}

		/** Apply switch function.
		 * Switches the elements in references according to current boolean
		 * value.  Swaps the elements if boolean is true, otherwise does nothing.
		 * It is templatized by the element type to be swapped.
		 * @return bool \c true if swapped, \c false otherwise
		 * @param F
		 * @param x reference to first element to be switched
		 * @param y reference to second element to be switched
		 */
		template <class Field>
		bool apply (const Field             &F,
			    typename Field::Element &x,
			    typename Field::Element &y) const;

		/** Apply switch transpose function.
		 * Switches the elements in references according to current boolean
		 * value.  Swaps the elements if boolean is true, otherwise does nothing.
		 * It is templatized by the element type to be swapped.
		 * @return bool \c true if swapped, \c false otherwise
		 * @param F
		 * @param x reference to first element to be switched
		 * @param y reference to second element to be switched
		 */
		template <class Field>
		bool applyTranspose (const Field             &F,
				     typename Field::Element &x,
				     typename Field::Element &y) const;

		template<typename _Tp1>
		struct rebind {
			typedef BooleanSwitch other;

		};

	protected:

		bool _s;

	}; // class boolean_switch

	/** Boolean switch factory
	 *
	 * This class facilitates construction of boolean switch objects by the
	 * butterfly matrix.
	 */
}// linbox
namespace LinBox {

	class BooleanSwitchFactory {
	public:
		/** Constructor from an STL vector of bools
		*/
		BooleanSwitchFactory (const std::vector<bool> &switches) :
		       	_switches (switches), _iter (switches.begin ())
		{}

		/** Construct and return a boolean switch object
		 *
		 * This function walks through the switches object given in the
		 * constructor, advancing on each invocation. It wraps around to the
		 * beginning of the vector when it reaches the end.
		 */
		BooleanSwitch makeSwitch ()
		{
			if (_iter == _switches.end ())
				_iter = _switches.begin ();

			return BooleanSwitch (*_iter++);
		}

	private:

		const std::vector<bool> &_switches;
		std::vector<bool>::const_iterator _iter;
	};

	template <class Field>
	inline bool BooleanSwitch::apply (const Field             &F,
					  typename Field::Element &x,
					  typename Field::Element &y) const
	{
		if (_s)
			std::swap (x, y);

		return _s;
	}

	template <class Field>
	inline bool BooleanSwitch::applyTranspose (const Field             &F,
						   typename Field::Element &x,
						   typename Field::Element &y) const
	{
		if (_s)
			std::swap (x, y);

		return _s;
	}

// End boolean switch
}// namespace LinBox

#endif // __LINBOX_butterfly_INL

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
