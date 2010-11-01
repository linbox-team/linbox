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
 * See COPYING for license information
 */

#ifndef __LINBOX_butterfly_H
#define __LINBOX_butterfly_H

#include <vector>
#include <linbox/blackbox/blackbox-interface.h>


// Namespace in which all LinBox library code resides
namespace LinBox
{

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
 * See \ref{BlackboxArchetype} for the specification of these methods.
 *
 * This matrix requires a dense vector to be used.  Sparse vectors must
 * somehow be converted to dense vectors before this matrix may
 * be applied to them.
 *
 * @param Vector LinBox dense vector type
 * @param Switch switch object type
\ingroup blackbox
 */
template <class _Field, class Switch>
class Butterfly : public BlackboxInterface
{
    public:
	typedef _Field Field;
    	typedef Butterfly<_Field, Switch> Self_t;
	typedef typename Field::Element Element;

	/** No-Op Constructor 
         */
    	Butterfly (const Field &F, size_t n) : _F (F), _VD (F), _n (n) {}
    
  

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
	 * @param S switch predicate object object
	 */
	Butterfly (const Field &F, size_t n, typename Switch::Factory &factory);

	/* Destructor. */
	~Butterfly () {}


	/*- Application of BlackBox matrix.
	 * y = A*x.
	 * Requires one vector conforming to the \ref{LinBox}
	 * vector {@link Archetypes archetype}.
	 * Required by abstract base class.
	 * For this matrix, this involves applying each switch in order to the 
	 * input vector.
	 * @return reference to vector y containing output (after switching).
	 * @param  x constant reference to vector to contain input 
	 * 			(before switching)       
	*/

	template<class OutVector, class InVector>
	OutVector& apply (OutVector& y, const InVector& x) const;

	/*- Application of BlackBox matrix transpose.
	 * y = transpose (A)*x.
	 * Requires one vector conforming to the \ref{LinBox}
	 * vector {@link Archetypes archetype}.
	 * Required by abstract base class.
	 * For this matrix, this involves applying the transpose of each switch 
	 * to the input vector in the reverse order of the apply function.
	 * @return reference to vector y containing output (after switching).
	 * @param  x constant reference to vector to contain input 
	 * 			(before switching)
	 */
	template<class OutVector, class InVector>
	OutVector& applyTranspose (OutVector& y, const InVector& x) const;

    template<typename _Tp1, typename _Sw1 = typename Switch::template rebind<_Tp1>::other>
    struct rebind
    { 
        typedef Butterfly<_Tp1, _Sw1> other;

        void operator() (other & Ap, const Self_t& A, const _Tp1& F) {
//             other LAp(F,A._n);
            Ap.n_vec() = A.n_vec();
            Ap.l_vec() = A.l_vec();
            Ap.indices() = A.indices();
            
            typename std::vector<Switch>::const_iterator sit = A.switchesBegin();
            
            for( ; sit != A.switchesEnd(); ++sit) {
                _Sw1 newsw;
                typename Switch::template rebind<_Tp1>() (newsw, *sit, F, A._F);
                Ap.switches().push_back( newsw );
            }
//             Ap = new other(LAp);
        }  
    };
      
    template<typename _Tp1, typename _Sw1>
    Butterfly (const Butterfly<_Tp1,_Sw1>& B, const Field &F) : _F (F), _VD (F), _n (B.rowdim()) {
        typename Butterfly<_Tp1,_Sw1>::template rebind<Field>() (*this, B, F);
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

	const Field& field() const {return _F;}


        // Required for rebind
        // Don't know how to tell that rebind should be friend ...
    	std::vector<size_t> n_vec() const { return this->_n_vec; }
        std::vector<size_t> l_vec() const { return this->_l_vec; }
    	std::vector< std::pair< size_t, size_t > > indices() const { return this->_indices; }
    	std::vector<size_t>& n_vec() { return this->_n_vec; }
        std::vector<size_t>& l_vec() { return this->_l_vec; }
    	std::vector< std::pair< size_t, size_t > >& indices() { return this->_indices; }
    	typename std::vector<Switch>::const_iterator switchesBegin() const { return this->_switches.begin();}
     	typename std::vector<Switch>::const_iterator switchesEnd() const { return this->_switches.end(); }
    	std::vector<Switch>& switches() { return _switches; }
    

    private:


	// Field over which we are working
	const Field _F;
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

// Implementation of methods

template <class Field, class Switch>
inline Butterfly<Field, Switch>::Butterfly (const Field &F, size_t n, typename Switch::Factory &factory)
	: _F (F), _VD (F), _n (n)
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
		switch_iter->apply (_F, y[idx_iter->first], y[idx_iter->second]);

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
		switch_iter->applyTranspose (_F, y[idx_iter->first], y[idx_iter->second]);

	return y;
}

template <class Field, class Switch>
void Butterfly<Field, Switch>::buildIndices () 
{
	for (size_t value (_n), l_p (0), n_p (1); 
	     n_p != 0; 
	     value >>= 1, l_p++, n_p <<= 1)
	{
		if (value & 1) {
			_l_vec.push_back (l_p);
			_n_vec.push_back (n_p);      
		}
	}

	// Create vector of indices to switch
	size_t n_p, l_p;   	// size of group and number of levels in group
	size_t level (0), difference (1);	// track levels done for powers of 2

	// Vector containing indices for last level of last power of 2.
	std::vector< std::pair< size_t, size_t > > p_ind;

	// Vector and iterator used for computing p_ind.
	std::vector< std::pair< size_t, size_t > > temp_ind;
	std::vector< std::pair< size_t, size_t > >::iterator iter;

	// Loop over sub-groups of powers of two
	for (size_t p (0), start_index (0); 
	     p < _n_vec.size (); 
	     p++, start_index += n_p)
	{
		// update size
		n_p = _n_vec[p];
		l_p = _l_vec[p];

		// loop over levels of sub-group network
		for ( ; level < l_p; level++, difference <<= 1) {
			// Create 
			temp_ind = p_ind;

			// the second sub group is a shift of the first
			for (iter = temp_ind.begin (); iter != temp_ind.end (); iter++) {
				iter->first += difference;
				iter->second += difference;
			}

			// add the second group to the first
			p_ind.insert (p_ind.end (), temp_ind.begin (), temp_ind.end ());

			// add switches to mix the two sub groups
			temp_ind = std::vector< std::pair<size_t, size_t> >
				(difference, std::pair<size_t, size_t> (0, 0));

			size_t i = 0;
			for (iter = temp_ind.begin (); iter != temp_ind.end (); i++, iter++) {
				iter->first += i;
				iter->second += i + difference;
			}

			// add the combining group to the first and second
			p_ind.insert (p_ind.end (), temp_ind.begin (), temp_ind.end ());
		}

		// Add this level to total list of indices and correct starting point
		temp_ind = p_ind;

		for (iter = temp_ind.begin (); iter != temp_ind.end (); iter++) {
			iter->first += start_index;
			iter->second += start_index;
		}

		_indices.insert (_indices.end (), temp_ind.begin (), temp_ind.end ());

		// Combine everything so far
		temp_ind = std::vector< std::pair<size_t, size_t> > (start_index, std::pair<size_t, size_t> (0, 0));

		iter = temp_ind.begin ();
		for (size_t index = 0; index < start_index; index++, iter++) {
			iter->first = index;
			iter->second += index + n_p;
		}

		_indices.insert (_indices.end (), temp_ind.begin (), temp_ind.end ());
	}
}

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
 * @param log reference to ostream for logging
 */
inline std::vector<bool> setButterfly (const std::vector<bool>& x, 
				size_t j = 0)
{
	size_t n = x.size ();
 
	commentator.start ("Setting butterfly switches", "setButterfly");

	std::ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);

	report << "Called set switches with vector of size " << n
	       << " and offset " << j << std::endl;

	// return empty vector if zero or one elements in x because
	// no switching will be done.
	if (x.size () <= 1) {
		commentator.indent (report);
		report << "No switches needed. Returning with empty vector." << std::endl;

		commentator.stop ("done");
		return std::vector<bool> ();
	}

	commentator.indent (report);
	report << "Counting the number of switches that exist." << std::endl;
 
	// break inputs into groups of size powers of 2.
	// calculate size of groups, and powers of 2 that give sizes
	// store these values in vectors n and l, respectively
	std::vector<size_t> l_vec, n_vec;

	for (size_t value (n), l_p (0), n_p (1);
	     n_p != 0;
	     value >>= 1, l_p++, n_p <<= 1)
	{
		commentator.indent (report);
		report << "  looping at value = " << value
		       << ", l_p = " << l_p
		       << ", n_p = " << n_p << std::endl;
 
		if (value & 1) {
			l_vec.push_back (l_p);
			n_vec.push_back (n_p);

			commentator.indent (report);
			report << "    inserted value = " << value
			       << ", l_p = " << l_p
			       << ", n_p = " << n_p << std::endl;
		}
	}
 
	// Calculate total number of switches required
	size_t s (0);
 
	for (size_t ii = 0; ii < n_vec.size (); ii++)
		s += n_vec[ii] * l_vec[ii] / 2;

	for (size_t ii = 0; ii < n_vec.size () - 1; ii++)
		for (size_t jj = 0; jj <= ii; jj++)
			s += n_vec[jj];
 
	commentator.indent (report);
	report << "There are a total of " << s << " switches" << std::endl;
 
	// Set largest power of 2 in decomposition of n = x.size ()
	size_t n_p (*n_vec.rbegin ());

	commentator.indent (report);
	report << "Found largest power of 2 in decomposition of " << n
	       << " as n_p = " << n_p << std::endl;

	if ( (n != n_p) && (j != 0) ) {
		commentator.indent (report);
		report << "Non-zero offset " << j
		       << " used with non-power size."
		       << "Offset reset to zero." << std::endl;

		j = 0;
	} else
		j %= n;

	if (n == n_p) {
		n_p /= 2;	  // >> is not portable!

		commentator.indent (report);
		report << "n = " << n << " is a power of two.  "
		       << "Resetting n_p to be half of n: n_p = " << n_p << std::endl;
	}

	// count true elements not in largest power of 2 block
	size_t r_1(0);
 
	for (std::vector<bool>::const_iterator iter = x.begin ();
	     iter != x.begin () + (n - n_p);
	     iter++)
		if (*iter) r_1++;

	// count total number of true elements in x.
	size_t r (r_1);
 
	for (std::vector<bool>::const_iterator iter = x.begin () + (n - n_p);
	     iter != x.end ();
	     iter++)
		if (*iter) r++;

	commentator.indent (report);
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
		commentator.indent (report);
		report << "There are no true Elements in x, so the recursion is"
		       << "being broken and a vector of false flags returned." << std::endl;

		commentator.stop ("done");
		return std::vector<bool> (s, false);
	}
	else if (r == n) {
		commentator.indent (report);
		report << "There are no false Elements in x, so the recursion is"
		       << "being broken and a vector of false flags returned." << std::endl;

		commentator.stop ("done");
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
	} else
		s_1 = 0;

	size_t s_2 = 0;

	if (commentator.isPrinted (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)) {
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

	commentator.indent (report);
	report << "The number of Elements in each of the three blocks of "
	       << "true Elements in the end result are"
	       << "s_1 = " << s_1
	       << ", s_2 = " << s_2
	       << ", and s_3 = " << s_3 << "." << std::endl;

	// Create empty vector for output. y_temp is used to retrieve output
	// from recursion before inserting into output.
	std::vector<bool> y_1, y_2, y_3 = std::vector<bool> (n - n_p, false);

	if ((s_1 + s_3) == r_1) {
		commentator.indent (report);
		report << "Case I: s_1 + s_3 == r_1 and s_2 == r - r_1."
		       << "No Elements are moved between the two sub-vectors." << std::endl;

		if (j < (n - n_p)) {
			commentator.indent (report);
			report << "  A: j < (n - n_p).  j_1 = j = " << j << ", j_2 = 0";

			y_1 = setButterfly (std::vector<bool>(x.begin (), x.begin () + (n - n_p)), j);
			y_2 = setButterfly (std::vector<bool>(x.begin () + (n - n_p), x.end ()), 0);
 
		} else {
			commentator.indent (report);
			report << "  A: j >= (n - n_p).  j_1 = 0, j_2 = j - (n - n_p) = "
			       << j - (n - n_p) << std::endl;

			// This case cannot occur for n != 2*n_p because j != 0

			y_1 = setButterfly (std::vector<bool>(x.begin (), x.begin () + (n - n_p)), 0);
			y_2 = setButterfly (std::vector<bool>(x.begin () + (n - n_p), x.end ()), j - (n - n_p));
		}
	}
	else if ((s_1 + s_3) > r_1) {
		commentator.indent (report);
		report << "Case II: s_1 + s_3 > r_1 and s_2 < r - r_1."
		       << "Elements are moved from the right sub-vector to the left." << std::endl;

		// This means that s_2 < n_p, so either s_1 = 0 or s_3 = 0 (or both).
 
		if (j < (n - n_p)) {
			commentator.indent (report);
			report << "  A: j < (n - n_p).  j_1 = j, j_2 = 2*n_p + j + r_1 - n = "
			       << 2*n_p + j + r_1 - n << std::endl;

			// In this case, s_1 > 0, so s_3 = 0, and wrap-around cannot occur.

			y_1 = setButterfly (std::vector<bool>(x.begin (), x.begin () + (n - n_p)), j);
			y_2 = setButterfly (std::vector<bool>(x.begin () + (n - n_p), x.end ()), 2*n_p + j + r_1 - n);

			for (std::vector<bool>::iterator iter = (y_3.begin () + (j + r_1));
			     iter != (y_3.begin () + (n - n_p));
			     iter++)
				*iter = true;
		} else {
			commentator.indent (report);
			report << "  A: j >= (n - n_p).  j_1 = j + r - n - r_1 = "
			       << j + r - n - r_1 << ", j_2 = j - (n - n_p) = "
			       << j - (n - n_p) << std::endl;

			// In this case, s_1 = 0, so s_3 >= 0, and wrap-around may occur.
			// This case cannot occur for n != 2*n_p because j != 0.

			y_1 = setButterfly (std::vector<bool>(x.begin (), x.begin () + (n - n_p)), j + r - n - r_1);
			y_2 = setButterfly (std::vector<bool>(x.begin () + (n - n_p), x.end ()), j - (n - n_p));
 
			for (std::vector<bool>::iterator iter = y_3.begin ();
			     iter != (y_3.begin () + (j + r - n - r_1));
			     iter++)
				*iter = true;
		}
	}
	else if ((s_1 + s_3) < r_1) {
		commentator.indent (report);
		report << "Case III: s_1 + s_3 < r_1 and s_2 > r - r_1."
		       << "Elements are moved from the left sub-vector to the right." << std::endl;

		// This case also means that s_1 + s_3 < n - n_p, or the contiguous 
		// block cannot encompass the entire first sub-vector.  For this 
		// reason, this case is not considered when n != 2*n_p (when j = 0).

		if (j < (n - n_p)) {
			commentator.indent (report);
			report << "  A: j < (n - n_p).  j_1 = j = " << j
			       << ", j_2 = j + r_1 - n + n_p = " << j + r_1 - n + n_p << std::endl;
			// In this case, s_1 > 0, so s_3 = 0, and wrap-around cannot occur.

			y_1 = setButterfly (std::vector<bool>(x.begin (), x.begin () + (n - n_p)), j);
			y_2 = setButterfly (std::vector<bool>(x.begin () + (n - n_p), x.end ()), j + r_1 - n + n_p);

			for (std::vector<bool>::iterator iter = (y_3.begin () + s_3);
			     iter != (y_3.begin () + (j + r_1 - n + n_p));
			     iter++)
				*iter = true;
		} else {
			commentator.indent (report);
			report << "  A: j >= (n - n_p).  j_1 = j + r - n_p - r_1 = "
			       << j + r - n_p - r_1 << ", j_2 = j - (n - n_p) = "
			       << j - (n - n_p) << std::endl;

			// In this case, s_1 = 0, so s_3 >= 0, and wrap-around may occur.
			// This case cannot occur for n != 2*n_p because j != 0.

			y_1 = setButterfly (std::vector<bool>(x.begin (), x.begin () + (n - n_p)), j + r - n_p - r_1);
			y_2 = setButterfly (std::vector<bool>(x.begin () + (n - n_p), x.end ()), j - (n - n_p));
 
			for (std::vector<bool>::iterator iter (y_3.begin () + (j + r - n_p - r_1));
			     iter != (y_3.begin () + (n - n_p));
			     iter++)
				*iter = true;
		}
	}

	// Create output vector.
	std::vector<bool> y (y_1);
	y.insert (y.end (), y_2.begin (), y_2.end ());
	y.insert (y.end (), y_3.begin (), y_3.end ());

	commentator.indent (report);
	report << "The output vector for n = " << n << " has " << y.size ()
	       << " entries."
	       << "  " << y_1.size () << " from the first sub-vector"
	       << "  " << y_2.size () << " from the second sub-vector"
	       << "  " << y_3.size () << " from recombining the two"
	       << "And the output vector y is:"
	       << "-------------------------- " << std::endl;

	for (size_t i = 0; i < y.size (); i++) {
		commentator.indent (report);
		report << "  " << i << ": " << y[i] << std::endl;
	}

	commentator.indent (report);
	report << "-------------------------- " << std::endl;

	commentator.stop ("done");
 
	return y;

} // std::vector<bool> setButterfly (const std::vector<bool>& x, size_t j)

//@}
} // namespace LinBox

#endif // __LINBOX_butterfly_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
