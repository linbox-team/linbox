/* File: src/library/objects/blackbox/butterfly.h
 * Author: William J Turner for the LinBox group
 */

#ifndef _BUTTERFLY_
#define _BUTTERFLY_

#ifdef TRACE
#include <iostream>
#endif // TRACE

#include <vector>

#include "LinBox/blackbox_archetype.h"

// Namespace in which all LinBox library code resides
namespace LinBox
{

  /** Butterfly Switching Network BlackBox Matrix Object
   * Implements butterfly switching network on a LinBox vector
   * as a black box matrix through the use of a switch object.
   *
   * This is a blackbox matrix object, and it implements all
   * purely virtual methods of the abstract base class 
   * Blackbox_archetype.
   *
   * This matrix requires a dense vector to be used.  Sparse vectors must
   * somehow be converted to dense vectors before this matrix may
   * be applied to them.
   *
   * @param Vector LinBox dense vector type
   * @param Switch switch object type
   */
  template <class Vector, class Switch>
  class butterfly : public Blackbox_archetype<Vector>
  {
  public:

    /** Constructor from an integer and a switch object.
     * The switch object is an object that is applied
     * to two references to elements to switch them.  It must have both
     * an apply and an applyTranspose method.
     * It must contain all information needed by the switch other 
     * than the elements themsleves.  This includes any random
     * numbers or sequences of values.  It must also be able to 
     * be applied as many times as needed.  In particular, it must be able
     * to create new random elements or repeat a stored sequence
     * of values.
     * This is not required by the abstract base class.
     * @param n integer size of vectors to be applied to
     * @param S switch predicate object object
     */
    butterfly(size_t n, const Switch& S);

    /** Destructor. */
    ~butterfly() {}

    /** Virtual constructor.
     * Required because constructors cannot be virtual.
     * Make a copy of the Blackbox_archetype object.
     * Required by abstract base class.
     * @return pointer to new blackbox object
     */
   Blackbox_archetype<Vector>* clone() const 
   { return new butterfly(*this); }

    /** Application of BlackBox matrix.
     * y = A*x.
     * Requires one vector conforming to the \Ref{LinBox}
     * vector {@link Archetypes archetype}.
     * Required by abstract base class.
     * For this matrix, this involves applying each switch in order to the 
     * input vector.
     * @return reference to vector y containing output (after switching).
     * @param  x constant reference to vector to contain input 
     * 			(before switching)
     */
    Vector& apply(const Vector& x) const;

    /** Application of BlackBox matrix transpose.
     * y = transpose(A)*x.
     * Requires one vector conforming to the \Ref{LinBox}
     * vector {@link Archetypes archetype}.
     * Required by abstract base class.
     * For this matrix, this involves applying the transpose of each switch 
     * to the input vector in the reverse order of the apply function.
     * @return reference to vector y containing output (after switching).
     * @param  x constant reference to vector to contain input 
     * 			(before switching)
     */
    Vector& applyTranspose(const Vector& x) const;

    /** Retreive row dimensions of BlackBox matrix.
     * This may be needed for applying preconditioners.
     * Required by abstract base class.
     * @return integer number of rows of black box matrix.
     */
    size_t rowdim(void) const { return _n; }
    
    /** Retreive column dimensions of BlackBox matrix.
     * Required by abstract base class.
     * @return integer number of columns of black box matrix.
     */
    size_t coldim(void) const { return _n; }

  private:

    // Number of rows and columns of square matrix.
    size_t _n;

    // Switch object to use
    Switch _switch;

    // Vectors of sizes of sub-groups and number of levels in each
    // These may not need to be stored in general.
    // They may only be used in the constructor
    std::vector<size_t> _n_vec, _l_vec;
   
    // Vector of index pairs.  These are the indices to be switched with
    // a given switch.
    std::vector< std::pair< size_t, size_t > > _indices;
    
  }; // template <class Field, class Vector> class butterfly

  // Implementation of methods

  template <class Vector, class Switch>
  inline butterfly<Vector, Switch>::butterfly(size_t n, 
					      const Switch& S)
  : _n(n), _switch(S)
  {
    // Ensure n is non-negative
    if (n < 0) n = 0;

#ifdef TRACE
    clog
      << "Called butterfly constructor from size_t " << n 
      << " and switch object" << endl
      << "Constructed empty vectors _l_vec and _n_vec:" << endl
      << "    _l_vec    _n_vec" << endl
      << "    ----------------" << endl;

    for (size_t i = 0; i != _n_vec.size(); i++)
      clog << "    " << _l_vec[i] << "        " << _n_vec[i] << endl;
    
    clog << "    ----------------" << endl;
#endif // TRACE

    for (size_t value(_n), l_p(0), n_p(1); 
	 n_p != 0; 
	 value >>= 1, l_p++, n_p <<= 1)
    {
#ifdef TRACE_LOOP
      clog
	<< "  looping at value = " << value 
	<< ", l_p = " << l_p 
	<< ", n_p = " << n_p << endl;
#endif // TRACE_LOOP
      
      if (value & 1)
      {
	_l_vec.push_back(l_p);
	_n_vec.push_back(n_p);      
#ifdef TRACE
      clog 
	<< "    inserted value = " << value 
	<< ", l_p = " << l_p 
	<< ", n_p = " << n_p << endl;
#endif // TRACE
      
      } // if (value & 1)

    } // for (size_t value(_n), l_p(0), n_p(1); n_p != 0; ...)
   
#ifdef TRACE
    clog
      << "Constructed vectors _l_vec and _n_vec:" << endl
      << "    _l_vec    _n_vec" << endl
      << "    ----------------" << endl;
    
    for (size_t i = 0; i != _n_vec.size(); i++)
      clog << "    " << _l_vec[i] << "        " << _n_vec[i] << endl;
    
    clog
      << "    ----------------" << endl
      << "Constructed empty vector of indices:" << endl
      << "    i        index 1        index 2" << endl
      << "    -------------------------------" << endl;
    
    for (size_t i = 0; i != _indices.size(); i++)

      clog 
	<< "    " << i << "        " << _indices[i].first
	<< "        " << _indices[i].second << endl;
    
    clog << "    -------------------------------" << endl;
#endif // TRACE

    // Create vector of indices to switch
    size_t n_p, l_p;   	// size of group and number of levels in group
    size_t level(0), difference(1);	// track levels done for powers of 2

    // Vector containing indices for last level of last power of 2.
    std::vector< std::pair< size_t, size_t > > p_ind;
    
    // Vector and iterator used for computing p_ind.
    std::vector< std::pair< size_t, size_t > > temp_ind;
    std::vector< std::pair< size_t, size_t > >::iterator iter;
    
    // Loop over sub-groups of powers of two
    for (size_t p(0), start_index(0); 
	 p < _n_vec.size(); 
	 p++, start_index += n_p)
    {
      // update size
      n_p = _n_vec[p];
      l_p = _l_vec[p];

#ifdef TRACE_LOOP
      clog 
	<< "Sub-group p = " << p 
	<< ", size n_p = " << n_p 
	<< ", levels l_p = " << l_p
	<< ", starting at index " << start_index 
	<< ", switches " << n_p*l_p/2 << endl;
#endif // TRACE_LOOP

      // loop over levels of sub-group network
      for ( ; level < l_p; level++, difference <<= 1)
      {
#ifdef TRACE_LOOP
	clog
	  << "  level " << level
	  << ", number of nodes " << 2*difference
	  << ", number of switches " << (level + 1) * difference << endl;
#endif // TRACE_LOOP

	// Create 
	temp_ind = p_ind;

	// the second sub group is a shift of the first
	for (iter = temp_ind.begin(); iter != temp_ind.end(); iter++)
	{
	  iter->first += difference;
	  iter->second += difference;
	} // for (size_t i = 0; i < temp_ind.size(); i++, iter++)

	// add the second group to the first
	p_ind.insert(p_ind.end(), temp_ind.begin(), temp_ind.end());

	// add switches to mix the two sub groups
	temp_ind = std::vector< pair<size_t, size_t> >(difference, 
						       make_pair(0, 0));

	size_t i = 0;
	for (iter = temp_ind.begin(); iter != temp_ind.end(); i++, iter++)
	{
	  iter->first += i;
	  iter->second += i + difference;
	} // for (iter = temp_ind.begin(); iter != temp_ind.end(); i++, iter++)

	// add the combining group to the first and second
	p_ind.insert(p_ind.end(), temp_ind.begin(), temp_ind.end());

#ifdef TRACE_LOOP
	clog
	  << "      i        x        y" << endl
	  << "      -------------------" << endl;
	
	for (size_t j = 0; j < p_ind.size(); j++)
	  clog
	    << "      " << j << "       " << p_ind[j].first << "       " 
	    << p_ind[j].second << endl;
	
	clog << "      -------------------" << endl;
#endif // TRACE_LOOP

      } // for (size_t level(0), difference(1); level < l_p; ...)
	
      // Add this level to total list of indices and correct starting point
      temp_ind = p_ind;
      
      for (iter = temp_ind.begin(); iter != temp_ind.end(); iter++)
      {
	iter->first += start_index;
	iter->second += start_index;
      } // for (iter = temp_ind.begin(); iter != temp_ind.end(); iter++)
      
      _indices.insert(_indices.end(), temp_ind.begin(), temp_ind.end());
      
#ifdef TRACE_LOOP
      clog
	<< "combining sub-group " << p 
	<< ", starting at index " << 0
	<< ", differences in indices " << start_index
	<< ", switches " << start_index << endl;
#endif // TRACE_LOOP

      // Combine everything so far
      temp_ind = std::vector< pair<size_t, size_t> >(start_index,
						     make_pair(0, 0));

      iter = temp_ind.begin();
      for (size_t index = 0; index < start_index; index++, iter++)
      {
	iter->first = index;
	iter->second += index + n_p;
      } // for (size_t index = 0; index < start_index; index++, iter++)
      
#ifdef TRACE_LOOP
	clog
	  << "      i        x        y" << endl
	  << "      -------------------" << endl;
	
	for (size_t j = 0; j < temp_ind.size(); j++)
	  clog
	    << "      " << j << "       " << temp_ind[j].first << "       " 
	    << temp_ind[j].second << endl;
	
	clog << "      -------------------" << endl;
#endif // TRACE_LOOP

      _indices.insert(_indices.end(), temp_ind.begin(), temp_ind.end());
    
  } // for (size_t p(0), start_index(0); p < _n_vec.size(); ...)

#ifdef TRACE
    clog
      << "Constructed vector of indices:" << endl
      << "    i        index 1        index 2" << endl
      << "    -------------------------------" << endl;
    
    for (size_t i = 0; i != _indices.size(); i++)
      clog
	<< "    " << i << "        " << _indices[i].first
	<< "        " << _indices[i].second << endl;
    
    clog << "    -------------------------------" << endl;
    
    // Calculate total number of switches required
    size_t s(0);
    
    for (size_t i = 0; i < _n_vec.size(); i++)
      s += _n_vec[i]*_l_vec[i]/2;
    
    if (_n_vec.size() > 0)
      for (size_t i = 0; i < _n_vec.size() - 1; i++)
	for (size_t j = 0; j <= i; j++)
	  s += _n_vec[j];

    clog
      << "The total number of switches and index pairs needed is " << s 
      << endl << "and we have " << _indices.size() 
      << " pairs of indices." << endl;
#endif // TRACE

  } // butterfly<>::butterfly(size_t, const Switch&)
  
  template <class Vector, class Switch>
  inline Vector& butterfly<Vector, Switch>::apply(const Vector& x) const
  {
#ifdef TRACE
    clog << "Called butterfly.apply(x)" << endl;
#endif // TRACE
    
    Vector* y_ptr(new Vector(x));
    std::vector< pair<size_t, size_t> >::const_iterator iter;
    Switch temp_switch(_switch);
    
    for (iter = _indices.begin(); iter != _indices.end(); iter++)
    {
#ifdef TRACE
      clog
	<< "  Switching x[" << iter->first << "] and x[" << iter->second 
	<< "]: ";
#endif // TRACE

      temp_switch.apply((*y_ptr)[iter->first], (*y_ptr)[iter->second]);

#ifdef TRACE
      clog << endl;
#endif // TRACE

    } // for (iter = _indices.begin(); iter != _indices.end(); iter ++)

    return *y_ptr;
  } // Vector& butterfly<Vector, Switch>::apply(const Vector& x) const

  template <class Vector, class Switch>
  inline Vector& 
  butterfly<Vector, Switch>::applyTranspose(const Vector& x) const
  {
#ifdef TRACE
    clog << "Called butterfly.applyTranspose(x)" << endl;
#endif // TRACE
    
    Vector* y_ptr(new Vector(x));
    std::vector< pair<size_t, size_t> >::const_reverse_iterator iter;
    Switch temp_switch(_switch);
    
    for (iter = _indices.rbegin(); iter != _indices.rend(); iter++)
    {
#ifdef TRACE
      clog
	<< "  Switching x[" << iter->first << "] and x[" << iter->second 
	<< "]: ";
#endif // TRACE

      temp_switch.applyTranspose((*y_ptr)[iter->first], (*y_ptr)[iter->second]);

#ifdef TRACE
      clog << endl;
#endif // TRACE

    } // for (iter = _indices.begin(); iter != _indices.end(); iter ++)

    return *y_ptr;
  } // Vector& butterfly<Vector, Switch>::applyTranspose(const Vector& x) const

  /** Set switches function.
   * This function takes an STL vector x of booleans, and returns
   * a vector y of booleans such that setting the switches marked
   * by true flags in y to be on (or to swap elements) the true
   * elements x will be switched to a given contiguous block
   * through the use of a butterfly switching network.
   * The integer parameter j marks where this block is to begin.
   * If x has r true elements, the butterfly switching network will place
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
  std::vector<bool> set_butterfly(const std::vector<bool>& x, 
				  size_t j = 0,
				  ostream& log = clog)
  {
    size_t n = x.size();
 
#ifdef TRACE
    log
      << "Called set switches with vector of size " << n
      << " and offset " << j << endl;
#endif // TRACE

    // return empty vector if zero or one elements in x because
    // no switching will be done.
    if (x.size() <= 1)
    {
#ifdef TRACE
    log << "  No switches needed.  Returning with empty vector." << endl;
#endif // TRACE

      return std::vector<bool>();
    } // if (x.size() <= 1)

#ifdef TRACE
    log << "Counting the number of switches that exist." << endl;
#endif // TRACE
 
    // break inputs into groups of size powers of 2.
    // calculate size of groups, and powers of 2 that give sizes
    // store these values in vectors n and l, respectively
    vector<size_t> l_vec, n_vec;
    for (size_t value(n), l_p(0), n_p(1);
  	 n_p != 0;
  	 value >>= 1, l_p++, n_p <<= 1)
    {
#ifdef TRACE_LOOP
      log << "  looping at value = " << value
          << ", l_p = " << l_p
	  << ", n_p = " << n_p << endl;
#endif // TRACE_LOOP
 
      if (value & 1)
      {
  	l_vec.push_back(l_p);
  	n_vec.push_back(n_p);
#ifdef TRACE_LOOP
  	log << "    inserted value = " << value
	    << ", l_p = " << l_p
	    << ", n_p = " << n_p << endl;
#endif // TRACE_LOOP

      } // if (value & 1)
 
    } //     for (size_t value(n), l_p(0), n_p(1); n_p != 0; ...)
 
    // Calculate total number of switches required
    size_t s(0);
 
    for (size_t ii = 0; ii < n_vec.size(); ii++)
      s += n_vec[ii]*l_vec[ii]/2;
 
    for (size_t ii = 0; ii < n_vec.size() - 1; ii++)
      for (size_t jj = 0; jj <= ii; jj++)
  	s += n_vec[jj];
 
#ifdef TRACE
    log << "There are a total of " << s << " switches" << endl;
#endif // TRACE
 
    // Set largest power of 2 in decomposition of n = x.size()
    size_t n_p(*n_vec.rbegin());
 
#ifdef TRACE
      log << "Found largest power of 2 in decomposition of " << n
	  << " as n_p = " << n_p << endl;
#endif // TRACE

      if ( (n != n_p) && (j != 0) )
      {
  	log << "Non-zero offset " << j
	    << " used with non-power size." << endl
	    << "Offset reset to zero." << endl;

  	j = 0;
      } // if ( (n != n_p) && (j != 0) )
      else // if !( (n != n_p) && (j != 0) )
      {
  	j %= n;
  	if (j < 0) j += n;
      } // else // if !( (n != n_p) && (j != 0) )

      if (n == n_p)
      {
  	n_p /= 2;	  // >> is not portable!
#ifdef TRACE
  	log << "n = " << n << " is a power of two.  "
	    << "Resetting n_p to be half of n: n_p = " << n_p << endl;
#endif // TRACE
      } // if (n = n_p)

      // count true elements not in largest power of 2 block
      size_t r_1(0);
 
    for (std::vector<bool>::const_iterator iter = x.begin();
  	 iter != x.begin() + (n - n_p);
  	 iter++)
      if (*iter) r_1++;

    // count total number of true elements in x.
    size_t r(r_1);
 
    for (std::vector<bool>::const_iterator iter = x.begin() + (n - n_p);
  	 iter != x.end();
  	 iter++)
      if (*iter) r++;
 
#ifdef TRACE
      log << "The vector x will be broken into two sub-vectors," << endl
	  << "x_1 = x[0,...," << n - n_p - 1 << "] and x_2 = x["
	  << n - n_p << ",...," << n - 1 << "]." << endl
	  << "There are a total of " << r << " true elements in x, " << endl
	  << r_1 << " of which occured in the first sub-vector." << endl
	  << "The output vector will have " << s << " entries and will" << endl
	  << "switch the true elements of x into a contiguous block" << endl
	  << "[" << j << "," << j + r
	  << ") = [" << j << "," << j + r - 1<< "]."<< endl;
#endif // TRACE

      if (r == 0)
      {
#ifdef TRACE
  	log << "There are no true elements in x, so the recursion is" << endl
	    << "being broken and a vector of false flags returned." << endl;
#endif // TRACE

  	return std::vector<bool>(s, false);
      } // if (r == 0)
      else if (r == n)
      {
#ifdef TRACE
  	log << "There are no false elements in x, so the recursion is" << endl
	    << "being broken and a vector of false flags returned." << endl;
#endif // TRACE

  	return std::vector<bool>(s, false);
      } // else if (r == n)

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

      if (j < n - n_p)
      {
  	if (j + r < n - n_p)
  	  s_1 = r;
  	else
  	  s_1 = n - n_p - j;
      }
      else
  	s_1 = 0;

#ifdef TRACE
      size_t s_2;

      if (j + r < n - n_p)
  	s_2 = 0;
      else
      {
  	if (j + r < n)
  	  s_2 = j + r;
  	else
  	  s_2 = n;

  	if (j < n - n_p)
  	  s_2 -= (n - n_p);
  	else
  	  s_2 -= j;
      }
#endif // TRACE
 
      size_t s_3 = ((j + r) > n) ? j + r - n : 0;

#ifdef TRACE
      log << "The number of elements in each of the three blocks of " << endl
	  << "true elements in the end result are" << endl
	  << "s_1 = " << s_1
	  << ", s_2 = " << s_2
	  << ", and s_3 = " << s_3 << "." << endl;
#endif // TRACE

      // Create empty vector for output. y_temp is used to retrieve output
      // from recursion before inserting into output.
      std::vector<bool> y_1, y_2, y_3 = std::vector<bool>(n - n_p, false);

      if ( (s_1 + s_3) == r_1 )
      {
#ifdef TRACE
  	log << "Case I: s_1 + s_3 == r_1 and s_2 == r - r_1." << endl
  	    << "No elements are moved between the two sub-vectors." << endl;
#endif // TRACE

  	if (j < (n - n_p))
  	{
#ifdef TRACE
  	  log << "  A: j < (n - n_p).  j_1 = j = " << j << ", j_2 = 0" << endl;
#endif // TRACE

  	  y_1 = set_butterfly(std::vector<bool>(x.begin(), 
						x.begin() + (n - n_p)), j);
	  
  	  y_2 = set_butterfly(std::vector<bool>(x.begin() + (n - n_p), 
						x.end()), 0);
 
  	} // if (j < (n - n_p))
  	else // if !(j < (n - n_p))
  	{
#ifdef TRACE
  	  log << "  A: j >= (n - n_p).  j_1 = 0, j_2 = j - (n - n_p) = "
	      << j - (n - n_p) << endl;
  	  // This case cannot occur for n != 2*n_p because j != 0
#endif // TRACE

  	  y_1 = set_butterfly(std::vector<bool>(x.begin(), 
						x.begin() + (n - n_p)), 0);

  	  y_2 = set_butterfly(std::vector<bool>(x.begin() + (n - n_p), x.end()),
  			     j - (n - n_p));
 
  	} // else // if !(j < (n - n_p))

      } // if ( (s_1 + s_3) == r_1 )
      else if ( (s_1 + s_3) > r_1 )
      {
#ifdef TRACE
  	log << "Case II: s_1 + s_3 > r_1 and s_2 < r - r_1." << endl
	    << "Elements are moved from the right sub-vector to the left." 
	    << endl;
  	// This means that s_2 < n_p, so either s_1 = 0 or s_3 = 0 (or both).
#endif // TRACE
 
  	if (j < (n - n_p))
  	{
#ifdef TRACE
  	  log << "  A: j < (n - n_p).  j_1 = j, j_2 = 2*n_p + j + r_1 - n = "
	      << 2*n_p + j + r_1 - n << endl;
  	  // In this case, s_1 > 0, so s_3 = 0, and wrap-around cannot occur.
#endif // TRACE

  	  y_1 = set_butterfly(std::vector<bool>(x.begin(), 
						x.begin() + (n - n_p)), j);

  	  y_2 = set_butterfly(std::vector<bool>(x.begin() + (n - n_p), x.end()),
			      2*n_p + j + r_1 - n);

  	  for (std::vector<bool>::iterator iter = (y_3.begin() + (j + r_1));
  	       iter != (y_3.begin() + (n - n_p));
  	       iter++)
  	    *iter = true;
 
  	} // if (j < (n - n_p))
  	else // if !(j < (n - n_p))
  	{
#ifdef TRACE
  	  log << "  A: j >= (n - n_p).  j_1 = j + r - n - r_1 = "
	      << j + r - n - r_1 << ", j_2 = j - (n - n_p) = "
	      << j - (n - n_p) << endl;
  	  // In this case, s_1 = 0, so s_3 >= 0, and wrap-around may occur.
  	  // This case cannot occur for n != 2*n_p because j != 0.
#endif // TRACE

  	  y_1 = set_butterfly(std::vector<bool>(x.begin(), 
						x.begin() + (n - n_p)),
			      j + r - n - r_1);

  	  y_2 = set_butterfly(std::vector<bool>(x.begin() + (n - n_p), x.end()),
  			     j - (n - n_p));
 
  	  for (std::vector<bool>::iterator iter = y_3.begin();
  	       iter != (y_3.begin() + (j + r - n - r_1));
  	       iter++)
  	    *iter = true;
 
  	} // else // if !(j < (n - n_p))

      } // else if ( (s_1 + s_3) > r_1 )
      else if ( (s_1 + s_3) < r_1 )
      {
#ifdef TRACE
  	log << "Case III: s_1 + s_3 < r_1 and s_2 > r - r_1." << endl
	    << "Elements are moved from the left sub-vector to the right." 
	    << endl;
  	// This case also means that s_1 + s_3 < n - n_p, or the contiguous 
	// block cannot encompass the entire first sub-vector.  For this 
	// reason, this case is not considered when n != 2*n_p (when j = 0).
#endif // TRACE

  	if (j < (n - n_p))
  	{
#ifdef TRACE
  	  log << "  A: j < (n - n_p).  j_1 = j = " << j
	      << ", j_2 = j + r_1 - n + n_p = " << j + r_1 - n + n_p << endl;
  	  // In this case, s_1 > 0, so s_3 = 0, and wrap-around cannot occur.
#endif // TRACE

  	  y_1 = set_butterfly(std::vector<bool>(x.begin(), 
						x.begin() + (n - n_p)), j);
	  
	  y_2 = set_butterfly(std::vector<bool>(x.begin() + (n - n_p), x.end()),
			      j + r_1 - n + n_p);

  	  for (std::vector<bool>::iterator iter = (y_3.begin() + s_3);
  	       iter != (y_3.begin() + (j + r_1 - n + n_p));
  	       iter++)
  	    *iter = true;

  	} // if (j < (n - n_p))
  	else // if !(j < (n - n_p))
  	{
#ifdef TRACE
  	  log << "  A: j >= (n - n_p).  j_1 = j + r - n_p - r_1 = "
	      << j + r - n_p - r_1 << ", j_2 = j - (n - n_p) = "
	      << j - (n - n_p) << endl;
  	  // In this case, s_1 = 0, so s_3 >= 0, and wrap-around may occur.
  	  // This case cannot occur for n != 2*n_p because j != 0.
#endif // TRACE

  	  y_1 = set_butterfly(std::vector<bool>(x.begin(), 
						x.begin() + (n - n_p)),
			      j + r - n_p - r_1);
	  
  	  y_2 = set_butterfly(std::vector<bool>(x.begin() + (n - n_p), x.end()),
  			     j - (n - n_p));
 
  	  for (std::vector<bool>::iterator iter(y_3.begin() + (j + r - n_p - r_1));
  	       iter != (y_3.begin() + (n - n_p));
  	       iter++)
  	    *iter = true;
 
  	} // else // if !(j < (n - n_p))

      } // else if ( (s_1 + s_3) < r_1 )


      // Calculate offsets for recursion on each of the two sub-vectors
      size_t j_1, j_2;
      if (j < (n - n_p))
      {
  	j_1 = j;
  	j_2 = 0;
      }
      else
      {
  	j_1 = 0;
  	j_2 = j - (n - n_p);
      }

      // Create output vector.
      std::vector<bool> y(y_1);
      y.insert(y.end(), y_2.begin(), y_2.end());
      y.insert(y.end(), y_3.begin(), y_3.end());

#ifdef TRACE
      log << "The output vector for n = " << n << " has " << y.size()
	  << " entries." << endl
	  << "  " << y_1.size() << " from the first sub-vector" << endl
	  << "  " << y_2.size() << " from the second sub-vector" << endl
	  << "  " << y_3.size() << " from recombining the two" << endl
	  << "And the output vector y is:" << endl
	  << "-------------------------- " << endl;

      for (size_t i = 0; i < y.size(); i++)
  	log << "  " << i << ": " << y[i] << endl;
      
      log << "-------------------------- " << endl;
#endif // TRACE
 
      return y;
 
  } // std::vector<bool> set_butterfly(const std::vector<bool>& x, size_t j)

} // namespace LinBox

#endif // _BUTTERFLY_
