/* File: src/examples/util/vector_utility.h
 * Author: William J Turner for the Linbox group
 */

#ifndef _VECTOR_UTILITY_
#define _VECTOR_UTILITY_

#include <iostream>
#include <string>
#include <fstream>
#include <vector>

#include "../test_base.h"

#include "LinBox/vector_traits.h"

/** Vector utility object.
 * Contains methods for input and output of LinBox vectors to be used in 
 * examples.
 * Templatized by field type, vector type, and vector trait type.
 * The vector trait object has a default value taken from the vector_traits
 * object.  This class is then specialized for dense and sparse \Ref{LinBox} 
 * vectors.
 */
template < class Field, 
           class Vector, 
	   class Trait = LinBox::vector_traits<Vector>::vector_category >
class vector_utility
{
public:
  /** Constructor from a field object.
   * @param F \Ref{LinBox} field object for field in which arithmetic is to
   *          be done
   * @param  prompt  boolean for whether to prompt user for input
   *                 (default = true)
   * @param  in  istream from which to read input (default = cin)
   * @param  out output stream to which to print output (default = cout)
   */
  vector_utility(const Field& F);

  /* Creates LinBox vector.
   * Reads LinBox vector from input stream, first the index and then the 
   * element, ending when gets index of -1.
   * If TRACE is defined, output includes all inputs.
   * @return pointer to new STL vector in dynamic memory
   * @param  m row dimensions of vector (defualt = 0)
   * @param  prompt  boolean for whether to prompt user for input
   *                 (default = true)
   * @param  in  istream from which to read input (default = cin)
   * @param  out output stream to which to print output (default = cout)
   */
  Vector* new_vector (size_t m=0, 
		      bool prompt = true, 
		      istream& in = cin,
		      ostream& out = cout) const;

  /** Write vector to output stream.
   * Writes "Element i: element" where i is the index and element is the
   * field element stored in the given index.  Writes one element per line.
   * @param out output stream reference
   * @param b constant reference to vector object.
   */
  ostream& write(ostream& out, const Vector& b) const;
  
  /** Read vector from input stream.
   * Reads first index and then field element from input stream.
   * Stops when receives index value of -1.
   * @param in reference to input stream
   * @param b reference to vector
   */
  istream& read(istream& is, Vector& b) const;
  
}; // template <class Field, class Vector, class Trait> class vector_utility

// Specialization for dense LinBox vectors
template <class Field, class Vector>
class vector_utility<Field, Vector, LinBox::vector_categories::dense_vector_tag>
{
public:
  vector_utility(const Field& F) : _F(F) {}

  Vector* new_vector (size_t m=0, 
		      bool prompt = true, 
		      istream& in = cin,
		      ostream& out = cout) const;
  ostream& write(ostream& out, const Vector& b) const;
  istream& read(istream& is, Vector& b) const;
  
private:
  // Field object
  Field _F;

}; // template <> class vector_utility<dense_vector_tag>

// Specialization for sparse sequence LinBox vectors
template <class Field, class Vector>
class vector_utility<Field, 
                     Vector, 
		     LinBox::vector_categories::sparse_sequence_vector_tag>
{
public:
  vector_utility(const Field& F) : _F(F) {}

  Vector* new_vector (size_t m=0, 
		      bool prompt = true, 
		      istream& in = cin,
		      ostream& out = cout) const;
  ostream& write(ostream& out, const Vector& b) const;
  istream& read(istream& is, Vector& b) const;
   
private:
  // Field object
  Field _F;

  // used in lower_bound as function object
  struct comp_w_ind; 
/*
  {
    bool operator() (const pair< size_t, typename Field::element >& entry, 
		     size_t col_in)
    { return entry.first < col_in; }
  }; // struct comp_w_ind
*/

}; // template <> class vector_utility<sparse_sequence_vector_tag>

// Specialization for sparse associative LinBox vectors
template <class Field, class Vector>
class vector_utility<Field, 
                     Vector, 
		     LinBox::vector_categories::sparse_associative_vector_tag>
{
public:
  vector_utility(const Field& F) : _F(F) {}

  Vector* new_vector (size_t m=0, 
		      bool prompt = true, 
		      istream& in = cin,
		      ostream& out = cout) const;
  ostream& write(ostream& out, const Vector& b) const;
  istream& read(istream& is, Vector& b) const;
   
private:
  // Field object
  Field _F;

}; // template <> class vector_utility<sparse_associative_vector_tag>


// Implementation of dense vector methods

template <class Field, class Vector>
Vector*
vector_utility<Field, Vector, LinBox::vector_categories::dense_vector_tag>
::new_vector(size_t m=0, 
	     bool prompt = true, 
	     istream& in = cin, 
	     ostream& out = cout) const
{
  if (m==0)
  {
    if (prompt) cout << "What is the length of the vector? ";
    in >> m ;

#ifdef TRACE
    out << endl << "The length of the vector is " << m << endl;
#endif
  } // if (m==0)
  
  typename Field::element zero;
  _F.init(zero, 0);
  Vector* b_ptr = new Vector(m, zero);
  
  if (prompt)
    cout << endl << "Input vector by entering index and value." << endl
                 << "Remember vector is indexed starting at 0." << endl
                 << "End with a index of -1." << endl;
 
  read(in, *b_ptr);

#ifdef TRACE
  out << endl 
      << "The right hand side contains the following elements: " << endl;
  write(out, *b_ptr);
  out << endl;
#endif

  return b_ptr;
} // template <> Vector* vector_utility<dense_vector_tag>::new_vector(...) const

template <class Field, class Vector>
ostream&
vector_utility<Field, Vector, LinBox::vector_categories::dense_vector_tag>
::write(ostream& out, const Vector& b) const
{
  typename Vector::const_iterator iter;
  size_t i = 0;
  
  for (iter = b.begin(); iter != b.end(); iter++, i++)
  {
    out << "Element " << i << ": ";
    _F.write(out, *iter); // can't do this!!!
    out << endl;
  }

  return out;
} // template <> ostream& vector_utility<dense_vector_tag>::write(...) const

template <class Field, class Vector>
istream&
vector_utility<Field, Vector, LinBox::vector_categories::dense_vector_tag>
::read(istream& in, Vector& b) const
{
  size_t i;
  typename Field::element elem;
  _F.init(elem, 0);
  
  while (in >> i)
  {
    if(i == size_t(-1)) break; // return also if row index is not positive
    _F.read(in, elem);
    b[i] = elem;
  } // while (is >> i)

  return in;
  
} // template <> istream& vector_utility<dense_vector_tag>::read(...) const

// Implementation of sparse sequence vector methods

template <class Field, class Vector>
Vector*
vector_utility<Field, Vector, LinBox::vector_categories::sparse_sequence_vector_tag>
::new_vector(size_t m=0, 
	     bool prompt = true, 
	     istream& in = cin, 
	     ostream& out = cout) const
{
  if (m==0)
  {
    if (prompt) cout << "What is the length of the vector? ";
    in >> m ;

#ifdef TRACE
    out << endl << "The length of the vector is " << m << endl;
#endif
  } // if (m==0)
  
  typename Field::element zero;
  _F.init(zero, 0);
  Vector* b_ptr = new Vector();
  
  if (prompt)
    cout << endl << "Input vector by entering index and value." << endl
                 << "Remember vector is indexed starting at 0." << endl
                 << "End with a index of -1." << endl;
 
  read(in, *b_ptr);

  // Check to make sure no extra elements were added.

#ifdef TRACE
  out << endl 
      << "The right hand side contains the following elements: " << endl;
  write(out, *b_ptr);
  out << endl;
#endif

  return b_ptr;
} // template <> Vector* vector_utility<sparse_sequence_vector_tag>::new_vector(..) const

template <class Field, class Vector>
ostream&
vector_utility<Field, Vector, LinBox::vector_categories::sparse_sequence_vector_tag>
::write(ostream& out, const Vector& b) const
{
  typename Vector::const_iterator iter;
  
  for (iter = b.begin(); iter != b.end(); iter++)
  {
//    out << "Element " << iter->first << ": "; // problems with deque
//    _F.write(out, iter->second);              //    "      "     "
    out << "Element " << (*iter).first << ": ";
    _F.write(out, (*iter).second);
    out << endl;
  } // for (iter = b.begin(); iter != b.end(); iter++)
  
  return out;
} // template <> Vector* vector_utility<sparse_sequence_vector_tag>::write(...) const

template <class Field, class Vector>
istream&
vector_utility<Field, Vector, LinBox::vector_categories::sparse_sequence_vector_tag>
::read(istream& in, Vector& b) const
{
  size_t i;
  typename Field::element elem;
  _F.init(elem, 0);

  typename Vector::iterator iter;
  bool found;
  
  while (in >> i)
  {
    if(i == size_t(-1)) break; // return also if row index is not positive
    _F.read(in, elem);
    
    // find appropriate location of element
    if( b.begin() == b.end() )
      iter = b.end();
    else
      iter = lower_bound( b.begin(), b.end(), i, comp_w_ind() );

    // Check to see if element already exists.
    if ( b.end() == iter )
      found = false;
    else if ( iter->first != i )
      found = false;
    else 
      found = true;

    // If element is already in row, replace old value with new.
    // Otherwise, insert the element in the row.
    if (found) 
    {
      if (_F.isZero(elem))
  	b.erase(iter);
      else
  	iter->second = elem;
    } // if (found)
    else if (!_F.isZero(elem))
      b.insert(iter, make_pair(i,elem));

  } // while (is >> i)

  return in;
} // template <> Vector* vector_utility<sparse_sequence_vector_tag>::read(...) const

template <class Field, class Vector>
struct
vector_utility<Field, Vector, LinBox::vector_categories::sparse_sequence_vector_tag>
::comp_w_ind
{
  bool operator() (const pair< size_t, typename Field::element >& entry, 
		   size_t col_in)
  { return entry.first < col_in; }
}; // struct comp_w_ind

// Implementation of sparse associative vector methods

template <class Field, class Vector>
Vector*
vector_utility<Field, Vector, LinBox::vector_categories::sparse_associative_vector_tag>
::new_vector(size_t m=0, 
	     bool prompt = true, 
	     istream& in = cin, 
	     ostream& out = cout) const
{
  if (m==0)
  {
    if (prompt) cout << "What is the length of the vector? ";
    in >> m ;

#ifdef TRACE
    out << endl << "The length of the vector is " << m << endl;
#endif
  } // if (m==0)
  
  Vector* b_ptr = new Vector();
  
  if (prompt)
    cout << endl << "Input vector by entering index and value." << endl
                 << "Remember vector is indexed starting at 0." << endl
                 << "End with a index of -1." << endl;
 
  read(in, *b_ptr);

  // Check to make sure no extra elements were added.

#ifdef TRACE
  out << endl 
      << "The right hand side contains the following elements: " << endl;
  write(out, *b_ptr);
  out << endl;
#endif

  return b_ptr;
} // template <> Vector* vector_utility<sparse_associative_vector_tag>::new_vector(..) const

template <class Field, class Vector>
ostream&
vector_utility<Field, Vector, LinBox::vector_categories::sparse_associative_vector_tag>
::write(ostream& out, const Vector& b) const
{
  typename Vector::const_iterator iter;
  
  for (iter = b.begin(); iter != b.end(); iter++)
  {
//    out << "Element " << iter->first << ": "; // problems with deque
//    _F.write(out, iter->second);              //    "      "     "
    out << "Element " << (*iter).first << ": ";
    _F.write(out, (*iter).second);
    out << endl;
  } // for (iter = b.begin(); iter != b.end(); iter++)
  
  return out;
} // template <> Vector* vector_utility<sparse_associative_vector_tag>::write(...) const

template <class Field, class Vector>
istream&
vector_utility<Field, Vector, LinBox::vector_categories::sparse_associative_vector_tag>
::read(istream& in, Vector& b) const
{
  size_t i;
  typename Field::element elem;
  _F.init(elem, 0);

  typename Vector::iterator iter;
  
  while (in >> i)
  {
    if(i == size_t(-1)) break; // return also if row index is not positive
    _F.read(in, elem);
   
    // Find element in associative container.  
    // If exists, replace value if not zero, or remove if value is zero.
    // If not found, insert non-zero element
    if ( (iter = b.find(i)) != b.end() )
    {
      if (_F.isZero(elem))
        b.erase(iter);
      else
        iter->second = elem;
    } // ( (iter = b.find(i)) != b.end() )
    else
    {
      if (!_F.isZero(elem))
        b.insert(make_pair(i, elem));
    }

  } // while (is >> i)

  return in;
} // template <> Vector* vector_utility<sparse_associative_vector_tag>::read(...) const

/** This class contains code for creating a random vector of integers.
 * It is derived from test_base to include functions and files for
 * creating input and output streams.
 * It prompts for the size of the matrix.  It then creates output that can be
 * read into the \Ref{vector_utility} programs.
 */
class random_vector : public test_base
{
public:
  /** Constructor from int and array of C-style strings.
   * Creates input and output streams.
   * @param n size of vector (default = 0)
   * @param  in  istream from which input is read (default = cin)
   * @param  out ostream to which output is written (default = cout)
   * @param  log ostream to which messages are logged (default = clog)
   */
  random_vector(size_t n = 0, 
		istream& in = cin, 
		ostream& out = cout, 
		ostream& log = clog);

  /// Destructor
  ~random_vector(void) {}

  /** Create random dense matrix.
   * Prompts for the size of the vector and then creates output that
   * can be read into the \Ref{vector_utility} programs.
   * If the dimension is zero, it will prompt for it as well.
   * @return boolean true if successfull, false otherwise
   */
  bool create(void) const;

private:

  // size of vector
  size_t _n;
  
}; // class random_vector : public test_base

// Implementation of methods

random_vector::random_vector(size_t n = 0, 
			     istream& in = cin, 
			     ostream& out = cout, 
			     ostream& log = clog)
: test_base(in, out, log), _n(n)
{
  while (_n <= 0)
  {
    if (prompt) cout << "Enter the dimension of the vector: ";
    *in_ptr >> _n;
  } // while (_n <= 0)

} // random_vector::random_vector(size_t n = 0, ...)


bool random_vector::create(void) const
{
#ifdef TRACE
  *log_ptr << "Creating a vector of size " << _n << endl;
#endif // TRACE
 
  if (prompt) 
    cout << "Enter an integer seed for the random number generator" << endl
         << "or enter -1 to obtain seed from process time." << endl;

  size_t seed;
  *in_ptr >> seed;

  if (seed == size_t(-1)) seed = time(NULL);

#ifdef TRACE
  *log_ptr << "Using " << seed << " to seed random number generator." << endl;
#endif // TRACE
 
  srand(seed);

  if (prompt) cout << "Enter a maximum value for the random integers: ";
  size_t max;
  *in_ptr >> max;

#ifdef TRACE
  *log_ptr << "Using " << max << " as maximum random integer." << endl;
#endif // TRACE

  *out_ptr << _n << endl;
 
  for (size_t i = 0; i < _n; i++)
    *out_ptr << i << " " 
      << static_cast<long>((double(rand())/RAND_MAX)*max) << endl;

  *out_ptr << -1 << endl;

  return true;
  
} // bool random_vector::create(void) const

#endif // _VECTOR_UTILITY_

