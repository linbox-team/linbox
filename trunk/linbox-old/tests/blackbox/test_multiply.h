/* File: src/examples/blackbox/test_multiply.h
 * Author: William J Turner for the LinBox group
 */

#ifndef _TEST_MULTIPLY_
#define _TEST_MULTIPLY_

#include "LinBox/sparsemat.h"
#include "LinBox/multiply.h"

#include "Examples/test_base.h"
#include "Examples/vector_utility.h"

/** Class to test multiply blackbox matrix.
 * Templatized by field and vector types, as well as vector trait.
 * The vector trait object has a default value taken from the vector_traits
 * object.  This class is then specialized for dense and sparse \Ref{LinBox} 
 * vectors.
 * @see multiply
 */
template < class Field, 
           class Vector, 
           class Trait = LinBox::vector_traits<Vector>::vector_category >
class test_multiply : public test_base, public vector_utility<Field, Vector>
{
public:

  // Field element type
  typedef typename Field::element Element;
  
  // Sparsemat matrix type to test multiply
  typedef LinBox::sparsemat<Field, std::map<size_t, Element>, Vector> Matrix;

  /** Constructor from Field object
   * @param  F field in which arithmetic is done
   * @param  mode blackbox matrix apply mode
   * @param  bbtime boolean flag on whether to use blackbox timer
   * @param  givtime boolean flag on whether to use Givaro timer
   * @param  in  istream from which input is read (default = cin)
   * @param  out ostream to which output is written (default = cout)
   * @param  log ostream to which messages are logged (default = clog)
   */
  test_multiply(const Field& F,
	       int mode,
	       bool bbtime,
	       bool givtime,
	       istream& in = cin,
	       ostream& out = cout,
	       ostream& log = clog);

  /** Run tests on LinBox multiply blackbox matrices.
   * Tests the template created in multiply.h using the field F.
   * Applies multiply matrix to elementary vectors to pull off
   * each of the columns of the product.
   * @return boolean true if ended successfully, false if not
   */
  bool test(void) const;

}; // class test_multiply : public test_base

// Specializationg for dense vectors
template <class Field, class Vector>
class test_multiply<Field, Vector, LinBox::vector_categories::dense_vector_tag>
: public test_base, public vector_utility<Field, Vector>
{
public:

  typedef typename Field::element Element;
  typedef LinBox::sparsemat<Field, std::map<size_t, Element>, Vector> Matrix;
  test_multiply(const Field& F,
	       int mode,
	       bool bbtime,
	       bool givtime,
	       istream& in = cin,
	       ostream& out = cout,
	       ostream& log = clog);
  bool test(void) const;

private:

  // Field in which arithmetic is done
  Field _F;

  // Blackbox matrix apply mode
  int _mode;

  // Boolean flags on whether to use timers.
  bool _bbtimer, _givtimer;

}; // template<> class test_multiply<dense_vector_tag>

// Specializationg for sparse sequence vectors
template <class Field, class Vector>
class test_multiply<Field, 
                    Vector, 
		    LinBox::vector_categories::sparse_sequence_vector_tag>
: public test_base, public vector_utility<Field, Vector>
{
public:

  typedef typename Field::element Element;
  typedef LinBox::sparsemat<Field, std::map<size_t, Element>, Vector> Matrix;
  test_multiply(const Field& F,
	       int mode,
	       bool bbtime,
	       bool givtime,
	       istream& in = cin,
	       ostream& out = cout,
	       ostream& log = clog);
  bool test(void) const;

private:

  // Field in which arithmetic is done
  Field _F;

  // Blackbox matrix apply mode
  int _mode;

  // Boolean flags on whether to use timers.
  bool _bbtimer, _givtimer;

}; // template<> class test_multiply<sparse_sequence_vector_tag>

// Specializationg for sparse associative vectors
template <class Field, class Vector>
class test_multiply<Field, 
                    Vector, 
		    LinBox::vector_categories::sparse_associative_vector_tag>
: public test_base, public vector_utility<Field, Vector>
{
public:

  typedef typename Field::element Element;
  typedef LinBox::sparsemat<Field, std::map<size_t, Element>, Vector> Matrix;
  test_multiply(const Field& F,
	       int mode,
	       bool bbtime,
	       bool givtime,
	       istream& in = cin,
	       ostream& out = cout,
	       ostream& log = clog);
  bool test(void) const;

private:

  // Field in which arithmetic is done
  Field _F;

  // Blackbox matrix apply mode
  int _mode;

  // Boolean flags on whether to use timers.
  bool _bbtimer, _givtimer;

}; // template<> class test_multiply<sparse_associative_vector_tag>


// Implementation of methods for dense vectors

template <class Field, class Vector>
test_multiply<Field, Vector, LinBox::vector_categories::dense_vector_tag>
::test_multiply(const Field& F,
		int mode,
		bool bbtime,
		bool givtime,
		istream& in = cin,
		ostream& out = cout,
		ostream& log = clog)
: test_base(in, out, log), vector_utility<Field, Vector>(F),
  _F(F), _mode(mode), _bbtimer(bbtime), _givtimer(givtime)
{
  // prompt for input of matrix
  if (prompt)
    cout << "Enter the name of the file from which" << endl
         << "to read the multiply matrix input." << endl
	 << "Enter 'cin' to read from standard input" << endl
	 << "or 'same' to continue reading from the" << endl
	 << "current input." << endl;

  string filename;
  *in_ptr >> filename;
  if (filename != string("same")) open_input(filename);

  if (in_ptr == &cin) 
    prompt = true;
  else 
    prompt = false;

} // test_multiply<dense_vector_tag>::test_multiply(...)

template <class Field, class Vector>
bool test_multiply<Field, Vector, LinBox::vector_categories::dense_vector_tag>
::test(void) const
{
  Matrix A(*LinBox::newSparsemat<Field, std::map<size_t, Element> >(_F,
								    0,
								    0,
								    prompt,
								    *in_ptr,
								    *out_ptr));

  Matrix B(*LinBox::newSparsemat<Field, std::map<size_t, Element> >(_F,
								    0,
								    0,
								    prompt,
								    *in_ptr,
								    *out_ptr));

  LinBox::multiply<Vector> C(&A, &B);

  Element zero, one;
  _F.init(zero, 0);
  _F.init(one, 1);

  Vector x, y;

  size_t cols = C.coldim();
 
  *out_ptr 
    << "The product of these two matrices contains the following elements: "
    <<  endl;
  
  for (size_t i = 0; i < cols; i++)
  {
    x = Vector(cols, zero);
    x[i] = one;
    y = C.apply(x);

    *out_ptr << "** Column " << i << ":" << endl;
    write(*out_ptr, y);
    *out_ptr << endl;

  } // for (size_t i = 0; i < cols; i++)
  
  return true;
} // bool test_multiply<dense_vector_tag>::test(void) const

// Implementation of methods for sparse sequence vectors

template <class Field, class Vector>
test_multiply<Field, 
              Vector, 
	      LinBox::vector_categories::sparse_sequence_vector_tag>
::test_multiply(const Field& F,
		int mode,
		bool bbtime,
		bool givtime,
		istream& in = cin,
		ostream& out = cout,
		ostream& log = clog)
: test_base(in, out, log), vector_utility<Field, Vector>(F),
  _F(F), _mode(mode), _bbtimer(bbtime), _givtimer(givtime)
{
  // prompt for input of matrix
  if (prompt)
    cout << "Enter the name of the file from which" << endl
         << "to read the multiply matrix input." << endl
	 << "Enter 'cin' to read from standard input" << endl
	 << "or 'same' to continue reading from the" << endl
	 << "current input." << endl;

  string filename;
  *in_ptr >> filename;
  if (filename != string("same")) open_input(filename);

  if (in_ptr == &cin) 
    prompt = true;
  else 
    prompt = false;

} // test_multiply<sparse_sequence_vector_tag>::test_multiply(...)

template <class Field, class Vector>
bool test_multiply<Field, 
                   Vector, 
		   LinBox::vector_categories::sparse_sequence_vector_tag>
::test(void) const
{
  Matrix A(*LinBox::newSparsemat<Field, std::map<size_t, Element> >(_F,
								    0,
								    0,
								    prompt,
								    *in_ptr,
								    *out_ptr));

  Matrix B(*LinBox::newSparsemat<Field, std::map<size_t, Element> >(_F,
								    0,
								    0,
								    prompt,
								    *in_ptr,
								    *out_ptr));

  LinBox::multiply<Vector> C(&A, &B);

  Element zero, one;
  _F.init(zero, 0);
  _F.init(one, 1);

  Vector x, y;

  size_t cols = C.coldim();
 
  *out_ptr 
    << "The product of these two matrices contains the following elements: "
    <<  endl;
  
  for (size_t i = 0; i < cols; i++)
  {
    x = Vector();
    x.push_back(make_pair(i, one));
    y = C.apply(x);

    *out_ptr << "** Column " << i << ":" << endl;
    write(*out_ptr, y);
    *out_ptr << endl;

  } // for (size_t i = 0; i < cols; i++)
  
  return true;
} // bool test_multiply<sparse_sequence_vector_tag>::test(void) const

// Implementation of methods for sparse associative vectors

template <class Field, class Vector>
test_multiply<Field, 
              Vector, 
	      LinBox::vector_categories::sparse_associative_vector_tag>
::test_multiply(const Field& F,
		int mode,
		bool bbtime,
		bool givtime,
		istream& in = cin,
		ostream& out = cout,
		ostream& log = clog)
: test_base(in, out, log), vector_utility<Field, Vector>(F),
  _F(F), _mode(mode), _bbtimer(bbtime), _givtimer(givtime)
{
  // prompt for input of matrix
  if (prompt)
    cout << "Enter the name of the file from which" << endl
         << "to read the multiply matrix input." << endl
	 << "Enter 'cin' to read from standard input" << endl
	 << "or 'same' to continue reading from the" << endl
	 << "current input." << endl;

  string filename;
  *in_ptr >> filename;
  if (filename != string("same")) open_input(filename);

  if (in_ptr == &cin) 
    prompt = true;
  else 
    prompt = false;

} // test_multiply<sparse_associative_vector_tag>::test_multiply(...)

template <class Field, class Vector>
bool test_multiply<Field, 
                   Vector, 
		   LinBox::vector_categories::sparse_associative_vector_tag>
::test(void) const
{
  Matrix A(*LinBox::newSparsemat<Field, std::map<size_t, Element> >(_F,
								    0,
								    0,
								    prompt,
								    *in_ptr,
								    *out_ptr));

  Matrix B(*LinBox::newSparsemat<Field, std::map<size_t, Element> >(_F,
								    0,
								    0,
								    prompt,
								    *in_ptr,
								    *out_ptr));

  LinBox::multiply<Vector> C(&A, &B);

  Element zero, one;
  _F.init(zero, 0);
  _F.init(one, 1);

  Vector x, y;

  size_t cols = C.coldim();
 
  *out_ptr 
    << "The product of these two matrices contains the following elements: "
    <<  endl;
  
  for (size_t i = 0; i < cols; i++)
  {
    x = Vector();
    x.insert(make_pair(i, one));
    y = C.apply(x);

    *out_ptr << "** Column " << i << ":" << endl;
    write(*out_ptr, y);
    *out_ptr << endl;

  } // for (size_t i = 0; i < cols; i++)
  
  return true;
} // bool test_multiply<sparse_associative_vector_tag>::test(void) const

#endif // _TEST_MULTIPLY_
