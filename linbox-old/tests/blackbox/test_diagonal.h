/* File: src/examples/blackbox/test_diagonal.h
 * Author: William J Turner for the LinBox group
 */

#ifndef _TEST_DIAGONAL_
#define _TEST_DIAGONAL_

#include "LinBox/diagonal.h"

#include "../utils/vector_utility.h"

/** Class to test diagonal blackbox matrix.
 * Templatized by field and vector types.
 * @see diagonal
 */
template <class Field, class Vector>
class test_diagonal : public test_base, public vector_utility<Field, Vector>
{
public:

  /// Diagonal matrix type
  typedef LinBox::Diagonal<Field, Vector> matrix_type;

  /** Constructor from Field object
   * @param  F field in which arithmetic is done
   * @param  mode blackbox matrix apply mode
   * @param  bbtime boolean flag on whether to use blackbox timer
   * @param  givtime boolean flag on whether to use Givaro timer
   * @param  in  istream from which input is read (default = cin)
   * @param  out ostream to which output is written (default = cout)
   * @param  log ostream to which messages are logged (default = clog)
   */
  test_diagonal(const Field& F,
	       int mode,
	       bool bbtime,
	       bool givtime,
	       istream& in = cin,
	       ostream& out = cout,
	       ostream& log = clog);

  /** Run tests on LinBox diagonal blackbox matrices.
   * Tests the template created in diagonal.h using the field F.
   * Applies diagonal matrix to vector x.
   * @return boolean true if ended successfully, false if not
   */
  bool test(void) const;

private:

  // Field in which arithmetic is done
  Field _F;

  // Blackbox matrix apply mode
  int _mode;

  // Boolean flags on whether to use timers.
  bool _bbtimer, _givtimer;

}; // class test_diagonal : public test_base

// Implementation of methods

template <class Field, class Vector>
test_diagonal<Field, Vector>::test_diagonal(const Field& F,
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
         << "to read the diagonal vector input." << endl
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

} // test_diagonal<Field, Vector>::test_diagonal(...)

template <class Field, class Vector>
bool test_diagonal<Field, Vector>::test(void) const
{
  // read diagonal vector from input
  if (prompt)
    cout << "What is the size of the square diagonal matrix to create? ";
  
  size_t n;
  *in_ptr >> n;

#ifdef TRACE
  *out_ptr << endl << "The size of the diagonal matrix is " << n << endl;
#endif

  typedef typename Field::element Element;
  Element elem;
  _F.init(elem, 0);
  std::vector<Element> d(n, elem);
  
  if (prompt)
    cout << endl
         << "Input diagonal entries by entering index and value." << endl
	 << "Remember matriices and vectors are indexed starting at 0." << endl
	 << "End with a index of -1." << endl;
 
  size_t i;
  
  while (*in_ptr >> i)
  {
    if(i == size_t(-1)) break; // return also if row index is not positive
    _F.read(*in_ptr, elem);
    d[i] = elem;
  } // while (is >> i)

  *out_ptr << endl
       << "The diagonal vector in the diagonal test is the following:" << endl;

  typename std::vector<Element>::iterator iter;
  i = 0;
  
  for (iter = d.begin(); iter != d.end(); iter++, i++)
  {
    *out_ptr << "Element " << i << ": ";
    _F.write(*out_ptr, *iter); // can't do this!!!
    *out_ptr << endl;
  }

  if (prompt)
    cout << "Enter a vector to be multiplied by the diagonal matrix:"
         << endl;

  Vector x = *new_vector(n,prompt,*in_ptr,*out_ptr);

  *out_ptr << endl
       << "The input vector to the diagonal test is the following:" << endl;

  write(*out_ptr, x);

  typename Field::element zero;
  _F.init(zero, 0);

  Vector b;

  LinBox::Diagonal<Field, Vector> D(_F, d);

  apply_blackbox(b, D, x, _mode);

  *out_ptr << endl
       << "Applying the diagonal matrix to the vector gives " << endl
       << "the following vector" << endl;

  write(*out_ptr, b);

  return true;

} // bool test_diagonal<Field, Vector>::test(void) const

#endif // _TEST_DIAGONAL_
