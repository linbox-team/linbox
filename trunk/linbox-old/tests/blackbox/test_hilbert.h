/* File: src/examples/blackbox/test_hilbert.h
 * Author: William J Turner for the LinBox group
 */

#ifndef _TEST_HILBERT_
#define _TEST_HILBERT_

#include "LinBox/hilbert.h"

#include "Examples/vector_utility.h"

/** Class to test hilbert blackbox matrix.
 * Templatized by field and vector types.
 * @see hilbert
 */
template <class Field, class Vector>
class test_hilbert : public test_base, public vector_utility<Field, Vector>
{
public:

  /// Hilbert matrix type
  typedef LinBox::hilbert<Field, Vector> matrix_type;

  /** Constructor from Field object
   * @param  F field in which arithmetic is done
   * @param  mode blackbox matrix apply mode
   * @param  bbtime boolean flag on whether to use blackbox timer
   * @param  givtime boolean flag on whether to use Givaro timer
   * @param  in  istream from which input is read (default = cin)
   * @param  out ostream to which output is written (default = cout)
   * @param  log ostream to which messages are logged (default = clog)
   */
  test_hilbert(const Field& F,
	       int mode,
	       bool bbtime,
	       bool givtime,
	       istream& in = cin,
	       ostream& out = cout,
	       ostream& log = clog);

  /** Run tests on LinBox hilbert blackbox matrices.
   * Tests the template created in hilbert.h using the field F.
   * Applies hilbert matrix to vector x.
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

}; // class test_hilbert : public test_base

// Implementation of methods

template <class Field, class Vector>
test_hilbert<Field, Vector>::test_hilbert(const Field& F,
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
         << "to read the Hilbert vector input." << endl
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

} // test_hilbert<Field, Vector>::test_hilbert(...)

template <class Field, class Vector>
bool test_hilbert<Field, Vector>::test(void) const
{
  // read vector from input
  Vector x = *new_vector(0,prompt,*in_ptr,*out_ptr);

  *out_ptr << endl
       << "The input vector to the Hilbert test is the following:" << endl;

  write(*out_ptr, x);

  size_t n = x.size();

  typename Field::element zero;
  _F.init(zero, 0);

  Vector b;

  LinBox::hilbert<Field, Vector> H(_F,n);

  apply_blackbox(b, H, x, _mode);

  *out_ptr << endl
       << "Applying the Hilbert matrix to the vector gives " << endl
       << "the following vector" << endl;

  write(*out_ptr, b);

  return true;

} // bool test_hilbert<Field, Vector>::test(void) const

#endif // _TEST_HILBERT_
