/* File: src/examples/blackbox/test_wiedemann_det.h
 * Author: William J Turner for the LinBox group
 */

#ifndef _TEST_WIEDEMANN_DET_
#define _TEST_WIEDEMANN_DET_

#include "LinBox/wiedemann_det.h"

#include "../utils/vector_utility.h"

/** Class to test the Wiedemann determinant algorithm for a generic
 * blackbox matrix.
 * Templatized by field and vector types.
 * @see wiedemann_det
 */
template <class Field, class Vector>
class test_wiedemann_det : public test_base, public vector_utility<Field, Vector>
{
public:

  // Field element type
  typedef typename Field::element Element;
  
  // Sparsemat matrix type to test compose
  typedef LinBox::Sparsemat<Field, std::map<size_t, Element>, Vector> Matrix;

  /** Constructor from Field object
   * @param  F field in which arithmetic is done
   * @param  mode blackbox matrix apply mode
   * @param  bbtime boolean flag on whether to use blackbox timer
   * @param  givtime boolean flag on whether to use Givaro timer
   * @param  in  istream from which input is read (default = cin)
   * @param  out ostream to which output is written (default = cout)
   * @param  log ostream to which messages are logged (default = clog)
   */
  test_wiedemann_det(
			const Field& F,
			int mode,
			bool bbtime,
			bool givtime,
			istream& in = cin,
			ostream& out = cout,
			ostream& log = clog
			);

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

}; // class test_wiedemann_det : public test_base

// Implementation of methods

template <class Field, class Vector>
test_wiedemann_det<Field, Vector>::test_wiedemann_det(
		const Field& F,
		int mode,
		bool bbtime,
		bool givtime,
		istream& in = cin,
		ostream& out = cout,
		ostream& log = clog
		)
	: test_base(in, out, log), vector_utility<Field, Vector>(F),
		_F(F), _mode(mode), _bbtimer(bbtime), _givtimer(givtime)
{
  // prompt for input of matrix
  if (prompt)
    cout << "Enter the name of the file from which" << endl
         << "to read the wiedemann_det matrix input." << endl
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

} // test_wiedemann_det<Field, Vector>::test_wiedemann_det(...)

template <class Field, class Vector>
bool test_wiedemann_det<Field, Vector>::test(void) const
{
	Matrix A(*LinBox::newSparsemat<Field, std::map<size_t, Element> >
			(_F, 0, 0, prompt, *in_ptr, *out_ptr));

	if (A.rowdim() != A.coldim())
	{
		*out_ptr << "The matrix is not square; no determinant found." << endl;
		return false;
	}

	typename Field::randIter r(_F);
	
	*out_ptr << "The determinant is ";
	_F.write(*out_ptr, LinBox::wiedemann_det(_F, A, r));
	*out_ptr << endl;

  return true;

} // bool test_wiedemann_det<Field, Vector>::test(void) const

#endif // _TEST_WIEDEMANN_DET_
