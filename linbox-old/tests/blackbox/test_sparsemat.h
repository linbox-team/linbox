/* File: src/examples/blackbox/test_sparsemat.h
 * Author: William J Turner for the LinBox group
 */

#ifndef _TEST_SPARSEMAT_
#define _TEST_SPARSEMAT_

#include "LinBox/sparsemat.h"
#include "LinBox/unparam_field.h"

#include "Examples/bbtimer.h"
#include "Examples/givtimer.h"
#include "Examples/test_base.h"
#include "Examples/vector_utility.h"

/** Class to test sparsemat blackbox matrix.
 * Templatized by field, vector, and row types.
 * @see sparsemat
 */
template <class Field, class Vector, class Row>
class test_sparsemat : public test_base, public vector_utility<Field, Vector>
{
public:

  /// Sparsemat matrix type
  typedef LinBox::sparsemat<Field, Row, Vector> matrix_type;

  /** Constructor from Field object
   * @param  F field in which arithmetic is done
   * @param  mode blackbox matrix apply mode
   * @param  bbtime boolean flag on whether to use blackbox timer
   * @param  givtime boolean flag on whether to use Givaro timer
   * @param  in  istream from which input is read (default = cin)
   * @param  out ostream to which output is written (default = cout)
   * @param  log ostream to which messages are logged (default = clog)
   */
  test_sparsemat(const Field& F,
		 int mode,
		 bool bbtime,
		 bool givtime,
		 istream& in = cin,
		 ostream& out = cout,
		 ostream& log = clog);

  /** Run tests on LinBox sparsemat blackbox matrices.
   * Tests the template created in sparsemat.h using the field F.
   * Creates sparse matrix A and right side vector b from input
   * stream by calling newMatrix() and newVector().
   * If not timing, copies matrix, swaps rows 3 and 4 and adds 4 times
   * row 1 to row 2 in the new matrix.
   * Then, when timing and when not, performs Gaussian elimination on
   * the matrix and right side.
   * Finally, if not timing, computes solution to linear equation
   * Ax = b.
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

}; // class test_sparsemat : public test_base

// Implementation of methods

template <class Field, class Vector, class Row>
test_sparsemat<Field, Vector, Row>::test_sparsemat(const Field& F,
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
         << "to read the Sparsemat matrix and vector input." << endl
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

} // test_sparsemat<Field, Vector, Row>::test_sparsemat(...)

template <class Field, class Vector, class Row>
bool test_sparsemat<Field, Vector, Row>::test(void) const
{
  // read matrix and vector from input
  matrix_type 
    A(*LinBox::newSparsemat<Field, Row>(_F,0,0,prompt,*in_ptr,*out_ptr));
  Vector b(*new_vector(A.get_rowdim(),prompt,*in_ptr,*out_ptr));
  
  // make copies of input to change
  matrix_type A1(A);
  Vector b1(b);

  // set boolean for overall timer
  // and create objects for timer
  // and create objects for timers
  bool timer(_bbtimer || _givtimer);
  bbtimer B;
  time_t time0, time1;
  clock_t clock0, clock1;
  Timer T;

  if (timer)
  {
    *log_ptr << "Starting test_sparsemat timers." << endl;

    if (_givtimer)
    {
      T.start();
      time( &time0 );
      clock0 = clock();
    } // if (_givtimer)

    if (_bbtimer) B = bbtimer();
  } // if (timer)
  else // if (!timer)
  {
    *out_ptr << "The matrix A contains the following elements: " << endl;
    A1.write(*out_ptr);
    *out_ptr << endl 
             << "The vector b contains the following elements: " << endl;
    write(*out_ptr, b1);

    if (A1.rowdim() >= 3)
    {
      A1.swaprow(2,3);
      *out_ptr << endl << "Swapping rows 2 and 3: " << endl << A1 << endl;
    }

    if (A1.rowdim() >= 1)
    {
      typename Field::element a;
      _F.init(a, 4);
      
      A1.addrow(0,1,a);
      *out_ptr << endl << "Adding 4*(row 0) to row 1:" << endl << A1;
    }

    A1 = A;

  } // if (!timer)

  *out_ptr << endl << "Performing Gaussian elimination" << endl;
  A1.gauss(b1);

  if (timer)
  {
    if (_givtimer)
    {
      time( &time1 );
      clock1 = clock();
      
      T.stop();
      
/**/
      *log_ptr << endl << "Time\n" << T << endl;
      
      *log_ptr << endl
               << "Processor time: \t"
	       << ( (double) (clock1 - clock0) ) / CLOCKS_PER_SEC
	       << " sec.  \tCalendar time: "
	       << difftime( time1, time0 ) << " sec." << endl;
/**/
      
    } // if (_givtimer)

    if (_bbtimer) *log_ptr << "Sparsetest BBTimer: " << B;

    *out_ptr << endl << "A(m,n) = ";
    _F.write(*out_ptr, A1[make_pair(A1.rowdim()-1,A1.coldim()-1)]);
    *out_ptr << endl;

  } // if (timer)
  else // if (!timer)
  {
    *out_ptr << "The matrix A contains the following elements: " << endl;
    A1.write(*out_ptr);
    *out_ptr << endl << "The vector b contains the following elements: " << endl;
    write(*out_ptr, b1);

    Vector x(A1.linsolve(b1));
    
    *out_ptr << endl 
        << "The solution to the system A*x = b contains the folowing elements:"
	<< endl;
    write(*out_ptr, x);
    
    Vector b2;
    apply_blackbox(b2, A, x, _mode);
    
    *out_ptr << endl 
         << "Applying the solution to the original matrix gives the following:"
 	 << endl;
    write(*out_ptr, b2);

  } // if (!timer)

  if (in_ptr != &cin) delete in_ptr;

  return true;

} // bool test_sparsemat<Field, Vector, Row>::test(void) const

/** This class contains code for creating a random matrix of integers.
 * It is derived from test_base to include functions and files for
 * creating input and output streams.
 * It prompts for the size of the matrix.  It then creates output that can be
 * read into the \Ref{test_sparsemat} programs.
 */
class random_sparsemat : public test_base
{
public:
  /** Constructor from references to input and output streams.
   * @param m row dimensions of matrix (default = 0)
   * @param n column dimensions of matrix (default = 0)
   * @param  in  istream from which input is read (default = cin)
   * @param  out ostream to which output is written (default = cout)
   * @param  log ostream to which messages are logged (default = clog)
   */
  random_sparsemat(size_t m = 0, 
		   size_t n = 0,
		   istream& in = cin, 
		   ostream& out = cout, 
		   ostream& log = clog);

  /// Destructor
  ~random_sparsemat(void) {}

  /** Create random dense matrix.
   * Prompts for the size of the matrix and then creates output that
   * can be read into the \Ref{test_sparsemat} programs.
   * It will prompt for where to place the output.
   * If the row and column dimensions are zero, it will prompt for those 
   * as well.
   * It will also ask whether to create a corresponding vector b for solving
   * the system A*x = b.  If it creates the vector, it will also ask for
   * the desired solution vector x.
   * @return boolean true if successfull, false otherwise
   */
  bool create(void) const;
 
private:

  // row and column dimensions
  size_t _m, _n;

}; // class random_sparsemat : public test_base

// Implementation of methods

random_sparsemat::random_sparsemat(size_t m = 0, 
				   size_t n = 0,
				   istream& in = cin, 
				   ostream& out = cout, 
				   ostream& log = clog)
: test_base(in, out, log), _m(m), _n(n)
{
  while ( (_m <= 0) || (_n <= 0) )
  {
    if (prompt) cout << "Enter the row and column dimensions of the matrix: ";
    *in_ptr >> _m >> _n;
  } // while ( (_m <= 0) || (_n <= 0) )

  if ( (_n + _m > 10) && (&out == &cout) )
  {
    if (prompt)
      cout << "Output will be too large for standard out." << endl
	   << "Enter the name of the file to which" << endl
	   << "to write output." << endl;
    string filename;
    *in_ptr >> filename;
    open_input(filename);

    if (&out == &cout) exit(1);

  } // if ( (_n + _m > 10) && (&out == &cout) )

} // random_sparsemat::random_sparsemat(size_t n = 0, ...)

bool random_sparsemat::create(void) const
{
#ifdef TRACE
  *log_ptr << "Creating a matrix of " << _m << " rows and " << _n 
           << " columns." << endl;
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
 
  string temp_filename = "/tmp/sparsemat";
  
  if (prompt) 
    cout << "Enter the number corresponding to the vector b to create:" << endl
         << "  0: no vector b" << endl
	 << "  1: random vector b" << endl
	 << "  2: b = zeros(n)" << endl
	 << "  3: b = ones(n)" << endl
	 << "  4: b = A*ones(n)" << endl
	 << "  5: b = A*(user-defined vector)" << endl;

  size_t value;
  *in_ptr >> value;

  if ( (value < 0) || (value > 4) )
  {
    cout << "Invalid response: " << value << ". Please try again: ";
    *in_ptr >> value;
  } // if ( (value < 0) || (value > 4) )

  typedef LinBox::sparsemat<LinBox::unparam_field<long>, 
                            std::map<size_t, long>, 
	   	            std::vector<long> > matrix_type;
  LinBox::unparam_field<long> F;
  matrix_type A(F, _m, _n);

  *out_ptr << _m << " " << _n << endl;
  
  int a;
  
  for (size_t i = 0; i < _m; i++)
    for (size_t j = 0; j < _n; j++)
    {
      a = static_cast<long>((double(rand())/RAND_MAX)*max);      
      *out_ptr << i << " " << j << " " << a << endl;
      if (value >= 4) A.put_value(make_pair(i, j), a);
    }

  *out_ptr << -1 << endl;

  if (value == 0) 
    return true;
  else if (value == 1)
  {
    for (size_t i = 0; i < _n; i++)
      *out_ptr << i << " " 
	       << static_cast<long>((double(rand())/RAND_MAX)*max) << endl;
    
    *out_ptr << -1 << endl;

    return true;
  } // else if (value == 1)
  else if (value == 2)
  {
    for (size_t i = 0; i < _n; i++)
      *out_ptr << i << " " << 0 << endl;
    
    *out_ptr << -1 << endl;

    return true;
  } // else if (value == 2)
  else if (value == 3)
  {
    for (size_t i = 0; i < _n; i++)
      *out_ptr << i << " " << 1 << endl;
    
    *out_ptr << -1 << endl;

    return true;
  } // else if (value == 3)
  
  std::vector<long> x(_m, 1);

  if (value == 5)
  {
    if (prompt)
      cout << endl << "Input vector by entering index and value." << endl
	           << "Remember vector is indexed starting at 0." << endl
		   << "End with a index of -1." << endl;

    size_t i;
    while (*in_ptr >> i)
    {
      if (i == size_t(-1)) break;
      *in_ptr >> x[i];
    } // while (*in_ptr >> i)

#ifdef TRACE
    *log_ptr << "Input vector:" << endl
             << "  i        x[i]" << endl
	     << "  -------------" << endl;
    
    for (size_t i = 0; i < _n; i++)
      *log_ptr << "  " << i << "    " << x[i] << endl;
#endif // TRACE

  } // else if (value == 5)

#ifdef TRACE
  *log_ptr << "Input vector:" << endl
           << "  i        x[i]" << endl
           << "  -------------" << endl;
    
  for (size_t i = 0; i < x.size(); i++)
    *log_ptr << "  " << i << "    " << x[i] << endl;
  
  *log_ptr << "Using matrix:" << endl << A << endl;
#endif // TRACE

  std::vector<long> y = A.apply(x);

#ifdef TRACE
  *log_ptr << "Output vector:" << endl
    << "  i        y[i]" << endl
    << "  -------------" << endl;
  
  for (size_t i = 0; i < y.size(); i++)
    *log_ptr << "  " << i << "    " << y[i] << endl;
#endif // TRACE

  for (size_t i = 0; i < y.size(); i++)
    *out_ptr << i << " " << y[i] << endl;

  *out_ptr << -1 << endl;

  return true;

} // bool random_sparsemat::create(void) const

#endif // _TEST_SPARSEMAT_
