/* File: src/examples/test_base.h
 * Author: William J Turner for the LinBox group
 */

#ifndef _TEST_BASE_
#define _TEST_BASE_

#include <iostream>
#include <fstream>
#include <string>

#include <vector>
#include <list>
#include <deque> // Has problems with iterator->()
#include <utility>

#include "LinBox/blackbox_archetype.h"
#include "Examples/fileutils.h"

/** Base class for running tests on LinBox code.
 * This class contains code for creating input and output streams
 * and also contains pointers to these streams.  Other LinBox test objects
 * may be derived from this to incorporate these streams and functions.
 */
class test_base
{
public:

  /** Constructor from int and array of C-style strings.
   * Creates input and output streams.
   * @param argc number of arguments
   * @param argv array of input arguments:
   *         argv[0]  program name, 
   *         argv[1]  file from which to read input (default = cin), 
   *         argv[2]  file to which to print output (default = cout), 
   *         argv[3]  file to which to log messages (default = clog)
   */
  test_base(int argc, char* argv[]);

  /** Constructor from references to input and output streams.
   * Sets pointers to correct streams.
   * @param  in  istream from which input is read (default = cin)
   * @param  out ostream to which output is written (default = cout)
   * @param  log ostream to which messages are logged (default = clog)
   */
  test_base(istream& in = cin, ostream& out = cout, ostream& log = clog);

  /// Destructor
  ~test_base(void);

protected:
  
  /** Print error and exit.
   * Prints error message
   * "ERROR: err err2"
   * to cerr and exits the program with return code 1.
   * @param  err   string containing error message.
   * @param  err2  string containing second error message.
   */
  void error(const string& err, const string& err2);

  /** Open input stream.
   * Opens input stream for file with name filename.
   * If filename is cin, sets pointer to cin.
   * @return boolean true if successful, false otherwise
   * @param  filename  name of input file to open
   */
  bool open_input(const string& filename);

  /** Open output stream.
   * Opens output stream for file with name filename.
   * If filename is cout, sets pointer to cout.
   * @return boolean true if successful, false otherwise
   * @param  filename  name of output file to open
   */
  bool open_output(const string& filename);

  /** Open output stream for logging.
   * Opens output stream for file with name filename.
   * If filename is clog, returns pointer to clog.
   * @return boolean true if successful, false otherwise
   * @param  filename  name of output file to open
   */
  bool open_log(const string& filename);

  /** Apply blackbox matrix archetype object to vector.
   * Templatized on \Ref{LinBox} vector type.
   * @return boolean true if ended successfully, false if not
   * @param  y reference to \Ref{LinBox} vector object to contain output
   * @param  A constant reference to blackbox matrix object
   * @param  x constant reference to \Ref{LinBox} vector object containing input
   * @param  mode integer marking which mode of apply to use
   */
  template <class Vector>
  bool apply_blackbox(Vector& y, 
		      const LinBox::Blackbox_archetype<Vector>& A, 
		      const Vector& x,
		      int mode) const;
 
  /** Apply transpose of blackbox matrix archetype object to vector.
   * Templatized on \Ref{LinBox} vector type.
   * @return boolean true if ended successfully, false if not
   * @param  y reference to \Ref{LinBox} vector object to contain output
   * @param  A constant reference to blackbox matrix object
   * @param  x constant reference to \Ref{LinBox} vector object containing input
   * @param  mode integer marking which mode of apply to use
   */
  template <class Vector>
  bool applyTranspose_blackbox(Vector& y, 
			       const LinBox::Blackbox_archetype<Vector>& A, 
			       const Vector& x,
			       int mode) const;
  
  // Pointers to input and output streams
  istream* in_ptr;
  ostream* out_ptr;
  ostream* log_ptr;

  // Boolean for whether to prompt for input
  bool prompt;
  
}; // class test_base

// Implementation of methods

test_base::test_base(int argc, char* argv[])
  : in_ptr(&cin), out_ptr(&cout), log_ptr(&clog), prompt(true)
{
#ifdef TRACE
  *log_ptr << argc << endl;

  for (int i  = 0; i < argc; i++) *log_ptr << argv[i] << endl;
#endif // TRACE

  if (argc > 3) open_log(argv[3]);
  if (argc > 1) open_input(argv[1]);
  if (argc > 2) open_output(argv[2]);
  
  if (in_ptr != &cin) prompt = false;

} // test_base::test_base(int argc, char* argv[])

test_base::test_base(istream& in = cin, 
		     ostream& out = cout, 
		     ostream& log = clog)
  : in_ptr(&in), out_ptr(&out), log_ptr(&log), prompt(true)
{ if (in_ptr != &cin) prompt = false; }

test_base::~test_base(void) 
{
/*
  if (in_ptr != &cin) delete in_ptr;
  if (out_ptr != &cout) delete out_ptr;
  if (log_ptr != &clog) delete log_ptr;
*/
} // test_base::~test_base(void)

void test_base::error(const string& err, const string& err2)
{
  cerr << "ERROR: " << err << err2 << endl;
  exit(1);
} // void test_base::error(const string& err, const string& err2)

bool test_base::open_input(const string& filename) 
{
  if (filename == string("cin")) 
  {
#ifdef TRACE
    *log_ptr << "Using standard input cin" << endl;
#endif // TRACE
    
    in_ptr = &cin;
    return true;
  } // if (filename == string("cin"))

#ifdef TRACE
  *log_ptr << "Opening input file " << filename << endl;
#endif // TRACE
  
  ifstream* input = new ifstream(filename.c_str());

  if (!input->is_open())
  {
    error("cannot open input file ",filename);
    return false;
  } // if (!input->is_open())
  else
  {
    in_ptr = input;
    return true;
  }
  
} // bool test_base::open_input(const string& filename)
  
bool test_base::open_output(const string& filename) 
{
  if (filename == string("cout")) 
  {
#ifdef TRACE
    *log_ptr << "Using standard output cout" << endl;
#endif // TRACE
    
    out_ptr = &cout;
    return true;
  } // if (filename == string("cout"))

#ifdef TRACE
  *log_ptr << "Opening output file " << filename << endl;
#endif // TRACE
  
  ofstream* output = new ofstream(filename.c_str());

  if (!output->is_open())
  {
    error("cannot open output file ",filename);
    return false;
  } // if (!output->is_open())
  else
  {
    out_ptr = output;
    return true;
  }
  
} // bool test_base::open_output(const string& filename) 

bool test_base::open_log(const string& filename) 
{
  if (filename == string("clog")) 
  {
#ifdef TRACE
    clog << "Using standard log clog" << endl;
#endif // TRACE
    
    log_ptr = &clog;
    return true;
  } // if (filename == string("clog")) 

#ifdef TRACE
  clog << "Opening log file " << filename << endl;
#endif // TRACE

  ofstream* log = new ofstream(filename.c_str());

  if (!log->is_open())
  {
    error("cannot open log file ", filename);
    return false;
  } // if (!log->is_open())
  else
  {
    log_ptr = log;
    return true;
  }
    
} // bool test_base::open_log(const string& filename) 

template <class Vector>
bool test_base::apply_blackbox(Vector& y, 
			       const LinBox::Blackbox_archetype<Vector>& A, 
			       const Vector& x,
			       int mode) const
{
  if (mode == 1)
  {
    y = A.apply(x);
  }
  else if (mode == 2)
    A.apply(y, x);
  else if (mode == 3)
  {
    y = x;
    A.applyin(y);
  }
  else
    return false;  // return failure flag if mode not known.

  return true;
  
} // bool test_base::apply_blackbox<Vector>(...) const

template <class Vector>
bool 
test_base::applyTranspose_blackbox(Vector& y, 
				   const LinBox::Blackbox_archetype<Vector>& A, 
				   const Vector& x,
				   int mode) const
{
  if (mode == 1)
  {
    y = A.applyTranspose(x);
  }
  else if (mode == 2)
    A.applyTranspose(y, x);
  else if (mode == 3)
  {
    y = x;
    A.applyinTranspose(y);
  }
  else
    return false;  // return failure flag if mode not known.

  return true;
  
} // bool test_base::applyTranspose_blackbox<Vector>(...) const

#endif // _TEST_BASE_
