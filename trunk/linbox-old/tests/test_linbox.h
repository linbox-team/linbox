/* File: src/examples/test_linbox.h
 * Author: William J Turner for the LinBox group
 */

#ifndef _TEST_LINBOX_
#define _TEST_LINBOX_

#include <iostream>
#include <vector>
#include <map>
#include <list>
#include <deque>
#include <utility>

#include "test_base.h"
#include "field/test_field.h"
#include "blackbox/test_sparsemat.h"
#include "blackbox/test_hilbert.h"
#include "blackbox/test_butterfly.h"
#include "blackbox/test_compose.h"

/** Class to test LinBox code.
 * This class is derived from \Ref{test_base} which contains code to 
 * open input and output streams and contains pointers to these 
 * streams.
 */
class test_linbox : public test_base
{
public:
  /** Constructor from int and array of C-style strings.
   * Creates input and output streams and calls \Ref{test_linbox}
   * to do actual test of LinBox.
   * @param argc number of arguments
   * @param argv array of input arguments:
   *         argv[0]  program name, 
   *         argv[1]  file from which to read input (default = cin), 
   *         argv[2]  file to which to print output (default = cout), 
   *         argv[3]  file to which to log messages (default = clog)
   */
  test_linbox(int argc, char* argv[]) : test_base(argc, argv) {}

  /// Destructor
  ~test_linbox(void) {}

  /** Test LinBox code.
   * Tests LinBox code for the field corresponding to the \Ref{field_categories}
   * structure tag.
   * Specialization for all_tag prompts user for which field arithmetic to use, 
   * and then calls specialization for that field tag.
   * Each specialization for a specific field tag then creates the field and 
   * calls the private method run_tests(const Field& F) which prompts for 
   * input and then creates and uses test objects for each blackbox matrix 
   * tested.
   * @see test_sparsemat
   * @see test_hilbert
   * @see test_butterfly
   * @return boolean true if succesfull, false otherwise
   */
  template <class Trait> bool test(void) const;

  /** List of vector categories.
   * This structure contains nineteen structures: one relating each field type
   * to be tested plus one for all fields combined.  
   * These tags are only needed by this class, so they are private members.
   */
  struct field_categories
  {
    struct all_tag {};
    struct double_tag {};
    struct float_tag {};
    struct modular_tag {};
    struct fuzzy_tag {};
    struct ntl_rr_tag{};
    struct ntl_ZZ_p_tag{};
    struct ntl_zz_p_tag{};
    struct param_modular_tag {};
    struct param_fuzzy_tag {};
    struct double_envelope_tag {};
    struct float_envelope_tag {};
    struct modular_envelope_tag {};
    struct fuzzy_envelope_tag {};
    struct ntl_rr_envelope_tag{};
    struct ntl_ZZ_p_envelope_tag{};
    struct ntl_zz_p_envelope_tag{};
    struct param_modular_envelope_tag {};
    struct param_fuzzy_envelope_tag {};
    struct abstract_double_tag {};
    struct abstract_float_tag {};
    struct abstract_modular_tag {};
    struct abstract_fuzzy_tag {};
//    struct abstract_ntl_rr_tag{};
//    struct abstract_ntl_ZZ_p_tag{};
//    struct abstract_ntl_zz_p_tag{};
    struct abstract_param_modular_tag{};
    struct abstract_param_fuzzy_tag{};
  }; // struct field_categories

private:

  // References to input and output streams

  /** Choose vector type.
   * Prompts user whether to test the field, 
   * and then for which LinBox vector type to use, 
   * and finally runs tests.
   * @param F Field object to pass to test_matrices
   */
  template <class Field> bool run_tests(const Field& F) const;

}; // class test_linbox : public test_base

// Implementation of methods

template <class Field> bool test_linbox::run_tests(const Field& F) const
{
  // Define element type
  typedef typename Field::element Element;

  // Which tests to run?

  if (prompt) 
    cout << "Enter the numbers corresponding to the tests you want to run." 
         << endl
         << "End with a test number of 0." << endl
         << "  1: field, elements, and random generators" << endl
         << "  2: sparsemat" << endl
         << "  3: hilbert" << endl
         << "  4: butterfly switching network" << endl
	 << "  5: multiplication of two sparsemat matrices" << endl;

  bool field(false), sparsemat(false), hilbert(false), butterfly(false), multiply(false);

  int value;
  while (*in_ptr >> value)
  {
    if (value == 0) break;

    if (value == 1) field = true;
    else if (value == 2) sparsemat = true;
    else if (value == 3) hilbert = true;
    else if (value == 4) butterfly = true;
    else if (value == 5) multiply = true;

  } // while (*in_ptr >> value)

#ifdef TRACE
  *log_ptr << "You have selected the following matrices to test:" << endl;
  if (sparsemat) *log_ptr << "    sparsemat" << endl;
  if (hilbert) *log_ptr << "    hilbert" << endl;
  if (butterfly) *log_ptr << "    butterfly switching network" << endl;
  if (multiply) *log_ptr << "    multiplication of two matrices" << endl;
#endif // TRACE

  // Test field
  if (field)  
  {
    test_field<Field> T(F, *in_ptr, *out_ptr, *log_ptr);
    T.test();
  }

  if (!(sparsemat || hilbert || butterfly || multiply))
  {
#ifdef TRACE
    *log_ptr << "You have not selected any matrices to test." << endl;
#endif // TRACE
    return true;
  }

  // Which vector type to use?

  if (prompt) 
    cout << "Enter the number corresponding to the vector type you want to use:"
         << endl
	 << "  1 : Dense vector implemented as STL vector" << endl
	 << "  2 : Sparse vector implemented as STL map" << endl
	 << "  3 : Sparse vector implemented as STL list" << endl
	 << "  4 : Sparse vector implemented as STL vector" << endl
	 << "  5 : Sparse vector implemented as STL deque" << endl;

  // retrieve input
  int vector;
  *in_ptr >> vector;
  
  // Print error and get new input if number is out of range
  while ( (vector < 1) || (vector > 5) )
  {
    *out_ptr << "Invalid response: " << vector << ".  Please try again: ";
    *in_ptr >> vector;
  }

#ifdef TRACE
  *log_ptr << "You have entered vector type " << vector << endl;
#endif // TRACE

  // Which blackbox matrix apply to use?

  // prompt for input
  if (prompt) 
    cout << "Enter the number corresponding to the black box matrix" << endl
         << "apply you want to use." << endl
         << "  1 : y = A.apply(x)" << endl
         << "  2 : A.apply(y, x)" << endl
         << "  3 : A.applyin(x)" << endl;
 
  int mode;
  *in_ptr >> mode;

  // Print error and get new input if number is out of range
  while ( (mode < 1) || (mode > 3) )
  {
    *out_ptr << "Invalid response: " << mode << ".  Please try again: ";
    *in_ptr >> mode;
  }

#ifdef TRACE
  *log_ptr << "You have entered mode " << mode << endl;
#endif // TRACE

  if (prompt) 
    cout << "Enter the numbers corresponding to the timers" << endl
         << "you want to use.  End with a timer number of 0." << endl
         << "  1: blackbox timer" << endl
         << "  2: Givaro timer" << endl;

  bool bbtimer(false), givtimer(false);

  while (*in_ptr >> value)
  {
    if (value == 0) break;

    if (value == 1) bbtimer = true;
    else if (value == 2) givtimer = true;

  } // while (*in_ptr >> value)

#ifdef TRACE
  *log_ptr << "You have selected the following timers to use:" << endl;
  if (bbtimer) *log_ptr << "    blackbox timer" << endl;
  if (givtimer) *log_ptr << "    Givaro timer" << endl;
#endif // TRACE

  if (sparsemat) 
  {
    if (prompt)
      cout << "Enter the number corresponding to the sparse vector type" << endl
	   << "you want to use for the rows of sparsemat." << endl
	   << "  1 : Sparse vector implemented as STL map" << endl
	   << "  2 : Sparse vector implemented as STL list" << endl
	   << "  3 : Sparse vector implemented as STL vector" << endl
	   << "  4 : Sparse vector implemented as STL deque" << endl;

    // retrieve input
    *in_ptr >> value;
  
    // Print error and get new input if number is out of range
    while ( (value < 1) || (value > 4) )
    {
      *out_ptr << "Invalid response: " << value << ".  Please try again: ";
      *in_ptr >> value;
    }

#ifdef TRACE
    *log_ptr << "You have entered row type " << value << endl;
#endif // TRACE

    if ( (vector == 1) && (value == 1) )
    {
#ifdef TRACE
      *log_ptr	<< "Testing sparsemat with :" << endl
	        << "  vector type: dense STL vectors" << endl
		<< "  row type: sparse STL maps" << endl
		<< "  apply mode: " << mode << endl
		<< "  timers: ";
      if (bbtimer) *log_ptr << "blackbox timer, ";
      if (givtimer) *log_ptr << "Givaro timer, ";
      *log_ptr << endl;
#endif TRACE

      test_sparsemat<Field, std::vector<Element>, std::map<size_t, Element> >
	T(F, mode, bbtimer, givtimer, *in_ptr, *out_ptr, *log_ptr);
      T.test();

    } // if ( (vector == 1) && (value == 1) )
    else if ( (vector == 1) && (value == 2) )
    {
#ifdef TRACE
      *log_ptr	<< "Testing sparsemat with :" << endl
	        << "  vector type: dense STL vectors" << endl
		<< "  row type: sparse STL lists" << endl
		<< "  apply mode: " << mode << endl
		<< "  timers: ";
      if (bbtimer) *log_ptr << "blackbox timer, ";
      if (givtimer) *log_ptr << "Givaro timer, ";
      *log_ptr << endl;
#endif TRACE

      test_sparsemat<Field, 
                     std::vector<Element>, 
                     std::list< pair<size_t, Element> > >
	T(F, mode, bbtimer, givtimer, *in_ptr, *out_ptr, *log_ptr);
      T.test();

    } // else if ( (vector == 1) && (value == 2) )
    else if ( (vector == 1) && (value == 3) )
    {
#ifdef TRACE
      *log_ptr	<< "Testing sparsemat with :" << endl
	        << "  vector type: dense STL vectors" << endl
		<< "  row type: sparse STL vectors" << endl
		<< "  apply mode: " << mode << endl
		<< "  timers: ";
      if (bbtimer) *log_ptr << "blackbox timer, ";
      if (givtimer) *log_ptr << "Givaro timer, ";
      *out_ptr << endl;
#endif TRACE

      test_sparsemat<Field, 
                     std::vector<Element>, 
		     std::vector< pair<size_t, Element> > >
	T(F, mode, bbtimer, givtimer, *in_ptr, *out_ptr, *log_ptr);
      T.test();

    } // else if ( (vector == 1) && (value == 3) )
    else if ( (vector == 1) && (value == 4) )
    {
#ifdef TRACE
      *log_ptr	<< "Testing sparsemat with :" << endl
	        << "  vector type: dense STL vectors" << endl
		<< "  row type: sparse STL deques" << endl
		<< "  apply mode: " << mode << endl
		<< "  timers: ";
      if (bbtimer) *log_ptr << "blackbox timer, ";
      if (givtimer) *log_ptr << "Givaro timer, ";
      *out_ptr << endl;
#endif TRACE

      test_sparsemat<Field, 
                     std::vector<Element>, 
		     std::deque< pair<size_t, Element> > >
	T(F, mode, bbtimer, givtimer, *in_ptr, *out_ptr, *log_ptr);
      T.test();

    } // else if ( (vector == 1) && (value == 4) )
    else if ( (vector == 2) && (value == 1) )
    {
#ifdef TRACE
      *log_ptr	<< "Testing sparsemat with :" << endl
	        << "  vector type: sparse STL maps" << endl
		<< "  row type: sparse STL maps" << endl
		<< "  apply mode: " << mode << endl
		<< "  timers: ";
      if (bbtimer) *log_ptr << "blackbox timer, ";
      if (givtimer) *log_ptr << "Givaro timer, ";
      *out_ptr << endl;
#endif TRACE

      test_sparsemat<Field, 
                     std::map<size_t, Element>, 
		     std::map<size_t, Element> >
	T(F, mode, bbtimer, givtimer, *in_ptr, *out_ptr, *log_ptr);
      T.test();

    } // else if ( (vector == 2) && (value == 1) )
    else if ( (vector == 2) && (value == 2) )
    {
#ifdef TRACE
      *log_ptr	<< "Testing sparsemat with :" << endl
	        << "  vector type: sparse STL maps" << endl
		<< "  row type: sparse STL lists" << endl
		<< "  apply mode: " << mode << endl
		<< "  timers: ";
      if (bbtimer) *log_ptr << "blackbox timer, ";
      if (givtimer) *log_ptr << "Givaro timer, ";
      *out_ptr << endl;
#endif TRACE

      test_sparsemat<Field, 
                     std::map<size_t, Element>, 
		     std::list< pair<size_t, Element> > >
	T(F, mode, bbtimer, givtimer, *in_ptr, *out_ptr, *log_ptr);
      T.test();

    } // else if ( (vector == 2) && (value == 2) )
    else if ( (vector == 2) && (value == 3) )
    {
#ifdef TRACE
      *log_ptr	<< "Testing sparsemat with :" << endl
	        << "  vector type: sparse STL maps" << endl
		<< "  row type: sparse STL vectors" << endl
		<< "  apply mode: " << mode << endl
		<< "  timers: ";
      if (bbtimer) *log_ptr << "blackbox timer, ";
      if (givtimer) *log_ptr << "Givaro timer, ";
      *out_ptr << endl;
#endif TRACE

      test_sparsemat<Field, 
                     std::map<size_t, Element>, 
		     std::vector< pair<size_t, Element> > >
	T(F, mode, bbtimer, givtimer, *in_ptr, *out_ptr, *log_ptr);
      T.test();

    } // else if ( (vector == 2) && (value == 3) )
    else if ( (vector == 2) && (value == 4) )
    {
#ifdef TRACE
      *log_ptr	<< "Testing sparsemat with :" << endl
	        << "  vector type: sparse STL maps" << endl
		<< "  row type: sparse STL deques" << endl
		<< "  apply mode: " << mode << endl
		<< "  timers: ";
      if (bbtimer) *log_ptr << "blackbox timer, ";
      if (givtimer) *log_ptr << "Givaro timer, ";
      *out_ptr << endl;
#endif TRACE

      test_sparsemat<Field, 
                     std::map<size_t, Element>, 
		     std::deque< pair<size_t, Element> > >
	T(F, mode, bbtimer, givtimer, *in_ptr, *out_ptr, *log_ptr);
      T.test();

    } // else if ( (vector == 2) && (value == 4) )
    else if ( (vector == 3) && (value == 1) )
    {
#ifdef TRACE
      *log_ptr	<< "Testing sparsemat with :" << endl
	        << "  vector type: sparse STL lists" << endl
		<< "  row type: sparse STL maps" << endl
		<< "  apply mode: " << mode << endl
		<< "  timers: ";
      if (bbtimer) *log_ptr << "blackbox timer, ";
      if (givtimer) *log_ptr << "Givaro timer, ";
      *out_ptr << endl;
#endif TRACE

      test_sparsemat<Field, 
                     std::list< pair<size_t, Element> >,
		     std::map<size_t, Element> >
	T(F, mode, bbtimer, givtimer, *in_ptr, *out_ptr, *log_ptr);
      T.test();

    } // else if ( (vector == 3) && (value == 1) )
    else if ( (vector == 3) && (value == 2) )
    {
#ifdef TRACE
      *log_ptr	<< "Testing sparsemat with :" << endl
	        << "  vector type: sparse STL lists" << endl
		<< "  row type: sparse STL lists" << endl
		<< "  apply mode: " << mode << endl
		<< "  timers: ";
      if (bbtimer) *log_ptr << "blackbox timer, ";
      if (givtimer) *log_ptr << "Givaro timer, ";
      *out_ptr << endl;
#endif TRACE

      test_sparsemat<Field, 
                     std::list< pair<size_t, Element> >,
                     std::list< pair<size_t, Element> > >
	T(F, mode, bbtimer, givtimer, *in_ptr, *out_ptr, *log_ptr);
      T.test();

    } // else if ( (vector == 3) && (value == 2) )
    else if ( (vector == 3) && (value == 3) )
    {
#ifdef TRACE
      *log_ptr	<< "Testing sparsemat with :" << endl
	        << "  vector type: sparse STL lists" << endl
		<< "  row type: sparse STL vectors" << endl
		<< "  apply mode: " << mode << endl
		<< "  timers: ";
      if (bbtimer) *log_ptr << "blackbox timer, ";
      if (givtimer) *log_ptr << "Givaro timer, ";
      *out_ptr << endl;
#endif TRACE

      test_sparsemat<Field, 
                     std::list< pair<size_t, Element> >,
                     std::vector< pair<size_t, Element> > >
	T(F, mode, bbtimer, givtimer, *in_ptr, *out_ptr, *log_ptr);
      T.test();

    } // else if ( (vector == 3) && (value == 3) )
    else if ( (vector == 3) && (value == 4) )
    {
#ifdef TRACE
      *log_ptr	<< "Testing sparsemat with :" << endl
	        << "  vector type: sparse STL lists" << endl
		<< "  row type: sparse STL deques" << endl
		<< "  apply mode: " << mode << endl
		<< "  timers: ";
      if (bbtimer) *log_ptr << "blackbox timer, ";
      if (givtimer) *log_ptr << "Givaro timer, ";
      *out_ptr << endl;
#endif TRACE

      test_sparsemat<Field,
                     std::list< pair<size_t, Element> >,
                     std::deque< pair<size_t, Element> > >
	T(F, mode, bbtimer, givtimer, *in_ptr, *out_ptr, *log_ptr);
      T.test();

    } // else if ( (vector == 3) && (value == 4) )
    else if ( (vector == 4) && (value == 1) )
    {
#ifdef TRACE
      *log_ptr	<< "Testing sparsemat with :" << endl
	        << "  vector type: sparse STL vectors" << endl
		<< "  row type: sparse STL maps" << endl
		<< "  apply mode: " << mode << endl
		<< "  timers: ";
      if (bbtimer) *log_ptr << "blackbox timer, ";
      if (givtimer) *log_ptr << "Givaro timer, ";
      *out_ptr << endl;
#endif TRACE

      test_sparsemat<Field, 
                     std::vector< pair<size_t, Element> >,
                     std::map<size_t, Element> >
	T(F, mode, bbtimer, givtimer, *in_ptr, *out_ptr, *log_ptr);
      T.test();

    } // else if ( (vector == 4) && (value == 1) )
    else if ( (vector == 4) && (value == 2) )
    {
#ifdef TRACE
      *log_ptr	<< "Testing sparsemat with :" << endl
	        << "  vector type: sparse STL vectors" << endl
		<< "  row type: sparse STL lists" << endl
		<< "  apply mode: " << mode << endl
		<< "  timers: ";
      if (bbtimer) *log_ptr << "blackbox timer, ";
      if (givtimer) *log_ptr << "Givaro timer, ";
      *out_ptr << endl;
#endif TRACE

      test_sparsemat<Field, 
                     std::vector< pair<size_t, Element> >,
                     std::list< pair<size_t, Element> > >
	T(F, mode, bbtimer, givtimer, *in_ptr, *out_ptr, *log_ptr);
      T.test();

    } // else if ( (vector == 4) && (value == 2) )
    else if ( (vector == 4) && (value == 3) )
    {
#ifdef TRACE
      *log_ptr	<< "Testing sparsemat with :" << endl
	        << "  vector type: sparse STL vectors" << endl
		<< "  row type: sparse STL vectors" << endl
		<< "  apply mode: " << mode << endl
		<< "  timers: ";
      if (bbtimer) *log_ptr << "blackbox timer, ";
      if (givtimer) *log_ptr << "Givaro timer, ";
      *out_ptr << endl;
#endif TRACE

      test_sparsemat<Field, 
                     std::vector< pair<size_t, Element> >,
                     std::vector< pair<size_t, Element> > >
	T(F, mode, bbtimer, givtimer, *in_ptr, *out_ptr, *log_ptr);
      T.test();

    } // else if ( (vector == 4) && (value == 3) )
    else if ( (vector == 4) && (value == 4) )
    {
#ifdef TRACE
      *log_ptr	<< "Testing sparsemat with :" << endl
	        << "  vector type: sparse STL vectors" << endl
		<< "  row type: sparse STL deques" << endl
		<< "  apply mode: " << mode << endl
		<< "  timers: ";
      if (bbtimer) *log_ptr << "blackbox timer, ";
      if (givtimer) *log_ptr << "Givaro timer, ";
      *out_ptr << endl;
#endif TRACE

      test_sparsemat<Field, 
                     std::vector< pair<size_t, Element> >,
                     std::deque< pair<size_t, Element> > >
	T(F, mode, bbtimer, givtimer, *in_ptr, *out_ptr, *log_ptr);
      T.test();

    } // else if ( (vector == 4) && (value == 4) )
    else if ( (vector == 5) && (value == 1) )
    {
#ifdef TRACE
      *log_ptr	<< "Testing sparsemat with :" << endl
	        << "  vector type: sparse STL deques" << endl
		<< "  row type: sparse STL maps" << endl
		<< "  apply mode: " << mode << endl
		<< "  timers: ";
      if (bbtimer) *log_ptr << "blackbox timer, ";
      if (givtimer) *log_ptr << "Givaro timer, ";
      *out_ptr << endl;
#endif TRACE

      test_sparsemat<Field,
                     std::deque< pair<size_t, Element> >,
                     std::map<size_t, Element> >
	T(F, mode, bbtimer, givtimer, *in_ptr, *out_ptr, *log_ptr);
      T.test();

    } // else if ( (vector == 5) && (value == 1) )
    else if ( (vector == 5) && (value == 2) )
    {
#ifdef TRACE
      *log_ptr	<< "Testing sparsemat with :" << endl
	        << "  vector type: sparse STL deques" << endl
		<< "  row type: sparse STL lists" << endl
		<< "  apply mode: " << mode << endl
		<< "  timers: ";
      if (bbtimer) *log_ptr << "blackbox timer, ";
      if (givtimer) *log_ptr << "Givaro timer, ";
      *out_ptr << endl;
#endif TRACE

      test_sparsemat<Field,
                     std::deque< pair<size_t, Element> >,
                     std::list< pair<size_t, Element> > >
	T(F, mode, bbtimer, givtimer, *in_ptr, *out_ptr, *log_ptr);
      T.test();

    } // else if ( (vector == 5) && (value == 2) )
    else if ( (vector == 5) && (value == 3) )
    {
#ifdef TRACE
      *log_ptr	<< "Testing sparsemat with :" << endl
	        << "  vector type: sparse STL deques" << endl
		<< "  row type: sparse STL vectors" << endl
		<< "  apply mode: " << mode << endl
		<< "  timers: ";
      if (bbtimer) *log_ptr << "blackbox timer, ";
      if (givtimer) *log_ptr << "Givaro timer, ";
      *out_ptr << endl;
#endif TRACE

      test_sparsemat<Field,
                     std::deque< pair<size_t, Element> >,
                     std::vector< pair<size_t, Element> > >
	T(F, mode, bbtimer, givtimer, *in_ptr, *out_ptr, *log_ptr);
      T.test();

    } // else if ( (vector == 5) && (value == 3) )
    else if ( (vector == 5) && (value == 4) )
    {
#ifdef TRACE
      *log_ptr	<< "Testing sparsemat with :" << endl
	        << "  vector type: sparse STL deques" << endl
		<< "  row type: sparse STL deques" << endl
		<< "  apply mode: " << mode << endl
		<< "  timers: ";
      if (bbtimer) *log_ptr << "blackbox timer, ";
      if (givtimer) *log_ptr << "Givaro timer, ";
      *out_ptr << endl;
#endif TRACE

      test_sparsemat<Field,
                     std::deque< pair<size_t, Element> >,
                     std::deque< pair<size_t, Element> > >
	T(F, mode, bbtimer, givtimer, *in_ptr, *out_ptr, *log_ptr);
      T.test();

    } // else if ( (vector == 5) && (value == 4) )

  } // if (sparsemat)

  if (multiply)
  {
    if (vector == 1)
    {
#ifdef TRACE
      *log_ptr	<< "Testing multiplication of two sparsemat matrices with:" 
	        << endl
	        << "  vector type: dense STL vectors" << endl
		<< "  row type: sparse STL maps" << endl
		<< "  apply mode: " << mode << endl
		<< "  timers: ";
      if (bbtimer) *log_ptr << "blackbox timer, ";
      if (givtimer) *log_ptr << "Givaro timer, ";
      *out_ptr << endl;
#endif TRACE

      test_multiply<Field, std::vector<Element> >
	T(F, mode, bbtimer, givtimer, *in_ptr, *out_ptr, *log_ptr);
      T.test();

    } // if (vector == 1)
    else if (vector == 2)
    {
#ifdef TRACE
      *log_ptr	<< "Testing multiplication of two sparsemat matrices with:" 
	        << endl
	        << "  vector type: sparse STL maps" << endl
		<< "  row type: sparse STL maps" << endl
		<< "  apply mode: " << mode << endl
		<< "  timers: ";
      if (bbtimer) *log_ptr << "blackbox timer, ";
      if (givtimer) *log_ptr << "Givaro timer, ";
      *out_ptr << endl;
#endif TRACE

      test_multiply<Field, std::map<size_t, Element> >
	T(F, mode, bbtimer, givtimer, *in_ptr, *out_ptr, *log_ptr);
      T.test();

    } // else if (vector == 2)
    else if (vector == 3)
    {
#ifdef TRACE
      *log_ptr	<< "Testing multiplication of two sparsemat matrices with:" 
	        << endl
	        << "  vector type: sparse STL lists" << endl
		<< "  row type: sparse STL maps" << endl
		<< "  apply mode: " << mode << endl
		<< "  timers: ";
      if (bbtimer) *log_ptr << "blackbox timer, ";
      if (givtimer) *log_ptr << "Givaro timer, ";
      *out_ptr << endl;
#endif TRACE

      test_multiply<Field, std::list< std::pair<size_t, Element> > >
	T(F, mode, bbtimer, givtimer, *in_ptr, *out_ptr, *log_ptr);
      T.test();
      
    } // else if (vector == 3)
    else if (vector == 4)
    {
#ifdef TRACE
      *log_ptr	<< "Testing multiplication of two sparsemat matrices with:" 
	        << endl
	        << "  vector type: sparse STL vectors" << endl
		<< "  row type: sparse STL maps" << endl
		<< "  apply mode: " << mode << endl
		<< "  timers: ";
      if (bbtimer) *log_ptr << "blackbox timer, ";
      if (givtimer) *log_ptr << "Givaro timer, ";
      *out_ptr << endl;
#endif TRACE

      test_multiply<Field, std::vector< std::pair<size_t, Element> > >
	T(F, mode, bbtimer, givtimer, *in_ptr, *out_ptr, *log_ptr);
      T.test();
      
    } // else if (vector == 4)
    else if (vector == 5)
    {
#ifdef TRACE
      *log_ptr	<< "Testing multiplication of two sparsemat matrices with:" 
	        << endl
	        << "  vector type: sparse STL deques" << endl
		<< "  row type: sparse STL maps" << endl
		<< "  apply mode: " << mode << endl
		<< "  timers: ";
      if (bbtimer) *log_ptr << "blackbox timer, ";
      if (givtimer) *log_ptr << "Givaro timer, ";
      *out_ptr << endl;
#endif TRACE

      test_multiply<Field, std::deque< std::pair<size_t, Element> > >
	T(F, mode, bbtimer, givtimer, *in_ptr, *out_ptr, *log_ptr);
      T.test();
      
    } // else if (vector == 5)
       
  } // if (multiply)

  if (hilbert) 
  {

    if (vector == 1)
    {
#ifdef TRACE
      *log_ptr	<< "Testing hilbert with :" << endl
	        << "  vector type: dense STL vectors" << endl
		<< "  apply mode: " << mode << endl
		<< "  timers: ";
      if (bbtimer) *log_ptr << "blackbox timer, ";
      if (givtimer) *log_ptr << "Givaro timer, ";
      *out_ptr << endl;
#endif TRACE

      test_hilbert<Field, std::vector<Element> >
	T(F, mode, bbtimer, givtimer, *in_ptr, *out_ptr, *log_ptr);
      T.test();

    } // if (vector == 1)
    else if (vector == 2)
    {
#ifdef TRACE
      *log_ptr	<< "Testing hilbert with :" << endl
	        << "  vector type: sparse STL maps" << endl
		<< "  apply mode: " << mode << endl
		<< "  timers: ";
      if (bbtimer) *log_ptr << "blackbox timer, ";
      if (givtimer) *log_ptr << "Givaro timer, ";
      *out_ptr << endl;
#endif TRACE

      test_hilbert<Field, std::map<size_t, Element> >
	T(F, mode, bbtimer, givtimer, *in_ptr, *out_ptr, *log_ptr);
      T.test();

    } // else if (vector == 2)
    else if (vector == 3)
    {
#ifdef TRACE
      *log_ptr	<< "Testing hilbert with :" << endl
	        << "  vector type: sparse STL lists" << endl
		<< "  apply mode: " << mode << endl
		<< "  timers: ";
      if (bbtimer) *log_ptr << "blackbox timer, ";
      if (givtimer) *log_ptr << "Givaro timer, ";
      *out_ptr << endl;
#endif TRACE

      test_hilbert<Field, std::list< pair<size_t, Element> > >
	T(F, mode, bbtimer, givtimer, *in_ptr, *out_ptr, *log_ptr);
      T.test();

    } // else if (vector == 3)
    else if (vector == 4)
    {
#ifdef TRACE
      *log_ptr	<< "Testing hilbert with :" << endl
	        << "  vector type: sparse STL vectors" << endl
		<< "  apply mode: " << mode << endl
		<< "  timers: ";
      if (bbtimer) *log_ptr << "blackbox timer, ";
      if (givtimer) *log_ptr << "Givaro timer, ";
      *out_ptr << endl;
#endif TRACE

      test_hilbert<Field, std::vector< pair<size_t, Element> > >
	T(F, mode, bbtimer, givtimer, *in_ptr, *out_ptr, *log_ptr);
      T.test();

    } // else if (vector == 4)
    else if (vector == 5)
    {
#ifdef TRACE
      *log_ptr	<< "Testing hilbert with :" << endl
	        << "  vector type: sparse STL deques" << endl
		<< "  apply mode: " << mode << endl
		<< "  timers: ";
      if (bbtimer) *log_ptr << "blackbox timer, ";
      if (givtimer) *log_ptr << "Givaro timer, ";
      *out_ptr << endl;
#endif TRACE

      test_hilbert<Field, std::deque< pair<size_t, Element> > >
	T(F, mode, bbtimer, givtimer, *in_ptr, *out_ptr, *log_ptr);
      T.test();

    } // else if (vector == 5)
 
  } // if (hilbert)
 
  if (butterfly) 
  {

    if (vector == 1)
    {
#ifdef TRACE
      *log_ptr	<< "Testing butterfly switching network matrix with :" << endl
	        << "  vector type: dense STL vectors" << endl
		<< "  apply mode: " << mode << endl
		<< "  timers: ";
      if (bbtimer) *log_ptr << "blackbox timer, ";
      if (givtimer) *log_ptr << "Givaro timer, ";
      *out_ptr << endl;
#endif TRACE

      test_butterfly<Field, std::vector<Element> >
	T(F, mode, bbtimer, givtimer, *in_ptr, *out_ptr, *log_ptr);
      T.test();

    } // if (vector == 1)
    else if (vector == 2)
    {
#ifdef TRACE
      *log_ptr	<< "Testing butterfly switching network matrix with :" << endl
	        << "  vector type: sparse STL maps" << endl
		<< "  apply mode: " << mode << endl
		<< "  timers: ";
      if (bbtimer) *log_ptr << "blackbox timer, ";
      if (givtimer) *log_ptr << "Givaro timer, ";
      *out_ptr << endl;
#endif TRACE

      *out_ptr << "Butterfly switching network blackbox matrix is only "
	       << "implemented for dense vectors." << endl;

    } // else if (vector == 2)
    else if (vector == 3)
    {
#ifdef TRACE
      *log_ptr	<< "Testing butterfly switching network matrix with :" << endl
	        << "  vector type: sparse STL lists" << endl
		<< "  apply mode: " << mode << endl
		<< "  timers: ";
      if (bbtimer) *log_ptr << "blackbox timer, ";
      if (givtimer) *log_ptr << "Givaro timer, ";
      *out_ptr << endl;
#endif TRACE

      *out_ptr << "Butterfly switching network blackbox matrix is only "
	       << "implemented for dense vectors." << endl;

    } // else if (vector == 3)
    else if (vector == 4)
    {
#ifdef TRACE
      *log_ptr	<< "Testing butterfly switching network matrix with :" << endl
	        << "  vector type: sparse STL vectors" << endl
		<< "  apply mode: " << mode << endl
		<< "  timers: ";
      if (bbtimer) *log_ptr << "blackbox timer, ";
      if (givtimer) *log_ptr << "Givaro timer, ";
      *out_ptr << endl;
#endif TRACE

      *out_ptr << "Butterfly switching network blackbox matrix is only "
	       << "implemented for dense vectors." << endl;

    } // else if (vector == 4)
    else if (vector == 5)
    {
#ifdef TRACE
      *log_ptr	<< "Testing butterfly switching network matrix with :" << endl
	        << "  vector type: sparse STL deques" << endl
		<< "  apply mode: " << mode << endl
		<< "  timers: ";
      if (bbtimer) *log_ptr << "blackbox timer, ";
      if (givtimer) *log_ptr << "Givaro timer, ";
      *out_ptr << endl;
#endif TRACE

      *out_ptr << "Butterfly switching network blackbox matrix is only "
	       << "implemented for dense vectors." << endl;

    } // else if (vector == 5)
 
  } // if (butterfly)

  return true;

} // template <class Field> bool test_linbox::run_tests(...)

#endif // _TEST_LINBOX_
