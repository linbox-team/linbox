/* File: src/examples/blackbox/test_butterfly.h
 * Author: William J Turner for the LinBox group
 */

#ifndef _TEST_BUTTERFLY_
#define _TEST_BUTTERFLY_

#include "LinBox/butterfly.h"
#include "LinBox/boolean_switch.h"
#include "LinBox/cekstv_switch.h"

#include "Examples/vector_utility.h"

/** Class to test butterfly switching network blackbox matrix.
 * Templatized by field and (dense) vector types.
 * @see butterfly
 */
template <class Field, class Vector>
class test_butterfly : public test_base, public vector_utility<Field, Vector>
{
public:

  /** Constructor from Field object
   * @param  F field in which arithmetic is done
   * @param  mode blackbox matrix apply mode
   * @param  bbtime boolean flag on whether to use blackbox timer
   * @param  givtime boolean flag on whether to use Givaro timer
   * @param  in  istream from which input is read (default = cin)
   * @param  out ostream to which output is written (default = cout)
   * @param  log ostream to which messages are logged (default = clog)
   */
  test_butterfly(const Field& F,
		 int mode,
		 bool bbtime,
		 bool givtime,
		 istream& in = cin,
		 ostream& out = cout,
		 ostream& log = clog);

  /** Run tests on LinBox butterfly switching network blackbox matrices.
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

}; // class test_butterfly : public test_base

// Implementation of methods

template <class Field, class Vector>
test_butterfly<Field, Vector>::test_butterfly(const Field& F,
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
         << "to read the Butterfly input." << endl
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

} // test_butterfly<Field, Vector>::test_butterfly(...)

template <class Field, class Vector>
bool test_butterfly<Field, Vector>::test(void) const
{
  if (prompt) cout << "Enter the size of the vector to precondition: ";

  size_t value;
  *in_ptr >> value;

  std::vector<bool> x(value, false);

  if (prompt)
    cout << "Enter the numbers corresponding to the linearly independent rows" 
         << endl
	 << "with the rows numbered from 0 to " << value - 1 << "." << endl
	 << "End with a row number of '-1'." << endl;

  while (*in_ptr >> value)
  {
    if (value == size_t(-1)) break;
    if (value < x.size()) x[value] = true;
  } // while (*in_ptr >> value)

  *out_ptr << "The test vector is" << endl
           << "  i    x[i]" << endl
	   << "  ---------" << endl;
  for (size_t i = 0; i < x.size(); i++)
    *out_ptr << "  " << i << "    " << x[i] << endl;

  if (prompt)
    cout << "Enter non-negative offset for contiguous block: ";

  *in_ptr >> value;

#ifdef TRACE
  *log_ptr << "The offset is " << value << endl;
#endif // TRACE

  std::vector<bool> y1 = LinBox::set_butterfly(x, value, *log_ptr);

  // Create vector to switch.
  typedef typename Field::element Element;
  Element zero, one, elem;
  _F.init(zero, 0);
  _F.init(one, 1);
  _F.init(elem, 0);
  
  Vector x1;

  for (std::vector<bool>::iterator iter = x.begin(); iter != x.end(); iter++)
    x1.push_back( (*iter) ? one : zero );

  *out_ptr << "before the switch the vector is" << endl
           << "  i    x[i]" << endl
	   << "  ---------" << endl;
  for (size_t i = 0; i < x1.size(); i++)
  {
    *out_ptr << "  " << i << "    ";
    _F.write(*out_ptr, x1[i]);
    *out_ptr << endl;
  } // for (size_t i = 0; i < x1.size(); i++)

  LinBox::boolean_switch switch1(y1);

  LinBox::butterfly < Vector, LinBox::boolean_switch > 
    butterfly1(x1.size(), switch1);

  x1 = butterfly1.apply(x1);
  
  *out_ptr << "after the switch the vector is" << endl
           << "  i    x[i]" << endl
	   << "  ---------" << endl;
  for (size_t i = 0; i < x1.size(); i++)
  {
    *out_ptr << "  " << i << "    ";
    _F.write(*out_ptr, x1[i]);
    *out_ptr << endl;
  } // for (size_t i = 0; i < x1.size(); i++)

  *out_ptr << "Now testing switches set by hand..." << endl;
    
  vector<Element> xx;

  for (size_t i = 0; i < x.size(); i++, _F.addin(elem, one))
    xx.push_back(elem);

  Vector x2(xx);
  
  // break inputs into groups of size powers of 2.
  // calculate size of groups, and powers of 2 that give sizes
  // store these values in vectors n and l, respectively
  value = xx.size();
  vector<size_t> l, n;
  for (size_t l_p(0), n_p(1); n_p != 0; value >>= 1, l_p++, n_p <<= 1)
  {
#ifdef TRACE_LOOP
    *log_ptr 
      << "  looping at value = " << value
      << ", l_p = " << l_p 
      << ", n_p = " << n_p << endl;
#endif // TRACE_LOOP
    
    if (value & 1)
    {
      l.push_back(l_p);
      n.push_back(n_p);      
#ifdef TRACE_LOOP
      *log_ptr 
	<< "    inserted value = " << value 
	<< ", l_p = " << l_p 
	<< ", n_p = " << n_p << endl;
#endif // TRACE_LOOP
      
    } // if (value & 1)
    
  } //     for (size_t value(_n), l_p(0), n_p(1); n_p != 0; ...)
  
  // Calculate total number of switches required
  size_t s(0);
 
  for (size_t i = 0; i < n.size(); i++)
    s += n[i]*l[i]/2;
 
  if (n.size() != 0)
    for (size_t i = 0; i < n.size() - 1; i++)
      for (size_t j = 0; j <= i; j++)
	s += n[j];

  *out_ptr << "Now testing boolean switches..." << endl << endl;

  std::vector<bool> y2(s, false);

  // Prompt for switches to set
  if (prompt)
    cout << "There are a total of " << s << " switches that can be set." << endl
         << "Enter the numbers corresponding to the switches, " << endl
	 << "numbered from 0 to " << s - 1 << " you want set." << endl
	 << "End with a switch number of '-1'." << endl;

  while (*in_ptr >> value)
  {
    if (value == size_t(-1)) break;
    if (value < s) y2[value] = true;
  } // while (*in_ptr >> value)

  *out_ptr << "before the switch the vector is" << endl
           << "  i    x[i]" << endl
	   << "  ---------" << endl;
  for (size_t i = 0; i < x2.size(); i++)
  {
    *out_ptr << "  " << i << "    ";
    _F.write(*out_ptr, x2[i]);
    *out_ptr << endl;
  } // for (size_t i = 0; i < x2.size(); i++)

  LinBox::boolean_switch switch2(y2);

  LinBox::butterfly < Vector, LinBox::boolean_switch > 
    butterfly2(x2.size(), switch2);

  x2 = butterfly2.apply(x2);
  
  *out_ptr << "after the switch the vector is" << endl
           << "  i    x[i]" << endl
	   << "  ---------" << endl;
  for (size_t i = 0; i < x2.size(); i++)
  {
    *out_ptr << "  " << i << "    ";
    _F.write(*out_ptr, x2[i]);
    *out_ptr << endl;
  } // for (size_t i = 0; i < x2.size(); i++)

  x2 = butterfly2.applyTranspose(x2);

  *out_ptr << "after the switch transpose the vector is" << endl
           << "  i    x[i]" << endl
	   << "  ---------" << endl;
  for (size_t i = 0; i < x2.size(); i++)
  {
    *out_ptr << "  " << i << "    ";
    _F.write(*out_ptr, x2[i]);
    *out_ptr << endl;
  } // for (size_t i = 0; i < x2.size(); i++)

  *out_ptr << "Now testing multiplication switches..." << endl << endl;

  Vector y3(s, zero);
  Vector x3(xx);
  
  // Prompt for switches to set
  if (prompt)
    cout << "There are a total of " << s << " switches that can be set." << endl
         << "Enter the numbers corresponding to the switches, " << endl
	 << "numbered from 0 to " << s - 1 << " you want set" << endl
	 << "and then the value of the switch." << endl
	 << "End with a switch number of '-1'." << endl;

  while (*in_ptr >> value)
  {
    if (value == size_t(-1)) break;
    _F.read(*in_ptr, elem);
    if (value < s) y3[value] = elem;
  } // while (*in_ptr >> value)

  *out_ptr << "before the switch the vector is" << endl
           << "  i    x[i]" << endl
	   << "  ---------" << endl;
  for (size_t i = 0; i < x3.size(); i++)
  {
    *out_ptr << "  " << i << "    ";
    _F.write(*out_ptr, x3[i]);
    *out_ptr << endl;
  } // for (size_t i = 0; i < x3.size(); i++)

  LinBox::cekstv_switch<Field> 
    switch3(_F, y3);

  LinBox::butterfly<Vector, LinBox::cekstv_switch<Field> >
    butterfly3(x3.size(), switch3);

  x3 = butterfly3.apply(x3);
  
  *out_ptr << "after the switch the vector is" << endl
           << "  i    x[i]" << endl
	   << "  ---------" << endl;
  for (size_t i = 0; i < x3.size(); i++)
  {
    *out_ptr << "  " << i << "    ";
    _F.write(*out_ptr, x3[i]);
    *out_ptr << endl;
  } // for (size_t i = 0; i < x3.size(); i++)

  *out_ptr << "Now testing random multiplication switches." << endl;

  if (prompt)
    cout << "Enter a sampling size and seed: ";

  size_t size, seed;
  *in_ptr >> size >> seed;

#ifdef TRACE
  *log_ptr << "You entered size " << size << " and seed " << seed << endl;
#endif // TRACE
  
  Vector x4(xx);
 
  *out_ptr << "before the switch the vector is" << endl
           << "  i    x[i]" << endl
	   << "  ---------" << endl;
  for (size_t i = 0; i < x4.size(); i++)
  {
    *out_ptr << "  " << i << "    ";
    _F.write(*out_ptr, x4[i]);
    *out_ptr << endl;
  } // for (size_t i = 0; i < x4.size(); i++)

  typename Field::randIter r(_F, size, seed);
  LinBox::cekstv_switch<Field> switch4(_F, r);

  LinBox::butterfly<Vector, LinBox::cekstv_switch<Field> >
    butterfly4(x4.size(), switch4);

  x4 = butterfly4.apply(x4);
  
  *out_ptr << "after the switch the vector is" << endl
           << "  i    x[i]" << endl
	   << "  ---------" << endl;
  for (size_t i = 0; i < x4.size(); i++)
  {
    *out_ptr << "  " << i << "    ";
    _F.write(*out_ptr, x4[i]);
    *out_ptr << endl;
  } // for (size_t i = 0; i < x4.size(); i++)

  return true;

} // bool test_butterfly<Field, Vector>::test(void) const

#endif // _TEST_BUTTERFLY_
