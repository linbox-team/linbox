/* File: src/examples/blackbox/test_field.h
 * Author: William J Turner for the LinBox group
 */

#ifndef _TEST_FIELD_
#define _TEST_FIELD_

#include "LinBox/integer.h"
#include "LinBox/faxpy.h"

#include "../test_base.h"

/** Class to test sparsemat blackbox matrix.
 * Templatized by field, vector, and row types.
 * @see sparsemat
 */
template <class Field>
class test_field : public test_base
{
public:

  /** Constructor from Field object
   * @param  F field in which arithmetic is done
   * @param  in  istream from which input is read (default = cin)
   * @param  out ostream to which output is written (default = cout)
   * @param  log ostream to which messages are logged (default = clog)
   */
  test_field(const Field& F, 
	     istream& in = cin, 
	     ostream& out = cout,
	     ostream& log = clog);

  /** Run tests on LinBox field.
   * Tests the LinBox field arithmetic.
   * @return boolean true if ended successfully, false if not
   */
  bool test(void) const;

private:

  // Field in which arithmetic is done
  Field _F;

}; // class test_field : public test_base

// Implementation of methods

template <class Field>
test_field<Field>::test_field(const Field& F,
			      istream& in = cin,
			      ostream& out = cout,
			      ostream& log = clog)
  : test_base(in, out, log), _F(F) 
{
  // prompt for input of matrix
  if (prompt)
    cout << "Enter the name of the file from which" << endl
         << "to read the field input." << endl
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

} // test_field<Field>::test_field(...)

template <class Field>
bool test_field<Field>::test(void) const
{
  *out_ptr << "Using the field ";
  _F.write(*out_ptr);
  *out_ptr << ", which has cardinality " << _F.cardinality() 
           << " and charactersitic " << _F.characteristic() << endl;

  typename Field::element temp, one, zero;
  LinBox::integer temp_int;
  
#ifdef TRACE
  *log_ptr << "Initializing temporary and constants elements" << endl
           << "and creating other elements through copy constructor." << endl;
#endif // TRACE

  _F.init(temp, temp_int);
  _F.init(one, 1);
  _F.init(zero, 0);
  typename Field::element a(temp), b(temp), c(temp);
  

  if (prompt) 
    cout << "Enter two integers to which to initialize three field elements: ";

  _F.read(*in_ptr, a);
  _F.read(*in_ptr, b);

  *out_ptr << "The two integers entered were ";
  _F.write(*out_ptr, a);
  *out_ptr << " and ";
  _F.write(*out_ptr, b);
  *out_ptr << endl;

#ifdef TRACE
  *log_ptr << "Assigning first element to third and then converting to integer."
           << endl;
#endif // TRACE

  _F.assign(c, a);
  _F.convert(temp_int, c);

  *out_ptr << "Testing assign and convert: ";
  _F.write(*out_ptr, a);
  *out_ptr << " = ";
  _F.write(*out_ptr, c);
  *out_ptr << " = " << temp_int << endl;

  *out_ptr << "Testing arithmetic functions." << endl;

  if (_F.areEqual(zero, zero))
    *out_ptr << "    0 == 0" << endl;
  else
    *out_ptr << "    0 != 0" << endl;

  if (_F.areEqual(zero, one))
    *out_ptr << "    0 == 1" << endl;
  else
    *out_ptr << "    0 != 1" << endl;

  *out_ptr << "    ";
  _F.write(*out_ptr, a);
  *out_ptr << " + ";
  _F.write(*out_ptr, b);
  *out_ptr << " = ";
  _F.write(*out_ptr, _F.add(c, a, b));
  *out_ptr << " = ";
  _F.write(*out_ptr, c);
  *out_ptr << endl;
  
  *out_ptr << "    ";
  _F.write(*out_ptr, a);
  *out_ptr << " - ";
  _F.write(*out_ptr, b);
  *out_ptr << " = ";
  _F.write(*out_ptr, _F.sub(c, a, b));
  *out_ptr << " = ";
  _F.write(*out_ptr, c);
  *out_ptr << endl;
  
  *out_ptr << "    ";
  _F.write(*out_ptr, a);
  *out_ptr << " * ";
  _F.write(*out_ptr, b);
  *out_ptr << " = ";
  _F.write(*out_ptr, _F.mul(c, a, b));
  *out_ptr << " = ";
  _F.write(*out_ptr, c);
  *out_ptr << endl;
  
  *out_ptr << "    ";
  _F.write(*out_ptr, a);
  *out_ptr << " / ";
  _F.write(*out_ptr, b);
  *out_ptr << " = ";
  _F.write(*out_ptr, _F.div(c, a, b));
  *out_ptr << " = ";
  _F.write(*out_ptr, c);
  *out_ptr << endl;
  
  *out_ptr << "    - ";
  _F.write(*out_ptr, a);
  *out_ptr << " = ";
  _F.write(*out_ptr, _F.neg(c, a));
  *out_ptr << " = ";
  _F.write(*out_ptr, c);
  *out_ptr << endl;
  
  *out_ptr << "    1 / ";
  _F.write(*out_ptr, a);
  *out_ptr << " = ";
  _F.write(*out_ptr, _F.inv(c, a));
  *out_ptr << " = ";
  _F.write(*out_ptr, c);
  *out_ptr << endl;
  
  *out_ptr << "Testing in-place arithmetic functions." << endl;

  if (_F.isZero(zero))
    *out_ptr << "    0 == 0" << endl;
  else
    *out_ptr << "    0 != 0" << endl;

  if (_F.isZero(one))
    *out_ptr << "    1 == 0" << endl;
  else
    *out_ptr << "    1 != 0" << endl;

  if (_F.isOne(one))
    *out_ptr << "    1 == 1" << endl;
  else
    *out_ptr << "    1 != 1" << endl;

  if (_F.isOne(zero))
    *out_ptr << "    0 == 1" << endl;
  else
    *out_ptr << "    0 != 1" << endl;

  _F.assign(c, a);
  *out_ptr << "    ";
  _F.write(*out_ptr, c);
  *out_ptr << " += ";
  _F.write(*out_ptr, b);
  *out_ptr << " = ";
  _F.write(*out_ptr, _F.addin(c, b));
  *out_ptr << " = ";
  _F.write(*out_ptr, c);
  *out_ptr << endl;
  
  _F.assign(c, a);
  *out_ptr << "    ";
  _F.write(*out_ptr, c);
  *out_ptr << " -= ";
  _F.write(*out_ptr, b);
  *out_ptr << " = ";
  _F.write(*out_ptr, _F.subin(c, b));
  *out_ptr << " = ";
  _F.write(*out_ptr, c);
  *out_ptr << endl;
  
  _F.assign(c, a);
  *out_ptr << "    ";
  _F.write(*out_ptr, c);
  *out_ptr << " *= ";
  _F.write(*out_ptr, b);
  *out_ptr << " = ";
  _F.write(*out_ptr, _F.mulin(c, b));
  *out_ptr << " = ";
  _F.write(*out_ptr, c);
  *out_ptr << endl;
  
  _F.assign(c, a);
  *out_ptr << "    ";
  _F.write(*out_ptr, c);
  *out_ptr << " /= ";
  _F.write(*out_ptr, b);
  *out_ptr << " = ";
  _F.write(*out_ptr, _F.divin(c, b));
  *out_ptr << " = ";
  _F.write(*out_ptr, c);
  *out_ptr << endl;
  
  _F.assign(c, a);
  *out_ptr << "    - ";
  _F.write(*out_ptr, c);
  *out_ptr << " = ";
  _F.write(*out_ptr, _F.negin(c));
  *out_ptr << " = ";
  _F.write(*out_ptr, c);
  *out_ptr << endl;
 
  _F.assign(c, a);
  *out_ptr << "    1 / ";
  _F.write(*out_ptr, c);
  *out_ptr << " = ";
  _F.write(*out_ptr, _F.invin(c));
  *out_ptr << " = ";
  _F.write(*out_ptr, c);
  *out_ptr << endl;

  *out_ptr << "Testing faxpy." << endl;

  _F.init(temp, 2);
  LinBox::faxpy<Field> Faxpy(_F, temp);

  *out_ptr << "    ";
  _F.write(*out_ptr, temp);
  *out_ptr << " * ";
  _F.write(*out_ptr, a);
  *out_ptr << " + ";
  _F.write(*out_ptr, b);
  *out_ptr << " = ";
  _F.write(*out_ptr, Faxpy.apply(c, a, b));
  *out_ptr << " = ";
  _F.write(*out_ptr, c);
  *out_ptr << endl;
 
  _F.assign(c, b);
  *out_ptr << "    ";
  _F.write(*out_ptr, c);
  *out_ptr << " += ";
  _F.write(*out_ptr, temp);
  *out_ptr << " * ";
  _F.write(*out_ptr, a);
  *out_ptr << " = ";
  _F.write(*out_ptr, Faxpy.applyin(c, a));
  *out_ptr << " = ";
  _F.write(*out_ptr, c);
  *out_ptr << endl;

  _F.init(temp, 5);
  Faxpy.assign(temp);

  *out_ptr << "    ";
  _F.write(*out_ptr, temp);
  *out_ptr << " * ";
  _F.write(*out_ptr, a);
  *out_ptr << " + ";
  _F.write(*out_ptr, b);
  *out_ptr << " = ";
  _F.write(*out_ptr, Faxpy.apply(c, a, b));
  *out_ptr << " = ";
  _F.write(*out_ptr, c);
  *out_ptr << endl;
 
  _F.assign(c, b);
  *out_ptr << "    ";
  _F.write(*out_ptr, c);
  *out_ptr << " += ";
  _F.write(*out_ptr, temp);
  *out_ptr << " * ";
  _F.write(*out_ptr, a);
  *out_ptr << " = ";
  _F.write(*out_ptr, Faxpy.applyin(c, a));
  *out_ptr << " = ";
  _F.write(*out_ptr, c);
  *out_ptr << endl;
 
  return true;

} // bool test_field<Field>::test(void) const

#endif // _TEST_FIELD_
