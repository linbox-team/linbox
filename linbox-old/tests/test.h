/* File: src/examples/field/test.h
 * Author: William J Turner for the LinBox group
 *
 * Include file to test all fields.
 */

#ifndef _TEST_
#define _TEST_

#include <iostream>
#include "Examples/fileutils.h"
#include "Examples/field_categories.h"

#include "Examples/test_double.h"
#include "Examples/test_float.h"
#include "Examples/test_modular.h"
#include "Examples/test_fuzzy.h"
//#include "Examples/test_param_modular.h"
//#include "Examples/test_param_fuzzy.h"

#include "Examples/test_double_envelope.h"
#include "Examples/test_float_envelope.h"
#include "Examples/test_modular_envelope.h"
#include "Examples/test_fuzzy_envelope.h"
//#include "Examples/test_param_modular_envelope.h"
//#include "Examples/test_param_fuzzy_envelope.h"

//#include "Examples/test_abstract_double.h"
//#include "Examples/test_abstract_float.h"
//#include "Examples/test_abstract_modular.h"
//#include "Examples/test_abstract_fuzzy.h"
//#include "Examples/test_param_modular.h"
//#include "Examples/test_param_fuzzy.h"

#include "Examples/test_linbox.h"

// Specialization of test_linbox for all_tag
template <> 
bool test_linbox::test<test_linbox::field_categories::all_tag>(void) const
{
  // prompt for input
  if (prompt) 
    *out_ptr << "Enter the number corresponding to the field you want to use:" 
         << endl
 	 << "  1  : doubles" << endl
	 << "  2  : floats" << endl
	 << "  3  : unparamteric finite field with prime modulus" << endl
	 << "  4  : unparamteric \"fuzzy\" doubles" << endl
	 << "  5  : unparamteric NTL finite field with prime modulus" << endl
	 << "  6  : parameteric finite field with prime modulus" << endl
	 << "  7  : parameteric \"fuzzy\" doubles" << endl
	 << "  8  : Number one wrapped with archetype envelope" << endl
	 << "  9  : Number two wrapped with archetype envelope" << endl
	 << "  10 : Number three wrapped with archetype envelope" << endl
	 << "  11 : Number four wrapped with archetype envelope" << endl
	 << "  12 : Number five wrapped with archetype envelope" << endl
	 << "  13 : Number six wrapped with archetype envelope" << endl
	 << "  14 : Number seven wrapped with archetype envelope" << endl
	 << "  15 : Number one derived from abstract base class" << endl
	 << "  16 : Number two derived from abstract base class" << endl
	 << "  17 : Number three derived from abstract base class" << endl
	 << "  18 : Number four derived from abstract base class" << endl
	 << "  19 : Number five derived from abstract base class" << endl
	 << "  20 : Number six derived from abstract base class" << endl
	 << "  21 : Number seven derived from abstract base class" << endl;
    
  // retrieve input
  int field;
  *in_ptr >> field;
  
  // Print error and get new input if number is out of range
  while ( (field < 1) || (field > 21) )
  {
    *out_ptr << "Invalid response: " << field << ".  Please try again: ";
    *in_ptr >> field;
  }
  
  switch (field)
  {
    case 1:
      return test<field_categories::double_tag>();
    case 2:
      return test<field_categories::float_tag>();
    case 3:
      return test<field_categories::modular_tag>();
    case 4:
      return test<field_categories::fuzzy_tag>();
/*
    case 5:
      return test<field_categories::ntl_modular_tag>();
    case 6:
      return test<field_categories::param_modular_tag>();
    case 7:
      return test<field_categories::param_fuzzy_tag>();
*/
    case 8:
      return test<field_categories::double_envelope_tag>();
    case 9:
      return test<field_categories::float_envelope_tag>();
    case 10:
      return test<field_categories::modular_envelope_tag>();
    case 11:
      return test<field_categories::fuzzy_envelope_tag>();
/*    
    case 12:
      return test<field_categories::ntl_modular_envelope_tag>();
    case 13:
      return test<field_categories::param_modular_envelope_tag>();
    case 14:
      return test<field_categories::param_fuzzy_envelope_tag>();
*/
/*
    case 15:
      return test<field_categories::abstract_double_tag>();
    case 16:
      return test<field_categories::abstract_float_tag>();
    case 17:
      return test<field_categories::abstract_modular_tag>();
    case 18:
      return test<field_categories::abstract_fuzzy_tag>();
    case 19:
      return test<field_categories::abstract_ntl_modular_tag>();
    case 20:
      return test<field_categories::abstract_param_modular_tag>();
    case 21:
      return test<field_categories::abstract_param_fuzzy_tag>();
*/
    default:
      return false;
  } // switch (field)
  
} // template <> bool test_linbox<all_tag> (...)

#endif // _TEST_

