/* File: src/examples/field/test.h
 * Author: William J Turner for the LinBox group
 *
 * Include file to test all fields.
 */

#ifndef _TEST_
#define _TEST_

#include <iostream>
#include "utils/fileutils.h"

#include "field/unparametric/test_double.h"
#include "field/unparametric/test_float.h"
#include "field/unparametric/test_modular.h"
#include "field/unparametric/test_fuzzy.h"
#include "field/unparametric/test_ntl_rr.h"
#include "field/unparametric/test_ntl_ZZ_p.h"
#include "field/unparametric/test_ntl_zz_p.h"
#include "field/parametric/test_param_modular.h"
#include "field/parametric/test_param_fuzzy.h"

#include "field/unparametric/test_double_envelope.h"
#include "field/unparametric/test_float_envelope.h"
#include "field/unparametric/test_modular_envelope.h"
#include "field/unparametric/test_fuzzy_envelope.h"
#include "field/unparametric/test_ntl_rr_envelope.h"
#include "field/unparametric/test_ntl_ZZ_p_envelope.h"
#include "field/unparametric/test_ntl_zz_p_envelope.h"
#include "field/parametric/test_param_modular_envelope.h"
#include "field/parametric/test_param_fuzzy_envelope.h"

#include "field/abstract/test_abstract_double.h"
#include "field/abstract/test_abstract_float.h"
#include "field/abstract/test_abstract_modular.h"
#include "field/abstract/test_abstract_fuzzy.h"
#include "field/abstract/test_abstract_param_modular.h"
#include "field/abstract/test_abstract_param_fuzzy.h"

#include "test_linbox.h"

// Specialization of test_linbox for all_tag
template <> 
bool test_linbox::test<test_linbox::field_categories::all_tag>(void) const
{
  // prompt for input
  if (prompt) 
    *out_ptr
      << "Enter the number corresponding to the field you want to use:" 
      << endl
      << "  1 : doubles" << endl
      << "  2 : floats" << endl
      << "  3 : unparametric finite field with prime modulus" << endl
      << "  4 : unparametric \"fuzzy\" doubles" << endl
      << "  5 : parametric finite field with prime modulus" << endl
      << "  6 : parametric \"fuzzy\" doubles" << endl
      << "  7 : unparametric NTL RR field" << endl
      << "  8 : unparametric NTL ZZ_p field with prime modulus" << endl
      << "  9 : unparametric NTL zz_p field with prime modulus" << endl;
    
  // retrieve input
  int field;
  *in_ptr >> field;
  
  // Print error and get new input if number is out of range
  while ( (field < 1) || (field > 9) )
  {
    *out_ptr << "Invalid response: " << field << ".  Please try again: ";
    *in_ptr >> field;
  }

  if (prompt)
  {
    *out_ptr
      << "Enter the number corresponding to the type of field you want to use:"
      << endl
      << "  1 : wrapped with unparam_field template" << endl 
      << "  2 : wrapped with field archetype envelope" << endl;
    if (field < 7)
      *out_ptr << "  3 : derived from abstract base class" << endl;
  } // if (prompt)
  
  int field_type;
  *in_ptr >> field_type;
	
  // Print error and get new input if number is out of range
  while ( (field_type < 1) || (field_type > 3) )
  {
    *out_ptr << "Invalid response: " << field_type << ".  Please try again: ";
    *in_ptr >> field_type;
  }

  if (field == 1)
  {
    if (field_type == 1) 
      return test<test_linbox::field_categories::double_tag>();
    if (field_type == 2)   
      return test<test_linbox::field_categories::double_envelope_tag>();
    if (field_type == 3)  
      return test<test_linbox::field_categories::abstract_double_tag>();
  }
  else if (field == 2)
  {
    if (field_type == 1)  
      return test<test_linbox::field_categories::float_tag>();
    if (field_type == 2)  
      return test<test_linbox::field_categories::float_envelope_tag>();
    if (field_type == 3)  
      return test<test_linbox::field_categories::abstract_float_tag>();
  }
  else if (field == 3)
  {
    if (field_type == 1)  
      return test<test_linbox::field_categories::modular_tag>();
    if (field_type == 2)  
      return test<test_linbox::field_categories::modular_envelope_tag>();
    if (field_type == 3)  
      return test<test_linbox::field_categories::abstract_modular_tag>();
  }
  else if (field == 4)
  {
    if (field_type == 1)  
      return test<test_linbox::field_categories::fuzzy_tag>();
    if (field_type == 2)  
      return test<test_linbox::field_categories::fuzzy_envelope_tag>();
    if (field_type == 3)  
      return test<test_linbox::field_categories::abstract_fuzzy_tag>();
  }
  else if (field == 5)
  {
    if (field_type == 1) 
      return test<test_linbox::field_categories::param_modular_tag>();
    if (field_type == 2) 
      return test<test_linbox::field_categories::param_modular_envelope_tag>();
    if (field_type == 3) 
      return test<test_linbox::field_categories::abstract_param_modular_tag>();
  }
  else if (field == 6)
  {
    if (field_type == 1) 
      return test<test_linbox::field_categories::param_fuzzy_tag>();
    if (field_type == 2) 
      return test<test_linbox::field_categories::param_fuzzy_envelope_tag>();
    if (field_type == 3) 
      return test<test_linbox::field_categories::abstract_param_fuzzy_tag>();
  }
  else if (field == 7)
  {
    if (field_type == 1)  
      return test<test_linbox::field_categories::ntl_rr_tag>();
    if (field_type == 2)  
      return test<test_linbox::field_categories::ntl_rr_envelope_tag>();
  }
  else if (field == 8)
  {
    if (field_type == 1)  
      return test<test_linbox::field_categories::ntl_ZZ_p_tag>();
    if (field_type == 2)  
      return test<test_linbox::field_categories::ntl_ZZ_p_envelope_tag>();
  }
  else if (field == 9)
  {
    if (field_type == 1)  
      return test<test_linbox::field_categories::ntl_zz_p_tag>();
    if (field_type == 2) 
      return test<test_linbox::field_categories::ntl_zz_p_envelope_tag>();
  }

  return false;
  
} // template <> bool test_linbox<all_tag> (...)

#endif // _TEST_

