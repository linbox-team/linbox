/* File: src/examples/field/unparametric/test_param_modular.h
 * Author: William J Turner for the LinBox group
 */

#ifndef _TEST_PARAM_MODULAR_
#define _TEST_PARAM_MODULAR_

#include "Examples/test_linbox.h"
#include "LinBox/param_modular.h"

// Specialization of setup_field for param_modular
template <> 
bool test_linbox::test<test_linbox::field_categories::param_modular_tag>(void) const
{
  long modulus; // prime modulus
  if (prompt)
    cout << endl << "Enter a prime number for the modulus of the field: ";
  *in_ptr >> modulus;
  LinBox::param_modular F(modulus);
  return run_tests(F);
} // template <> bool test_linbox<param_modular_tag>(void)

#endif // _TEST_PARAM_MODULAR_
