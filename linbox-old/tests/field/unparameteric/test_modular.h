/* File: src/examples/field/unparameteric/test_modular.h
 * Author: William J Turner for the LinBox group
 */

#ifndef _TEST_MODULAR_
#define _TEST_MODULAR_

#include <iostream>
#include "Examples/test_linbox.h"
#include "LinBox/unparam_field.h"
#include "LinBox/modular.h"

// Specialization of setup_field for unparam<modular>
template <> 
bool test_linbox::test<test_linbox::field_categories::modular_tag>(void) const
{
  long modulus; // prime modulus for the mathematical field
  if (prompt)
    cout << endl << "Enter a prime number for the modulus of the field: ";
  *in_ptr >> modulus;
  LinBox::modular::put_modulus(modulus);
  LinBox::unparam_field<LinBox::modular> F;
  return run_tests(F);
} // template <> bool test_linbox<modular_tag>(void)

#endif // _TEST_MODULAR_
