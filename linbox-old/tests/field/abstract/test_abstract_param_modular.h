/* File: src/examples/field/unparametric/test_abstract_param_modular.h
 * Author: William J Turner for the LinBox group
 */

#ifndef _TEST_ABSTRACT_PARAM_MODULAR_
#define _TEST_ABSTRACT_PARAM_MODULAR_

#include "Examples/test_linbox.h"
#include "LinBox/field_archetype.h"
#include "LinBox/abstract_param_modular.h"

// Specialization of setup_field for abstract_param_modular
template <> 
bool test_linbox::test<test_linbox::field_categories::abstract_param_modular_tag>(void) const
{
  long modulus; // prime modulus
  if (prompt)
    cout << endl << "Enter a prime number for the modulus of the field: ";
  *in_ptr >> modulus;
  LinBox::abstract_param_modular F(modulus);
  LinBox::abstract_param_modular::element e;
  LinBox::abstract_param_modular::randIter r(F);
  LinBox::Field_archetype A(&F, &e, &r);
  return run_tests(A);
} // template <> bool test_linbox<param_modular_tag>(void)

#endif // _TEST_ABSTRACT_PARAM_MODULAR_
