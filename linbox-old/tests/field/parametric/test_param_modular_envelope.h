/* File: src/examples/field/unparametric/test_param_modular_envelope.h
 * Author: William J Turner for the LinBox group
 */

#ifndef _TEST_PARAM_MODULAR_ENVELOPE_
#define _TEST_PARAM_MODULAR_ENVELOPE_

#include <iostream>
#include "Examples/test_linbox.h"
#include "LinBox/param_modular.h"
#include "LinBox/field_archetype.h"
#include "LinBox/field_envelope.h"

// Specialization of setup_field for modular_envelope
template <> 
bool test_linbox::test<test_linbox::field_categories::param_modular_envelope_tag>
(void) const
{
  long modulus; // prime modulus for the mathematical field
  if (prompt)
    cout << endl << "Enter a prime number for the modulus of the field: ";
  *in_ptr >> modulus;
  LinBox::param_modular F(modulus);
  LinBox::Field_envelope< LinBox::param_modular> E(F);
  LinBox::Field_envelope< LinBox::param_modular>::element e;
  LinBox::Field_envelope< LinBox::param_modular>::randIter r(E);
  LinBox::Field_archetype A(&E, &e, &r);
  return run_tests(A);
} // template <> bool test_linbox<param_modular_envelope_tag>(void)

#endif // _TEST_PARAM_MODULAR_ENVELOPE_
