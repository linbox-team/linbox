/* File: src/examples/field/unparametric/test_modular.h
 * Author: William J Turner for the LinBox group
 */

#ifndef _TEST_MODULAR_ENVELOPE_
#define _TEST_MODULAR_ENVELOPE_

#include <iostream>
#include "../../test_linbox.h"
#include "LinBox/unparam_field.h"
#include "LinBox/modular.h"
#include "LinBox/field_archetype.h"
#include "LinBox/field_envelope.h"

// Specialization of setup_field for modular_envelope
template <> 
bool test_linbox::test<test_linbox::field_categories::modular_envelope_tag>
(void) const
{
  long modulus; // prime modulus for the mathematical field
  if (prompt)
    cout << endl << "Enter a prime number for the modulus of the field: ";
  *in_ptr >> modulus;
  LinBox::modular::put_modulus(modulus);
  LinBox::unparam_field<LinBox::modular> F;
  LinBox::Field_envelope< LinBox::unparam_field<LinBox::modular> > E(F);
  LinBox::Field_envelope< LinBox::unparam_field<LinBox::modular> >::element e;
  LinBox::Field_envelope< LinBox::unparam_field<LinBox::modular> >::randIter 
    r(E);
  LinBox::Field_archetype A(&E, &e, &r);
  return run_tests(A);
} // template <> bool test_linbox<modular_tag>(void)

#endif // _TEST_MODULAR_ENVELOPE_
