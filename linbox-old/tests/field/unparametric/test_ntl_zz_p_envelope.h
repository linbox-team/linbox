/* File: src/examples/field/unparametric/test_ntl_zz_p_envelope.h
 * Author: William J Turner for the LinBox group
 */

#ifndef _TEST_NTL_zz_P_ENVELOPE_
#define _TEST_NTL_zz_P_ENVELOPE_

#include <iostream>
#include "Examples/test_linbox.h"
#include "LinBox/unparam_field.h"
#include "LinBox/ntl.h"
#include "LinBox/field_archetype.h"
#include "LinBox/field_envelope.h"
#include "LinBox/randiter_envelope.h"
#include "LinBox/randiter_archetype.h"

// Specialization of setup_field for unparam<ntl_zz_p_envelope>
template <> 
bool test_linbox::test<test_linbox::field_categories::ntl_zz_p_envelope_tag>(void) const
{
  long modulus; // prime modulus for the mathematical field
  if (prompt)
    cout << endl << "Enter a prime integer for the modulus of the field: ";
  *in_ptr >> modulus;
  NTL::zz_p::init(modulus);
  LinBox::unparam_field<NTL::zz_p> F;
  LinBox::Field_envelope< LinBox::unparam_field<NTL::zz_p> > E(F);
  LinBox::Field_envelope< LinBox::unparam_field<NTL::zz_p> >::element e;
  LinBox::Field_envelope< LinBox::unparam_field<NTL::zz_p> >::randIter r(E);
  LinBox::Field_archetype A(&E, &e, &r);
  return run_tests(A);
} // template <> bool test_linbox<ntl_zz_p_envelope_tag>(void)

#endif // _TEST_NTL_zz_P_ENVELOPE_
