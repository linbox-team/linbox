/* File: src/examples/field/unparametric/test_ntl_ZZ_p_envelope.h
 * Author: William J Turner for the LinBox group
 */

#ifndef _TEST_NTL_ZZ_P_ENVELOPE_
#define _TEST_NTL_ZZ_P_ENVELOPE_

#include <iostream>
#include "../../test_linbox.h"
#include "LinBox/unparam_field.h"
#include "LinBox/ntl.h"
#include "LinBox/field_archetype.h"
#include "LinBox/field_envelope.h"
#include "LinBox/randiter_envelope.h"
#include "LinBox/randiter_archetype.h"

// Specialization of setup_field for unparam<ntl_ZZ_p_envelope>
template <> 
bool test_linbox::test<test_linbox::field_categories::ntl_ZZ_p_envelope_tag>(void) const
{
  long modulus; // prime modulus for the mathematical field
  if (prompt)
    cout << endl << "Enter a prime integer for the modulus of the field: ";
  *in_ptr >> modulus;
  NTL::ZZ_p::init(NTL::to_ZZ(modulus));
  LinBox::unparam_field<NTL::ZZ_p> F;
  LinBox::Field_archetype A(&F);
  return run_tests(A);
} // template <> bool test_linbox<ntl_ZZ_p_envelope_tag>(void)

#endif // _TEST_NTL_ZZ_P_ENVELOPE_
