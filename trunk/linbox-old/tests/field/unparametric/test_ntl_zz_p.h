/* File: src/examples/field/unparametric/test_ntl_zz_p.h
 * Author: William J Turner for the LinBox group
 */

#ifndef _TEST_NTL_zz_P_
#define _TEST_NTL_zz_P_

#include <iostream>
#include "Examples/test_linbox.h"
#include "LinBox/unparam_field.h"
#include "LinBox/ntl.h"

// Specialization of setup_field for unparam<ntl_zz_p>
template <> 
bool test_linbox::test<test_linbox::field_categories::ntl_zz_p_tag>(void) const
{
  long modulus; // prime modulus for the mathematical field
  if (prompt)
    cout << endl << "Enter a prime integer for the modulus of the field: ";
  *in_ptr >> modulus;
  NTL::zz_p::init(modulus);
  LinBox::unparam_field<NTL::zz_p> F;
  return run_tests(F);
} // template <> bool test_linbox<ntl_zz_p_tag>(void)

#endif // _TEST_NTL_zz_P_
