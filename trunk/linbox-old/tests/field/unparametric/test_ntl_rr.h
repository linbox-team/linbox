/* File: src/examples/field/unparametric/test_ntl_rr.h
 * Author: William J Turner for the LinBox group
 */

#ifndef _TEST_NTL_RR_
#define _TEST_NTL_RR_

#include "Examples/test_linbox.h"
#include "LinBox/unparam_field.h"
#include "LinBox/ntl.h"

// Specialization of setup_field for unparam<NTL::RR>
template <> 
bool test_linbox::test<test_linbox::field_categories::ntl_rr_tag>(void) const
{
  LinBox::unparam_field<NTL::RR> F;
  return run_tests(F);
} // template <> bool test_linbox<wntl_rr_tag>(void)

#endif // _TEST_NTL_RR_
