/* File: src/examples/field/unparametric/test_double.h
 * Author: William J Turner for the LinBox group
 */

#ifndef _TEST_DOUBLE_
#define _TEST_DOUBLE_

#include "Examples/test_linbox.h"
#include "LinBox/unparam_field.h"

// Specialization of setup_field for unparam<double>
template <> 
bool test_linbox::test<test_linbox::field_categories::double_tag>(void) const
{
  LinBox::unparam_field<double> F;
  return run_tests(F);
} // template <> bool test_linbox<double_tag>(void)

#endif // _TEST_DOUBLE_
