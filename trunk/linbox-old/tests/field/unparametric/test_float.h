/* File: src/examples/field/unparametric/test_float.h
 * Author: William J Turner for the LinBox group
 */

#ifndef _TEST_FLOAT_
#define _TEST_FLOAT_

#include "../../test_linbox.h"
#include "LinBox/unparam_field.h"

// Specialization of setup_field for unparam<float>
template <> 
bool test_linbox::test<test_linbox::field_categories::float_tag>(void) const
{
  LinBox::unparam_field<float> F;
  return run_tests(F);
} // template <> bool test_linbox<float_tag>(void)

#endif // _TEST_FLOAT_
