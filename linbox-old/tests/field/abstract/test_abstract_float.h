/* File: src/examples/field/unparametric/test_float.h
 * Author: William J Turner for the LinBox group
 */

#ifndef _TEST_ABSTRACT_FLOAT_
#define _TEST_ABSTRACT_FLOAT_

#include "../../test_linbox.h"
#include "LinBox/unparam_field.h"
#include "LinBox/field_archetype.h"
#include "LinBox/abstract_float.h"

// Specialization of setup_field for abstract_float
template <> 
bool test_linbox::test<test_linbox::field_categories::abstract_float_tag>(void) const
{
  LinBox::abstract_float F;
  LinBox::abstract_float::element e;
  LinBox::abstract_float::randIter r(F);
  LinBox::Field_archetype A(&F, &e, &r);
  return run_tests(A);
} // template <> bool test_linbox<float_tag>(void)

#endif // _TEST_ABSTRACT_FLOAT_
