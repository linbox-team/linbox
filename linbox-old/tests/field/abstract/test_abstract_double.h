/* File: src/examples/field/unparametric/test_double.h
 * Author: William J Turner for the LinBox group
 */

#ifndef _TEST_ABSTRACT_DOUBLE_
#define _TEST_ABSTRACT_DOUBLE_

#include "../../test_linbox.h"
#include "LinBox/unparam_field.h"
#include "LinBox/field_archetype.h"
#include "LinBox/abstract_double.h"

// Specialization of setup_field for abstractdouble
template <> 
bool test_linbox::test<test_linbox::field_categories::abstract_double_tag>(void) const
{
  LinBox::abstract_double F;
  LinBox::abstract_double::element e;
  LinBox::abstract_double::randIter r(F);
  LinBox::Field_archetype A(&F, &e, &r);
  return run_tests(A);
} // template <> bool test_linbox<double_tag>(void)

#endif // _TEST_ABSTRACT_DOUBLE_
