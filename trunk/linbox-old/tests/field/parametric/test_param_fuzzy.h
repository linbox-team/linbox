/* File: src/examples/field/unparametric/test_fuzzy.h
 * Author: William J Turner for the LinBox group
 */

#ifndef _TEST_PARAM_FUZZY_
#define _TEST_PARAM_FUZZY_

#include "../../test_linbox.h"
#include "LinBox/unparam_field.h"
#include "LinBox/field_archetype.h"
#include "LinBox/param_fuzzy.h"

// Specialization of setup_field for param_fuzzy
template <> 
bool test_linbox::test<test_linbox::field_categories::param_fuzzy_tag>(void) const
{
  double fuzz;
  if (prompt)
    cout << endl << "Enter a fuzz value: ";
  *in_ptr >> fuzz;
  LinBox::param_fuzzy F(fuzz);
  return run_tests(F);
} // template <> bool test_linbox<param_fuzzy_tag>(void)

#endif // _TEST_PARAM_FUZZY_
