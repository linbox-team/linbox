/* File: src/examples/field/unparametric/test_fuzzy.h
 * Author: William J Turner for the LinBox group
 */

#ifndef _TEST_FUZZY_
#define _TEST_FUZZY_

#include <iostream>
#include "../../test_linbox.h"
#include "LinBox/unparam_field.h"
#include "LinBox/fuzzy.h"

// Specialization of setup_field for unparam<fuzzy>
template <> 
bool test_linbox::test<test_linbox::field_categories::fuzzy_tag>(void) const
{
  double fuzz;
  if (prompt) cout << endl << "Enter a fuzz value: ";
  *in_ptr >> fuzz;
  LinBox::fuzzy::put_fuzz(fuzz);
  LinBox::unparam_field<LinBox::fuzzy> F;
  return run_tests(F);
} // template <> bool test_linbox<fuzzy_tag>(void)

#endif // _TEST_FUZZY_
