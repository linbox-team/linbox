/* File: src/examples/field/unparametric/test_fuzzy.h
 * Author: William J Turner for the LinBox group
 */

#ifndef _TEST_ABSTRACT_PARAM_FUZZY_
#define _TEST_ABSTRACT_PARAM_FUZZY_

#include "Examples/test_linbox.h"
#include "LinBox/unparam_field.h"
#include "LinBox/field_archetype.h"
#include "LinBox/abstract_param_fuzzy.h"

// Specialization of setup_field for abstract_param_fuzzy
template <> 
bool test_linbox::test<test_linbox::field_categories::abstract_param_fuzzy_tag>(void) const
{
  double fuzz;
  if (prompt)
    cout << endl << "Enter a fuzz value: ";
  *in_ptr >> fuzz;
  LinBox::abstract_param_fuzzy F(fuzz);
  LinBox::abstract_param_fuzzy::element e;
  LinBox::abstract_param_fuzzy::randIter r(F);
  LinBox::Field_archetype A(&F, &e, &r);
  return run_tests(A);
} // template <> bool test_linbox<fuzzy_tag>(void)

#endif // _TEST_ABSTRACT_PARAM_FUZZY_
