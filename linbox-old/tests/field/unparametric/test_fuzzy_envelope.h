/* File: src/examples/field/unparametric/test_fuzzy.h
 * Author: William J Turner for the LinBox group
 */

#ifndef _TEST_FUZZY_ENVELOPE_
#define _TEST_FUZZY_ENVELOPE_

#include <iostream>
#include "../../test_linbox.h"
#include "LinBox/unparam_field.h"
#include "LinBox/fuzzy.h"
#include "LinBox/field_archetype.h"
#include "LinBox/field_envelope.h"

// Specialization of setup_field for fuzzy_envelope
template <> 
bool test_linbox::test<test_linbox::field_categories::fuzzy_envelope_tag>
(void) const
{
  double fuzz;
  if (prompt) cout << endl << "Enter a fuzz value: ";
  *in_ptr >> fuzz;
  LinBox::fuzzy::put_fuzz(fuzz);
  LinBox::unparam_field<LinBox::fuzzy> F;
  LinBox::Field_archetype A(&F);
  return run_tests(A);
} // template <> bool test_linbox<fuzzy_tag>(void)

#endif // _TEST_FUZZY_ENVELOPE_
