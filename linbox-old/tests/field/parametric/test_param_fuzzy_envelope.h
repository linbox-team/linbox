/* File: src/examples/field/unparametric/test_param_fuzzy_envelope.h
 * Author: William J Turner for the LinBox group
 */

#ifndef _TEST_PARAM_FUZZY_ENVELOPE_
#define _TEST_PARAM_FUZZY_ENVELOPE_

#include <iostream>
#include "Examples/test_linbox.h"
#include "LinBox/param_fuzzy.h"
#include "LinBox/field_archetype.h"
#include "LinBox/field_envelope.h"

// Specialization of setup_field for fuzzy_envelope
template <> 
bool test_linbox::test<test_linbox::field_categories::param_fuzzy_envelope_tag>
(void) const
{
  double fuzz;
  if (prompt)
    cout << endl << "Enter a fuzz value: ";
  *in_ptr >> fuzz;
  LinBox::param_fuzzy F(fuzz);
  LinBox::Field_envelope< LinBox::param_fuzzy> E(F);
  LinBox::Field_envelope< LinBox::param_fuzzy>::element e;
  LinBox::Field_envelope< LinBox::param_fuzzy>::randIter r(E);
  LinBox::Field_archetype A(&E, &e, &r);
  return run_tests(A);
} // template <> bool test_linbox<param_fuzzy_envelope_tag>(void)

#endif // _TEST_PARAM_FUZZY_ENVELOPE_
