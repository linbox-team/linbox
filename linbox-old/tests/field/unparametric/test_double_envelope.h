/* File: src/examples/field/unparametric/test_double_envelope.h
 * Author: William J Turner for the LinBox group
 */

#ifndef _TEST_DOUBLE_ENVELOPE_
#define _TEST_DOUBLE_ENVELOPE_

#include "../../test_linbox.h"
#include "LinBox/unparam_field.h"
#include "LinBox/field_archetype.h"
#include "LinBox/field_envelope.h"
#include "LinBox/randiter_envelope.h"
#include "LinBox/randiter_archetype.h"

// Specialization of setup_field for double_envelope
template <> 
bool test_linbox::test<test_linbox::field_categories::double_envelope_tag>(void) const
{
  LinBox::unparam_field<double> F;
  LinBox::Field_archetype A(&F);
  return run_tests(A);
} // template <> bool test_linbox<double_envelope_tag>(void)

#endif // _TEST_DOUBLE_ENVELOPE_
