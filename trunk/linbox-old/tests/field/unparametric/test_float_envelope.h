/* File: src/examples/field/unparametric/test_float_envelope.h
 * Author: William J Turner for the LinBox group
 */

#ifndef _TEST_FLOAT_ENVELOPE_
#define _TEST_FLOAT_ENVELOPE_

#include "../../test_linbox.h"
#include "LinBox/unparam_field.h"
#include "LinBox/field_archetype.h"
#include "LinBox/field_envelope.h"

// Specialization of setup_field for float_envelope
template <> 
bool test_linbox::test<test_linbox::field_categories::float_envelope_tag>(void) const
{
  LinBox::unparam_field<float> F;
  LinBox::Field_archetype A(&F);
  return run_tests(A);
} // template <> bool test_linbox<float_envelope_tag>(void)

#endif // _TEST_FLOAT_ENVELOPE_
