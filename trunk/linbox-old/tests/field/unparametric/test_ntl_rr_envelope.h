/* File: src/examples/field/unparametric/test_ntl_rr_envelope.h
 * Author: William J Turner for the LinBox group
 */

#ifndef _TEST_NTL_RR_ENVELOPE_
#define _TEST_NTL_RR_ENVELOPE_

#include "../../test_linbox.h"
#include "LinBox/unparam_field.h"
#include "LinBox/ntl.h"
#include "LinBox/field_archetype.h"
#include "LinBox/field_envelope.h"
#include "LinBox/randiter_envelope.h"
#include "LinBox/randiter_archetype.h"

// Specialization of setup_field for unparam<NTL::RR>
template <> 
bool test_linbox::test<test_linbox::field_categories::ntl_rr_envelope_tag>(void) const
{
  LinBox::unparam_field<NTL::RR> F;
  LinBox::Field_envelope< LinBox::unparam_field<NTL::RR> > E(F);
  LinBox::Field_envelope< LinBox::unparam_field<NTL::RR> >::element e;
  LinBox::Field_envelope< LinBox::unparam_field<NTL::RR> >::randIter r(E);
  LinBox::Field_archetype A(&E, &e, &r);
  return run_tests(A);
} // template <> bool test_linbox<ntl_rr_envelope_tag>(void)

#endif // _TEST_NTL_RR_ENVELOPE_
