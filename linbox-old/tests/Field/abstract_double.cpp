/* File:	tests/abstract_double.cpp
 * Author:	William J. Turner for the LinBox group.
 * 
 * Testing field abstract_double
 */

#include "LinBox/field_archetype.h"
#include "LinBox/abstract_double.h"
#include "field.h"

int main(void)
{
	LinBox::abstract_double F;
	LinBox::abstract_double::element e;
	LinBox::abstract_double::randIter r(F);
	LinBox::Field_archetype A(&F, &e, &r);
	test_field<LinBox::Field_archetype>(A);
}


