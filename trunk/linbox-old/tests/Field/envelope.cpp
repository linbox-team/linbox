/* File:	tests/envelope.cpp
 * Author:	William J. Turner for the LinBox group.
 * 
 * Testing field envelope on unparam_field<double>
 */

#include "LinBox/unparam_field.h"
#include "LinBox/field_archetype.h"
#include "field.h"

int main(void)
{
	LinBox::unparam_field<double> F;
	LinBox::Field_archetype A(&F);
	test_field<LinBox::Field_archetype>(A);
}

