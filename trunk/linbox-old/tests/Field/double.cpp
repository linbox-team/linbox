/* File:	tests/double.cpp
 * Author:	William J. Turner for the LinBox group.
 * 
 * Testing field unparam_field<double>
 */

#include "LinBox/unparam_field.h"
#include "field.h"

int main(void)
{
	typedef LinBox::unparam_field<double> Field;
	Field F;
	test_field<Field>(F);
}

