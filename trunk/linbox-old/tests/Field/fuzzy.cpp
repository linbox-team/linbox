/* File:	tests/fuzzy.cpp
 * Author:	William J. Turner for the LinBox group.
 * 
 * Testing field unparam_field<fuzzy>
 */

#include <iostream>
#include "LinBox/unparam_field.h"
#include "LinBox/fuzzy.h"
#include "field.h"

int main(void)
{
	typedef LinBox::unparam_field<LinBox::fuzzy> Field;

	double fuzz;
	cout << endl << "Enter a fuzz value: ";
	cin >> fuzz;
	LinBox::fuzzy::put_fuzz(fuzz);
	
	Field F;
	test_field<Field>(F);
}
