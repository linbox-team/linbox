/* File:	tests/param_fuzzy.cpp
 * Author:	William J. Turner for the LinBox group.
 * 
 * Testing field param_fuzzy
 */

#include <iostream>
#include "LinBox/param_fuzzy.h"
#include "field.h"

int main(void)
{
	typedef LinBox::param_fuzzy Field;

	double fuzz;
	cout << endl << "Enter a fuzz value: ";
	cin >> fuzz;
	
	Field F(fuzz);
	test_field<Field>(F);
}
