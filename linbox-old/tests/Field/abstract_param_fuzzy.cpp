/* File:	tests/abstract_param_fuzzy.cpp
 * Author:	William J. Turner for the LinBox group.
 * 
 * Testing field param_fuzzy
 */

#include <iostream>
#include "LinBox/field_archetype.h"
#include "LinBox/abstract_param_fuzzy.h"
#include "field.h"

int main(void)
{
	double fuzz;
	cout << endl << "Enter a fuzz value: ";
	cin >> fuzz;
	LinBox::abstract_param_fuzzy F(fuzz);
	LinBox::abstract_param_fuzzy::element e;
	LinBox::abstract_param_fuzzy::randIter r(F);
	LinBox::Field_archetype A(&F, &e, &r);
	test_field<LinBox::Field_archetype>(A);
}
