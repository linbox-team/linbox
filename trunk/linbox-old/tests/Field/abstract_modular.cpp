/* File:	tests/abstract_modular.cpp
 * Author:	William J. Turner for the LinBox group.
 * 
 * Testing field abstract_modular
 */

#include <iostream>
#include "LinBox/field_archetype.h"
#include "LinBox/abstract_modular.h"
#include "field.h"

int main(void)
{
	LinBox::integer modulus; 	// prime modulus
	cout << endl << "Enter a prime number for the modulus of the field: ";
	cin >> modulus;
	LinBox::abstract_modular F(modulus);
	LinBox::abstract_modular::element e;
	LinBox::abstract_modular::randIter r(F);
	LinBox::Field_archetype A(&F, &e, &r);
	test_field<LinBox::Field_archetype>(A);
}

