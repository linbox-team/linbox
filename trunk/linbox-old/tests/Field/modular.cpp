/* File:	tests/modular.cpp
 * Author:	William J. Turner for the LinBox group.
 * 
 * Testing field unparam_field<modular>
 */

#include <iostream>
#include "LinBox/unparam_field.h"
#include "LinBox/modular.h"
#include "field.h"

int main(void)
{
	typedef LinBox::unparam_field<LinBox::modular> Field;

	LinBox::integer modulus; 	// prime modulus
	cout << endl << "Enter a prime number for the modulus of the field: ";
	cin >> modulus;
	LinBox::modular::put_modulus(modulus);
	Field F;
	test_field<Field>(F);
}

