/* File:	tests/param_modular.cpp
 * Author:	William J. Turner for the LinBox group.
 * 
 * Testing field param_modular
 */

#include <iostream>
#include "LinBox/param_modular.h"
#include "field.h"

int main(void)
{
	typedef LinBox::param_modular Field;

	LinBox::integer modulus; 	// prime modulus
	cout << endl << "Enter a prime number for the modulus of the field: ";
	cin >> modulus;
	Field F(modulus);
	test_field<Field>(F);
}

