/* File:	tests/modular.cpp
 * Author:	William J. Turner for the LinBox group.
 * 
 * Testing field unparam_field<modular>
 */

#include <iostream>
#include "LinBox/ntl.h"
#include "field.h"

int main(void)
{
	typedef LinBox::unparam_field<NTL::zz_p> Field;

	long modulus; 	// prime modulus
	cout << endl << "Enter a prime number for the modulus of the field: ";
	cin >> modulus;
	NTL::zz_p::init(modulus);
	Field F;
	test_field<Field>(F);
}

