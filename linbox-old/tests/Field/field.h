/* File:	tests/field.h
 * Author:	William J. Turner for the LinBox group.
 * 
 * Testing function for LinBox fields.
 * For its use, please see double.cpp, modular.cpp, etc.
 */

#ifndef _FIELD_
#define _FIELD_

#include <iostream>
#include "LinBox/faxpy.h"

template<class Field> bool test_field(const Field& F)
{
	LinBox::integer temp_int;

	cout << "Using the field ";
	F.write(cout);
	cout << ", which has cardinality " << F.cardinality(temp_int) 
		<< " and charactersitic " << F.characteristic(temp_int) 
		<< endl;

	cout << "Initializing temporary and constants elements" << endl
		<< "and creating other elements through copy constructor." 
		<< endl;

	typename Field::element temp, one, zero;
	F.init(temp, temp_int);
	F.init(one, 1);
	F.init(zero, 0);

	typename Field::element a(temp), b(temp), c(temp);

	cout << "Enter two integers to which to initialize three field elements: ";

	F.read(cin, a);
	F.read(cin, b);

	cout << "The two integers entered were ";
	F.write(cout, a);
	cout << " and ";
	F.write(cout, b);
	cout << endl;

	cout << "Enter an integer size for the set of random numbers: ";
	LinBox::integer size;
	cin >> size;

	cout << "Enter integer seed for random number generator: ";
	LinBox::integer seed;
	cin >> seed;
	
	typename Field::randIter r(F, size, seed);
	
	cout << "Assigning first element to third and then converting to integer."
	   << endl;

	F.assign(c, a);
	F.convert(temp_int, c);

	cout << "Testing assign and convert: ";
	F.write(cout, a);
	cout << " = ";
	F.write(cout, c);
	cout << " = " << temp_int << endl;

	cout << "Testing arithmetic functions." << endl;

	if (F.areEqual(zero, zero))
	cout << "\t0 == 0" << endl;
	else
	cout << "\t0 != 0" << endl;

	if (F.areEqual(zero, one))
	cout << "\t0 == 1" << endl;
	else
	cout << "\t0 != 1" << endl;

	cout << "\t";
	F.write(cout, a);
	cout << " + ";
	F.write(cout, b);
	cout << " = ";
	F.write(cout, F.add(c, a, b));
	cout << " = ";
	F.write(cout, c);
	cout << endl;

	cout << "\t";
	F.write(cout, a);
	cout << " - ";
	F.write(cout, b);
	cout << " = ";
	F.write(cout, F.sub(c, a, b));
	cout << " = ";
	F.write(cout, c);
	cout << endl;

	cout << "\t";
	F.write(cout, a);
	cout << " * ";
	F.write(cout, b);
	cout << " = ";
	F.write(cout, F.mul(c, a, b));
	cout << " = ";
	F.write(cout, c);
	cout << endl;

	cout << "\t";
	F.write(cout, a);
	cout << " / ";
	F.write(cout, b);
	cout << " = ";
	F.write(cout, F.div(c, a, b));
	cout << " = ";
	F.write(cout, c);
	cout << endl;

	cout << "\t- ";
	F.write(cout, a);
	cout << " = ";
	F.write(cout, F.neg(c, a));
	cout << " = ";
	F.write(cout, c);
	cout << endl;

	cout << "\t1 / ";
	F.write(cout, a);
	cout << " = ";
	F.write(cout, F.inv(c, a));
	cout << " = ";
	F.write(cout, c);
	cout << endl;

	cout << "Testing in-place arithmetic functions." << endl;

	if (F.isZero(zero))
	cout << "\t0 == 0" << endl;
	else
	cout << "\t0 != 0" << endl;

	if (F.isZero(one))
	cout << "\t1 == 0" << endl;
	else
	cout << "\t1 != 0" << endl;

	if (F.isOne(one))
	cout << "\t1 == 1" << endl;
	else
	cout << "\t1 != 1" << endl;

	if (F.isOne(zero))
	cout << "\t0 == 1" << endl;
	else
	cout << "\t0 != 1" << endl;

	F.assign(c, a);
	cout << "\t";
	F.write(cout, c);
	cout << " += ";
	F.write(cout, b);
	cout << " = ";
	F.write(cout, F.addin(c, b));
	cout << " = ";
	F.write(cout, c);
	cout << endl;

	F.assign(c, a);
	cout << "\t";
	F.write(cout, c);
	cout << " -= ";
	F.write(cout, b);
	cout << " = ";
	F.write(cout, F.subin(c, b));
	cout << " = ";
	F.write(cout, c);
	cout << endl;

	F.assign(c, a);
	cout << "\t";
	F.write(cout, c);
	cout << " *= ";
	F.write(cout, b);
	cout << " = ";
	F.write(cout, F.mulin(c, b));
	cout << " = ";
	F.write(cout, c);
	cout << endl;

	F.assign(c, a);
	cout << "\t";
	F.write(cout, c);
	cout << " /= ";
	F.write(cout, b);
	cout << " = ";
	F.write(cout, F.divin(c, b));
	cout << " = ";
	F.write(cout, c);
	cout << endl;

	F.assign(c, a);
	cout << "\t- ";
	F.write(cout, c);
	cout << " = ";
	F.write(cout, F.negin(c));
	cout << " = ";
	F.write(cout, c);
	cout << endl;

	F.assign(c, a);
	cout << "\t1 / ";
	F.write(cout, c);
	cout << " = ";
	F.write(cout, F.invin(c));
	cout << " = ";
	F.write(cout, c);
	cout << endl;

	cout << "Testing faxpy." << endl;

	F.init(temp, 2);
	LinBox::faxpy<Field> Faxpy(F, temp);

	cout << "\t";
	F.write(cout, temp);
	cout << " * ";
	F.write(cout, a);
	cout << " + ";
	F.write(cout, b);
	cout << " = ";
	F.write(cout, Faxpy.apply(c, a, b));
	cout << " = ";
	F.write(cout, c);
	cout << endl;

	F.assign(c, b);
	cout << "\t";
	F.write(cout, c);
	cout << " += ";
	F.write(cout, temp);
	cout << " * ";
	F.write(cout, a);
	cout << " = ";
	F.write(cout, Faxpy.applyin(c, a));
	cout << " = ";
	F.write(cout, c);
	cout << endl;

	F.init(temp, 5);
	Faxpy.assign(temp);

	cout << "\t";
	F.write(cout, temp);
	cout << " * ";
	F.write(cout, a);
	cout << " + ";
	F.write(cout, b);
	cout << " = ";
	F.write(cout, Faxpy.apply(c, a, b));
	cout << " = ";
	F.write(cout, c);
	cout << endl;

	F.assign(c, b);
	cout << "\t";
	F.write(cout, c);
	cout << " += ";
	F.write(cout, temp);
	cout << " * ";
	F.write(cout, a);
	cout << " = ";
	F.write(cout, Faxpy.applyin(c, a));
	cout << " = ";
	F.write(cout, c);
	cout << endl;

	int n = 10;
	cout << "Generating " << n << " random field elements." << endl;

	for (int i = 0; i < 10; i++)
	{
		cout << "\t";
		F.write(cout, r());
		cout << endl;
	}

	return true;

} // template<class Field> bool test_field(const Field& F)

#endif // _FIELD_
