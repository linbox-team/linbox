/* File:	tests/BlackBox/butterfly.h
 * Author:	William J. Turner for the LinBox group
 *
 * Run tests on matrix Diagonal
 */

#include <iostream>
#include <vector>
#include "LinBox/giverror.C"
#include "LinBox/param_modular.h"
#include "LinBox/butterfly.h"
#include "LinBox/boolean_switch.h"
#include "LinBox/cekstv_switch.h"

int main(void)
{
	typedef LinBox::param_modular Field;
	typedef Field::randIter RandIter;
	typedef Field::element Element;
	typedef std::vector<Element> Vector;  // Only dense vectors allowed

	LinBox::integer modulus; 	// prime modulus
	cout << endl << "Enter a prime number for the modulus of the field: ";
	cin >> modulus;
	Field F(modulus);

	cout << "Enter an integer size for the set of random numbers: ";
	LinBox::integer size;
	cin >> size;

	cout << "Enter integer seed for random number generator: ";
	LinBox::integer seed;
	cin >> seed;
	
	Field::randIter r(F, size, seed);

	cout << "Enter the size of the vector to precondition: ";
	size_t n;
	cin >> n;

	cout << "Enter a vector to be switched by the butterfly matrix." << endl
		<< "Input the vector by entering index and value." << endl
		<< "Remember matrices and vectors are indexed starting at 0." << endl
		<< "End with a index of -1." << endl;

	Element zero;
	F.init(zero, 0);
	Element elem(zero);

	Vector x(n, zero), y(n, zero);
	Vector::iterator x_iter, y_iter;

	size_t i;
	bool found;

	while (cin >> i)
	{
		// return also if row index is not positive
		if(i == size_t(-1)) break; 
		
		F.read(cin, elem);

		// Record element in dense vector
 		x[i] = elem;

	}

	// Calculate number of butterfly switches needed
	size_t s(LinBox::count_butterfly(n));

	cout << "Now testing boolean switches..." << endl << endl;

	std::vector<bool> s1(s, false);

	cout << "There are a total of " << s << " switches that can be set." << endl
		<< "Enter the numbers corresponding to the switches, " << endl
		<< "numbered from 0 to " << s - 1 << " you want set." << endl
		<< "End with a switch number of '-1'." << endl;
	
	size_t value;
	while (cin >> value)
	{
		if (value == size_t(-1)) break;
		if (value < s) s1[value] = true;
	}

	LinBox::boolean_switch switch1(s1);
	
	cout << "*** Running tests of hand set boolean switches." << endl;

	cout << "Dense vector x:" << endl;
	i = 0;
	for (x_iter = x.begin(); x_iter != x.end(); x_iter++, i++)
	{
		cout << "\tx[" << i << "] = ";
		F.write(cout, *x_iter); 
		cout << endl;
	}


	LinBox::butterfly < Vector, LinBox::boolean_switch > 
		B1(x.size(), switch1);
	LinBox::Blackbox_archetype<Vector>& A1 = B1;

	A1.apply(y,x);

	cout << "Using A1.apply(y,x) gives vector y:" << endl;
	i = 0;
	for (y_iter = y.begin(); y_iter != y.end(); y_iter++, i++)
	{
		cout << "\ty[" << i << "] = ";
		F.write(cout, *y_iter); 
		cout << endl;
	}

	y = Vector(n, zero);
	y = A1.apply(x);

	cout << "Using y = A1.apply(x) gives vector y:" << endl;
	i = 0;
	for (y_iter = y.begin(); y_iter != y.end(); y_iter++, i++)
	{
		cout << "\ty[" << i << "] = ";
		F.write(cout, *y_iter); 
		cout << endl;
	}

	y = x;
	y = A1.applyin(y);

	cout << "Using y = x; y = A1.applyin(y) gives vector y:" << endl;
	i = 0;
	for (y_iter = y.begin(); y_iter != y.end(); y_iter++, i++)
	{
		cout << "\ty[" << i << "] = ";
		F.write(cout, *y_iter); 
		cout << endl;
	}

	y = Vector(n, zero);
	A1.applyTranspose(y,x);

	cout << "Using A1.applyTranspose(y,x) gives vector y:" << endl;
	i = 0;
	for (y_iter = y.begin(); y_iter != y.end(); y_iter++, i++)
	{
		cout << "\ty[" << i << "] = ";
		F.write(cout, *y_iter); 
		cout << endl;
	}

	y = Vector(n, zero);
	y = A1.applyTranspose(x);

	cout << "Using y = A1.applyTranspose(x) gives vector y:" << endl;
	i = 0;
	for (y_iter = y.begin(); y_iter != y.end(); y_iter++, i++)
	{
		cout << "\ty[" << i << "] = ";
		F.write(cout, *y_iter); 
		cout << endl;
	}

	y = x;
	y = A1.applyTransposein(y);

	cout << "Using y = x; y = A1.applyTransposein(y) gives vector y:" << endl;
	i = 0;
	for (y_iter = y.begin(); y_iter != y.end(); y_iter++, i++)
	{
		cout << "\ty[" << i << "] = ";
		F.write(cout, *y_iter); 
		cout << endl;
	}

	cout << "Now testing multiplication (CEKSTV) switches..." << endl << endl;

	Vector s2(s, zero);

	cout << "There are a total of " << s << " switches that can be set." << endl
		<< "Enter the numbers corresponding to the switches, " << endl
		<< "numbered from 0 to " << s - 1 << " you want set" << endl
		<< "and then the value of the switch." << endl
		<< "End with a switch number of '-1'." << endl;

	while (cin >> value)
	{
		if (value == size_t(-1)) break;
		F.read(cin, elem);
		if (value < s) s2[value] = elem;
	}

	LinBox::cekstv_switch<Field> switch2(F, s2);
	
	cout << "*** Running tests of hand set multiplication switches." << endl;

	cout << "Dense vector x:" << endl;
	i = 0;
	for (x_iter = x.begin(); x_iter != x.end(); x_iter++, i++)
	{
		cout << "\tx[" << i << "] = ";
		F.write(cout, *x_iter); 
		cout << endl;
	}

	LinBox::butterfly < Vector, LinBox::cekstv_switch<Field> > 
		B2(x.size(), switch2);

	LinBox::Blackbox_archetype<Vector>& A2 = B2;

	A2.apply(y,x);

	cout << "Using A2.apply(y,x) gives vector y:" << endl;
	i = 0;
	for (y_iter = y.begin(); y_iter != y.end(); y_iter++, i++)
	{
		cout << "\ty[" << i << "] = ";
		F.write(cout, *y_iter); 
		cout << endl;
	}

	y = Vector(n, zero);
	y = A2.apply(x);

	cout << "Using y = A2.apply(x) gives vector y:" << endl;
	i = 0;
	for (y_iter = y.begin(); y_iter != y.end(); y_iter++, i++)
	{
		cout << "\ty[" << i << "] = ";
		F.write(cout, *y_iter); 
		cout << endl;
	}

	y = x;
	y = A2.applyin(y);

	cout << "Using y = x; y = A2.applyin(y) gives vector y:" << endl;
	i = 0;
	for (y_iter = y.begin(); y_iter != y.end(); y_iter++, i++)
	{
		cout << "\ty[" << i << "] = ";
		F.write(cout, *y_iter); 
		cout << endl;
	}

	y = Vector(n, zero);
	A2.applyTranspose(y,x);

	cout << "Using A2.applyTranspose(y,x) gives vector y:" << endl;
	i = 0;
	for (y_iter = y.begin(); y_iter != y.end(); y_iter++, i++)
	{
		cout << "\ty[" << i << "] = ";
		F.write(cout, *y_iter); 
		cout << endl;
	}

	y = Vector(n, zero);
	y = A2.applyTranspose(x);

	cout << "Using y = A2.applyTranspose(x) gives vector y:" << endl;
	i = 0;
	for (y_iter = y.begin(); y_iter != y.end(); y_iter++, i++)
	{
		cout << "\ty[" << i << "] = ";
		F.write(cout, *y_iter); 
		cout << endl;
	}

	y = x;
	y = A2.applyTransposein(y);

	cout << "Using y = x; y = A2.applyTransposein(y) gives vector y:" << endl;
	i = 0;
	for (y_iter = y.begin(); y_iter != y.end(); y_iter++, i++)
	{
		cout << "\ty[" << i << "] = ";
		F.write(cout, *y_iter); 
		cout << endl;
	}

}
