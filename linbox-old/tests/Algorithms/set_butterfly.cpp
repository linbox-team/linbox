/* File:	tests/Algorithms/set_butterfly.cpp
 * Author:	William J. Turner for the LinBox group
 *
 * Run tests set_butterfly algorithm to set switches needed for
 * given contiguous block.
 */

#include <iostream>
#include <vector>
#include "LinBox/giverror.C"
#include "LinBox/unparam_field.h"
#include "LinBox/butterfly.h"
#include "LinBox/boolean_switch.h"

int main(void)
{
	typedef LinBox::unparam_field<double> Field;
	typedef Field::randIter RandIter;
	typedef Field::element Element;
	typedef std::vector<Element> Vector;

	Field F;

	cout << "Enter the number of inputs to consider: ";
	size_t n;
	cin >> n;

	std::vector<bool> x(n, false);

	cout << "Enter the numbers corresponding to the linearly independent rows" << endl
		<< "with the rows numbered from 0 to " << n - 1 << "." << endl
		<< "End with a row number of '-1'." << endl;

	size_t value;
	while (cin >> value)
	{
		if (value == size_t(-1)) break;
		if (value < x.size()) x[value] = true;
	}
	
	cout << "The test vector is:" << endl;
	for (size_t i = 0; i < x.size(); i++)
		cout << "\tx[" << i << "] = " << x[i] << endl;

	cout << "Enter non-negative offset for contiguous block: ";

	cin >> value;

	cout << "The offset is " << value << endl;

	std::vector<bool> s1 = LinBox::set_butterfly(x, value, cout);
	LinBox::boolean_switch switch1(s1);

	// Create vector to switch.
	Element zero, one;
	F.init(zero, 0);
	F.init(one, 1);
	Element elem(zero);

	Vector x1, x2(n, zero), y(n, zero);
	Vector::iterator x_iter, y_iter;

	for (std::vector<bool>::iterator iter = x.begin(); iter != x.end(); iter++)
		x1.push_back( (*iter) ? one : zero );

	cout << "*** Running tests." << endl;

	cout << "Dense vector x1:" << endl;
	size_t i = 0;
	for (x_iter = x1.begin(); x_iter != x1.end(); x_iter++, i++)
	{
		cout << "\tx1[" << i << "] = ";
		F.write(cout, *x_iter); 
		cout << endl;
	}

	LinBox::butterfly < Vector, LinBox::boolean_switch > 
		B1(x.size(), switch1);
	LinBox::Blackbox_archetype<Vector>& A1 = B1;

	A1.apply(y,x1);

	cout << "Using A1.apply(y,x1) gives vector y:" << endl;
	i = 0;
	for (y_iter = y.begin(); y_iter != y.end(); y_iter++, i++)
	{
		cout << "\ty[" << i << "] = ";
		F.write(cout, *y_iter); 
		cout << endl;
	}

	cout << "Dense vector x2:" << endl;
	i = 0;
	for (x_iter = x2.begin(); x_iter != x2.end(); x_iter++, i++)
	{
		F.assign(*x_iter, i);
		cout << "\tx2[" << i << "] = ";
		F.write(cout, *x_iter); 
		cout << endl;
	}

	A1.apply(y,x2);

	cout << "Using A1.apply(y,x2) gives vector y:" << endl;
	i = 0;
	for (y_iter = y.begin(); y_iter != y.end(); y_iter++, i++)
	{
		cout << "\ty[" << i << "] = ";
		F.write(cout, *y_iter); 
		cout << endl;
	}

}
