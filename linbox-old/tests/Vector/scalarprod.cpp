/* File:	tests/BlackBox/scalarprod.h
 * Author:	William J. Turner for the LinBox group
 *
 * Run tests on vector utility scalarprod
 */

#include <iostream>
#include <vector>
#include <list>
#include <utility>
#include <map>
#include <algorithm>
#include "LinBox/param_modular.h"
#include "LinBox/scalarprod.h"

template <class Field> struct comp_w_ind
{ 
	bool operator() (const std::pair< size_t, typename Field::element >& entry, 
			size_t col_in)
	{ return entry.first < col_in; }
};

int main(void)
{
	typedef LinBox::param_modular Field;
	typedef Field::element Element;

	LinBox::integer modulus; 	// prime modulus
	cout << endl << "Enter a prime number for the modulus of the field: ";
	cin >> modulus;
	Field F(modulus);

	cout << "Enter integer to be used as scalar multiple: ";
	Element a;
	LinBox::integer a_int;
	cin >> a_int;
	F.init(a, a_int);

	cout << "What is the size of the vector to create? ";

	size_t n;
	cin >> n;

	cout << "Input the vector by entering index and value." << endl
		<< "Remember matrices and vectors are indexed starting at 0." << endl
		<< "End with a index of -1." << endl;

	typedef std::vector<Element> Vector1;
	typedef std::list< std::pair<size_t, Element> > Vector2;
	typedef std::map<size_t, Element> Vector3;

	Element elem, zero;
	F.init(zero, 0);
	F.init(elem, 0);

	Vector1 x1(n, zero), y1(n, zero);
	Vector2 x2, y2;
	Vector3 x3, y3;

	Vector1::iterator x1_iter, y1_iter;
	Vector2::iterator x2_iter, y2_iter;
	Vector3::iterator x3_iter, y3_iter;

	bool found;

	size_t i;

	while (cin >> i)
	{
		// return also if row index is not positive
		if(i == size_t(-1)) break; 
		
		F.read(cin, elem);

		// Record element in dense vector
 		x1[i] = elem;

		// Record element in sparse sequence vector

		// find appropriate location of element
		if( x2.begin() == x2.end() )
			x2_iter = x2.end();
		else
			x2_iter = lower_bound( x2.begin(), x2.end(), i, comp_w_ind<Field>() );

		// Check to see if element already exists.
		if ( x2.end() == x2_iter )
			found = false;
		else if ( x2_iter->first != i )
			found = false;
		else 
			found = true;

		// If element is already in row, replace old value with new.
		// Otherwise, insert the element in the row.
		if (found) 
		{
			if (F.isZero(elem))
				x2.erase(x2_iter);
			else
				x2_iter->second = elem;
		} // if (found)
		else if (!F.isZero(elem))
			x2.insert(x2_iter, make_pair(i,elem));

		// Record element in sparse associative vector

		// Find element in associative container.  
		// If exists, replace value if not zero, or remove if value is zero.
		// If not found, insert non-zero element
		if ( (x3_iter = x3.find(i)) != x3.end() )
		{
			if (F.isZero(elem))
				x3.erase(x3_iter);
			else
				x3_iter->second = elem;
		}
		else
		{
			if (!F.isZero(elem))
			x3.insert(make_pair(i, elem));
		}

	} // while (cin >> i)

	cout << "*** Running tests with dense vector." << endl;

	cout << "Dense vector x1:" << endl;
	i = 0;
	for (x1_iter = x1.begin(); x1_iter != x1.end(); x1_iter++, i++)
	{
		cout << "\tx1[" << i << "] = ";
		F.write(cout, *x1_iter); 
		cout << endl;
	}

	LinBox::scalarprod(F, y1, a, x1);

	cout << "Using scalarprod(F, y1, a, x1) gives vector y1:" << endl;
	i = 0;
	for (y1_iter = y1.begin(); y1_iter != y1.end(); y1_iter++, i++)
	{
		cout << "\ty1[" << i << "] = ";
		F.write(cout, *y1_iter); 
		cout << endl;
	}

	y1 = Vector1(n, zero);
	y1 = LinBox::scalarprod(F, a, x1);

	cout << "Using y1 = scalarprod(F, a, x1) gives vector y1:" << endl;
	i = 0;
	for (y1_iter = y1.begin(); y1_iter != y1.end(); y1_iter++, i++)
	{
		cout << "\ty1[" << i << "] = ";
		F.write(cout, *y1_iter); 
		cout << endl;
	}

	y1 = x1;
	y1 = LinBox::scalarprodin(F, x1, a);

	cout << "Using y1 = x1; scalarprodin(F, x1, a) gives vector y1:" << endl;
	i = 0;
	for (y1_iter = y1.begin(); y1_iter != y1.end(); y1_iter++, i++)
	{
		cout << "\ty1[" << i << "] = ";
		F.write(cout, *y1_iter); 
		cout << endl;
	}

	cout << "*** Running tests with sparse sequence vector." << endl;

	cout << "Sparse sequence vector x2:" << endl;
	for (x2_iter = x2.begin(); x2_iter != x2.end(); x2_iter++)
	{
		cout << "\tx2[" << x2_iter->first << "] = ";
		F.write(cout, x2_iter->second); 
		cout << endl;
	}

	LinBox::scalarprod(F, y2, a, x2);

	cout << "Using scalarprod(F, y2, a, x2) gives vector y2:" << endl;
	for (y2_iter = y2.begin(); y2_iter != y2.end(); y2_iter++)
	{
		cout << "\ty2[" << y2_iter->first << "] = ";
		F.write(cout, y2_iter->second); 
		cout << endl;
	}

	y2 = Vector2();
	y2 = LinBox::scalarprod(F, a, x2);

	cout << "Using y2 = scalarprod(F, a, x2) gives vector y2:" << endl;
	for (y2_iter = y2.begin(); y2_iter != y2.end(); y2_iter++)
	{
		cout << "\ty2[" << y2_iter->first << "] = ";
		F.write(cout, y2_iter->second); 
		cout << endl;
	}

	y2 = x2;
	y2 = LinBox::scalarprodin(F, y2, a);

	cout << "Using y2 = x2; scalarprodin(F, y2, a) gives vector y2:" << endl;
	for (y2_iter = y2.begin(); y2_iter != y2.end(); y2_iter++)
	{
		cout << "\ty2[" << y2_iter->first << "] = ";
		F.write(cout, y2_iter->second); 
		cout << endl;
	}

	cout << "*** Running tests with sparse associative vector." << endl;

	cout << "Sparse sequence vector x3:" << endl;
	for (x3_iter = x3.begin(); x3_iter != x3.end(); x3_iter++)
	{
		cout << "\tx3[" << x3_iter->first << "] = ";
		F.write(cout, x3_iter->second); 
		cout << endl;
	}

	LinBox::scalarprod(F, y3, a, x3);

	cout << "Using scalarprod(F, y3, a, x3) gives vector y3:" << endl;
	for (y3_iter = y3.begin(); y3_iter != y3.end(); y3_iter++)
	{
		cout << "\ty3[" << y3_iter->first << "] = ";
		F.write(cout, y3_iter->second); 
		cout << endl;
	}

	y3 = Vector3();
	y3 = LinBox::scalarprod(F, a, x3);

	cout << "Using y3 = scalarprod(F, a, x3) gives vector y3:" << endl;
	for (y3_iter = y3.begin(); y3_iter != y3.end(); y3_iter++)
	{
		cout << "\ty3[" << y3_iter->first << "] = ";
		F.write(cout, y3_iter->second); 
		cout << endl;
	}

	y3 = x3;
	y3 = LinBox::scalarprodin(F, y3, a);

	cout << "Using y3 = x3; scalarprodin(F, y3, a) gives vector y3:" << endl;
	for (y3_iter = y3.begin(); y3_iter != y3.end(); y3_iter++)
	{
		cout << "\ty3[" << y3_iter->first << "] = ";
		F.write(cout, y3_iter->second); 
		cout << endl;
	}

}
