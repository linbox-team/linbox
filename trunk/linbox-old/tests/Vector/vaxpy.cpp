/* File:	tests/BlackBox/vaxpy.h
 * Author:	William J. Turner for the LinBox group
 *
 * Run tests on vector utility vaxpy
 */

#include <iostream>
#include <vector>
#include <list>
#include <utility>
#include <map>
#include <algorithm>
#include "LinBox/param_modular.h"
#include "LinBox/vaxpy.h"

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

	cout << "What is the size of the vectors to create? ";

	size_t n;
	cin >> n;

	cout << "Enter the first vector, x." << endl
		<< "Input the vector by entering index and value." << endl
		<< "Remember matrices and vectors are indexed starting at 0." << endl
		<< "End with a index of -1." << endl;

	typedef std::vector<Element> Vector1;
	typedef std::list< std::pair<size_t, Element> > Vector2;
	typedef std::map<size_t, Element> Vector3;

	Element elem, zero;
	F.init(zero, 0);
	F.init(elem, 0);

	Vector1 x1(n, zero), y1(n, zero), z1(n, zero);
	Vector2 x2, y2, z2;
	Vector3 x3, y3, z3;

	Vector1::iterator x1_iter, y1_iter, z1_iter;
	Vector2::iterator x2_iter, y2_iter, z2_iter;
	Vector3::iterator x3_iter, y3_iter, z3_iter;

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

	cout << "Enter the second vector, y." << endl
		<< "Input the vector by entering index and value." << endl
		<< "End with a index of -1." << endl;

	while (cin >> i)
	{
		// return also if row index is not positive
		if(i == size_t(-1)) break; 
		
		F.read(cin, elem);

		// Record element in dense vector
 		y1[i] = elem;

		// Record element in sparse sequence vector

		// find appropriate location of element
		if( y2.begin() == y2.end() )
			y2_iter = y2.end();
		else
			y2_iter = lower_bound( y2.begin(), y2.end(), i, comp_w_ind<Field>() );

		// Check to see if element already exists.
		if ( y2.end() == y2_iter )
			found = false;
		else if ( y2_iter->first != i )
			found = false;
		else 
			found = true;

		// If element is already in row, replace old value with new.
		// Otherwise, insert the element in the row.
		if (found) 
		{
			if (F.isZero(elem))
				y2.erase(y2_iter);
			else
				y2_iter->second = elem;
		} // if (found)
		else if (!F.isZero(elem))
			y2.insert(y2_iter, make_pair(i,elem));

		// Record element in sparse associative vector

		// Find element in associative container.  
		// If exists, replace value if not zero, or remove if value is zero.
		// If not found, insert non-zero element
		if ( (y3_iter = y3.find(i)) != y3.end() )
		{
			if (F.isZero(elem))
				y3.erase(y3_iter);
			else
				y3_iter->second = elem;
		}
		else
		{
			if (!F.isZero(elem))
			y3.insert(make_pair(i, elem));
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

	cout << "Dense vector y1:" << endl;
	i = 0;
	for (y1_iter = y1.begin(); y1_iter != y1.end(); y1_iter++, i++)
	{
		cout << "\ty1[" << i << "] = ";
		F.write(cout, *y1_iter); 
		cout << endl;
	}

	LinBox::vaxpy(F, z1, a, x1, y1);

	cout << "Using vaxpy(F, z1, a, x1, y1) gives vector z1:" << endl;
	i = 0;
	for (z1_iter = z1.begin(); z1_iter != z1.end(); z1_iter++, i++)
	{
		cout << "\tz1[" << i << "] = ";
		F.write(cout, *z1_iter); 
		cout << endl;
	}

	z1 = Vector1(n, zero);
	z1 = LinBox::vaxpy(F, a, x1, y1);

	cout << "Using z1 = vaxpy(F, a, x1, y1) gives vector z1:" << endl;
	i = 0;
	for (z1_iter = z1.begin(); z1_iter != z1.end(); z1_iter++, i++)
	{
		cout << "\tz1[" << i << "] = ";
		F.write(cout, *z1_iter); 
		cout << endl;
	}

	z1 = y1;
	z1 = LinBox::vaxpyin(F, z1, a, x1);

	cout << "Using z1 = y1; vaxpyin(F, z1, a, x1) gives vector z1:" << endl;
	i = 0;
	for (z1_iter = z1.begin(); z1_iter != z1.end(); z1_iter++, i++)
	{
		cout << "\tz1[" << i << "] = ";
		F.write(cout, *z1_iter); 
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

	cout << "Sparse sequence vector y2:" << endl;
	for (y2_iter = y2.begin(); y2_iter != y2.end(); y2_iter++)
	{
		cout << "\ty2[" << y2_iter->first << "] = ";
		F.write(cout, y2_iter->second); 
		cout << endl;
	}

	LinBox::vaxpy(F, z2, a, x2, y2);

	cout << "Using vaxpy(F, z2, a, x2, y2) gives vector z2:" << endl;
	for (z2_iter = z2.begin(); z2_iter != z2.end(); z2_iter++)
	{
		cout << "\tz2[" << z2_iter->first << "] = ";
		F.write(cout, z2_iter->second); 
		cout << endl;
	}

	z2 = Vector2();
	z2 = LinBox::vaxpy(F, a, x2, y2);

	cout << "Using z2 = vaxpy(F, a, x2, y2) gives vector z2:" << endl;
	for (z2_iter = z2.begin(); z2_iter != z2.end(); z2_iter++)
	{
		cout << "\tz2[" << z2_iter->first << "] = ";
		F.write(cout, z2_iter->second); 
		cout << endl;
	}

	z2 = y2;
	z2 = LinBox::vaxpyin(F, z2, a, x2);

	cout << "Using z2 = y2; vaxpyin(F, z2, a, x2) gives vector z2:" << endl;
	for (z2_iter = z2.begin(); z2_iter != z2.end(); z2_iter++)
	{
		cout << "\tz2[" << z2_iter->first << "] = ";
		F.write(cout, z2_iter->second); 
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

	cout << "Sparse sequence vector y3:" << endl;
	for (y3_iter = y3.begin(); y3_iter != y3.end(); y3_iter++)
	{
		cout << "\ty3[" << y3_iter->first << "] = ";
		F.write(cout, y3_iter->second); 
		cout << endl;
	}

	LinBox::vaxpy(F, z3, a, x3, y3);

	cout << "Using vaxpy(F, z3, a, x3, y3) gives vector z3:" << endl;
	for (z3_iter = z3.begin(); z3_iter != z3.end(); z3_iter++)
	{
		cout << "\tz3[" << z3_iter->first << "] = ";
		F.write(cout, z3_iter->second); 
		cout << endl;
	}

	z3 = Vector3();
	z3 = LinBox::vaxpy(F, a, x3, y3);

	cout << "Using z3 = vaxpy(F, a, x3, y3) gives vector z3:" << endl;
	for (z3_iter = z3.begin(); z3_iter != z3.end(); z3_iter++)
	{
		cout << "\tz3[" << z3_iter->first << "] = ";
		F.write(cout, z3_iter->second); 
		cout << endl;
	}

	z3 = y3;
	z3 = LinBox::vaxpyin(F, z3, a, x3);

	cout << "Using z3 = y3; vaxpyin(F, z3, a, x3) gives vector z1:" << endl;
	i = 0;
	for (z3_iter = z3.begin(); z3_iter != z3.end(); z3_iter++)
	{
		cout << "\tz3[" << z3_iter->first << "] = ";
		F.write(cout, z3_iter->second); 
		cout << endl;
	}


}
