/* File:	tests/BlackBox/diagonal.h
 * Author:	William J. Turner for the LinBox group
 *
 * Run tests on matrix Diagonal
 */

#include <iostream>
#include <vector>
#include <list>
#include <utility>
#include <map>
#include <algorithm>
#include "LinBox/param_modular.h"
//#include "LinBox/diagonal.h"

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

	cout << "What is the size of the square diagonal matrix to create? ";

	size_t n;
	cin >> n;

	cout << "The size of the diagonal matrix is " << n << endl;

	Element elem, zero;
	F.init(zero, 0);
	F.init(elem, 0);
	std::vector<Element> d(n, zero);

	cout << endl
		<< "Input diagonal entries by entering index and value." << endl
		<< "Remember matrices and vectors are indexed starting at 0." << endl
		<< "End with a index of -1." << endl;

	size_t i;

	while (cin >> i)
	{
		// return also if row index is not positive
		if(i == size_t(-1)) break; 

		F.read(cin, elem);
		d[i] = elem;
	} // while (is >> i)

	cout << "The diagonal vector in the diagonal test is the following:" << endl;

	std::vector<Element>::iterator iter;
	i = 0;

	for (iter = d.begin(); iter != d.end(); iter++, i++)
	{
		cout << "\td[" << i << "] = ";
		F.write(cout, *iter); // can't do this!!!
		cout << endl;
	}

	cout << "Enter a vector to be multiplied by the diagonal matrix." << endl
		<< "Input the vector by entering index and value." << endl
		<< "Remember matrices and vectors are indexed starting at 0." << endl
		<< "End with a index of -1." << endl;

	typedef std::vector<Element> Vector1;
	typedef std::list< std::pair<size_t, Element> > Vector2;
	typedef std::map<size_t, Element> Vector3;

	Vector1 x1(n, zero), y1(n, zero);
	Vector2 x2, y2;
	Vector3 x3, y3;

	Vector1::iterator x1_iter, y1_iter;
	Vector2::iterator x2_iter, y2_iter;
	Vector3::iterator x3_iter, y3_iter;

	bool found;

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

	} // while (cin >> i)

	cout << "*** Running tests with dense vector." << endl;

	cout << "Dense vector x1:" << endl;
	i = 0;
	for (x1_iter = x1.begin(); x1_iter != x1.end(); x1_iter++, i++)
	{
		cout << "\tx1[" << i << "] = ";
		F.write(cout, *x1_iter); // can't do this!!!
		cout << endl;
	}


//	LinBox::Diagonal<Field, Vector1> A1(F,d);
	
//	A1.apply(y1,x1);

	cout << "Using A1.apply(y1,x1) gives vector y1:" << endl;
	i = 0;
	for (y1_iter = y1.begin(); y1_iter != y1.end(); y1_iter++, i++)
	{
		cout << "\ty1[" << i << "] = ";
		F.write(cout, *y1_iter); // can't do this!!!
		cout << endl;
	}

	y1 = Vector1(n, zero);

//	y1 = A1.apply(x1);

	cout << "Using y1 = A1.apply(x1) gives vector y1:" << endl;
	i = 0;
	for (y1_iter = y1.begin(); y1_iter != y1.end(); y1_iter++, i++)
	{
		cout << "\ty1[" << i << "] = ";
		F.write(cout, *y1_iter); // can't do this!!!
		cout << endl;
	}

	y1 = x1;

//	y1 = A1.applyin(y1);

	cout << "Using y1 = x1; y1 = A1.applyin(y1) gives vector y1:" << endl;
	i = 0;
	for (y1_iter = y1.begin(); y1_iter != y1.end(); y1_iter++, i++)
	{
		cout << "\ty1[" << i << "] = ";
		F.write(cout, *y1_iter); // can't do this!!!
		cout << endl;
	}





/*
	Vector x = *new_vector(n,prompt,cin,cout);

	cout << endl
	<< "The input vector to the diagonal test is the following:" << endl;

	write(cout, x);

	typename Field::element zero;
	F.init(zero, 0);

	Vector b;

	LinBox::Diagonal<Field, Vector> D(F, d);

	apply_blackbox(b, D, x, _mode);

	cout << endl
	<< "Applying the diagonal matrix to the vector gives " << endl
	<< "the following vector" << endl;

	write(cout, b);

	return true;

*/

}
