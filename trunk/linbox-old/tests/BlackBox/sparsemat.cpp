/* File:	tests/BlackBox/sparsemat.h
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
#include "LinBox/giverror.C"
#include "LinBox/param_modular.h"
#include "LinBox/sparsemat.h"

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
	typedef std::list< pair<size_t, Element> > Row;

	LinBox::integer modulus; 	// prime modulus
	cout << endl << "Enter a prime number for the modulus of the field: ";
	cin >> modulus;
	Field F(modulus);

	LinBox::sparsemat_base<Field, Row> 
		B(*LinBox::newSparsemat<Field, Row>(F));

	size_t m, n;
	m = B.get_rowdim();
	n = B.get_coldim();

	cout << "The sparesemat matrix contains the following elements:" 
		<< endl << B;

	cout << "Enter a vector to be multiplied by the sparsemat matrix." << endl
		<< "Input the vector by entering index and value." << endl
		<< "Remember matrices and vectors are indexed starting at 0." << endl
		<< "End with a index of -1." << endl;

	typedef std::vector<Element> Vector1;
	typedef std::list< std::pair<size_t, Element> > Vector2;
	typedef std::map<size_t, Element> Vector3;

	Element zero;
	F.init(zero, 0);
	Element elem(zero);

	Vector1 x1(n, zero), y1(m, zero);
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


	LinBox::sparsemat<Field, Row, Vector1> S1(B);
	LinBox::Blackbox_archetype<Vector1>& A1 = S1;

	A1.apply(y1,x1);

	cout << "Using A1.apply(y1,x1) gives vector y1:" << endl;
	i = 0;
	for (y1_iter = y1.begin(); y1_iter != y1.end(); y1_iter++, i++)
	{
		cout << "\ty1[" << i << "] = ";
		F.write(cout, *y1_iter); 
		cout << endl;
	}

	y1 = Vector1(m, zero);
	y1 = A1.apply(x1);

	cout << "Using y1 = A1.apply(x1) gives vector y1:" << endl;
	i = 0;
	for (y1_iter = y1.begin(); y1_iter != y1.end(); y1_iter++, i++)
	{
		cout << "\ty1[" << i << "] = ";
		F.write(cout, *y1_iter); 
		cout << endl;
	}

	y1 = x1;
	y1 = A1.applyin(y1);

	cout << "Using y1 = x1; y1 = A1.applyin(y1) gives vector y1:" << endl;
	i = 0;
	for (y1_iter = y1.begin(); y1_iter != y1.end(); y1_iter++, i++)
	{
		cout << "\ty1[" << i << "] = ";
		F.write(cout, *y1_iter); 
		cout << endl;
	}

	y1 = Vector1(m, zero);
	A1.applyTranspose(y1,x1);

	cout << "Using A1.applyTranspose(y1,x1) gives vector y1:" << endl;
	i = 0;
	for (y1_iter = y1.begin(); y1_iter != y1.end(); y1_iter++, i++)
	{
		cout << "\ty1[" << i << "] = ";
		F.write(cout, *y1_iter); 
		cout << endl;
	}

	y1 = Vector1(m, zero);
	y1 = A1.applyTranspose(x1);

	cout << "Using y1 = A1.applyTranspose(x1) gives vector y1:" << endl;
	i = 0;
	for (y1_iter = y1.begin(); y1_iter != y1.end(); y1_iter++, i++)
	{
		cout << "\ty1[" << i << "] = ";
		F.write(cout, *y1_iter); 
		cout << endl;
	}

	y1 = x1;
	y1 = A1.applyTransposein(y1);

	cout << "Using y1 = x1; y1 = A1.applyTransposein(y1) gives vector y1:" << endl;
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

	LinBox::sparsemat<Field, Row, Vector2> S2(B);
	LinBox::Blackbox_archetype<Vector2>& A2 = S2;
	
	A2.apply(y2,x2);

	cout << "Using A2.apply(y2,x2) gives vector y2:" << endl;
	for (y2_iter = y2.begin(); y2_iter != y2.end(); y2_iter++)
	{
		cout << "\ty2[" << y2_iter->first << "] = ";
		F.write(cout, y2_iter->second); 
		cout << endl;
	}

	y2 = Vector2();
	y2 = A2.apply(x2);

	cout << "Using y2 = A2.apply(x2) gives vector y2:" << endl;
	for (y2_iter = y2.begin(); y2_iter != y2.end(); y2_iter++)
	{
		cout << "\ty2[" << y2_iter->first << "] = ";
		F.write(cout, y2_iter->second); 
		cout << endl;
	}

	y2 = x2;
	y2 = A2.applyin(y2);

	cout << "Using y2 = x2; y2 = A2.applyin(y2) gives vector y2:" << endl;
	for (y2_iter = y2.begin(); y2_iter != y2.end(); y2_iter++)
	{
		cout << "\ty2[" << y2_iter->first << "] = ";
		F.write(cout, y2_iter->second); 
		cout << endl;
	}

	y2 = Vector2();
	A2.applyTranspose(y2,x2);

	cout << "Using A2.applyTranspose(y2,x2) gives vector y2:" << endl;
	for (y2_iter = y2.begin(); y2_iter != y2.end(); y2_iter++)
	{
		cout << "\ty2[" << y2_iter->first << "] = ";
		F.write(cout, y2_iter->second); 
		cout << endl;
	}

	y2 = Vector2();
	y2 = A2.applyTranspose(x2);

	cout << "Using y2 = A2.applyTranspose(x2) gives vector y2:" << endl;
	for (y2_iter = y2.begin(); y2_iter != y2.end(); y2_iter++)
	{
		cout << "\ty2[" << y2_iter->first << "] = ";
		F.write(cout, y2_iter->second); 
		cout << endl;
	}

	y2 = x2;
	y2 = A2.applyTransposein(y2);

	cout << "Using y2 = x2; y2 = A2.applyTransposein(y2) gives vector y2:" << endl;
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

	LinBox::sparsemat<Field, Row, Vector3> S3(B);
	LinBox::Blackbox_archetype<Vector3>& A3 = S3;
	
	A3.apply(y3,x3);

	cout << "Using A3.apply(y3,x3) gives vector y3:" << endl;
	for (y3_iter = y3.begin(); y3_iter != y3.end(); y3_iter++)
	{
		cout << "\ty3[" << y3_iter->first << "] = ";
		F.write(cout, y3_iter->second); 
		cout << endl;
	}

	y3 = Vector3();
	y3 = A3.apply(x3);

	cout << "Using y3 = A3.apply(x3) gives vector y3:" << endl;
	for (y3_iter = y3.begin(); y3_iter != y3.end(); y3_iter++)
	{
		cout << "\ty3[" << y3_iter->first << "] = ";
		F.write(cout, y3_iter->second); 
		cout << endl;
	}

	y3 = x3;
	y3 = A3.applyin(y3);

	cout << "Using y3 = x3; y3 = A3.applyin(y3) gives vector y3:" << endl;
	for (y3_iter = y3.begin(); y3_iter != y3.end(); y3_iter++)
	{
		cout << "\ty3[" << y3_iter->first << "] = ";
		F.write(cout, y3_iter->second); 
		cout << endl;
	}

	y3 = Vector3();
	A3.applyTranspose(y3,x3);

	cout << "Using A3.applyTranspose(y3,x3) gives vector y3:" << endl;
	for (y3_iter = y3.begin(); y3_iter != y3.end(); y3_iter++)
	{
		cout << "\ty3[" << y3_iter->first << "] = ";
		F.write(cout, y3_iter->second); 
		cout << endl;
	}

	y3 = Vector3();
	y3 = A3.applyTranspose(x3);

	cout << "Using y3 = A3.applyTranspose(x3) gives vector y3:" << endl;
	for (y3_iter = y3.begin(); y3_iter != y3.end(); y3_iter++)
	{
		cout << "\ty3[" << y3_iter->first << "] = ";
		F.write(cout, y3_iter->second); 
		cout << endl;
	}

	y3 = x3;
	y3 = A3.applyTransposein(y3);

	cout << "Using y3 = x3; y3 = A3.applyTransposein(y3) gives vector y3:" << endl;
	for (y3_iter = y3.begin(); y3_iter != y3.end(); y3_iter++)
	{
		cout << "\ty3[" << y3_iter->first << "] = ";
		F.write(cout, y3_iter->second); 
		cout << endl;
	}

}
