/* File:	tests/Algorithms/wiedemann_minpoly.h
 * Author:	William J. Turner for the LinBox group
 *
 * Run tests on Wiedemann determinant algorithm
 */

#include <iostream>
#include <vector>
#include <list>
#include <utility>
#include "LinBox/giverror.C"
#include "LinBox/param_modular.h"
#include "LinBox/sparsemat.h"
#include "LinBox/wiedemann_minpoly.h"

int main(void)
{
	typedef LinBox::param_modular Field;
	typedef Field::element Element;
	typedef Field::randIter RandIter;
	typedef std::list< pair<size_t, Element> > Row;
	typedef std::vector<Element> Vector;

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

	cout << "Input a square matrix of which to find the minimum polynomial:" << endl;
	LinBox::sparsemat<Field, Row, Vector>
		A(*LinBox::newSparsemat<Field, Row>(F));

	if (A.rowdim() != A.coldim())
	{
		cout << "The matrix is not square; no minimum polynomial found." << endl;
		return false;
	}

	cout << "The sparesemat matrix contains the following elements:" 
		<< endl << A;

	std::vector<Element> minpoly = LinBox::wiedemann_minpoly(F, A, r);
	
	cout << "The minimum polynomial is" << endl << "\t";
	F.write(cout, minpoly[0]);
	for (size_t i = 1; i < minpoly.size(); i++)
	{
		cout << " + ";
		F.write(cout, minpoly[i]);
		cout << " * x^" << i;
	}
	cout << endl;

}
