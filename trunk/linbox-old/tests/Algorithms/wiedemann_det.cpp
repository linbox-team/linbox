/* File:	tests/Algorithms/wiedemann_det.h
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
#include "LinBox/wiedemann_det.h"

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

	cout << "Input a square matrix of which to find the determinant:" << endl;
	LinBox::sparsemat<Field, Row, Vector>
		A(*LinBox::newSparsemat<Field, Row>(F));

	if (A.rowdim() != A.coldim())
	{
		cout << "The matrix is not square; no determinant found." << endl;
		return false;
	}

	cout << "The sparesemat matrix contains the following elements:" 
		<< endl << A;


	cout << "The determinant is ";
	F.write(cout, LinBox::wiedemann_det(F, A, r));
	cout << endl;

}
