/* File:	tests/Algorithms/wiedemann_linsolve1.h
 * Author:	William J. Turner for the LinBox group
 *
 * Run tests on Wiedemann algorithm for solving nonhomogeneous linear
 * equations
 */

#include <iostream>
#include <vector>
#include <list>
#include <utility>
#include "LinBox/giverror.C"
#include "LinBox/param_modular.h"
#include "LinBox/sparsemat.h"
#include "LinBox/wiedemann_linsolve1.h"

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

	cout << "Input a square matrix:" << endl;
	LinBox::sparsemat<Field, Row, Vector>
		S(*LinBox::newSparsemat<Field, Row>(F));

	if (S.rowdim() != S.coldim())
	{
		cout << "The matrix is not square; no solution found." << endl;
		return false;
	}

	size_t n = S.rowdim();

	cout << "Enter the right hand side vector, b." << endl
		<< "Input the vector by entering index and value." << endl
		<< "Remember matrices and vectors are indexed starting at 0." << endl
		<< "End with a index of -1." << endl;

	Element elem, zero;
	F.init(zero, 0);
	F.init(elem, 0);

	Vector b(n, zero), x(n, zero);

	Vector::iterator b_iter, x_iter;

	bool found;

	size_t i;

	while (cin >> i)
	{
		// return also if row index is not positive
		if(i == size_t(-1)) break; 
		
		F.read(cin, elem);

		if (i < n) b[i] = elem;

	}

	cout << "The sparesemat matrix contains the following elements:" 
		<< endl << S;

	cout << "Vector b:" << endl;
	i = 0;
	for (b_iter = b.begin(); b_iter != b.end(); b_iter++, i++)
	{
		cout << "\tb[" << i << "] = ";
		F.write(cout, *b_iter); 
		cout << endl;
	}

	LinBox::Blackbox_archetype<Vector>& A = S;

	LinBox::wiedemann_linsolve1(F, x, A, b, r);
	
	cout << "Using wiedemann_linsolve1(F, x, A, b, r) gives vector x:" << endl;
	i = 0;
	for (x_iter = x.begin(); x_iter != x.end(); x_iter++, i++)
	{
		cout << "\tx[" << i << "] = ";
		F.write(cout, *x_iter); 
		cout << endl;
	}

}
