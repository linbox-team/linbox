/* File:	src/library/algorithms/BlackBox/wiedemann_det.h
 * Author:	William J. Turner for the LinBox group
 *
 * This file contains an algorithm to compute the determinant of a
 * generic black box matrix using the Wiedemann algorithm.
 */

#ifndef _BERLEKAMP_MASSEY_
#define _BERLEKAMP_MASSEY_

#include <vector>

// Namespace in which all LinBox library code resides
namespace LinBox
{
	/** Berlekamp-Massey algorithm
	 * No early termination; obtains limit from length of input sequence.
	 * @return	STL vector of coefficients of minimal polynomial
	 * @param	F	Field for arithmetic
	 * @param	a	STL vector of sequence elements
	 */
	template<class Field>
	std::vector<typename Field::element>
	berlekamp_massey(Field& F, std::vector<typename Field::element>& a);
}

template<class Field>
std::vector<typename Field::element>
LinBox::berlekamp_massey(Field& F, std::vector<typename Field::element>& a)
{
	
	// Types

	typedef typename Field::element Element;

	// Constants and initialization

	Element zero, one;
	F.init(zero, 0);
	F.init(one, 1);
	Element temp(zero), scale(zero);

	long n = a.size();		// length of sequence
	std::vector<Element> Lambda(1, one), B, temp_Lambda;
	long L(0), n_r, s;
	Element Delta(one), Delta_r;

	// Loop over sequence

	for (long r = 1; r <= n; r++)
	{
		// Degree of polynomial Lambda
		n_r = Lambda.size() - 1;

		// Compute Delta_r
		F.assign(Delta_r, zero);
		for (long i = 0; i <= n_r; i++)
		F.addin(Delta_r, F.mul(temp, Lambda[i], a[r - i - 1]));

#ifdef TRACE
		cout << "*****" << endl;
		cout << "r = " << r << ", n_r = " << n_r << ", Delta_r = ";
		F.write(cout, Delta_r);
		cout << ", L = " << L << endl;
#endif // TRACE

		if (F.isZero(Delta_r)) 
			B.insert(B.begin(), zero);
		else if (2 * L < r)
		{
			temp_Lambda = Lambda;
			// Lambda -= Delta_r/Delta * x * B
			F.div(scale, Delta_r, Delta);
			s = B.size();
			for (long i = n_r+1; i <= s; i++)
				Lambda.push_back(zero);
			for (long i = 0; i < s; i++)
				F.subin(Lambda[i+1], F.mul(temp, scale, B[i]));
			B = temp_Lambda;
			L = r - L;
			F.assign(Delta, Delta_r);
		}
		else
		{
			// Lambda -= Delta_r/Delta * x * B
			F.div(scale, Delta_r, Delta);
			s = B.size();
			for (long i = n_r+1; i <= s; i++)
				Lambda.push_back(zero);
			for (long i = 0; i < s; i++)
				F.subin(Lambda[i+1], F.mul(temp, scale, B[i]));
			B.insert(B.begin(), zero);
		}

#ifdef TRACE
		cout << "Lambda = " << endl;
		for (unsigned long i = 0; i < Lambda.size(); i++)
		{
			cout << "  Lambda[" << i << "] = ";
			F.write(cout, Lambda[i]);
			cout << endl;
		}

		cout << "B = " << endl;
		for (unsigned long i = 0; i < B.size(); i++)
		{
			cout << "  B[" << i << "] = ";
			F.write(cout, B[i]);
			cout << endl;
		}

		cout << "L = " << L << endl;
		cout << "Delta = ";
		F.write(cout, Delta);
		cout << endl;

#endif // TRACE

	} // for (unsigned long i = 0; i < n; i++)

	// Return minimal polynomial = x^L * Lambda(1/x)

	std::vector<Element> minpoly(L+1, zero);
	s = Lambda.size();
	for (long i = 0; i < s; i++)
		F.assign(minpoly[L-i], Lambda[i]);

#ifdef TRACE
	cout << "minpoly = " << endl;
	for (unsigned long i = 0; i < minpoly.size(); i++)
	{
		cout << "  minpoly[" << i << "] = ";
		F.write(cout, minpoly[i]);
		cout << endl;
	}

#endif // TRACE

	return minpoly;

} // template<class Field> std::vector<...> berlekamp-massey(...)

#endif // _BERLEKAMP_MASSEY_
