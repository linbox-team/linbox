/* File:	src/library/algorithms/BlackBox/wiedemann_det.h
 * Author:	William J. Turner for the LinBox group
 *
 * This file contains an algorithm to compute the determinant of a
 * generic black box matrix using the Wiedemann algorithm.
 */

#ifndef _WIEDEMANN_DET_
#define _WIEDEMANN_DET_

#include <vector>

#include "LinBox/blackbox_archetype.h"
#include "LinBox/diagonal.h"
#include "LinBox/compose.h"
#include "LinBox/vector_traits.h"
#include "LinBox/vector_traits.h"

#include "dotprod.h"
#include "random_vector.h"

// Namespace in which all LinBox library code resides
namespace LinBox
{
  /** Wiedemann black box determinant algorithm for a generic field.
   * This function is templatized by the field used.
	 * @return	field element containing determinant of matrix
	 * @param	F	Field in which arithmetic is done
	 * @param	BB	Black box matrix
	 * @param	r	random field element generator
   */
  template<class Field, class Vector>
  typename Field::element 
  wiedemann_det(Field& F,
		Blackbox_archetype<Vector>& BB,
		typename Field::randIter& r
		)
  {
    // Types

    typedef typename Field::element Element;
    typedef typename Field::randIter RandIter;

    Element zero, one;
    F.init(zero, 0);
    F.init(one, 1);

    // Dimensions

    integer n(BB.rowdim());
    
    // Precondition matrix with random diagonal matrix
    
    std::vector<Element> d(random_vector<Field,std::vector<Element> >(F, n, r)); 
		Element prod(one);
    for (long i = 0; i < n; i++)
			F.mulin(prod, d[i]);

		for (long k = 0; (k < 64) && (F.isZero(prod)); k++)
		{
			F.assign(prod, one);
			for (long i = 0; i < n; i++)
				F.mulin(prod, d[i]);
		}

    Diagonal<Field, Vector> D(F, d);	// preconditioning matrix
    Compose< Vector > DA(&D, &BB);	// preconditioned matrix
    
#ifdef TRACE
    cout << "The random elements of the diagonal matrix are:" << endl;
    for (long i = 0; i < n; i++)
    {
      cout << "i = " << i << ", \t D[i,i] = ";
      F.write(cout, d[i]);
      cout << endl;
    }
		cout << "Their product is ";
		F.write(cout, prod);
		cout << endl;
#endif // TRACE
    
    // Random vectors u and v

//random_vector<Field, Vector>(F, n, r);

    Vector u(random_vector<Field, Vector>(F, n, r)); 
    Vector v(random_vector<Field, Vector>(F, n, r)); 
		
    integer N = 2*n;
    std::vector<Element> a(N);

		F.init(a[0], 0);
		F.assign(a[0], dotprod(F, u, v));

    // Compute sequence u^T A^i v and minimal polynomial

    for (long i = 1; i < N; i++)
    {
      DA.applyin(v);		// (DA)^i * v
			F.init(a[i], 0);
			F.assign(a[i], dotprod(F, u, v));
    } // for (long i = 1; i < 2*n; i++)

#ifdef TRACE
    cout << "The scalar sequence is:" << endl;
    for (unsigned long i = 0; i < a.size(); i++)
    {
      cout << "i = " << i << ", \t a[i] = ";
      F.write(cout, a[i]);
      cout << endl;
    }
#endif // TRACE

    std::vector<Element> minpoly(berlekamp_massey(F, a));

		if (F.isZero(minpoly[0]))
			return zero;
		
		if (n != integer(minpoly.size() - 1))
		{
			cerr << "Bad projection in Wiedemann algorithm" << endl;
			return zero;
		}

		Element det(minpoly[0]);
		F.divin(det, prod);
		if (0 != (n % 2)) F.negin(det);

#ifdef TRACE
		cout << "The determinant is ";
		F.write(cout, det);
		cout << endl;
#endif // TRACE

    return det;

  } // wiedemann_det(F, BB, r)

  /* Berlekamp-Massey algorithm
   * No early termination; obtains limit from length of input sequence.
   * @return	STL vector of coefficients of minimal polynomial
   * @param	F	Field for arithmetic
   * @param	a	STL vector of sequence elements
   */
  template<class Field>
  std::vector<typename Field::element>
  berlekamp_massey(Field& F, std::vector<typename Field::element>& a)
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
      n_r = Lambda.size() - 1;	// degree of polynomial Lambda

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

} // namespace LinBox

#endif // _WIEDEMANN_DET_
