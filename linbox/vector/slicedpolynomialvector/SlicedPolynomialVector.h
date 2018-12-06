#ifndef __LINBOX_matrix_SlicedPolynomialVector_SlicedPolynomialVector_H
#define __LINBOX_matrix_SlicedPolynomialVector_SlicedPolynomialVector_H

#include <vector>
#include "linbox/vector/blas-vector.h"
#include <givaro/modular.h>
#include <givaro/givpoly1dense.h>
#include <givaro/givpoly1denseops.inl>

namespace Linbox
{
	template <class _Field, class _Storage, class _VectorElement = double>
	class SlicedPolynomialVector
	/*Requirement: variable class _Field must be a representation of a field of order p^n
	 *and have functions that return p (characteristic()) and p^n(cardinality()).
	 *For example, GivaroGfq or LidiaGfq can be used.
	 */

        /*Description: this class implements representation of a vector over GF
	 *as a vector of length n of vectors over F (BlasVectors),
	 *which corresponds to representation of this vector as a polynomial of degree n-1 with vector-coefficients.
	 */
	{
	public:
		typedef _Field Field;
		typedef typename Field::Element Element;
		typedef _Storage Rep;
		typedef _VectorElement VectorElement;
		typedef Givaro::Modular<VectorElement, VectorElement> IntField;
		typedef SlicedPolynomialVector<Field, Rep, VectorElement> Self_t;
		typedef Givaro::Poly1Dom<IntField, Dense>::Rep polynomial;
	private:
		Field GF;
		IntField F;
	private:
		int n; //GF.cardinality() == p^n
		std::vector<BlasVector<IntField>> V;
	public:
		polynomial irreducible;
	private:
		polynomial& modulo(polynomial& g, polynomial&h, polynomial& f);
		bool rabin(int n, int q, int a, int b);
		void setIrreduciblePlynomial(int max_steps = 1000000);

						////////////////
		        			//Constructors//
						////////////////

	public:
		/*! Allocates a vector of new zero vectors of size 0 (shaped and ready). Irreducible polynomial is chosen randomly.
		 */
		SlicedPolynomialVector (const Field &BF);

		/*Allocates a vector of new vectors of size m (shaped and ready). Irreducible polynomial is chosen randomly.
		 */
		SlicedPolynomialVector (const Field &BF, const size_t &m);
		
		/*! Allocates a vector of new zero vectors of size 0 (shaped and ready).
		 */
		SlicedPolynomialVector (const Field &BF, polynomial& pp);

		/*Allocates a vector of new vectors of size m (shaped and ready).
		 */
		SlicedPolynomialVector (const Field &BF, const size_t &m, polynomial& pp);

						///////////////
						// Destructor//
						///////////////

	public:
		~SlicedPolynomialVector();

						////////////////////////
	          				//dimensions of vector//
						////////////////////////

	public:
		/*Get length of V.
		 * @returns length of V
		 */
		size_t length() const;
		
		/*Get the number of rows in a vector.
		 * @returns Number of rows in a vector
		 */
		size_t rowdim() const;
		
	                    			/////////////////
	                    			//return fields//
	                    			/////////////////
	
	public:
		const Field& fieldGF() const;

		const IntField& fieldF() const;

						/////////////////////////
		        			//functions for entries//
						/////////////////////////

	public:
		/* Set the entry of the m-th matrix-coefficient at the (k) position to a_mk.
		 * @param m matrix-coefficient number, 0...length() - 1
		 * @param k Row number 0...rowdim () - 1
		 * @param a_mk Element to set
		 */
		const VectorElement& setEntry (size_t m, size_t k, const VectorElement &a_mk);
		
	private:
		/* Get a writeable reference to the m-th matrix-coefficient at the (k) position.
		 * @param m matrix-coefficient number, 0...length() - 1
		 * @param k Row number 0...rowdim () - 1
		 * @returns Reference to matrix entry
		 */
		VectorElement &refEntry (size_t m, size_t k);

	public:
		/* Get a read-only reference to the m-th matrix-coefficient at the (k) position.
		 * @param m matrix-coefficient number, 0...length() - 1
		 * @param k Row number 0...rowdim () - 1
		 * @returns Const reference to matrix entry
		 */
		VectorElement &getEntry (size_t m, size_t k);
		
						/////////////////////////////////////
		                		//functions for matrix-coefficients//
						/////////////////////////////////////
						
	public:
		/* Set the m-th vector-coefficient to V_m.
		 * @param m vector-coefficient number, 0...length() - 1
		 * @param V_m matrix to set
		 */
		void setVectorCoefficient (size_t m, const BlasVector<IntField> &V_m) ;

	private:
		/* Get a writeable reference to the m-th matrix-coefficient.
		 * @param m matrix-coefficient number, 0...length() - 1
		 * @returns Reference to matrix-coefficent
		 */
		BlasVector<IntField> &refVectorCoefficient (size_t m) ;

	public:
		/** Get a read-only reference to the m-th matrix-coefficient
		 * @param m matrix-coefficient number, 0...length() - 1
		 * @returns Const reference to matrix-coefficent
		 */
		const BlasVector<IntField> &getVectorCoefficient (size_t m) const ;

						/////////
		                		//swaps//
						/////////

	public:
		/* Swap i1-th and i2-th rows of matrices.
		 * This is done inplace.
		 */
		void swapRows(size_t k1, size_t k2);
		
								//////////////////
		                				//input / output//
								//////////////////
	
	public:
		std::istream &read (std::istream &file);

		std::ostream &write (std::ostream &os) const;
	};
}

#include "SlicedPolynomialVector.inl"

#endif

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
