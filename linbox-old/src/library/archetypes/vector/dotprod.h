/* File:	src/library/archetypes/vector/dotprod.h
 * Author:	William J. Turner for the LinBox group
 */

#ifndef _DOTPROD_
#define _DOTPROD_

#include "LinBox/vector_traits.h"

// Namespace in which all LinBox library code lives
namespace LinBox
{

	/** Vector Dot Product
	 * This template function takes the dot, or inner, product of two LinBox 
	 * vectors: transpose(u) * v.
	 * It is templatized by the field and the vector types being used.
	 * This function calls another function by the same name with an additional
	 * parameter of the vector category of the vector it is called with.
	 * This mechanism is used because functions cannot have partial template
	 * specializations like classes can.
	 * This new, extended function can be specialized for specific fields
	 * and vectors to allow for better performance.  For example, the function 
	 * may be specialized for a field of integers modulo a prime number so that
	 * the modular reduction is only done once instead of at every multiplication
	 * and addition as in the default implementation.
	 * @return	field element containing dot, or inner, product of two vectors
	 * @param	F	Field in which arithmetic is done
	 * @param	u	First of two vectors
	 * @param	v	Second vector
	 */
	template <class Field, class Vector>
		inline
		typename Field::element
		dotprod(const Field& F, const Vector& u, const Vector& v)
		{ return dotprod(F, u, v, vector_traits<Vector>::vector_category()); }

	/* Vector Dot Product for dense vectors
	 * This template function takes the dot, or inner, product of two LinBox 
	 * dense vectors: transpose(u) * v.
	 * It is templatized by the field and the vector types being used.
	 * This function can be specialized for specific fields
	 * and vectors to allow for better performance.  For example, the function 
	 * may be specialized for a field of integers modulo a prime number so that
	 * the modular reduction is only done once instead of at every multiplication
	 * and addition as in the default implementation.
	 * @return	field element containing dot, or inner, product of two vectors
	 * @param	F	Field in which arithmetic is done
	 * @param	u	First of two dense vectors
	 * @param	v	Second dense vector
	 * @param tag	category of vector obtained from vector trait
	 */
	template <class Field, class Vector>
		typename Field::element
		dotprod(
				const Field& F, 
				const Vector& u, 
				const Vector& v,
				vector_categories::dense_vector_tag tag
				)
		{
#ifdef TRACE
			cout << "Called dense vector dot product" << endl;
#endif // TRACE

			typename Field::element prod;
			F.init(prod, 0);

			if (u.size() != v.size())
			{
				cerr 
					<< "Dimensions of dense vectors are not compatible in dot product"
					<< endl
					<< "   u.size() = " << u.size() << " and v.size() = " << v.size()
					<< endl;
				return prod;
			}
			
			typename Field::element temp(prod);
			typename Vector::const_iterator u_iter, v_iter;
			
      v_iter = v.begin();
      for (u_iter = u.begin(); u_iter != u.end(); u_iter++, v_iter++)
        F.addin(prod, F.mul(temp, *u_iter, *v_iter));

			return prod;
		} // dotprod(F, u, v, vector_categories::dense_vector_tag& tag)

	/* Vector Dot Product for sparse sequence vectors
	 * This template function takes the dot, or inner, product of two LinBox 
	 * sparse sequence vectors: transpose(u) * v.
	 * It is templatized by the field and the vector types being used.
	 * This function can be specialized for specific fields
	 * and vectors to allow for better performance.  For example, the function 
	 * may be specialized for a field of integers modulo a prime number so that
	 * the modular reduction is only done once instead of at every multiplication
	 * and addition as in the default implementation.
	 * @return	field element containing dot, or inner, product of two vectors
	 * @param	F	Field in which arithmetic is done
	 * @param	u	First of two dense vectors
	 * @param	v	Second dense vector
	 * @param tag	category of vector obtained from vector trait
	 */
	template <class Field, class Vector>
		typename Field::element
		dotprod(
				const Field& F, 
				const Vector& u, 
				const Vector& v,
				vector_categories::sparse_sequence_vector_tag tag
				)
		{
#ifdef TRACE
			cout << "Called sparse sequence vector dot product" << endl;
#endif // TRACE

			typename Field::element prod;
			F.init(prod, 0);
			
			// Check if either input vector is empty.  If so, return zero.
			if ( ( u.begin() == u.end() ) || ( v.begin() == v.end() ) )
				return prod;

			typename Field::element temp(prod);
			typename Vector::const_iterator u_iter, v_iter(v.begin());

			// Iterate through first vector's nonzero elements
			for (u_iter = u.begin(); u_iter != u.end(); u_iter++)
			{
				// Find same element in second vector
				// All other elements are ignored as they would be multiplied by zero
				// Note: (*iter).first is used instead of iter->first because
				//       deque has problems with the latter.
				while (	(v_iter != v.end()) && ((*v_iter).first < (*u_iter).first) )
				{ v_iter++; }

				// Multiply corresponding elements and add to dot product
				if ( (v_iter != v.end()) && ((*v_iter).first == (*u_iter).first) )
					F.addin(prod, F.mul(temp, (*u_iter).second, (*v_iter).second));

			} // for (u_iter = u.begin(); u_iter != u.end(); u_iter++)

			return prod;

		} // dotprod(F, u, v, vector_categories::sparse_sequence_vector_tag& tag)
		  
	/* Vector Dot Product for sparse associative vectors
	 * This template function takes the dot, or inner, product of two LinBox 
	 * sparse associative vectors: transpose(u) * v.
	 * It is templatized by the field and the vector types being used.
	 * This function can be specialized for specific fields
	 * and vectors to allow for better performance.  For example, the function 
	 * may be specialized for a field of integers modulo a prime number so that
	 * the modular reduction is only done once instead of at every multiplication
	 * and addition as in the default implementation.
	 * @return	field element containing dot, or inner, product of two vectors
	 * @param	F	Field in which arithmetic is done
	 * @param	u	First of two dense vectors
	 * @param	v	Second dense vector
	 * @param tag	category of vector obtained from vector trait
	 */
	template <class Field, class Vector>
		typename Field::element
		dotprod(
				Field& F, 
				Vector& u, 
				Vector& v,
				vector_categories::sparse_associative_vector_tag tag
				)
		{
#ifdef TRACE
			cout << "Called sparse associative vector dot product" << endl;
#endif // TRACE

			typename Field::element prod;
			F.init(prod, 0);
			
			// Check if either input vector is empty.  If so, return zero.
			if ( u.empty() || v.empty() ) return prod;

			typename Field::element temp(prod);
			typename Vector::const_iterator u_iter, v_iter(v.begin());

			// Iterate through first vector's nonzero elements
			for (u_iter = u.begin(); u_iter != u.end(); u_iter++)
			{
				// Find same element in second vector
				// All other elements are ignored as they would be multiplied by zero
				while (	(v_iter != v.end()) && (v_iter->first < u_iter->first) )
				{ v_iter++; }

				// Multiply corresponding elements and add to dot product
				if ( (v_iter != v.end()) && (v_iter->first == u_iter->first) )
					F.addin(prod, F.mul(temp, u_iter->second, v_iter->second));

			} // for (u_iter = u.begin(); u_iter != u.end(); u_iter++)

			return prod;

		} // dotprod(F, u, v, vector_categories::sparse_associative_vector_tag& tag)
		  
} // namespace LinBox

#endif // _DOTPROD_

