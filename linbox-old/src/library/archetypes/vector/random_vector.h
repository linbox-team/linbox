/* File:	src/library/archetypes/vector/random_vector.h
 * Author:	William J. Turner for the LinBox group
 */

#ifndef _RANDOM_VECTOR_
#define _RANDOM_VECTOR_

#include <utility>

#include "LinBox/vector_traits.h"

// Namespace in which all LinBox library code lives
namespace LinBox
{
	/** Random vector generator
	 * This templated function takes a field and a random field element
	 * generator and returns a vector of random field elements.
	 * The vector is dense in the field elements, even if the vector is
	 * a sparse LinBox vector.
	 * The funtion is templatized by the field and the vector types being used.
	 * This function calls another function by the same name with an additional
	 * parameter of the vector category of the vector it is called with.
	 * This mechanism is used because functions cannot have partial template
	 * specializations like classes can.
	 * This new, extended function can be specialized for specific fields
	 * and vectors to allow for better performance.
	 * @return	v	vector of random field elements
	 * @param	F	Field in which arithmetic is done
	 * @param	n	integer number of elements in vector
	 * @param	r	Random field element generator
	 */
	template <class Field, class Vector>
		inline
		Vector
		random_vector(Field& F, integer& n, typename Field::randIter& r)
		{ 
			return 
				random_vector<Field, Vector>(
						F, 
						n, 
						r, 
						vector_traits<Vector>::vector_category()
						); 
		}

	/* Random dense vector generator
	 * This templated function takes a field and a random field element
	 * generator and returns a vector of random field elements.
	 * The funtion is templatized by the field and the vector types being used.
	 * This function can be specialized for specific fields
	 * and vectors to allow for better performance.
	 * @return	v	vector of random field elements
	 * @param	F	Field in which arithmetic is done
	 * @param	n	integer number of elements in vector
	 * @param	r	Random field element generator
	 * @param tag	category of vector obtained from vector trait
	 */
	template <class Field, class Vector>
		inline
		Vector
		random_vector(
				Field& F, 
				integer& n,
				typename Field::randIter& r,
				vector_categories::dense_vector_tag tag
				)
		{
#ifdef TRACE
			cout << "Called dense random vector" << endl;
#endif // TRACE

			Vector v(n);
			typename Vector::iterator iter;

			for (iter = v.begin(); iter != v.end(); iter++)
				*iter = r();

			return v;
  
		} // random_vector(F, r, vector_categories::dense_vector_tag& tag)

	/* Random sparse sequence vector generator
	 * This templated function takes a field and a random field element
	 * generator and returns a vector of random field elements.
	 * The funtion is templatized by the field and the vector types being used.
	 * This function can be specialized for specific fields
	 * and vectors to allow for better performance.
	 * @return	v	vector of random field elements
	 * @param	F	Field in which arithmetic is done
	 * @param	n	integer number of elements in vector
	 * @param	r	Random field element generator
	 * @param tag	category of vector obtained from vector trait
	 */
	template <class Field, class Vector>
		inline
		Vector
		random_vector(
				Field& F, 
				integer& n,
				typename Field::randIter& r,
				vector_categories::sparse_sequence_vector_tag tag
				)
		{
#ifdef TRACE
			cout << "Called sparse sequence random vector" << endl
				<< "   return vector v = " << endl;
#endif // TRACE

			Vector v(n);
			typename Vector::iterator iter;
			size_t i = 0;
			for (iter = v.begin(); iter != v.end(); iter++, i++)
			{
				iter->first = i;
				iter->second = r();
#ifdef TRACE
				cout << "      v[" << iter->first << "] = ";
				F.write(cout, iter->second);
				cout << endl;
#endif // TRACE
			}

			return v;
  
		} // random_vector(F, r, vector_categories::sparse_sequence_vector_tag& tag)

	/* Random sparse associative vector generator
	 * This templated function takes a field and a random field element
	 * generator and returns a vector of random field elements.
	 * The funtion is templatized by the field and the vector types being used.
	 * This function can be specialized for specific fields
	 * and vectors to allow for better performance.
	 * @return	v	vector of random field elements
	 * @param	F	Field in which arithmetic is done
	 * @param	n	integer number of elements in vector
	 * @param	r	Random field element generator
	 * @param tag	category of vector obtained from vector trait
	 */
	template <class Field, class Vector>
		inline
		Vector
		random_vector(
				Field& F, 
				integer& n,
				typename Field::randIter& r,
				vector_categories::sparse_associative_vector_tag tag
				)
		{
#ifdef TRACE
			cout << "Called sparse associative random vector" << endl
				<< "   return vector v = " << endl;
#endif // TRACE

			Vector v;
			for (size_t i = 0; i < size_t(n); i++)
			{
				v.insert(make_pair(i, r()));
#ifdef TRACE
				cout << "      v[" << i << "] = ";
				F.write(cout, v[i]);
				cout << endl;
#endif // TRACE
			}

			return v;
  
		} // random_vector(F, r, vector_categories::sparse_associative_vector_tag& tag)

}

#endif // _RANDOM_VECTOR_
