/* File:	src/library/archetypes/vector/vaxpy.h
 * Author:	William J. Turner for the LinBox group
 */

#ifndef _VAXPY_
#define _VAXPY_

#include "LinBox/vector_traits.h"

// Namespace in which all LinBox library code lives
namespace LinBox
{

	/** Vector "axpy"
	 * a * x + y
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
	 * @return	vector containing sum of two vectors
	 * @param	F	Field in which arithmetic is done
	 * @param	a	Field element for multiplication
	 * @param	x	First of two vectors
	 * @param	y	Second vector
	 */
	template <class Field, class Vector>
	inline Vector vaxpy(
			const Field& F, 
			const typename Field::element& a,
			const Vector& x, 
			const Vector& y
			)
	{ 
		Vector z;
		return vaxpy(F, z, a, x, y); 
	}

	/** Vector "axpy"
	 * z = a * x + y.
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
	 * @return	vector containing sum of two vectors
	 * @param	F	Field in which arithmetic is done
	 * @param	z	vector to contain output
	 * @param	a	Field element for multiplication
	 * @param	x	First of two vectors
	 * @param	y	Second vector
	 */
	template <class Field, class Vector>
	inline Vector& vaxpy(
			const Field& F, 
			Vector& z, 
			const typename Field::element& a,
			const Vector& x, 
			const Vector& y
			)
	{ return vaxpy(F, z, a, x, y, vector_traits<Vector>::vector_category()); }

	/** Vector "axpyin"
	 * y += a * x.
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
	 * @return	vector containing sum of two vectors
	 * @param	F	Field in which arithmetic is done
	 * @param	y	vector to contain output
	 * @param	a	Field element for multiplication
	 * @param	x	vector
	 */
	template <class Field, class Vector>
	inline Vector& vaxpyin(
			const Field& F, 
			Vector& y, 
			const typename Field::element& a,
			const Vector& x
			)
	{ return vaxpyin(F, y, a, x, vector_traits<Vector>::vector_category()); }

	// Vector axpy for dense vectors
	template <class Field, class Vector>
	Vector& vaxpy(
			const Field& F, 
			Vector& z, 
			const typename Field::element& a,
			const Vector& x, 
			const Vector& y,
			vector_categories::dense_vector_tag tag
			);

	// Vector axpyin for dense vectors
	template <class Field, class Vector>
	Vector& vaxpyin(
			const Field& F, 
			Vector& y, 
			const typename Field::element& a,
			const Vector& x, 
			vector_categories::dense_vector_tag tag
			);

	// Vector axpy for sparse sequence vectors
	template <class Field, class Vector>
	Vector& vaxpy(
			const Field& F, 
			Vector& z, 
			const typename Field::element& a,
			const Vector& x, 
			const Vector& y,
			vector_categories::sparse_sequence_vector_tag tag
			);

	// Vector axpy for sparse sequence vectors
	template <class Field, class Vector>
	Vector& vaxpyin(
			const Field& F, 
			Vector& y, 
			const typename Field::element& a,
			const Vector& x, 
			vector_categories::sparse_sequence_vector_tag tag
			)
	{
		Vector z;
		return y = vaxpy(F, z, a, x, y, tag);
	}

	// Vector axpy for sparse associative vectors
	template <class Field, class Vector>
	inline Vector& vaxpy(
			const Field& F, 
			Vector& z, 
			const typename Field::element& a,
			const Vector& x, 
			const Vector& y,
			vector_categories::sparse_associative_vector_tag tag
			);

	// Vector axpy for sparse associative  vectors
	template <class Field, class Vector>
	inline Vector& vaxpyin(
			const Field& F, 
			Vector& y, 
			const typename Field::element& a,
			const Vector& x, 
			vector_categories::sparse_associative_vector_tag tag
			)
	{
		Vector z;
		return y = vaxpy(F, z, a, x, y, tag);
	}



} // namespace LinBox

// Vector axpy for dense vectors
template <class Field, class Vector>
Vector& LinBox::vaxpy(
		const Field& F, 
		Vector& z, 
		const typename Field::element& a,
		const Vector& x, 
		const Vector& y,
		vector_categories::dense_vector_tag tag
		)
{
#ifdef TRACE
	cout << "Called dense vector axpy" << endl;
#endif // TRACE

	typename Field::element zero;
	F.init(zero, 0);
	if (x.size() != y.size())
	{
		cerr << "Dimensions of dense vectors are not compatible in vaxpy"
			<< endl
			<< "   x.size() = " << x.size() << " and y.size() = " << y.size()
			<< endl;
		return z;
	}

	z = Vector(x.size(), zero);
	typename Vector::iterator z_iter;
	typename Vector::const_iterator x_iter, y_iter;
	
	z_iter = z.begin();
	y_iter = y.begin();
	for (x_iter = x.begin(); x_iter != x.end(); x_iter++, y_iter++, z_iter++)
		F.axpy(*z_iter, a, *x_iter, *y_iter);

	return z;

} // vaxpy(F, u, v, vector_categories::dense_vector_tag& tag)

// Vector axpyin for dense vectors
template <class Field, class Vector>
Vector& LinBox::vaxpyin(
		const Field& F, 
		Vector& y, 
		const typename Field::element& a,
		const Vector& x, 
		vector_categories::dense_vector_tag tag
		)
{
#ifdef TRACE
	cout << "Called dense vector axpyin" << endl;
#endif // TRACE

	typename Field::element zero;
	F.init(zero, 0);
	if (x.size() != y.size())
	{
		cerr << "Dimensions of dense vectors are not compatible in vaxpy"
			<< endl
			<< "   x.size() = " << x.size() << " and y.size() = " << y.size()
			<< endl;
		return y;
	}

	typename Vector::iterator y_iter;
	typename Vector::const_iterator x_iter;
	
	y_iter = y.begin();
	for (x_iter = x.begin(); x_iter != x.end(); x_iter++, y_iter++)
		F.axpyin(*y_iter, a, *x_iter);

	return y;

} // vaxpyin(vector_categories::dense_vector_tag& tag)

// Vector axpy for sparse sequence vectors
template <class Field, class Vector>
Vector& LinBox::vaxpy(
		const Field& F, 
		Vector& z, 
		const typename Field::element& a,
		const Vector& x, 
		const Vector& y,
		vector_categories::sparse_sequence_vector_tag tag
		)
{
#ifdef TRACE
	cout << "Called sparse sequence vector axpy" << endl;
#endif // TRACE

	// Check to see if a is the zero Field element.
	// If so, no addition is performed.
	if (F.isZero(a)) return z = y;

	// Check to see if x is empty.  If so, no addition is preformed.
	if ( x.begin() == x.end() ) return z = y;

	typename Field::element temp;
	F.init(temp, 0);
	z = Vector();
	typename Vector::const_iterator x_iter, y_iter(y.begin());

	// Iterate through first vector's nonzero elements
	for (x_iter = x.begin(); x_iter != x.end(); x_iter++)
	{
		// Find same element in second vector
		while (	(y_iter != y.end()) && (y_iter->first < x_iter->first) )
		{ 
			z.push_back(*y_iter);
			y_iter++; 
		}

		// Multiply corresponding elements and add
		if ( (y_iter != y.end()) && (y_iter->first == x_iter->first) )
		{
			if (!F.isZero(F.axpy(temp, a, x_iter->second, y_iter->second)))
				z.push_back(make_pair(x_iter->first, temp));
			y_iter++;
		}
		else
			z.push_back(make_pair(x_iter->first, 
						F.mul(temp, a, x_iter->second)));
      
	}

	while (y_iter != y.end())
	{ 
		z.push_back(*y_iter);
		y_iter++; 
	}

	return z;

} // vaxpy(vector_categories::sparse_sequence_vector_tag& tag)
 
// Vector axpy for sparse sequence vectors
/* The following code attempts to do the axpy inplace, however it is
 * not guaranteed to work because the C++ standard allows insertion 
 * and deletion to invalidate iterators.
 * 
template <class Field, class Vector>
Vector& LinBox::vaxpyin(
		const Field& F, 
		Vector& y, 
		const typename Field::element& a,
		const Vector& x, 
		vector_categories::sparse_sequence_vector_tag tag
		)
{
#ifdef TRACE
	cout << "Called sparse sequence vector axpyin" << endl;
#endif // TRACE

	// Check to see if a is the zero Field element.
	// If so, no addition is performed.
	if (F.isZero(a)) return y;

	// Check to see if x is empty.  If so, no addition is preformed.
	if ( x.begin() == x.end() ) return y;

//	LinBox::integer k;
	size_t k;
	typename Field::element temp;
	F.init(temp, 0);
	typename Vector::const_iterator x_iter, y_iter(y.begin()), iter;

	bool found(true);
	
	// Iterate through first vector's nonzero elements
	for (x_iter = x.begin(); x_iter != x.end(); x_iter++)
	{
		found = true;
		k = x_iter->first;	// mark current index

		// Find same element in second vector
		while (	(y_iter != y.end()) && (y_iter->first < x_iter->first) )
			y_iter++; 

		// Check if row j has element for column k.
		if ( ( y_iter == y.end() ) || ( y_iter->first != k ) ) 
			found = false;
		
		// If y contains element for index k, perform sum.
		// Otherwise, sum = a * x[k]
		if (found) 
		{
			if (F.isZero(F.axpyin(y_iter->second, a, x_iter->second)))
			{
				iter = y_iter++;
				y.erase(iter);
			}
			else
				y_iter++;
			
		}
		else
			y.insert(y_iter, make_pair(k, F.mul(temp,a,x_iter->second)));
		
	}

	return y;

} // vaxpyin(vector_categories::sparse_sequence_vector_tag& tag)

*/

// Vector axpy for sparse associative vectors
template <class Field, class Vector>
Vector& LinBox::vaxpy(
		const Field& F, 
		Vector& z, 
		const typename Field::element& a,
		const Vector& x, 
		const Vector& y,
		vector_categories::sparse_associative_vector_tag tag
		)
{
#ifdef TRACE
	cout << "Called sparse associative vector axpy" << endl;
#endif // TRACE

	// Check to see if a is the zero Field element.
	// If so, no addition is performed.
	if (F.isZero(a)) return z = y;

	// Check to see if x is empty.  If so, no addition is preformed.
	if ( x.begin() == x.end() ) return z = y;

	typename Field::element temp;
	F.init(temp, 0);
	z = Vector();
	typename Vector::const_iterator x_iter, y_iter(y.begin());

	// Iterate through first vector's nonzero elements
	for (x_iter = x.begin(); x_iter != x.end(); x_iter++)
	{
		// Find same element in second vector
		while (	(y_iter != y.end()) && (y_iter->first < x_iter->first) )
		{ 
			z.insert(*y_iter);
			y_iter++; 
		}

		// Multiply corresponding elements and add
		if ( (y_iter != y.end()) && (y_iter->first == x_iter->first) )
		{
			if (!F.isZero(F.axpy(temp, a, x_iter->second, y_iter->second)))
				z.insert(make_pair(x_iter->first, temp));
			y_iter++;
		}
		else
			z.insert(make_pair(x_iter->first, 
						F.mul(temp, a, x_iter->second)));
      
	}

	while (y_iter != y.end())
	{ 
		z.insert(*y_iter);
		y_iter++; 
	}

	return z;

} // vaxpy(vector_categories::sparse_associative_vector_tag& tag)
 
#endif // _VAXPY_

