/* File:	src/library/archetypes/vector/scalarprod.h
 * Author:	William J. Turner for the LinBox group
 */

#ifndef _SCALARPROD_
#define _SCALARPROD_

#include "LinBox/vector_traits.h"

// Namespace in which all LinBox library code lives
namespace LinBox
{

	/** Vector scalar product
	 * a * x
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
	 * @return	vector containing scalar product
	 * @param	F	Field in which arithmetic is done
	 * @param	a	Field element for multiplication
	 * @param	x	Vector
	 */
	template <class Field, class Vector>
	inline Vector scalarprod(
			const Field& F, 
			const typename Field::element& a,
			const Vector& x
			)
	{ 
		Vector y;
		return scalarprod(F, y, a, x); 
	}

	/** Vector scalar product
	 * y = a * x.
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
	 * @return	vector containing scalar product
	 * @param	F	Field in which arithmetic is done
	 * @param	y	vector to contain output
	 * @param	a	Field element for multiplication
	 * @param	x	Vector
	 */
	template <class Field, class Vector>
	inline Vector& scalarprod(
			const Field& F, 
			Vector& y, 
			const typename Field::element& a,
			const Vector& x
			)
	{ return scalarprod(F, y, a, x, vector_traits<Vector>::vector_category()); }

	/** Vector inplace scalar product
	 * x *= a.
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
	 * @return	vector containing scalar product
	 * @param	F	Field in which arithmetic is done
	 * @param	x	Vector
	 * @param	a	Field element for multiplication
	 */
	template <class Field, class Vector>
	inline Vector& scalarprodin(
			const Field& F, 
			Vector& x, 
			const typename Field::element& a
			)
	{ return scalarprodin(F, x, a, vector_traits<Vector>::vector_category()); }

	// Vector inplace scalar product for dense vectors
	template <class Field, class Vector>
	Vector& scalarprod(
			const Field& F, 
			Vector& y, 
			const typename Field::element& a,
			const Vector& x, 
			vector_categories::dense_vector_tag tag
			);

	// Vector scalar productin for dense vectors
	template <class Field, class Vector>
	Vector& scalarprodin(
			const Field& F, 
			Vector& x, 
			const typename Field::element& a,
			vector_categories::dense_vector_tag tag
			);

	// Vector inplace scalar product for sparse sequence vectors
	template <class Field, class Vector>
	Vector& scalarprod(
			const Field& F, 
			Vector& y, 
			const typename Field::element& a,
			const Vector& x, 
			vector_categories::sparse_sequence_vector_tag tag
			);

	// Vector scalar product for sparse sequence vectors
	template <class Field, class Vector>
	Vector& scalarprodin(
			const Field& F, 
			Vector& x, 
			const typename Field::element& a,
			vector_categories::sparse_sequence_vector_tag tag
			);

	// Vector scalar product for sparse associative vectors
	template <class Field, class Vector>
	inline Vector& scalarprod(
			const Field& F, 
			Vector& y, 
			const typename Field::element& a,
			const Vector& x, 
			vector_categories::sparse_associative_vector_tag tag
			)
	// same as sparse sequence vectors
	{ return scalarprod(F, y, a, x, vector_categories::sparse_sequence_vector_tag()); }

	// Vector inplace scalar product for sparse associative  vectors
	template <class Field, class Vector>
	inline Vector& scalarprodin(
			const Field& F, 
			Vector& x, 
			const typename Field::element& a,
			vector_categories::sparse_associative_vector_tag tag
			)
	// same as sparse sequence vectors
	{ return scalarprodin(F, x, a, vector_categories::sparse_sequence_vector_tag()); }

} // namespace LinBox

// Vector scalar product for dense vectors
template <class Field, class Vector>
Vector& LinBox::scalarprod(
		const Field& F, 
		Vector& y, 
		const typename Field::element& a,
		const Vector& x, 
		vector_categories::dense_vector_tag tag
		)
{
#ifdef TRACE
	cout << "Called dense vector scalar product" << endl;
#endif // TRACE

	if (F.isOne(a)) return y = x;

	typename Field::element zero;
	F.init(zero, 0);

	y = Vector(x.size(), zero);
	if (F.isZero(a)) return y;

	typename Vector::iterator y_iter;
	typename Vector::const_iterator x_iter;
	
	y_iter = y.begin();
	for (x_iter = x.begin(); x_iter != x.end(); x_iter++, y_iter++)
		F.mul(*y_iter, a, *x_iter);

	return y;

} // scalarprod(vector_categories::dense_vector_tag& tag)

// Vector inplace scalar productin for dense vectors
template <class Field, class Vector>
Vector& LinBox::scalarprodin(
		const Field& F, 
		Vector& x, 
		const typename Field::element& a,
		vector_categories::dense_vector_tag tag
		)
{
#ifdef TRACE
	cout << "Called dense vector scalar productin" << endl;
#endif // TRACE

	if (F.isOne(a)) return x;

	typename Field::element zero;
	F.init(zero, 0);

	if (F.isZero(a)) return x = Vector(x.size(), zero);

	typename Vector::iterator x_iter;
	
	for (x_iter = x.begin(); x_iter != x.end(); x_iter++)
		F.mulin(*x_iter, a);

	return x;

} // scalarprodin(vector_categories::dense_vector_tag& tag)

// Vector scalar product for sparse sequence vectors
template <class Field, class Vector>
Vector& LinBox::scalarprod(
		const Field& F, 
		Vector& y, 
		const typename Field::element& a,
		const Vector& x, 
		vector_categories::sparse_sequence_vector_tag tag
		)
{
#ifdef TRACE
	cout << "Called sparse sequence vector scalar product" << endl;
#endif // TRACE

	if (F.isOne(a)) return y = x;

	if (F.isZero(a)) return y = Vector();;

	y = x;
	
	typename Vector::iterator y_iter;
	
	for (y_iter = y.begin(); y_iter != y.end(); y_iter++)
		F.mulin(y_iter->second, a);

	return y;

} // scalarprod(vector_categories::sparse_sequence_vector_tag& tag)
 
// Vector inplace scalar product for sparse sequence vectors
template <class Field, class Vector>
Vector& LinBox::scalarprodin(
		const Field& F, 
		Vector& x, 
		const typename Field::element& a,
		vector_categories::sparse_sequence_vector_tag tag
		)
{
#ifdef TRACE
	cout << "Called sparse sequence vector scalar product" << endl;
#endif // TRACE

	if (F.isOne(a)) return x;

	if (F.isZero(a)) return x = Vector();;

	typename Vector::iterator x_iter;
	
	for (x_iter = x.begin(); x_iter != x.end(); x_iter++)
		F.mulin(x_iter->second, a);

	return x;

} // scalarprodin(vector_categories::sparse_sequence_vector_tag& tag)
 
#endif // _SCALARPROD_

