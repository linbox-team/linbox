/* File: src/library/objects/blackbox/diagonal.h
 * Authro: William J Turner for the LinBox group
 */

#ifndef __DIAGONAL_H
#define __DIAGONAL_H

#include "LinBox/blackbox_archetype.h"
#include "LinBox/vector_traits.h"

// Namespace in which all LinBox library code resides
namespace LinBox
{

	/** Blackbox Diagonal matrix.
	 * This is a class of n by n diagonal matrices templatized by the 
	 * {@link Fields field} in 
	 * which the elements reside.  The class conforms to the 
	 * {@link Archetypes archetype} for \Ref{BlackBox Matrices}.
	 *
	 * The matrix itself is not stored in memory.  Rather, its apply
	 * methods use a vector of {@link Fields field} elements, which are 
	 * used to "multiply" the matrix to a vector.
	 * 
	 * This class has three template parameters.  The first is the field in 
	 * which the arithmetic is to be done.  The second is the type of 
	 * \Ref{LinBox} vector to which to apply the matrix.  The 
	 * third is chosen be defualt to be the \Ref{LinBox} vector trait
	 * of the vector.  This class is then specialized for dense and sparse 
	 * vectors.
	 * 
	 * The default class is not implemented.  It's functions should never
	 * be called because partial template specialization should always be
	 * done on the vector traits.
	 * @param Field \Ref{LinBox} field
	 * @param Vector \Ref{LinBox} dense or sparse vector of field elements
	 * @param Trait  Marker whether to use dense or sparse LinBox vector 
	 *               implementation.  This is chosen by a default parameter 
	 *               and partial template specialization.
	 */
	template <class Field, class Vector, class Trait = vector_traits<Vector>::vector_category>
	class Diagonal : public Blackbox_archetype<Vector>
	{
		public:

		/// Element definition
		typedef typename Field::element element;

		/** Constructor from field and dense vector of field elements.
		 * @param F	LinBox field in which to do arithmetic
		 * @param v LinBox dense vector of field elements to be used 
		 * 		as the diagonal of the matrix.
		 */
		Diagonal (const Field F, const std::vector<typename Field::element>& v);

#if 0
		/** Constructor from field and random iterator over the field
		 * @param F    LinBox field in which to do arithmetic
		 * @param iter Random iterator from which to get the diagonal elements
		 */
		Diagonal (const Field F, const RandomIterator &iter);
#endif

		/** Virtual constructor.
		 * Required because constructors cannot be virtual.
		 * Make a copy of the Blackbox_archetype object.
		 * Required by abstract base class.
		 * @return pointer to new blackbox object
		 */
		Blackbox_archetype<Vector>* clone() const;

		/** Application of BlackBox matrix.
		 * y= A*x.
		 * Requires one vector conforming to the \Ref{LinBox}
		 * vector {@link Archetypes archetype}.
		 * Required by abstract base class.
		 * @return reference to vector y containing output.
		 * @param  x constant reference to vector to contain input
		 */
		Vector& apply(Vector& y, const Vector& x) const;

		/** Application of BlackBox matrix transpose.
		 * y= transpose(A)*x.
		 * Requires one vector conforming to the \Ref{LinBox}
		 * vector {@link Archetypes archetype}.
		 * Required by abstract base class.
		 * Because the diagonal matrix is symmetric, this is the same as calling 
		 * the apply function.
		 * @return reference to vector y containing output.
		 * @param  x constant reference to vector to contain input
		 */
		Vector& applyTranspose(Vector& y, const Vector& x) const;

		/** Retreive row dimensions of BlackBox matrix.
		 * This may be needed for applying preconditioners.
		 * Required by abstract base class.
		 * @return integer number of rows of black box matrix.
		 */
		size_t rowdim(void) const;
    
		/** Retreive column dimensions of BlackBox matrix.
		 * Required by abstract base class.
		 * @return integer number of columns of black box matrix.
		 */
		size_t coldim(void) const;

	}; // template <Field, Vector> class Diagonal
 
	// Specialization of diagonal for LinBox dense vectors
	template <class Field, class Vector>
	class Diagonal<Field, Vector, vector_categories::dense_vector_tag>
		: public Blackbox_archetype<Vector>
	{
	    public:

		typedef typename Field::element element;
		Diagonal(const Field F, const std::vector<typename Field::element>& v);
		Blackbox_archetype<Vector>* clone() const 
			{ return new Diagonal(*this); }
		Vector& apply(Vector& y, const Vector& x) const;
		Vector& applyTranspose(Vector& y, const Vector& x) const { return apply(y, x); }
		size_t rowdim(void) const { return _n; } 
		size_t coldim(void) const { return _n; } 

	    private:

		// Field for arithmetic
		Field _F;

		// Number of rows and columns of square matrix.
		size_t _n;

		// STL vector of field elements used in applying matrix.
		std::vector<element> _v;
    
	}; // template <Field, Vector> class Diagonal<dense_vector_tag>
   
	// Specialization of diagonal for LinBox sparse sequence vectors
	template <class Field, class Vector>
	class Diagonal<Field, Vector, vector_categories::sparse_sequence_vector_tag>
		: public Blackbox_archetype<Vector>
	{
	    public:

		typedef typename Field::element element;
		Diagonal(const Field F, const std::vector<typename Field::element>& v);
		Blackbox_archetype<Vector>* clone() const 
			{ return new Diagonal(*this); }
		Vector& apply(Vector& y, const Vector& x) const;
		Vector& applyTranspose(Vector& y, const Vector& x) const { return apply(y, x); }
		size_t rowdim(void) const { return _n; } 
		size_t coldim(void) const { return _n; } 

	    private:

		// Field for arithmetic
		Field _F;

		// Number of rows and columns of square matrix.
		size_t _n;

		// STL vector of field elements used in applying matrix.
		std::vector<element> _v;
    
	}; // template <Field, Vector> class Diagonal<sparse_sequence_vector_tag>

	// Specialization of diagonal for LinBox sparse associative vectors
	template <class Field, class Vector>
	class Diagonal<Field, Vector, vector_categories::sparse_associative_vector_tag>
		: public Blackbox_archetype<Vector>
	{
	    public:

		typedef typename Field::element element;
		Diagonal(const Field F, const std::vector<typename Field::element>& v);
		Blackbox_archetype<Vector>* clone() const 
			{ return new Diagonal(*this); }
		Vector& apply(Vector& y, const Vector& x) const;
		Vector& applyTranspose(Vector& y, const Vector& x) const { return apply(y, x); }
		size_t rowdim(void) const { return _n; } 
		size_t coldim(void) const { return _n; } 

	    private:

		// Field for arithmetic
		Field _F;

		// Number of rows and columns of square matrix.
		size_t _n;

		// STL vector of field elements used in applying matrix.
		std::vector<element> _v;
    
	}; // template <Field, Vector> class Diagonal<sparse_associative_vector_tag>

	// Method implementations for dense vectors
 
	template <class Field, class Vector>
	inline Diagonal<Field, Vector, vector_categories::dense_vector_tag>
		::Diagonal(const Field F, const std::vector<typename Field::element>& v)
		: _F(F), _n(v.size()), _v(v) {}

	template <class Field, class Vector>
	inline Vector& Diagonal<Field, Vector, vector_categories::dense_vector_tag>
		::apply(Vector& y, const Vector& x) const
	{
		// Create zero vector to hold output
		element temp;
		_F.init (temp, 0);
		//Vector* y_ptr = new Vector (_n, temp);
 
		if (_n != x.size ()) {
			cerr << endl << "ERROR:  Input vector not of right size." << endl
			     << endl;
			return y;
		}
 
		// Create iterators for input, output, and stored vectors
		std::vector<element>::const_iterator v_iter;
		typename Vector::const_iterator x_iter;
		typename Vector::iterator y_iter;
 
		// Start at beginning of _v and x vectors
		v_iter = _v.begin ();
		x_iter = x.begin ();

		// Iterate through all three vectors, multiplying input and stored
		// vector elements to create output vector element.
		for (y_iter = y.begin ();
		     y_iter != y.end ();
		     y_iter++, v_iter++, x_iter++)
			_F.mul (*y_iter, *v_iter, *x_iter);
 
		return y;
	} // Vector& Diagonal<dense_vector_tag>::apply(Vector& y, const Vector&) const
  
	// Method implementations for sparse sequence vectors
 
	template <class Field, class Vector>
	inline Diagonal<Field, Vector, vector_categories::sparse_sequence_vector_tag>
		::Diagonal(const Field F, const std::vector<typename Field::element>& v)
		: _F(F), _n(v.size()), _v(v) {}

	template <class Field, class Vector>
	inline Vector &Diagonal<Field, Vector, vector_categories::sparse_sequence_vector_tag>
		::apply(Vector& y, const Vector& x) const
	{
		// Create zero vector to hold output
		//Vector* y_ptr = new Vector ();
 
		if ((!x.empty ()) && (_n < x.back ().first) ) {
			cerr << endl << "ERROR:  Input vector not of right size." << endl
			     << endl;
			return y;
		}
		y.clear (); // we'll overwrite using push_backs.

		// create field elements and size_t to be used in calculations
		size_t i;
		element zero, entry;
		_F.init (zero, 0);
		_F.init (entry, 0);

		// Create iterators for input and stored vectors
		std::vector<element>::const_iterator v_iter;
		typename Vector::const_iterator x_iter;
 
		// Start at beginning of _v vector
		v_iter = _v.begin ();
 
		// Iterator over indices of input vector.
		// For each element, multiply input element with corresponding element
		// of stored vector and insert non-zero elements into output vector
		for (x_iter = x.begin (); x_iter != x.end (); x_iter++) {
			i = (*x_iter).first;
			_F.mul (entry, *(v_iter + i), (*x_iter).second);
			if (!_F.isZero (entry)) y.push_back (make_pair (i, entry));
		} // for (x_iter = x.begin (); x_iter != x.end (); x_iter++)

		return y;
	} // Vector& Diagonal<sparse_sequence_vector_tag>::apply(Vector& y, const Vector&) const

	// Method implementations for sparse associative vectors
 
	template <class Field, class Vector>
	inline Diagonal<Field, Vector, vector_categories::sparse_associative_vector_tag>
		::Diagonal(const Field F, const std::vector<typename Field::element>& v)
		: _F(F), _n(v.size()), _v(v) {}

	template <class Field, class Vector>
	inline Vector& Diagonal<Field, Vector, vector_categories::sparse_associative_vector_tag>
		::apply(Vector& y, const Vector& x) const
	{
		// Create zero vector to hold output
 
		if ((!x.empty ()) && (_n < x.rbegin ()->first)) {
			cerr << endl << "ERROR:  Input vector not of right size." << endl
			     << endl;
			return y;
		}

		y.clear (); // we'll overwrite using inserts

		// create field elements and size_t to be used in calculations
		size_t i;
		element zero, entry;
		_F.init (zero, 0);
		_F.init (entry, 0);

		// Create iterators for input and stored vectors
		std::vector<element>::const_iterator v_iter;
		typename Vector::const_iterator x_iter;
 
		// Start at beginning of _v vector
		v_iter = _v.begin ();
 
		// Iterator over indices of input vector.
		// For each element, multiply input element with corresponding element
		// of stored vector and insert non-zero elements into output vector
		for (x_iter = x.begin (); x_iter != x.end (); x_iter++)
		{
			i = x_iter->first;
			_F.mul (entry, *(v_iter + i), (*x_iter).second);
			if (!_F.isZero (entry)) y.insert (y.end (), make_pair (i, entry));
		}

		return y;
	} // Vector& Diagonal<sparse_associative_vector_tag>::apply(...) const

} // namespace LinBox

#endif // __DIAGONAL_H
