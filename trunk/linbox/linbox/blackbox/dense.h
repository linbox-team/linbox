/* -*- mode: c; style: linux -*- */

/* linbox/blackbox/dense-matrix1.h
 *
 * evolved from dense-matrix.h by -bds, Zhendong Wan
 */

#ifndef __DENSE_MATRIX_H
#define __DENSE_MATRIX_H

#include <iostream>
#include <vector>
#include <fstream>

#include "linbox/blackbox/archetype.h"
#include "linbox/field/dense-vector-domain.h"

namespace LinBox
{
	/** Blackbox dense matrix class.  This class is templatized by the
	 * {@link Fields field} in which the elements reside. 
	 * The matrix is stored as a one dimensional STL vector of
	 * the elements, by rows. The interface provides for iteration
	 * over rows and over columns.
	 *
	 * The class also conforms to the {@link Archetypes archetype} for
	 * \Ref{Blackbox Matrices}.
	 *
	 * @param Field \Ref{LinBox} field
	 */

	template <class Field>
	class DenseMatrix1
		: public BlackboxArchetype< std::vector<typename Field::Element> >
	{
	    public:
		typedef typename Field::Element        Element;
		typedef std::vector<Element>           Vector;
		typedef std::vector<Element>::iterator pointer;
      
		/** Constructor.
		 * @param  F the field of entries; passed so that arithmetic may be done on elements. 
		 * @param  m  row dimension
		 * @param  n  column dimension
		 */
		DenseMatrix1 (Field &F, size_t m, size_t n)
			: _F (F), _rep (m*n), _VD (F), _rows(m), _cols(n)
		{}

		/** Copy constructor
		 */
		DenseMatrix1 (const DenseMatrix1 &M)
			: _F (M._F), _rep (M._rep), _VD (M._F), _rows(M._rows), _cols(M._cols)
		{}

		/// Blackbox interface

		BlackboxArchetype<Vector> *clone () const 
		  { return new DenseMatrix1 (*this);}
		
		/* try later
		template<class Vect1, class Vect2>
		Vect1& apply (Vect1& y, const Vect2& x) const;
		*/
		Vector& apply (Vector& y, const Vector& x) const;
		 
		//template<class Vect1, class Vect2>
		//Vect1& applyTranspose (Vect1& y, const Vect2& x) const;
		Vector& applyTranspose (Vector& y, const Vector& x) const;

		size_t rowdim (void) const; // aka m
		
		size_t coldim (void) const; // aka n
	

		/** entry access raw view.  
		 *  Size m*n vector in row major order
		 *  (sometimes called C order).
		 */
		   
		typedef Vector::iterator RawIterator;
		typedef Vector::const_iterator ConstRawIterator;

		RawIterator rawBegin();		  
		RawIterator rawEnd();

		ConstRawIterator rawBegin() const;
		 
		ConstRawIterator rawEnd() const;

		/// col sequence of rows view
		typedef Vector::iterator RowIterator;
		typedef Vector::const_iterator ConstRowIterator;

		/// A Row has the static vector interface.
		class Row;

		typedef const Row ConstRow;    
		
		class ColOfRowsIterator;

		typedef const ColOfRowsIterator ConstColOfRowsIterator;
		
		ColOfRowsIterator colOfRowsBegin();

		ConstColOfRowsIterator colOfRowsBegin() const;

		ColOfRowsIterator colOfRowsEnd();

		ConstColOfRowsIterator colOfRowsEnd() const;

		// row sequence of cols view
		class ColIterator;
	      
		typedef const ColIterator ConstColIterator;

		/** A Col has most of the static vector interface.
		 *  The contiguity promise of the vector interface
		 *  is not provided.  Thus &C[0] is not an array.
		 *  The Col iterators or &C[0] with stride _len must be used.
		 */

		class Col;
		typedef const Col ConstCol;

		class RowOfColsIterator;

		typedef const RowOfColsIterator ConstRowOfColsIterator;

		RowOfColsIterator rowOfColsBegin();

		RowOfColsIterator rowOfColsEnd();

		ConstRowOfColsIterator rowOfColsBegin() const;
		
		ConstRowOfColsIterator rowOfColsEnd() const;

		/** Set the entry at (i, j)
		 * @param i Row number, 0...rowdim () - 1
		 * @param j Column number 0...coldim () - 1
		 * @param a_ij Element to set
		 */
		void setEntry (size_t i, size_t j, Element& a_ij) ;

		Element& getEntry (size_t i, size_t j, Element& a_ij) ;

		/** Read the matrix from an input stream
		 * @param file Input stream from which to read
		 */
		void read (std::istream &file);

		/** Write the matrix to an output stream
		 * @param os Output stream to which to write
		 */
		std::ostream &write(std::ostream &os = std::cout);

	    protected:

		std::vector<Element>                 _rep;
		Field                                _F;
		DenseVectorDomain<Field>             _VD;
		size_t _rows, _cols;
	};
}
#include "dense-matrix1.C"

#endif
