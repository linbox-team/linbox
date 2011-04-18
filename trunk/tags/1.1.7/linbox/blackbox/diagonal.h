/* linbox/blackbox/diagonal.h
 * Copyright (C) 1999-2001 William J Turner,
 *               2001 Bradford Hovinen
 *
 * Written by William J Turner <wjturner@math.ncsu.edu>,
 *            Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * ------------------------------------
 * Modified by Dmitriy Morozov <linbox@foxcub.org>. May 28, 2002.
 *
 * Added parametrization of VectorCategory tags by VectorTraits. See 
 * vector-traits.h for more details.
 * 
 * ------------------------------------
 *
 * See COPYING for license information.
 */
#ifndef __LINBOX_diagonal_H
#define __LINBOX_diagonal_H

#include <vector>
#include "linbox/vector/vector-traits.h"
#include "linbox/util/debug.h"
#include "linbox/linbox-config.h"
#include "linbox/field/hom.h"

// Namespace in which all LinBox library code resides
namespace LinBox
{

	/** 
	 * \brief Random diagonal matrices are used heavily as preconditioners.

	 * This is a class of n by n diagonal matrices templatized by the 
	 * field in 
	 * which the elements reside.  The class conforms to the 
	 * BlackboxArchetype.
	 *
	 * The matrix itself is not stored in memory.  Rather, its <tt>apply</tt>
	 * methods use a vector of field elements, which are 
	 * used to "multiply" the matrix to a vector.
	 * 
	 * This class has two template parameters.  The first is the field in 
	 * which the arithmetic is to be done.  
	 * The second is the vector trait indicating dense or 
	 * sparse vector interface (dense by default).
	 * This class is then specialized for dense and sparse vectors.
	 * 
	 * The default class is not implemented.  It's functions should never
	 * be called because partial template specialization should always be
	 * done on the vector traits.
	 * \ingroup blackbox
	 * @param Field \ref LinBox field
	 * @param Trait  Marker whether to use dense or sparse LinBox vector 
	 *               implementation.  This is chosen by a default parameter 
	 *               and partial template specialization.
	 */
	template <class Field,
		  class Trait = typename VectorTraits<typename LinBox::Vector<Field>::Dense>::VectorCategory>
	class Diagonal
	{
	    private:
		Diagonal () {}
	};
	
 
	/** diagonal.h linbox/blackbox/diagonal.h 
	\brief Specialization of Diagonal for application to dense vectors
	 */
	template <class _Field>
	class Diagonal<_Field, VectorCategories::DenseVectorTag>
	{
            typedef Diagonal<_Field, VectorCategories::DenseVectorTag> Self_t;
	    public:

		typedef _Field Field;
		typedef typename Field::Element    Element;

		/// \brief cstor from vector of elements
		Diagonal(const Field F, const std::vector<typename Field::Element>& v);

        // construct random nonsingular n by n diagonal matrix.
		Diagonal(const Field F, const size_t n, bool nonsing=true);

		Diagonal(const Field F, const size_t n, typename Field::RandIter& iter);

		template <class OutVector, class InVector>
		OutVector &apply (OutVector &y, const InVector &x) const;

		template <class OutVector, class InVector>
		OutVector &applyTranspose (OutVector &y, const InVector &x) const { return apply (y, x); }

		size_t rowdim(void) const { return _n; } 

		size_t coldim(void) const { return _n; } 

		void random();
		void randomNonsingular();

		/// \brief the field of the entries
		const Field& field() const{ return _F; }

                /** Get an entry and store it in the given value
                 * This form is more in the Linbox style and is provided for interface
                 * compatibility with other parts of the library
                 * @param x Element in which to store result
                 * @param i Row index
                 * @param j Column index
                 * @return Reference to x
                 */
            Element &getEntry (Element &x, size_t i, size_t j) const {
                return (i==j?_F.assign(x,_v[i]):_F.init(x,0));
            }
                    

            template<typename _Tp1>
            struct rebind
            { typedef Diagonal<_Tp1, VectorCategories::DenseVectorTag> other; 

                void operator() (other & Ap, const Self_t& A, const _Tp1& F) 
                    {
                        
                        Hom<typename Self_t::Field, _Tp1> hom(A.field(), F);

                        typename std::vector<typename _Tp1::Element>::iterator nit = Ap.getData().begin();
                        typename std::vector<Element>::const_iterator oit = A.getData().begin();
                        for( ; oit != A.getData().end() ; ++nit, ++oit)
                            hom.image (*nit, *oit);
                    }
 
            };

            template<typename _Tp1, typename _Vc1>
            Diagonal(const Diagonal<_Tp1,_Vc1>& D, const Field& F)
                    : _F(F), _n(D.rowdim()), _v(D.rowdim()) {
                typename Diagonal<_Tp1,_Vc1>::template rebind<Field>() (*this, D, F);
            }
            
            


            std::ostream& write(std::ostream& out) {	
                out << "diag(";
                for (typename std::vector<Element>::iterator p = _v.begin(); p != _v.end(); ++p)
                    _F.write(out, *p) << ", ";
                return out << "\b\b)";
            }

            const std::vector<Element>& getData() const { return _v; }
            std::vector<Element>& getData() { return _v; }

	    private:

		// Field for arithmetic
		Field _F;

		// Number of rows and columns of square matrix.
		size_t _n;

		// STL vector of field elements used in applying matrix.
		std::vector<Element> _v;
    
	}; // template <Field, Vector> class Diagonal<DenseVectorTag>
   
	// Specialization of diagonal for LinBox sparse sequence vectors
	/** 
	\brief Specialization of Diagonal for application to sparse sequence vectors
	 */
	template <class _Field>
	class Diagonal<_Field, VectorCategories::SparseSequenceVectorTag >
	{
            typedef Diagonal<_Field, VectorCategories::SparseSequenceVectorTag > Self_t;
	    public:

		typedef _Field Field;
		typedef typename Field::Element    Element;

		Diagonal(const Field F, const std::vector<typename Field::Element>& v);

		Diagonal(const Field F, const size_t n, typename Field::RandIter& iter);
		
		template<class OutVector, class InVector>
		OutVector& apply(OutVector& y, const InVector& x) const;

		template<class OutVector, class InVector>
		OutVector& applyTranspose(OutVector& y, const InVector& x) const { return apply(y, x); }

		size_t rowdim(void) const { return _n; } 
		size_t coldim(void) const { return _n; } 
		const Field& field() const {return _F;}
                /** Get an entry and store it in the given value
                 * This form is more in the Linbox style and is provided for interface
                 * compatibility with other parts of the library
                 * @param x Element in which to store result
                 * @param i Row index
                 * @param j Column index
                 * @return Reference to x
                 */
            Element &getEntry (Element &x, size_t i, size_t j) const {
                return (i==j?_F.assign(x,_v[i]):_F.init(x));
            }
                    


            template<typename _Tp1>
            struct rebind
            { typedef Diagonal<_Tp1, VectorCategories::SparseSequenceVectorTag> other; 
                void operator() (other & Ap, const Self_t& A, const _Tp1& F) 
                    {

                        Hom<typename Self_t::Field, _Tp1> hom(A.field(), F);

                        typename std::vector<typename _Tp1::Element>::iterator nit = Ap.getData().begin();
                        typename std::vector<Element>::const_iterator oit = A.getData().begin();
                        for( ; oit != A.getData().end() ; ++nit, ++oit)
                            hom.image (*nit, *oit);
                    }
 
            };

            template<typename _Tp1, typename _Vc1>
            Diagonal(const Diagonal<_Tp1,_Vc1>& D, const Field& F)
                    : _F(F), _n(D.rowdim()), _v(D.rowdim()) {
                typename Diagonal<_Tp1,_Vc1>::template rebind<Field>() (*this, D, F);
            }
            
            

            const std::vector<Element>& getData() const { return _v; }
            std::vector<Element>& getData() { return _v; }


	    private:

		// Field for arithmetic
		Field _F;

		// Number of rows and columns of square matrix.
		size_t _n;

		// STL vector of field elements used in applying matrix.
		std::vector<Element> _v;
    
	}; // template <Field, Vector> class Diagonal<SparseSequenceVectorTag>

	// Specialization of diagonal for LinBox sparse associative vectors
	/** 
	\brief Specialization of Diagonal for application to sparse associative vectors
	 */
	template <class _Field>
	class Diagonal<_Field, VectorCategories::SparseAssociativeVectorTag >
	{
            typedef Diagonal<_Field, VectorCategories::SparseAssociativeVectorTag > Self_t;
	    public:


		typedef _Field Field;
		typedef typename Field::Element    Element;

		Diagonal(const Field F, const std::vector<typename Field::Element>& v);

		Diagonal(const Field F, const size_t n, typename Field::RandIter& iter);

		template<class OutVector, class InVector>
		OutVector& apply(OutVector& y, const InVector& x) const;

		template<class OutVector, class InVector>
		OutVector& applyTranspose(OutVector& y, const InVector& x) const { return apply(y, x); } 


		size_t rowdim(void) const { return _n; } 
		size_t coldim(void) const { return _n; } 
		const Field field() const { return _F; }
                /** Get an entry and store it in the given value
                 * This form is more in the Linbox style and is provided for interface
                 * compatibility with other parts of the library
                 * @param x Element in which to store result
                 * @param i Row index
                 * @param j Column index
                 * @return Reference to x
                 */
            Element &getEntry (Element &x, size_t i, size_t j) const {
                return (i==j?_F.assign(x,_v[i]):_F.init(x));
            }
                    


            template<typename _Tp1>
            struct rebind
            { typedef Diagonal<_Tp1, VectorCategories::SparseAssociativeVectorTag> other; 
                void operator() (other & Ap, const Self_t& A, const _Tp1& F) 
                    {
                        Hom<typename Self_t::Field, _Tp1> hom(A.field(), F);

                        typename std::vector<typename _Tp1::Element>::iterator nit = Ap.getData().begin();
                        typename std::vector<Element>::const_iterator oit = A.getData().begin();
                        for( ; oit != A.getData().end() ; ++nit, ++oit)
                            hom.image (*nit, *oit);
                    }
 
            };

            template<typename _Tp1, typename _Vc1>
            Diagonal(const Diagonal<_Tp1,_Vc1>& D, const Field& F)
                    : _F(F), _n(D.rowdim()), _v(D.rowdim()) {
                typename Diagonal<_Tp1,_Vc1>::template rebind<Field>() (*this, D, F);
            }
            
            const std::vector<Element>& getData() const { return _v; }
            std::vector<Element>& getData() { return _v; }

            

	    private:

		// Field for arithmetic
		Field _F;

		// Number of rows and columns of square matrix.
		size_t _n;

		// STL vector of field elements used in applying matrix.
		std::vector<Element> _v;
    
	}; // template <Field, Vector> class Diagonal<SparseAssociativeVectorTag>


	// Method implementations for dense vectors
 
	template <class Field>
	inline Diagonal<Field, VectorCategories::DenseVectorTag >
		::Diagonal(const Field F, const std::vector<typename Field::Element>& v)
		: _F(F), _n(v.size()), _v(v)
	{}


	template <class _Field>
	inline Diagonal<_Field, VectorCategories::DenseVectorTag>
		::Diagonal(const Field F, const size_t n, bool nonsing)
	: _F(F), _n(n), _v(n)
	{   typename Field::RandIter r(F);
		typedef typename std::vector<typename Field::Element>::iterator iter;
		if (nonsing)
		randomNonsingular();
		//for (iter i = _v.begin(); i < _v.end(); ++i) 
		//	while (_F.isZero(r.random(*i))) ;
		else
		random();
		//for (iter i = _v.begin(); i < _v.end(); ++i) 
		//	r.random(*i);
	}


	template <class Field>
	inline Diagonal<Field, VectorCategories::DenseVectorTag >
		::Diagonal(const Field F, const size_t n, typename Field::RandIter& iter)
		: _F(F), _n(n), _v(n)
	{	//for (typename std::vector<typename Field::Element>::iterator 
		//		i = _v.begin(); i != _v.end(); ++i) 
		//	iter.random(*i); 
		random();
	}

	template <class _Field>
	inline void Diagonal<_Field, VectorCategories::DenseVectorTag>
		::random()
	{   typename Field::RandIter r(_F);
		typedef typename std::vector<typename Field::Element>::iterator iter;
		for (iter i = _v.begin(); i < _v.end(); ++i) 
			r.random(*i);
	}

	template <class _Field>
	inline void Diagonal<_Field, VectorCategories::DenseVectorTag>
		::randomNonsingular()
	{   typename Field::RandIter r(_F);
		typedef typename std::vector<typename Field::Element>::iterator iter;
		for (iter i = _v.begin(); i < _v.end(); ++i) 
			while (_F.isZero(r.random(*i))) ;
	}

	template <class Field>
	template <class OutVector, class InVector>
	inline OutVector &Diagonal<Field, VectorCategories::DenseVectorTag >
		::apply (OutVector &y, const InVector &x) const
	{
		linbox_check (_n == x.size ());
 
		// Create iterators for input, output, and stored vectors
		typename std::vector<Element>::const_iterator v_iter;
		typename InVector::const_iterator x_iter;
		typename OutVector::iterator y_iter;
 
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
	} // Vector& Diagonal<DenseVectorTag>::apply(Vector& y, const Vector&) const
  
	// Method implementations for sparse sequence vectors
 
	template <class Field>
	inline Diagonal<Field, VectorCategories::SparseSequenceVectorTag >
		::Diagonal(const Field F, const std::vector<typename Field::Element>& v)
		: _F(F), _n(v.size()), _v(v)
	{}

	template <class Field>
	template<class OutVector, class InVector>
	inline OutVector &Diagonal<Field, VectorCategories::SparseSequenceVectorTag >
		::apply(OutVector& y, const InVector& x) const
	{
		linbox_check ((!x.empty ()) && (_n >= x.back ().first));

		y.clear (); // we'll overwrite using push_backs.

		// create field elements and size_t to be used in calculations
		size_t i;
		Element zero, entry;
		_F.init (zero, 0);
		_F.init (entry, 0);

		// Create iterators for input and stored vectors
		typename std::vector<Element>::const_iterator v_iter;
		typename InVector::const_iterator x_iter;
 
		// Start at beginning of _v vector
		v_iter = _v.begin ();
 
		// Iterator over indices of input vector.
		// For each element, multiply input element with corresponding element
		// of stored vector and insert non-zero elements into output vector
		for (x_iter = x.begin (); x_iter != x.end (); x_iter++) {
			i = (*x_iter).first;
			_F.mul (entry, *(v_iter + i), (*x_iter).second);
			if (!_F.isZero (entry)) y.push_back ( std::pair<size_t, Element>(i, entry));
		} // for (x_iter = x.begin (); x_iter != x.end (); x_iter++)

		return y;
	} // Vector& Diagonal<SparseSequenceVectorTag>::apply(Vector& y, const Vector&) const

	// Method implementations for sparse associative vectors
 
	template <class Field>
	inline Diagonal<Field, VectorCategories::SparseAssociativeVectorTag >
		::Diagonal(const Field F, const std::vector<typename Field::Element>& v)
		: _F(F), _n(v.size()), _v(v)
	{}

	template <class Field>
	template<class OutVector, class InVector>
	inline OutVector& Diagonal<Field, VectorCategories::SparseAssociativeVectorTag >
		::apply(OutVector& y, const InVector& x) const
	{
		linbox_check ((!x.empty ()) && (_n >= x.rbegin ()->first));

		y.clear (); // we'll overwrite using inserts

		// create field elements and size_t to be used in calculations
		size_t i;
		Element zero, entry;
		_F.init (zero, 0);
		_F.init (entry, 0);

		// Create iterators for input and stored vectors
		typename std::vector<Element>::const_iterator v_iter;
		typename InVector::const_iterator x_iter;
 
		// Start at beginning of _v vector
		v_iter = _v.begin ();
 
		// Iterator over indices of input vector.
		// For each element, multiply input element with corresponding element
		// of stored vector and insert non-zero elements into output vector
		for (x_iter = x.begin (); x_iter != x.end (); x_iter++)
		{
			i = x_iter->first;
			_F.mul (entry, *(v_iter + i), (*x_iter).second);
			if (!_F.isZero (entry)) y.insert (y.end (), std::pair<size_t, Element>(i, entry));
		}

		return y;
	} // Vector& Diagonal<SparseAssociativeVectorTag>::apply(...) const

	//@}



} // namespace LinBox

#endif // __LINBOX_diagonal_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
