/* linbox/blackbox/triplesbb.h
 * Copyright (C) 2002 Rich Seagraves,  see COPYING for details.
 *
 * Written by Rich Seagraves <seagrave@cis.udel.edu>
 * with mods by bds
 */

#ifndef __LINBOX_triplesbb_H
#define __LINBOX_triplesbb_H

#include <algorithm>
using std::max;
#include "linbox/linbox-config.h"
#include "linbox/util/debug.h"
#include "linbox/util/field-axpy.h"
#include <linbox/blackbox/blackbox-interface.h>
#include <linbox/field/hom.h>

#include <vector>

namespace LinBox 
{

	/** \brief wrapper for NAG Sparse Matrix format.
	 *
\ingroup blackbox
	 * This class acts as a wrapper for a pre-existing NAGSparse Matrix.
	 * To be used for interface between LinBox and computer algebra systems such
	 * as Maple that can encode sparse matrices in the NAGSparse format
	 */ 

	template<class _Field>
		class TriplesBB : public BlackboxInterface{	 
	 
		public:
		typedef _Field Field;
		typedef typename Field::Element Element;
	    typedef TriplesBB<Field> Self_t;


		// Default constructor.
		TriplesBB() {}

		// Takes 3 vectors and copies(bad) them.
		TriplesBB(	Field F, 
					std::vector<Element> values, 
					std::vector<size_t> rowP, 
					std::vector<size_t> colP, 
					size_t rows, 
					size_t cols, 
					bool RowSortFlag = false, 
					bool ColSortFlag = false);

		// Alternate constructor.  Allows for use of addEntry operation.
		TriplesBB(Field F, size_t rows, size_t cols, size_t reserve = 0);

		~TriplesBB() {};

		// Copy Constructor
		TriplesBB(const TriplesBB<Field> &);

		// Assignment operator for use in STL map
		const TriplesBB<Field> & operator=(const TriplesBB<Field> & );

		template<class OutVector, class InVector>
		OutVector & apply(OutVector &, const InVector &) const; // y = Ax;

		template<class OutVector, class InVector>
		OutVector & applyTranspose(OutVector &, const InVector &) const; // y = ATx
		size_t rowdim() const { return _rows; }

		size_t coldim() const { return _cols; }

        template<typename _Tp1> 
        struct rebind { 
            typedef TriplesBB<_Tp1> other; 
            void operator() (other & Ap, const Self_t& A, const _Tp1& F)
                {
                    Hom <typename Self_t::Field, _Tp1> hom( A.field(), F);
                    
                    typedef typename _Tp1::Element otherElt;
                    typedef typename std::vector<otherElt> othervec;
                    typedef typename std::vector<Element> selfvec;
                    typedef typename othervec::iterator otheriter;
                    typedef typename selfvec::const_iterator selfiter;
                    otheriter vp_p; selfiter v_p;

                    Ap._values.resize(A._values.size());
                    for (v_p = A._values.begin(), vp_p = Ap._values.begin();
                         v_p != A._values.end(); ++ v_p, ++ vp_p)
                        hom.image (*vp_p, *v_p);
                }
        };

        template<typename _Tp1> 
        TriplesBB(const TriplesBB<_Tp1>& T, const Field& F) 
                : _F(F), _values(T.size()), _RowV(T.getRows()), _ColV(T.getCols()), _rows(T.rowdim()), _cols(T.coldim()), _faxpy(max(T.getRows(),T.getCols()), FieldAXPY<Field>(F)), _RowSortFlag(T.isRowSorted()), _ColSortFlag(T.isColSorted())
		{}



		/* Returns number of non-zero entries */
		size_t size() const { return _values.size(); }

		// Add entry function, element e is added in the i,j position.  Zero based?
		void addEntry(const Element & e, const size_t i, const size_t j);

		const Field & field() const { return _F; }

		/* Data accessors.  Used to access the 3 vectors containing Matrix data
		 */
		const std::vector<Element> & getData() const { return _values; }
		const std::vector<size_t> & getRows() const { return _RowV; }
		const std::vector<size_t> & getCols() const { return _ColV; }
                    bool isRowSorted() { return _RowSortFlag; }
                    bool isColSorted() { return _ColSortFlag; }
                    

		protected:
		Field _F; // The field used by this class

		/// _values contains the nonzero elements of the BlackBox
		std::vector<Element> _values;

		/// _RowV & _ColV are vectors containing the row & column indices 
		std::vector<size_t> _RowV, _ColV;

		/// The number of rows, columns 
		size_t _rows, _cols;

		/* _apply is the generic apply utility funtion called by apply() and
		 * applyTranspose().  Walks through the non-zero elements of the Matrix and
		 * performs the proper calculation using a vector of FieldAxpy's
		 */
		template<class OutVector, class InVector>
		void _apply(OutVector &, const InVector &, std::vector<size_t>::const_iterator, std::vector<size_t>::const_iterator) const;

		/* STL vector of FieldAXPY objects.  Supports delayed modding out, a feature
		 * which contributes a significant speed boost when performing apply &
		 * applyTranspose calculations over a field of multi-precision integers
		 */

		mutable std::vector<FieldAXPY<Field> > _faxpy;

		/* Sort flag.  Used by the sort function to determine whether a sort
		 * operation is needed.  Also used by _apply for a slightly optimized
		 * apply operation.  Note, "sorted" is considered sorted if the
		 * matrix is row-sorted, IE if the entries go row 1, row 1, row 1, row 2, row 2, etc
		 */
		bool _RowSortFlag, _ColSortFlag;


	};

	/*  Constructor for the TriplesBB class.  This is the constructor that is
	 * expected to be used.  To use it, you must pass in a field element that
	 * will work over the data (F), pointers to the 3 arrays used by the NAGSparse
	 * format (values, rowP, colP), the number of rows and columns (rows and
	 * cols), the number of non-zero elements (NNz) and the ordering, which
	 * defaults to 0 (no ordering implied).
	 */
	template<class Field>
		TriplesBB<Field>::TriplesBB(Field F, 
								    std::vector<Element> values, 
								    std::vector<size_t> RowV, 
									std::vector<size_t> ColV, 
									size_t rows, 
									size_t cols, 
									bool RowSortFlag, 
									bool ColSortFlag) :
		_F(F), _values(values), _RowV(RowV), _ColV(ColV), _rows(rows), _cols(cols), _faxpy(max(rows,cols), FieldAXPY<Field>(F)), _RowSortFlag(RowSortFlag), _ColSortFlag(ColSortFlag)
		{}

	/* Better constructor that only takes the field, m, n and recommended
	 * reserve (optional arguement) for use with STL vector reserve option
	 * (reduce the number of memory management moves).  Meant to be used in
	 * conjuction with the addEntry() method
	 */
	template<class Field>
		TriplesBB<Field>::TriplesBB( Field F, size_t rows, size_t cols, size_t res):
		_F(F), _rows(rows), _cols(cols), _faxpy( max(rows, cols), FieldAXPY<Field>(F)), _RowSortFlag(false), _ColSortFlag(false)
		{
			if(res != 0) {
				_values.reserve(res);
				_RowV.reserve(res);
				_ColV.reserve(res);
			}
		}



	template<class Field>
		TriplesBB<Field>::TriplesBB(const TriplesBB<Field> &In) :
		_faxpy( max(In._rows, In._cols), FieldAXPY<Field>(In._F)),
			_F ( In._F ),
			_values ( In._values ),
			_RowV ( In._RowV ),
			_ColV ( In._ColV ),
			_rows ( In._rows ), 
			_cols ( In._cols ),
			_RowSortFlag ( In._RowSortFlag ),
			_ColSortFlag ( In._ColSortFlag )
		{ }


	template<class Field>
		const TriplesBB<Field> & TriplesBB<Field>::operator=(const TriplesBB<Field> & rhs)
		{
			_F = rhs._F;
			_values = rhs._values;
			_RowV = rhs._RowV;
			_ColV = rhs._ColV;
			_rows = rhs._rows; _cols = rhs._cols;
			_RowSortFlag = rhs._RowSortFlag;
			_ColSortFlag  = rhs._ColSortFlag;

			_faxpy.resize(rhs._faxpy.size(), FieldAXPY<Field>(_F));

			return *this;
		}


	template<class Field>
		template<class OutVector, class InVector>
		OutVector & TriplesBB<Field>::apply(OutVector & y, const InVector & x) const
		{

			_apply( y, x, _RowV.begin(), _ColV.begin() );
			return y;
		}

	/* BlackBoxArchetype applyTranspose function.  Performs the y = ATx, where
	 * y and x are vectors passed in applyTranspose(y,x), and A is the present
	 * Matrix.  Returns a reference to y.  As this is a tranpose calculation,
	 * the indexing is reversed, so y is indexed by the columns, while x is indexed
	 * by the rows.  Thus, as in apply above, takes advantage of this fact by
	 * switching on the ordering.
	 */

	template<class Field>
		template<class OutVector, class InVector>
		OutVector & TriplesBB<Field>::applyTranspose(OutVector & y, const InVector & x) const
		{
			_apply( y, x, _ColV.begin(), _RowV.begin() );
			return y;
		}


	template<class Field>
		template<class OutVector, class InVector>
		void TriplesBB<Field>::_apply(OutVector & y, const InVector & x, std::vector<size_t>::const_iterator i, std::vector<size_t>::const_iterator j) const
		{
			typename OutVector::iterator yp;
			typename InVector::const_iterator xp;
			typename Field::Element zero;
			typename std::vector<Element>::const_iterator v;
			typename std::vector<FieldAXPY<Field> >::iterator fa_i;

			_F.init(zero,0);

			for(fa_i = _faxpy.begin(); fa_i != _faxpy.end(); ++fa_i) 
				fa_i->assign(zero);

			for( v = _values.begin(), fa_i = _faxpy.begin() - 1, xp = x.begin() - 1; v != _values.end(); ++i, ++j, ++v) 
				(fa_i + *i)->mulacc(*v,  *(xp + *j));



			for(fa_i = _faxpy.begin(), yp = y.begin(); yp != y.end(); ++yp, ++fa_i) 
				fa_i->get(*yp);

  

		}


	/* addEntry method.  Allows user to add entries on the fly.  Meant to be used
	 * with the "copyless" constructor above.  Note, will automatically set the
	 * _sortFlag to false, as you can't be sure the entries are still sorted afterwards
	 */

	template<class Field>
		void TriplesBB<Field>::addEntry(const Element &Elem, const size_t i, const size_t j) {
		_RowSortFlag = _ColSortFlag = false;
		_values.push_back(Elem);
		_RowV.push_back(i);
		_ColV.push_back(j);
	}

	/*
	  template<class Field, class Vector>
	  void TriplesBB<Field, Vector>::SortByRow()
	  {
	  RowWiseLessThan<Field,Vector> rwlt;
	  if(_RowSortFlag) return; // If already sorted, bail

	  std::sort( rawIndexedBegin(), rawIndexedEnd(), rwlt  );
	  _RowSortFlag = true;     // Sets the row sort flag
	  _ColSortFlag = false;    // Unset the col sort flag

	  }

	  template<class Field, class Vector>
	  void TriplesBB<Field, Vector>::SortByCol()
	  {

	  ColWiseLessThan<Field,Vector> cwlt;
	  if(_ColSortFlag) return;  // If already sorted, bail

	  std::sort( rawIndexedBegin(), rawIndexedEnd(), cwlt );
	  _ColSortFlag = true;     // Sets the Col sort flag
	  _RowSortFlag = false;    // Unset the Row sort flag
	  }
	*/

} // namespace LinBox

#endif // __LINBOX_triplesbb_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
