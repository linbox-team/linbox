/* -*- mode: C++; style: linux -*- */

/* linbox/blackbox/nag-sparse.h
 * Copyright (C) 2002 Rich Seagraves
 *
 * Written by Rich Seagraves <seagrave@cis.udel.edu>
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __NAG_SPARSE_H
#define __NAG_SPARSE_H

#include "linbox-config.h"
#include "linbox/vector/vector-traits.h"
#include "linbox/util/debug.h"


// if __LINBOX_XML support in LinBox
#ifdef __LINBOX_XMLENABLED

#include "linbox/util/xml/linbox-reader.h"
#include "linbox/util/xml/linbox-writer.h"

using LinBox::Reader;
using LinBox::Writer;

#include <iostream>
#include <string>
#include <vector>

using std::istream;
using std::ostream;
using std::string;
using std::vector;

using std::cout;
using std::endl;

#endif




// For STL pair in RawIndexIterator
#include <utility>

namespace LinBox
{

        /** BlackBox wrapper for NAG Sparse Matrix format.
         *
         * This class acts as a wrapper for a pre-existing NAGSparse Matrix.
         * To be used for interface between LinBox and computer algebra systems such
         * as Maple that can encode sparse matrices in the NAGSparse format
         */

        template<class _Field>
        class NAGSparse{

	public:
		typedef _Field Field;
                typedef typename Field::Element Element;
                typedef size_t Index;
	
		
		// Default constructor, do nothing.
		NAGSparse();
		// The real constructor
		NAGSparse(Field F, Element* values, Index* rowP, Index* colP, Index rows, Index cols, Index NNz);
		// Destructor, once again do nothing
		~NAGSparse();
		
#ifdef __LINBOX_XMLENABLED
		NAGSparse(Reader &R) : _F(R.Down(1))
		{
			typedef typename Field::Element Element;
			vector<size_t> rows, cols;
			vector<Element> elem;
			Element e;
			size_t i;

			R.Up(1);
			if(!R.expectTagName("MatrixOver") ) return;
			if(!R.expectAttributeNum("rows", _rows) || !R.expectAttributeNum("cols", _cols)) return;
		
			if(!R.expectChildTag()) return;

			R.traverseChild();
			if(!R.expectTagName("field")) return;
			R.upToParent();

			if(!R.getNextChild() || !R.expectChildTag()) return;
			R.traverseChild();

			if(R.checkTagName("diag")) {
				if(!R.expectChildTag()) return;
				R.traverseChild();
				if(!R.expectTagName("entry") || !R.expectTagNumVector(elem)) return;
		       
				dynamic = true;
				_rowP = new size_t[elem.size()];
				_colP = new size_t[elem.size()];
				_values = new Element[elem.size()];
				_nnz = elem.size();
				
				for(i = 0; i < _nnz; ++i) {
					_rowP[i] = _colP[i] = i;
					_values[i] = elem[i];
				}
				
				R.upToParent();
				R.upToParent();
				R.getPrevChild();

				return;
			}
			else if(R.checkTagName("scalar")) {

				if(!R.expectChildTag()) return;
				R.traverseChild();
				if(!R.expectTagNum(e)) return;
				R.upToParent();

				dynamic = true;
				_rowP = new size_t[_rows];
				_colP = new size_t[_rows];
				_values = new Element[_rows];

				_nnz = _rows;
			
				for(i = 0; i < _nnz; ++i) {
					_rowP[i] = _colP[i] = i;
					_values[i] = e;
				}
				
				R.upToParent();
				R.getPrevChild();
			
				return;
			}
			else if(R.checkTagName("zero-one")) {
				if(!R.expectChildTag()) return;
				R.traverseChild();

				if(!R.expectTagName("index") || !R.expectTagNumVector(rows)) return;

				R.upToParent();
				if(!R.getNextChild() || !R.expectChildTag()) return;
				R.traverseChild();
				if(!R.expectTagName("index") || !R.expectTagNumVector(cols)) return;

				dynamic = true;
				_rowP = new size_t[rows.size()];
				_colP = new size_t[rows.size()];
				_values = new Element[rows.size()];

				_nnz = rows.size();
				
				_F.init(e, 1);
				for(i = 0; i < _nnz; ++i) {
					_rowP[i] = rows[i];
					_colP[i] = cols[i];
					_values[i] = e;
				}
				
				R.upToParent();
				R.upToParent();
				R.getPrevChild();
				
				return;
			}
			else if(!R.expectTagName("sparseMatrix") )
				return;
			else {

				if(!R.expectChildTag()) return;
				R.traverseChild();
				if(!R.expectTagName("index") || !R.expectTagNumVector(rows)) return;
			
				R.upToParent();
				if(!R.getNextChild() || !R.expectChildTag()) return;
				R.traverseChild();

				if(!R.expectTagName("index") || !R.expectTagNumVector(cols)) return;
			
			
				R.upToParent();
				if(!R.getNextChild() || !R.expectChildTag()) return;
				R.traverseChild();
			
				if(!R.expectTagName("entry") || !R.expectTagNumVector(elem)) return;

				R.upToParent();
				
				dynamic = true;
				_rowP = new size_t[rows.size()];
				_colP = new size_t[rows.size()];
				_values = new Element[rows.size()];
				_nnz = rows.size();

				for(i = 0; i < _nnz; ++i) {
					_rowP[i] = rows[i];
					_colP[i] = cols[i];
					_values[i] = elem[i];
				}
		
				R.upToParent();
				R.upToParent();
				R.getPrevChild();
				
				return;
			}
		}


		NAGSparse(const NAGSparse<Field, Vector> &M) : _F(M._F)
		{
			size_t i;
			
			dynamic = true;
			_nnz = M._nnz;
			_rowP = new size_t[M._nnz];
			_colP = new size_t[M._nnz];
			_values = new Element[M._nnz];
			_rows = M._rows;
			_cols = M._cols;

			for(i = 0; i < _nnz; ++i) {
				_rowP[i] = M._rowP[i];
				_colP[i] = M._colP[i];
				_values = M._values[i];
			}

			return;
		}


#endif

		/** BlackBoxArchetype apply function.
		 * Take constant vector x and
		 * vector y, and perform the calculation y = Ax.  Uses one of the three
		 * private utility functions. It calls the generalized utility function
		 * _apply if there is no special ordering, _fyapply if there is C_ordering
		 * or _fxapply if there is fortran_ordering
		 */
		template <class OutVector, class InVector>
		OutVector & apply(OutVector &, const InVector &) const; // y = Ax;

		/* BlackBoxArchetype applyTranspose function. Take constant vector x and
		 * vector y, and perform the calculation y = ATx.  Uses one of the three
		 * private utility functions, in the manner described above.  Worthy of
		 * note is the fact that applyTranspose works by passing the column
		 * positions to the _apply functions as if they were rows, and row positions
		 * as if they were columns, as if the matrix had been transposed.
		 */
		template <class OutVector, class InVector>
		OutVector & applyTranspose(OutVector &, const InVector &) const; // y = ATx

		/* BlackBoxArchetype rowdim function.  Passes back the number of rows of
		 * the matrix represented.  Note that's the number of rows of the matrix
		 * as if it were in dense format, not in it's actual representation here.
		 */

		size_t rowdim() const;

		/* BlackBoxArchetype coldim function.  Passes back the number of columns
		 * of the matrix represented.  Not much more to say about this.
		 */

		size_t coldim() const;

		/* Passes the ordering of the data.  There are three options:  C,
		 * FORTRAN, and ARB, whose values are defined by
		 * the public static ints below. C implies that the data is sorted
		 * by row, FORTRAN implies that the data is sorted by column, and
		 * ARB implies that there is no sorting.
		 */

		/* Non blackbox function.  Tells the number of nonzero entries
		 */
		size_t nnz() const;

		const Field& field() const { return  _F; }

		/* RawIterator class.  Iterates straight through the values of the matrix
		 */

#ifdef __LINBOX_XMLENABLED

		// because the NagSparse blackbox has no internal state, you cannot initalize it w/
		// XML.  So the read functions just return false

		ostream &write(ostream &) const;
		bool toTag(Writer &) const;

#endif




		class RawIterator {
		public:
			typedef Element value_type;

			RawIterator() {}

			RawIterator(value_type* ptr) :
				_ptr(ptr) {}

			RawIterator(const RawIterator& In) :
				_ptr(In._ptr) {}

			const RawIterator& operator=(const RawIterator& rhs) {
				_ptr = rhs._ptr;
				return *this;
			}


			bool operator==(const RawIterator &rhs) {
				return _ptr == rhs._ptr; 
			}

			bool operator!=(const RawIterator &rhs) {
				return _ptr != rhs._ptr;
			}

			RawIterator& operator++() {
				++_ptr;	
				return *this;
			}

			RawIterator operator++(int) {
						
				RawIterator tmp = *this;
				++_ptr;
				return tmp;
			}

			value_type &operator*() {
				return *_ptr;
			}

			const value_type &operator*() const {
				return *_ptr;
			}

		private:
			value_type *_ptr;
		};

		/* STL standard Begin and End functions.  Used to get
		 * the beginning and end of the data.  So that RawIterator
		 * can be used in algorithms like a normal STL iterator.
		 */

		RawIterator rawBegin() { return RawIterator( _values); }
		RawIterator rawEnd() { return RawIterator( _values + _nnz); }
		const RawIterator rawBegin() const { return RawIterator(_values); }
		const RawIterator rawEnd() const { return RawIterator(_values + _nnz); }

		/* RawIndexIterator - Iterates through the i and j of the current element
		 * and when accessed returns an STL pair containing the coordinates
		 */
		class RawIndexIterator {
		public:
			typedef std::pair<size_t, size_t> value_type;

			RawIndexIterator(size_t* row, size_t* col) :
				_row(row), _col(col) {}

			RawIndexIterator(const RawIndexIterator &In) :
				_row(In._row), _col(In._col) {}

			const RawIndexIterator &operator=(const RawIndexIterator &rhs) {
				_row = rhs._row;
				_col = rhs._col;
				return *this;
			}

			bool operator==(const RawIndexIterator &rhs) {
				return (_row == rhs._row && _col == rhs._col );
			}

			bool operator!=(const RawIndexIterator &rhs) {
				return !( *this == rhs);
			}

			const RawIndexIterator& operator++() {
				++_row;	++_col;
				return *this;
			}

			const RawIndexIterator operator++(int) {
				RawIndexIterator tmp = *this;
				++(*this);
				return tmp;
			}

			value_type operator*() {
				return std::pair<size_t,size_t>(*_row, *_col);
			}

			const value_type operator*() const {
				return std::pair<size_t,size_t>(*_row, *_col);
			}
		private:
			size_t* _row, *_col;
			value_type _curr;
		};

		RawIndexIterator indexBegin() {
			return RawIndexIterator(_rowP, _colP);
		}

		const RawIndexIterator indexBegin() const {
			return RawIndexIterator(_rowP, _colP);
		}

		RawIndexIterator indexEnd() {
			return RawIndexIterator(_rowP + _nnz,  _colP + _nnz );
		}

		const RawIndexIterator indexEnd() const {
			return RawIndexIterator(_rowP + _nnz, _colP + _nnz);
		}

		//                        int order() const;

		/* Passes the indexing of the data.  There are two options:  C and
		 * FORTRAN.  C imlies that the data is indexed using the 0 based scheme,
		 * array[0] is the first element of array.  FORTRAN implies that the data
		 * is indexed using the 1 based scheme, that array[1] is the first element
		 * of array.
		 */

		//                       int indexing() const;

		/* constant Static integer values which enumberate what the values of the
		 * order function mean.  These values are publically available so the
		 * application programmer can test for their value.
		 */

		//                        static const int C;
		//                        static const int FORTRAN;
		//                        static const int ARB;


	private:
		Field _F; // The field used by this class

		/* A pointer to an array of elements of the field.  In this case this also
		 * happens to be the first of 3 arrays in NAGSparse format.  This is the
		 * values of the data in the NAGSparse Matrix.
		 */

		Element *_values;

		/* _rowP is a pointer to an array of row indexes.  _colP is a pointer
		 * to an array of column indexes. These two are the other arrays of a
		 * NAGSparse format Matrix.  _rows and _cols are the number of rows and
		 * columns of the Matrix if it were in dense format.  _nnz is the Number of
		 * Non-Zero elements in the Matrix.  It also happens to be the length of
		 * the three NAGSparse arrays.
		 */

		Index *_rowP, *_colP, _rows, _cols, _nnz;

		/* _apply is the generic apply utility funtion called by apply() and
		 * applyTranspose().  Walks through the non-zero elements of the Matrix and
		 * performs the proper calculation using the axpyin method defined by the
		 * field element above.
		 */

		bool dynamic;

		template<class OutVector, class InVector>
		void _apply(OutVector &, const InVector &, Index*, Index*) const;

        };

	/* Default constructor.  Not really useful, just kinda there.
	 */
	template<class Field>
        NAGSparse<Field>::NAGSparse() : dynamic(false) {}


        /*  Constructor for the NAGSparse class.  This is the constructor that is
         * expected to be used.  To use it, you must pass in a field element that
         * will work over the data (F), pointers to the 3 arrays used by the NAGSparse
         * format (values, rowP, colP), the number of rows and columns (rows and
         * cols), the number of non-zero elements (NNz) and the ordering, which
         * defaults to 0 (no ordering implied).
         */

        template<class Field>
        NAGSparse<Field>::NAGSparse(Field F, Element* values, Index* rowP, Index* colP, Index rows, Index cols, Index NNz):
		_F(F), _values(values), _rowP(rowP), _colP(colP), _rows(rows), _cols(cols),_nnz(NNz), dynamic(false) 
        {}

	/* NAGSparse destructor.  If dynamic is true, delete memory allocated
	 * dynamically
	 */
	template<class Field>
	NAGSparse<Field>::~NAGSparse() 
	{
		if(dynamic) {
			delete [] _rowP;
			delete [] _colP;
			delete [] _values;
		}
	}



      
        /* BlackBoxArchetype rowdim function.  Not much to say here. Returns the number
         * of rows of the Matrix were it in dense format.
         */

        template<class Field>
        size_t NAGSparse<Field>::rowdim() const
        {
                return _rows;
        }

        /* BlackBoxArchetype coldim function.  Not much to say here either. Returns the
         * number of columns of the Matrix were it in dense format.
         */

        template<class Field>
        size_t NAGSparse<Field>::coldim() const
        {
                return _cols;
        }

        /* BlackBoxArchetype apply function.  Performs the y = Ax calculation, where
         * y and x are vectors passed in apply(y,x), and A is the present matrix.
         * Switches on the ordering of the data to
         * see which apply function to call.  In this case, y will be indexed using the
         * row index array, while x will be indexed using the column index array. The
         * switch determines which index has been sorted and thus which variable, x
         * or y, could be accessed fastest using the optimization.  If no ordering is
         * known, just call the generic apply algorithm.
         */

        template<class Field>
	template<class OutVector, class InVector>
        OutVector & NAGSparse<Field>::apply(OutVector & y, const InVector & x) const
        {
                _apply(y,x,_rowP,_colP);
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
        OutVector& NAGSparse<Field>::applyTranspose(OutVector & y, const InVector & x) const
        {
                _apply(y,x,_colP,_rowP);
                return y;
        }

        // Simple generic apply algorithm.  Works for completely unsorted NAG-sparse
        // arrays.  Doesn't take advantage of speed gains from use of pointers,
        // sorting of rows( C_order) or columns (Fortran_order).  Written to produce
        // something that would work.  Used by both apply and apply transpose

        template<class Field>
	template<class OutVector, class InVector>
        void NAGSparse<Field>::_apply(OutVector & y, const InVector & x, Index* _i, Index* _j) const
        {
                typename OutVector::iterator yp;
                typename InVector::const_iterator xp;
                Index* ip, *jp;
                Element* vp;

                // 0 out y.  Note, this implementation assumes a dense vector.
                for(yp = y.begin(); yp != y.end(); ++yp)
			_F.init(*yp , 0);

		//		switch(_index) {
		//			case 1:	
		//        			      yp = y.begin() - 1;
		//                     		      xp = x.begin() - 1;
		//                       case 1:
		yp = y.begin();
		xp = x.begin();
		//                             break;
		//                }


                for(ip = _i, jp = _j, vp = _values; ip < _i + _nnz; ++ip, ++jp, ++vp) 
                        _F.axpyin( *(yp + *ip), *vp, *(xp + *jp) );



        }


        // This function is a version of the _apply utility function that takes
        // Advantage of C or Fortran style ordering.  There is a version that
        // Uses the ordering for the y variable, and a version that uses it for
        // the x variable (depending if you are calling from apply() or
        // applyTranspose().  This is the y version.  Takes two vector, y & x,
        // returns y as y <- Ax (or ATx)
	/*
	  template<class Field, class Vector>
	  void NAGSparse<Field,Vector>::_fyapply(Vector & y, const Vector & x, Index* _i, Index* _j) const
	  {
	  Index* ip = _i, *jp = _j, i_val;
	  typename Vector::iterator yp, ybegin;
	  typename Vector::const_iterator xp;
	  Element* vp = _values;

	  for(yp = y.begin(); yp != y.end(); ++yp)
	  _F.init(*yp, 0);
	*/

	/*- Switches on the indexing.  A workaround for the problem of attempting
	 * to work with row and column index arrays that use fortran style indexing
	 * Everything is 1 off.  This sets the starting point of the vector pointers
	 * 1 back, putting everything into place.
	 */
	/*
	  switch(_index) {

	  case 2:
	  ybegin = y.begin() - 1;
	  xp = x.begin() - 1;
	  break;
	  case 1:
	  ybegin = y.begin();
	  xp = x.begin();
	  break;
	  }

	  i_val = *ip;
	  yp = ybegin + i_val;
	  while(ip < _i + _nnz ) {
	  while( (*ip == i_val) && (ip < _i + _nnz) ) {
	  _F.axpyin(*yp, *vp, *(xp + *jp) );
	  ++vp; ++ip; ++jp;
	  }
	  i_val = *ip;
	  yp = ybegin + i_val;
	  }

	  return;
	  }
	*/
        // This function is a version of the _apply utility function that takes
        // Advantage of C or Fortran style ordering.  There is a version that
        // Uses the ordering for the x variable, and a version that uses it for
        // the y variable (depending if you are calling from apply() or
        // applyTranspose().  This is the x version.  Takes two vector, y & x,
        // returns y as y <- Ax (or ATx)
	/*       template<class Field, class Vector>
		 void NAGSparse<Field,Vector>::_fxapply(Vector & y, const Vector & x, Index* _i, Index* _j) const
		 {
		 Index* ip = _i, *jp = _j, j_val;
		 typename Vector::iterator yp;
		 typename Vector::const_iterator xp, xbegin;
		 Element* vp = _values;

		 // 0 out y
		 for(yp = y.begin(); yp != y.end(); yp++)
		 _F.init(*yp, 0);
	*/
	/*- Switches on the indexing.  A workaround for the problem of attempting
	 * to work with row and column index arrays that use fortran style indexing
	 * Everything is 1 off.  This sets the starting point of the vector pointers
	 * 1 back, evening everything else out.
	 */
	/*               switch(_index) {

	case 2:
	yp = y.begin() - 1;
	xbegin = x.begin() - 1;
	break;
	case 1:
	yp = y.begin();
	xbegin = x.begin();
	break;
	}

	j_val = *jp;
	xp = xbegin + j_val;
	while(jp < _j + _nnz) {
	while( (*jp == j_val) && (jp < _j + _nnz)) {
	_F.axpyin( *(yp + *ip), *vp, *xp);
	++vp; ++ip; ++jp;
	}
	j_val = *jp;
	xp = xbegin + j_val;
	}

	return;
        }
	*/
	/*       template<class Field, class Vector>
		 int NAGSparse<Field,Vector>::order() const
		 {
		 return _order;
		 }

		 template<class Field, class Vector>
		 int NAGSparse<Field,Vector>::indexing() const
		 {
		 return _index;
		 }
	*/

#ifdef __LINBOX_XMLENABLED

	// As the NAGSparse Blackbox is meant to be a BlackBox wrapper around
	// existing data in NagSparse foramt, you cannot initalize it from XML data.
	// as such, the read & from Tag methods return false
			
		
	template<class Field, class Vector>
	ostream &NAGSparse<Field, Vector>::write(ostream & out) const
	{
		Writer W;
		if( toTag(W) ) 
			W.write(out);

		return out;
	}

	template<class Field, class Vector>
	bool NAGSparse<Field, Vector>::toTag(Writer &W) const
	{
		vector<size_t> row, col;
		vector<Element> elem;
		size_t i;
		string s;

		W.setTagName("MatrixOver");
		W.setAttribute("rows", Writer::numToString(s, _rows));
		W.setAttribute("cols", Writer::numToString(s, _cols));
		W.setAttribute("implDetail", "nag-sparse");

		W.addTagChild();
		_F.toTag(W);
		W.upToParent();

		W.addTagChild();
		W.setTagName("sparseMatrix");
		

		for(i = 0; i < _nnz; ++i ) {
			row.push_back(_rowP[i]);
			col.push_back(_colP[i]);
			elem.push_back(_values[i]);
		}
		
		W.addTagChild();
		W.setTagName("index");
		W.addNumericalList(row);
		W.upToParent();

		W.addTagChild();
		W.setTagName("index");
		W.addNumericalList(col);
		W.upToParent();
		
		W.addTagChild();
		W.setTagName("entry");
		W.addNumericalList(elem);
		W.upToParent();

		W.upToParent();

		return true;
	}


#endif



        template<class Field>
	size_t NAGSparse<Field>::nnz() const
	{
	        return _nnz;
	}
		

} // namespace LinBox

#endif // __NAG_SPARSE_H
