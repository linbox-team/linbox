/* linbox/blackbox/mapleBB.h
 * Copyright (C) 2002 Rich Seagraves
 *
 *
 * Written by Rich Seagraves <seagrave@cis.udel.edu>
 *
 *  This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */




#ifndef __MAPLEBB_H
#define __MAPLEBB_H

#include "linbox/blackbox/archetype.h"
#include "linbox/vector/vector-traits.h"
#include "linbox/util/debug.h"
#include "linbox/util/field-axpy.h"
#include <vector>
#include <algorithm>
#include <iostream>

namespace LinBox {


/** BlackBox wrapper for NAG Sparse Matrix format.
 *
 * This class acts as a wrapper for a pre-existing NAGSparse Matrix.
 * To be used for interface between LinBox and computer algebra systems such
 * as Maple that can encode sparse matrices in the NAGSparse format
 */ 

 template<class Field, class Vector>
 class MapleBB: public BlackboxArchetype<Vector> {

   typedef typename Field::Element Element;
 public:

   // Default constructor, do nothing.
   MapleBB();
   // A very bad constructor.  Takes 3 vectors and copies them
   MapleBB(Field F, std::vector<Element> values, std::vector<size_t> rowP, std::vector<size_t> colP, size_t rows, size_t cols, bool RowSortFlag = false, bool ColSortFlag = false);
   // Alternate constructor.  Allows for use of addEntry operation
   MapleBB(Field F, size_t rows, size_t cols, size_t reserve = 0);
   ~MapleBB() {};
   // Copy Constructor
   MapleBB(const MapleBB<Field,Vector> &);

   // Assignment operator for use in STL map
   const MapleBB<Field,Vector> & operator=(const MapleBB<Field,Vector> & );

   /** BlackBoxArchetype clone function.
    * Creates a copy of the NAGSparse Matrix and passes a pointer to it.
    * In this case it isn't too helpful
    * as this clone will of course suffer from the "siamese twin" problem.
    * The clonse created willwill point to the same data as the parent.
    * @return pointer to a new NAG format blackbox
    */

   BlackboxArchetype<Vector>* clone() const;

   /** BlackBoxArchetype apply function.
    * Take constant vector x and
    * vector y, and perform the calculation y = Ax.  Uses one of the three
    * private utility functions. It calls the generalized utility function
    * _apply if there is no special ordering, _fyapply if there is C_ordering
    * or _fxapply if there is fortran_ordering
    */

   Vector & apply(Vector &, const Vector &) const; // y = Ax;

   /* BlackBoxArchetype applyTranspose function. Take constant vector x and
    * vector y, and perform the calculation y = ATx.  Uses one of the three
    * private utility functions, in the manner described above.  Worthy of
    * note is the fact that applyTranspose works by passing the column
    * positions to the _apply functions as if they were rows, and row positions
    * as if they were columns, as if the matrix had been transposed.
    */

   Vector & applyTranspose(Vector &, const Vector &) const; // y = ATx

   /* BlackBoxArchetype rowdim function.  Passes back the number of rows of
    * the matrix represented.  Note that's the number of rows of the matrix
    * as if it were in dense format, not in it's actual representation here.
    */

   size_t rowdim() const;

   /* BlackBoxArchetype coldim function.  Passes back the number of columns
    * of the matrix represented.  Not much more to say about this.
    */

   size_t coldim() const;

   /* length function.  Returns number of non-zero entries */
   size_t size() const;

   // Add entry function, for intialization of matrix
   void addEntry(const Element &, const size_t, const size_t);

   // Sort functions.  When invoked, calls the STL sort utility
   // using indexedIterator.  If the entries are already sorted,
   // it does nothing.  SortbyRow sorts with row indices taking
   // precidence, while SortbyCol sorts with col indices taking precidence
   // Use _RowSortFlag and _ColSortFlag to see if the sort call is needed
   void SortByRow();
   void SortByCol();




   /* entryRep class.  Used as a value_type for the 
      RawIndexedIterator class.  Supports two important
      features, operator< and operator=.  These two enable 
      blackboxes to use STL functions such as sort()
   */


   // Forward declaration
   class entryRep;


   class RawIndexedIterator {
     public:
     // Types

     typedef MapleBB<Field, Vector>::entryRep value_type;
     typedef MapleBB<Field, Vector>::entryRep & reference;
     typedef MapleBB<Field, Vector>::entryRep * pointer;
     typedef typename std::iterator_traits<std::vector<Element>::iterator>::iterator_category iterator_category;
     typedef typename std::iterator_traits<std::vector<Element>::iterator>::difference_type difference_type;
     
     // Constructors
     RawIndexedIterator() : _element(0), _row(0), _col(0) {}
     RawIndexedIterator(std::vector<Element>::iterator element, std::vector<size_t>::iterator row, std::vector<size_t>::iterator col): _element(element), _row(row), _col(col) {}
     
     RawIndexedIterator( const RawIndexedIterator &In) :
       _element(In._element), _row(In._row), _col(In._col) {}

     const RawIndexedIterator &operator=(const RawIndexedIterator& rhs) {
       _element = rhs._element;
       _row = rhs._row;
       _col = rhs._col;
       return *this;
     }

     bool operator==(const RawIndexedIterator &rhs) const {
       return _element == rhs._element && _row == rhs._row && _col == rhs._col;
     }

     bool operator!=(const RawIndexedIterator &rhs) const {
       return _element != rhs._element || _row != rhs._row || _col != rhs._col;
     }

     bool operator<(const RawIndexedIterator &rhs) const {
       return _element < rhs._element; 
     }



     MapleBB<Field, Vector>::entryRep operator*() {
       return MapleBB<Field,Vector>::entryRep(_element, _row, _col, true);
     }

     RawIndexedIterator &operator++() {
       ++_element; ++_row; ++_col;
       return *this;
     }

     RawIndexedIterator operator++(int) {
       RawIndexedIterator tmp = *this;
       ++_element; ++_row; ++_col;
       return tmp;
     }

     RawIndexedIterator &operator--() {
       --_element; --_row; --_col;
       return *this;
     }

     RawIndexedIterator operator--(int) {
       RawIndexedIterator tmp = *this;
       --_element; --_row; --_col;
       return tmp;
     }

     RawIndexedIterator operator+(difference_type rhs) {
       return RawIndexedIterator(_element + rhs, _row + rhs, _col + rhs);
     }

     difference_type operator+(const RawIndexedIterator &rhs) {
       return _element + rhs._element;
     }

     RawIndexedIterator operator-(difference_type rhs) {
       return RawIndexedIterator(_element - rhs, _row - rhs, _col - rhs);
     }

     difference_type operator-(const RawIndexedIterator &rhs) {
       return _element - rhs._element;
     }

     RawIndexedIterator &operator+=(difference_type rhs) {
       _element += rhs; _row += rhs; _col += rhs;
       return *this;
     }

     RawIndexedIterator &operator-=(difference_type rhs) {
       _element -= rhs; _row -= rhs; _col -= rhs;
       return *this;
     }

     reference operator[](difference_type index) {
       return *( *this + index);
     }


      private:
         std::vector<Element>::iterator _element;
	 std::vector<size_t>::iterator _row;
	 std::vector<size_t>::iterator _col;
     
     };


   // begin() and end() functions for the RawIndexed iterator.  Used for the 
   // SortByRow and SortByCol

   RawIndexedIterator rawIndexedBegin() {
     return RawIndexedIterator( _values.begin(), _RowV.begin(), _ColV.begin() );
   }

   const RawIndexedIterator rawIndexedBegin() const {
     return RawIndexedIterator( _values.begin(), _RowV.begin(), _ColV.begin() );
   }

   RawIndexedIterator rawIndexedEnd() {
     return RawIndexedIterator( _values.end(), _RowV.end(), _ColV.end() );
   }

   const RawIndexedIterator rawIndexedEnd() const {
     return RawIndexedIterator( _values.end(), _RowV.end(), _ColV.end() );
   }

   /* Field accessor.  Will be used in by several functions to get the
    * field used by the class, to be passed into the linbox solution functions
    * such as rank, det, minpoly, or ssolve
    * 
    * NOTE NOTE - This function was changed, originally it took a field as an arguement and
    * initalized it.  This could break things!!!!!!
    */

   const Field & getField() const;

   /* Data accessors.  Used to access the 3 vectors containing Matrix data
    */

   const std::vector<Element> & getData() const;
   const std::vector<size_t> & getRows() const;
   const std::vector<size_t> & getCols() const;


 private:
   Field _F; // The field used by this class

   /* _values is a vector containing the elements of the BlackBox
    */

   std::vector<Element> _values;

   /* _RowV & _ColV are vectors containing the row & column indeces of the
    * Matrix
    */

   std::vector<size_t> _RowV, _ColV;

   /* The number of rows, columns and nonzero entries
    */

   size_t _rows, _cols;

   /* _apply is the generic apply utility funtion called by apply() and
    * applyTranspose().  Walks through the non-zero elements of the Matrix and
    * performs the proper calculation using the a vector of FieldAxpy's
    */

   void _apply(Vector &, const Vector &, std::vector<size_t>::const_iterator, std::vector<size_t>::const_iterator) const;

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

 template<class Field, class Vector>
 class RowWiseLessThan {
   public:
   bool operator()(const LinBox::MapleBB<Field,Vector>::entryRep &lhs, const LinBox::MapleBB<Field,Vector>::entryRep &rhs) {
     size_t lrow, rrow, lcol, rcol;
     lrow = lhs.getRow();
     rrow = rhs.getRow();
     lcol = lhs.getCol();
     rcol = rhs.getCol();

     if( lrow < rrow) return true;
     else if( lrow == rrow) {
        return (lcol < rcol); 
     }
     else return false;
   }
 };

 template<class Field, class Vector>
 class ColWiseLessThan {
   public:
   bool operator()(const LinBox::MapleBB<Field, Vector>::entryRep &lhs, const LinBox::MapleBB<Field,Vector>::entryRep &rhs) {
     size_t lrow, rrow, lcol, rcol;
     lrow = lhs.getRow();
     lcol = lhs.getCol();
     rrow = rhs.getRow();
     rcol = rhs.getCol();

     if( lcol < rcol) return true;
     else if( lcol == rcol) {
       return (lrow < rrow);
     }
     else return false;
   }
 };

 template<class Field, class Vector>
 class MapleBB<Field,Vector>::entryRep {
   friend class RowWiseLessThan<Field,Vector>;
   friend class ColWiseLessThan<Field,Vector>;
     public:
     // Default constructor
        entryRep() : _oelement(), _orow(0), _ocol(0), _flag(false) { std::cout << "Default constructor" << std::endl;}
	// Main constructor
	entryRep( std::vector<Element>::iterator element, std::vector<size_t>::iterator row, std::vector<size_t>::iterator col, bool flag = false ) : _flag(flag) {
		std::cout << "Main Constructor" << std::endl;
		if(flag) {
			_element = element;
			_row = row;
			_col = col;
		}
		else {
			_oelement = *element;
			_orow = *row;
			_ocol = *col;
		}
	}

	// Copy constructor
	entryRep(const entryRep &In) : _flag(false) {
		std::cout << "Copy Constructor" << std::endl;
		if( In._flag) {
			_oelement = *(In._element);
			_orow = *(In._row);
			_ocol = *(In._col);

		}
		else {
			_oelement = In._oelement;
			_orow = In._orow;
			_ocol = In._ocol;
		}
	} 
	
	~entryRep() {}


	
	// Assignment operator
	const entryRep &operator=(const entryRep &rhs) {
		Element setE;
		size_t setR, setC;
		if( rhs._flag) {
			setE = *(rhs._element);
			setR = *(rhs._row);
			setC = *(rhs._col);
		}
		else {
			setE = rhs._oelement;
			setR = rhs._orow;
			setC = rhs._ocol;
		}
			
		if( _flag) {
			*_element = setE;
			*_row = setR;
			*_col = setC;
		}
		else {
			_oelement = setE;
			_orow = setR;
			_ocol = setC;

		}
	  return *this;
	}

	const Element &getElement() const {
	  if( _flag) 
		return *_element;
	  else 	
		return _oelement;
	}

	const size_t &getRow() const {
	  if( _flag)
	  	return *_row;
	  else
		return _orow;

	}

	const size_t &getCol() const {
	  if( _flag) 
		  return *_col;
	  else
		return _ocol;

	}
	bool printFlag() {
		return _flag;
	}


     private:
	std::vector<Element>::iterator _element;
	std::vector<size_t>::iterator _row;
	std::vector<size_t>::iterator _col;
	Element _oelement;
	size_t	_orow;
	size_t 	_ocol;
	bool _flag;	  
   };


 /*  Constructor for the MapleLB class.  This is the constructor that is
  * expected to be used.  To use it, you must pass in a field element that
  * will work over the data (F), pointers to the 3 arrays used by the NAGSparse
  * format (values, rowP, colP), the number of rows and columns (rows and
  * cols), the number of non-zero elements (NNz) and the ordering, which
  * defaults to 0 (no ordering implied).
  */

 template<class Field, class Vector>
 MapleBB<Field, Vector>::MapleBB() {}

 /* bleck, "copying" constructor that takes 3 vectors as input.  Not
  * recommended
  */


 template<class Field, class Vector>
 MapleBB<Field, Vector>::MapleBB(Field F, std::vector<Element> values, std::vector<size_t> RowV, std::vector<size_t> ColV, size_t rows, size_t cols, bool RowSortFlag, bool ColSortFlag) :
   _F(F), _values(values), _RowV(RowV), _ColV(ColV), _rows(rows), _cols(cols), _RowSortFlag(RowSortFlag), _ColSortFlag(ColSortFlag)
 {
   int i;
   if( _rows > _cols)

     for( i = 0; i < _rows; i++)
       _faxpy.push_back( FieldAXPY<Field>(_F));

   else
     for( i = 0; i < _cols; i++)
       _faxpy.push_back( FieldAXPY<Field>(_F));

 }

 /* Better constructor that only takes the field, m, n and recommended
  * reserve (optional arguement) for use with STL vector reserve option
  * (reduce the number of memory management moves).  Meant to be used in
  * conjuction with the addEntry() method
  */

 template<class Field, class Vector>
 MapleBB<Field, Vector>::MapleBB( Field F, size_t rows, size_t cols, size_t res):
 	_F(F), _rows(rows), _cols(cols)
 {
	if(res != 0) {
		_values.reserve(res);
		_RowV.reserve(res);
		_ColV.reserve(res);
   }

	int i;
	if( _rows > _cols)
		for( i = 0; i < _rows; i++)
			_faxpy.push_back( FieldAXPY<Field>(_F));

	else
		for(i = 0; i < _cols; i++)
			_faxpy.push_back( FieldAXPY<Field>(_F));

	_RowSortFlag = false;
	_ColSortFlag = false;
 }



 template<class Field, class Vector>
 MapleBB<Field,Vector>::MapleBB(const MapleBB<Field,Vector> &In)
 {
   int i;

   _F = In._F;
   _values = In._values;
   _RowV = In._RowV;
   _ColV = In._COLV;
   _rows = In.rows; _cols = In._cols;
   _RowSortFlag = In._RowSortFlag;
   _ColSortFlag = In._ColSortFlag;

   if( _rows > _cols)

     for( i = 0; i < _rows; i++)
       _faxpy.push_back( FieldAXPY<Field>(_F));

   else
     for( i = 0; i < _cols; i++)
       _faxpy.push_back( FieldAXPY<Field>(_F));
 }

 template<class Field, class Vector>
 const MapleBB<Field,Vector> & MapleBB<Field,Vector>::operator=(const MapleBB<Field,Vector> & rhs)
 {
   int i, j;
   _F = rhs._F;
   _values = rhs._values;
   _RowV = rhs.RowV;
   _ColV = rhs.ColV;
   _rows = rhs._rows; _cols = rhs._cols;
   _RowSortFlag = rhs._RowSortFlag;
   _ColSortFlag  = rhs._ColSortFlag;

   i = rhs._faxpy.size();
   j = _faxpy.size();

   if( i > j )

     for( ; j < i; ++j )
       _faxpy.push_back(FieldAXPY<Field>(_F));

   else if( i < j)

     for(; j > i; --j)
       _faxpy.pop_back();

   return *this;
 }

/* BlackBoxArchetype clone function.  Creates a another NAGSparse Matrix
 * and returns a pointer to it.  Very simple in construction, just uses the
 * new operator.  Of course needs to be deleted to prevent a memory leak.
 * Note, the BlackBox created by this clone function is not an independant
 * entity.  A NAGSparse is little more than a wrapper over a pre-created
 * NAG Sparse Matrix that allows the linbox algorithms to work on that matrix.
 * Thus, this constructor creates another wrapper for the same data pointed to
 * by this object.
 */


 template<class Field, class Vector>
 BlackboxArchetype<Vector>* MapleBB<Field,Vector>::clone() const
 {
   BlackboxArchetype<Vector>* p = new MapleBB<Field,Vector>(_F,_values,_RowV,_ColV,_rows,_cols);
   return p;
 }

/* BlackBoxArchetype rowdim function.  Not much to say here. Returns the number
 * of rows of the Matrix were it in dense format.
 */

 template<class Field, class Vector>
 size_t MapleBB<Field,Vector>::rowdim() const
 {
   return _rows;
 }

/* BlackBoxArchetype coldim function.  Not much to say here either. Returns the
 * number of columns of the Matrix were it in dense format.
 */

 template<class Field, class Vector>
 size_t MapleBB<Field,Vector>::coldim() const
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

 template<class Field, class Vector>
 Vector & MapleBB<Field,Vector>::apply(Vector & y, const Vector & x) const
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

 template<class Field, class Vector>
 Vector & MapleBB<Field,Vector>::applyTranspose(Vector & y, const Vector & x) const
 {
   _apply( y, x, _ColV.begin(), _RowV.begin() );
   return y;
 }


 template<class Field, class Vector>
 void MapleBB<Field,Vector>::_apply(Vector & y, const Vector & x, std::vector<size_t>::const_iterator i, std::vector<size_t>::const_iterator j) const
 {
   typename Vector::iterator yp;
   typename Vector::const_iterator xp;
   typename Field::Element zero;
   std::vector<Element>::const_iterator v;
   std::vector<FieldAXPY<Field> >::iterator fa_i;

   _F.init(zero,0);

   for(fa_i = _faxpy.begin(); fa_i != _faxpy.end(); ++fa_i)
     fa_i->assign(zero);


   for( v = _values.begin(), fa_i = _faxpy.begin() - 1, xp = x.begin() - 1; v != _values.end(); ++i, ++j, ++v)
        (fa_i + *i)->accumulate(*v,*(xp + *j));


   for(fa_i = _faxpy.begin(), yp = y.begin(); yp != y.end(); ++yp, ++fa_i)
     fa_i->get(*yp);

 }

 /* size function.  Returns the number of non-zero entries.
  */

 template<class Field, class Vector>
 size_t MapleBB<Field, Vector>::size() const {
	 return _values.size();
 }

 /* addEntry method.  Allows user to add entries on the fly.  Meant to be used
  * with the "copyless" constructor above.  Note, will automatically set the
  * _sortFlag to false, as you can't be sure the entries are still sorted afterwards
  */

 template<class Field, class Vector>
 void MapleBB<Field, Vector>::addEntry(const Element &Elem, const size_t i, const size_t j) {
	 _RowSortFlag = _ColSortFlag = false;
	 _values.push_back(Elem);
	 _RowV.push_back(i);
	 _ColV.push_back(j);
 }


template<class Field, class Vector>
void MapleBB<Field, Vector>::SortByRow()
{
  RowWiseLessThan<Field,Vector> rwlt;
  if(_RowSortFlag) return; // If already sorted, bail

  std::sort( rawIndexedBegin(), rawIndexedEnd(), rwlt  );
  _RowSortFlag = true;     // Sets the row sort flag
  _ColSortFlag = false;    // Unset the col sort flag

}

template<class Field, class Vector>
void MapleBB<Field, Vector>::SortByCol()
{

 ColWiseLessThan<Field,Vector> cwlt;
 if(_ColSortFlag) return;  // If already sorted, bail

 std::sort( rawIndexedBegin(), rawIndexedEnd(), cwlt );
 _ColSortFlag = true;     // Sets the Col sort flag
 _RowSortFlag = false;    // Unset the Row sort flag
}


 template<class Field, class Vector>
 const Field & MapleBB<Field, Vector>::getField() const
 {
   return _F;
 }

 template<class Field, class Vector>
 const std::vector<typename Field::Element> & MapleBB<Field, Vector>::getData() const
 {
   return _values;
 }

 template<class Field, class Vector>
 const std::vector<size_t> & MapleBB<Field, Vector>::getRows() const
 {
   return _RowV;
 }

 template<class Field, class Vector>
 const std::vector<size_t> & MapleBB<Field, Vector>::getCols() const
 {
   return _ColV;
 }


#ifdef check
#include <iostream>
#include <fstream>
#include <test_common.h>

/* Linbox testing routine.  Tests the class to ensure that it functions
 * properly.  Uses the standard diagnostic model used by all Linbox BlackBox's
 * Runs 3 tests, a form test to ensure that results produced are proper, and a
 * function test to ensure that all BlackBoxArchetype API works properly.
 */

template<class Field,class Vector>
bool test(ofstream & report, Field &F)
{
    typedef Field::Element Element;
    Element blank;
    bool res1 = true, res2 = true;
    cout << "MapleBB blackbox suite" << endl;
    report << "MapleBB blackbox suite" << endl;

    cout << "Form Test. . .";
    report << "Form Test. . .";
    report << "A:: Values: 1 3 5  7 9 11 13 15" << endl;
    report << "   columns: 0 0 0  1 2  3  3  3" << endl;
    report << "and   rows: 0 1 4  1 2  0  3  4" << endl;
    report << "A:: mxn = 4x5.  NNZ = 8.  " << endl;
    report << "Now creating Matrix." << endl;

    std::vector<Element> val;
    for(int k = 0; k < 8; k++) {
      F.init(blank, 2 * k - 1);
      val.push_back(blank);
    }
    std::vector<size_t> cols;
    std::vector<size_t> rows;
    cols.push_back(0); cols.push_back(0); cols.push_back(0); cols.push_back(1);
    cols.push_back(2); cols.push_back(3); cols.push_back(3); cols.push_back(3);
    rows.push_back(0); rows.push_back(1); rows.push_back(4); rows.push_back(1);
    rows.push_back(2); rows.push_back(0); rows.push_back(3); rows.push_back(4);

    MapleBB<Field, Vector> testMaple(F, val, rows, cols, 4, 5, 8);

    int y_check[4][5];
    Vector x[4], y[5];
    report<< "Now checking each row:" << endl;
    report << "Create x as a m element vector w/ a 1 to look at each row" << endl;

    for(int i = 0; i < 4; ++i) {
      for(int j = 0; j < 5; ++j) {
	y_check[i][j] = 0;
      }
    }
    y_check[0][0] = 1;  y_check[0][1] = 3; y_check[0][4] = 4; y_check[1][1] = 7; y_check[2][2] = 9; y_check[3][0] = 11; y_check[3][3] = 13; y_check[3][4] = 15;
    for(int k = 0; k < testMaple.coldim(); ++k) {

      report << "Checking row " << k + 1 << endl;
      if(k > 0) x[k-1] = 0;
      x[k] = 1;
      testMaple.apply(y,x);
      // Now checks y against y_check
      for(int i = 0; i < 5; ++i) {
	if(y[i] != y_check[k][i]) {
	   res = false;
	   report << "ERROR:  expected " << y_check[k][i] << " at position (" << k << "," << i << "), got " y[i] << ".  Not cool." << endl;
	}
      }
    }
    if(res == false) {
      cout << "FAILED." << endl;
      report << "Form Test. . .FAILED." << endl;
    }
    else {
      cout << "PASSED." << endl;
      report << "Form Test. . .PASSED." << endl;
    }
    // Now on to the functional test
    cout << "FUNCTIONAL TEST. . .";
    report << "FUNCTIONAL TEST. . ." << endl;
    report << "Uses BBGeneralTest(A,F,report)" << endl;

    res2 = BBGeneralTest(testMaple, F, report);
    if(res2 == false) {
      cout << "FAILED " << endl;
      report << "FUNCTIONAL TEST. . .FAILED" << endl;
    }
    else {
      cout << "PASSED." << endl;
      report << "FUNCTIONAL TEST. . .PASSED." << endl;
    }
    return res1 && res2;
  }
#endif // ifdef check

} // namespace LinBox

#endif // ifdef __MAPLEBB_H
