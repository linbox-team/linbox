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




#ifndef __TRIPLESBB_H
#define __TRIPLESBB_H

#include "linbox-config.h"
#include "linbox/blackbox/archetype.h"
//#include "linbox/vector/vector-traits.h"
#include "linbox/util/debug.h"
#include "linbox/util/field-axpy.h"

#ifdef XMLENABLED
// For LinBox XML support.  For more information, check
// linbox/util/xml/README

#include "linbox/util/xml/linbox-reader.h"
#include "linbox/util/xml/linbox-writer.h"

using LinBox::Reader;
using LinBox::Writer;

#endif

#include <vector>
#include <iostream>

using std::istream;
using std::ostream;


namespace LinBox {


/** BlackBox wrapper for NAG Sparse Matrix format.
 *
 * This class acts as a wrapper for a pre-existing NAGSparse Matrix.
 * To be used for interface between LinBox and computer algebra systems such
 * as Maple that can encode sparse matrices in the NAGSparse format
 */ 

 template<class Field, class Vector>
 class TriplesBB: public BlackboxArchetype<Vector> {


   typedef typename Field::Element Element;
 public:

   // Default constructor, do nothing.
   TriplesBB();
   // A very bad constructor.  Takes 3 vectors and copies them
   TriplesBB(Field F, std::vector<Element> values, std::vector<size_t> rowP, std::vector<size_t> colP, size_t rows, size_t cols, bool RowSortFlag = false, bool ColSortFlag = false);
   // Alternate constructor.  Allows for use of addEntry operation
   TriplesBB(Field F, size_t rows, size_t cols, size_t reserve = 0);
   ~TriplesBB() {};
   // Copy Constructor
   TriplesBB(const TriplesBB<Field,Vector> &);


#ifdef XMLENABLED
   TriplesBB(Reader &);
#endif

   // Assignment operator for use in STL map
   const TriplesBB<Field,Vector> & operator=(const TriplesBB<Field,Vector> & );

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

   // only if XML reading & writing are enabled 
#ifdef XMLENABLED
   // XML in & out functions
   ostream &write(ostream &) const;
   bool toTag(Writer &) const;

#endif


   // Add entry function, for intialization of matrix
   void addEntry(const Element &, const size_t, const size_t);

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


   // small util function that determines the larger of two input size_t's
   size_t _max(size_t,size_t) const; 

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

 template<class Field, class Vector>
 TriplesBB<Field, Vector>::TriplesBB() {}

 /* bleck, "copying" constructor that takes 3 vectors as input.  Not
  * recommended
  */


 template<class Field, class Vector>
 TriplesBB<Field, Vector>::TriplesBB(Field F, std::vector<Element> values, std::vector<size_t> RowV, std::vector<size_t> ColV, size_t rows, size_t cols, bool RowSortFlag, bool ColSortFlag) :
   _F(F), _values(values), _RowV(RowV), _ColV(ColV), _rows(rows), _cols(cols), _faxpy(_max(rows,cols), FieldAXPY<Field>(F)), _RowSortFlag(RowSortFlag), _ColSortFlag(ColSortFlag)
 {}

 /* Better constructor that only takes the field, m, n and recommended
  * reserve (optional arguement) for use with STL vector reserve option
  * (reduce the number of memory management moves).  Meant to be used in
  * conjuction with the addEntry() method
  */

 template<class Field, class Vector>
 TriplesBB<Field, Vector>::TriplesBB( Field F, size_t rows, size_t cols, size_t res):
   _F(F), _rows(rows), _cols(cols), _faxpy( _max(rows, cols), FieldAXPY<Field>(F)), _RowSortFlag(false), _ColSortFlag(false)
 {
   if(res != 0) {
     _values.reserve(res);
     _RowV.reserve(res);
     _ColV.reserve(res);
   }

 }



 template<class Field, class Vector>
 TriplesBB<Field,Vector>::TriplesBB(const TriplesBB<Field,Vector> &In) :
    _faxpy( _max(In._rows, In._cols), FieldAXPY<Field>(In._F))

 {
   _F = In._F;
   _values = In._values;
   _RowV = In._RowV;
   _ColV = In._COLV;
   _rows = In.rows; _cols = In._cols;
   _RowSortFlag = In._RowSortFlag;
   _ColSortFlag = In._ColSortFlag;

 }
#ifdef XMLENABLED

 template<class Field, class Vector>
 TriplesBB<Field, Vector>::TriplesBB(Reader &R) : _F(R.Down(1))
 {
	 size_t i;
	 Element e;
	 _RowSortFlag = _ColSortFlag = false;
	 R.Up(1);


	 if(!R.expectTagName("MatrixOver")) return;

	 if(!R.expectAttributeNum("rows", _rows) || !R.expectAttributeNum("cols", _cols)) return;

	 if(!R.expectChildTag()) return;

	 R.traverseChild();
	 if(!R.expectTagName("field")) return;
	 R.upToParent();

	 if(!R.getNextChild()) {
		 R.setErrorString("Couldn't find matrix description in <matrix> tag");
		 R.setErrorCode(Reader::OTHER);
		 return;
	 }
	 R.traverseChild();
	 

	 // three divergent cases.  By default, a TriplesBB
	 // can be a sparseMatrix, a diag, or a scalar
	 // This is the code for each
	 //
	 if(R.checkTagName("diag")) {
		 if(!R.expectChildTag()) return;
		 R.traverseChild();
		 if(!R.expectTagName("entry") 
		    || !R.expectTagNumVector(_values)) 
			 return;
		 
		 _RowV.clear();
		 _ColV.clear();
		 if( _rows <= _cols) 
			 for(i = 0; i < _rows; ++i) {
				 _RowV.push_back(i);
				 _ColV.push_back(i);
			 }
		 else
			 for(i = 0; i < _cols; ++i) {
				 _RowV.push_back(i);
				 _ColV.push_back(i);
			 }
	 }
	 else if(R.checkTagName("scalar")) {
		 if(!R.expectChildTag()) return;
		 R.traverseChild();
		 if(!R.expectTagNum(e)) return;
		 R.upToParent();

		 _RowV.clear();
		 _ColV.clear();
		 if(_rows <= _cols)
			 for(i = 0; i < _rows; ++i) {
				 _RowV.push_back(i);
				 _ColV.push_back(i);
				 _values.push_back(e);
			 }		
		 else
			 for(i = 0; i < _cols; ++i) {
				 _RowV.push_back(i);
				 _ColV.push_back(i);
				 _values.push_back(e);
			 }
	 }
	 else if(R.checkTagName("zero-one")) {
		 if(!R.expectChildTag()) return;
		 R.traverseChild();
		 if(!R.expectTagName("index") || !R.expectTagNumVector(_RowV)) return;
		 R.upToParent();

		 if(!R.getNextChild()) {
			 R.setErrorString("Couldn't find column indices for zero-one matrix");
			 R.setErrorCode(Reader::OTHER);
			 return;
		 }

		 if(!R.expectChildTag()) return;

		 R.traverseChild();
		 if(!R.expectTagName("index") || !R.expectTagNumVector(_ColV)) return;
		 R.upToParent();

		 R.upToParent();

		 _values.clear();
		 _F.init(e, 1);
		 for(i = 0; i < _RowV.size(); ++i) {
			 _values.push_back(e);
		 }
		 
		 if(_faxpy.size() != _max(_rows, _cols))
			 _faxpy.resize(_max(_rows, _cols), FieldAXPY<Field>(_F));

	 }


	 else if(!R.expectTagName("sparseMatrix"))
		 return;
	 else {
		 
		 if(!R.expectChildTag()) return;
		 R.traverseChild();
		 if(!R.expectTagName("index") 
		    || !R.expectTagNumVector(_RowV)) 
			 return;
		 R.upToParent();
		 
		 if(!R.getNextChild() ) {
			 R.setErrorString("Couldn't find columnar indices in sparse matrix");
			 R.setErrorCode(Reader::OTHER);
			 return;
		 }

		  if(!R.expectChildTag()) 
			 return;
		 R.traverseChild();
		 if(!R.expectTagName("index")
		    || !R.expectTagNumVector(_ColV))
			 return;
		 R.upToParent();

		 if(!R.getNextChild()) {
			 R.setErrorString("Couldn't find matrix entries in sparse matrix");
			 R.setErrorCode(Reader::OTHER);
			 return;
		 }
			 
		 if( !R.expectChildTag())
			 return;
		 R.traverseChild();
		 if(!R.expectTagName("entry") 
		    || !R.expectTagNumVector(_values))
			 return;

		 R.upToParent();
		 R.upToParent();

	 }
	 // now we have to reset the _faxby object
	 if(_faxpy.size() != _max(_rows, _cols)) 
		 _faxpy.resize(_max(_rows, _cols), FieldAXPY<Field>(_F));

	 return;
 }

#endif


 template<class Field, class Vector>
 const TriplesBB<Field,Vector> & TriplesBB<Field,Vector>::operator=(const TriplesBB<Field,Vector> & rhs)
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
 BlackboxArchetype<Vector>* TriplesBB<Field,Vector>::clone() const
 {
   BlackboxArchetype<Vector>* p = new TriplesBB<Field,Vector>(_F,_values,_RowV,_ColV,_rows,_cols);
   return p;
 }

/* BlackBoxArchetype rowdim function.  Not much to say here. Returns the number
 * of rows of the Matrix were it in dense format.
 */

 template<class Field, class Vector>
 size_t TriplesBB<Field,Vector>::rowdim() const
 {
   return _rows;
 }

/* BlackBoxArchetype coldim function.  Not much to say here either. Returns the
 * number of columns of the Matrix were it in dense format.
 */

 template<class Field, class Vector>
 size_t TriplesBB<Field,Vector>::coldim() const
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
 Vector & TriplesBB<Field,Vector>::apply(Vector & y, const Vector & x) const
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
 Vector & TriplesBB<Field,Vector>::applyTranspose(Vector & y, const Vector & x) const
 {
   _apply( y, x, _ColV.begin(), _RowV.begin() );
   return y;
 }


 template<class Field, class Vector>
 void TriplesBB<Field,Vector>::_apply(Vector & y, const Vector & x, std::vector<size_t>::const_iterator i, std::vector<size_t>::const_iterator j) const
 {
   typename Vector::iterator yp;
   typename Vector::const_iterator xp;
   typename Field::Element zero;
   typename std::vector<Element>::const_iterator v;
   typename std::vector<FieldAXPY<Field> >::iterator fa_i;

   _F.init(zero,0);

   for(fa_i = _faxpy.begin(); fa_i != _faxpy.end(); ++fa_i) 
     fa_i->assign(zero);

   for( v = _values.begin(), fa_i = _faxpy.begin() - 1, xp = x.begin() - 1; v != _values.end(); ++i, ++j, ++v) 
        (fa_i + *i)->accumulate(*v,  *(xp + *j));



   for(fa_i = _faxpy.begin(), yp = y.begin(); yp != y.end(); ++yp, ++fa_i) 
     fa_i->get(*yp);

  

 }

 /* size function.  Returns the number of non-zero entries.
  */

 template<class Field, class Vector>
 size_t TriplesBB<Field, Vector>::size() const {
	 return _values.size();
 }

#ifdef XMLENABLED

 // Takes in an ostream and tries to write to it
 template<class Field, class Vector>
 ostream &TriplesBB<Field, Vector>::write(ostream &o) const {
	 Writer W;
	 if( toTag(W) ) 
		 W.write(o);

	 return o;
 }

 template<class Field, class Vector>
 bool TriplesBB<Field, Vector>::toTag(LinBox::Writer &W) const {
	 string buffer;

	 W.setTagName("MatrixOver");
	 W.setAttribute("rows", LinBox::Writer::numToString(buffer, _rows));
	 W.setAttribute("cols", LinBox::Writer::numToString(buffer, _cols));
	 W.setAttribute("implDetail", "triples");

	 W.addTagChild();
	 _F.toTag(W);
	 W.upToParent();

	 W.addTagChild();
	 W.setTagName("sparseMatrix");
		
	 W.addTagChild();
	 W.setTagName("index");
	 W.addNumericalList(_RowV, true);
	 W.upToParent();

	 W.addTagChild();
	 W.setTagName("index");
	 W.addNumericalList(_ColV, true);
	 W.upToParent();
	 
	 W.addTagChild();
	 W.setTagName("entry");
	 W.addNumericalList(_values);
	 W.upToParent();

	 W.upToParent();

	 return true;
 }

#endif




 /* addEntry method.  Allows user to add entries on the fly.  Meant to be used
  * with the "copyless" constructor above.  Note, will automatically set the
  * _sortFlag to false, as you can't be sure the entries are still sorted afterwards
  */

 template<class Field, class Vector>
 void TriplesBB<Field, Vector>::addEntry(const Element &Elem, const size_t i, const size_t j) {
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

 template<class Field, class Vector>
 const Field & TriplesBB<Field, Vector>::getField() const
 {
   return _F;
 }

 template<class Field, class Vector>
 const std::vector<typename Field::Element> & TriplesBB<Field, Vector>::getData() const
 {
   return _values;
 }

 template<class Field, class Vector>
 const std::vector<size_t> & TriplesBB<Field, Vector>::getRows() const
 {
   return _RowV;
 }

 template<class Field, class Vector>
 const std::vector<size_t> & TriplesBB<Field, Vector>::getCols() const
 {
   return _ColV;
 }

template<class Field, class Vector>
size_t TriplesBB<Field, Vector>::_max(size_t s1, size_t s2) const
{
	return s1 > s2 ? s1 : s2;
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
    cout << "TriplesBB blackbox suite" << endl;
    report << "TriplesBB blackbox suite" << endl;

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

    TriplesBB<Field, Vector> testMaple(F, val, rows, cols, 4, 5, 8);

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
