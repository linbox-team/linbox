/* MapleBB.h
 * Copyright (C) 2002 Rich Seagraves
 *
 *---------------------------------------------------
 * See COPYING for license information 
 */


#ifndef __MAPLEBB_H
#define __MAPLEBB_H

#include "linbox/blackbox/archetype.h"
#include "linbox/vector/vector-traits.h"
#include "linbox/util/debug.h"
#include "linbox/util/field-axpy.h"
#include <vector>

namespace LinBox {

/** BlackBox wrapper for NAG Sparse Matrix format.
 * 
 * This class acts as a wrapper for a pre-existing MapleBB Matrix.
 * To be used for interface between LinBox and computer algebra systems such
 * as Maple that can encode sparse matrices in the MapleBB format 
 */


 template<class Field, class Vector>
 class MapleBB: public BlackboxArchetype<Vector> {

   typedef typename Field::Element Element;
 public:

   // Default constructor, do nothing.
   MapleBB();
   // The real constructor
   MapleBB(Field F, Element* values, size_t* rowP, size_t* colP, size_t rows, size_t cols, size_t NNz);
   // Destructor, once again do nothing
   ~MapleBB() {};


   /** BlackBoxArchetype clone function.  
    * Creates a copy of the MapleBB Matrix and passes a pointer to it.  
    * In this case it isn't too helpful
    * as this clone will of course suffer from the "siamese twin" problem.
    * The clonse created will point to the same data as the parent.
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

   /* Passes the ordering of the data.  There are three options:  C,
    * FORTRAN, and ARB, whose values are defined by
    * the public static ints below. C implies that the data is sorted
    * by row, FORTRAN implies that the data is sorted by column, and
    * ARB implies that there is no sorting.  
    */

 private:
   Field _F; // The field used by this class
   /* A pointer to an array of elements of the field.  In this case this also
    * happens to be the first of 3 arrays in MapleBB format.  This is the
    * values of the data in the MapleBB Matrix.
    */

   Element *_values; 

   /* _rowP is a pointer to an array of row indexes.  _colP is a pointer
    * to an array of column indexes. These two are the other arrays of a 
    * MapleBB format Matrix.  _rows and _cols are the number of rows and 
    * columns of the Matrix if it were in dense format.  _nnz is the Number of
    * Non-Zero elements in the Matrix.  It also happens to be the length of
    * the three MapleBB arrays.
    */

   size_t *_rowP, *_colP, _rows, _cols, _nnz;

   /* _apply is the generic apply utility funtion called by apply() and
    * applyTranspose().  Walks through the non-zero elements of the Matrix and
    * performs the proper calculation using the axpyin method defined by the
    * field element above.
    */

   void _apply(Vector &, const Vector &, size_t*, size_t*) const;

   mutable std::vector<FieldAXPY<Field> > _faxpy;

 };

 /*  Constructor for the MapleBB class.  This is the constructor that is
  * expected to be used.  To use it, you must pass in a field element that
  * will work over the data (F), pointers to the 3 arrays used by the MapleBB
  * format (values, rowP, colP), the number of rows and columns (rows and
  * cols), the number of non-zero elements (NNz) and the ordering, which
  * defaults to 0 (no ordering implied).
  */

 template<class Field, class Vector>
 MapleBB<Field, Vector>::MapleBB() {}

 template<class Field, class Vector>
 MapleBB<Field, Vector>::MapleBB(Field F, Element* values, size_t* rowP, size_t* colP, size_t rows, size_t cols, size_t NNz):
   _F(F), _values(values), _rowP(rowP), _colP(colP), _rows(rows), _cols(cols),_nnz(NNz) { }

/* BlackBoxArchetype clone function.  Creates a another MapleBB Matrix
 * and returns a pointer to it.  Very simple in construction, just uses the
 * new operator.  Of course needs to be deleted to prevent a memory leak.
 * Note, the BlackBox created by this clone function is not an independant 
 * entity.  A MapleBB is little more than a wrapper over a pre-created
 * NAG Sparse Matrix that allows the linbox algorithms to work on that matrix.
 * Thus, this constructor creates another wrapper for the same data pointed to
 * by this object.
 */


 template<class Field, class Vector>
 BlackboxArchetype<Vector>* MapleBB<Field,Vector>::clone() const 
 {
   MapleBB<Field,Vector>* p = new MapleBB<Field,Vector>(_F,_values,_rowP,_colP,_rows,_cols,_nnz);
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

 template<class Field, class Vector>
 Vector & MapleBB<Field,Vector>::applyTranspose(Vector & y, const Vector & x) const
 {
   _apply(y,x,_colP,_rowP);
   return y;
 }

// Simple generic apply algorithm.  Works for completely unsorted NAG-sparse
// arrays.  Doesn't take advantage of speed gains from use of pointers,
// sorting of rows( C_order) or columns (Fortran_order).  Written to produce
// something that would work.  Used by both apply and apply transpose


 template<class Field, class Vector>
 void MapleBB<Field,Vector>::_apply(Vector & y, const Vector & x, size_t* _i, size_t* _j) const
 {
   typename Vector::iterator yp;
   typename Vector::const_iterator xp;
   size_t* ip, *jp;
   Element* vp;
   std::vector<FieldAXPY<Field> >::iterator fa_i;
   
   if(_faxpy.size() == 0) {

     if(_cols > _rows) 
       for( int i = _cols; i--;)
	 _faxpy.push_back(FieldAXPY<Field>(_F));

     else
       for(int i = _rows; i--;)
	 _faxpy.push_back(FieldAXPY<Field>(_F));

   }
   else {
     typename Field::Element zero;
     _F.init(zero,0);
     
     for(fa_i = _faxpy.begin(); fa_i != _faxpy.end(); ++fa_i)
       fa_i->assign(zero);
   }

   fa_i = _faxpy.begin() - 1;
   xp = x.begin() - 1;
   
   for(ip = _i, jp = _j, vp = _values; ip < _i + _nnz; ++ip, ++jp, ++vp)   
     (fa_i + *ip)->accumulate(*vp,*(xp + *jp));


   for(fa_i = _faxpy.begin(), yp = y.begin(); yp != y.end(); ++yp, ++fa_i)
     fa_i->get(*yp);

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

template<class Vector,class Field>
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

    Element val[8];
    for(int k = 0; k < 8; k++) {
      val[k] = F.init(blank, 2 * k - 1);
    }
    Index cols[8] = {0,0,0,1,2,3,3,3};
    Index rows[8] = {0,1,4,1,2,0,3,4};
    MapleBB<Field, Vector> testNAG(F, val, cols, rows, 4, 5, 8, 2, 1);

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
    for(int k = 0; k < testNAG.coldim(); ++k) {
  
      report << "Checking row " << k + 1 << endl;
      if(k > 0) x[k-1] = 0;
      x[k] = 1;
      testNAG.apply(y,x);
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
    
    res2 = BBGeneralTest(testNAG, F, report);
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

#endif // ifdef __NAGSPARSE_H
