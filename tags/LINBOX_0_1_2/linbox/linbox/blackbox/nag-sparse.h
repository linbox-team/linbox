/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

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

#include "linbox/blackbox/archetype.h"
#include "linbox/vector/vector-traits.h"
#include "linbox/util/debug.h"

namespace LinBox 
{

        /** BlackBox wrapper for NAG Sparse Matrix format.
         * 
         * This class acts as a wrapper for a pre-existing NAGSparse Matrix.
         * To be used for interface between LinBox and computer algebra systems such
         * as Maple that can encode sparse matrices in the NAGSparse format 
         */

        template<class Field, class Vector>
        class NAGSparse: public BlackboxArchetype<Vector> {

                typedef typename Field::Element Element;
                typedef size_t Index;
                public:

                        // Default constructor, do nothing.
                        NAGSparse();
                        // The real constructor
                        NAGSparse(Field F, Element* values, Index* rowP, Index* colP, Index rows, Index cols, Index NNz, int order = 0, int index = 1);
                         // Destructor, once again do nothing
                         ~NAGSparse() {};


                        /** BlackBoxArchetype clone function.  
                         * Creates a copy of the NAGSparse Matrix and passes a pointer to it.  
                         * In this case it isn't too helpful
                         * as this clone will of course suffer from the "siamese twin" problem.
                         * The clonse created will point to the same data as the parent.
                         * @return pointer to a new NAG format blackbox
                         */

                        BlackBoxArchetype<Vector>* clone() const;

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

                        int order() const;

                        /* Passes the indexing of the data.  There are two options:  C and
                         * FORTRAN.  C imlies that the data is indexed using the 0 based scheme,
                         * array[0] is the first element of array.  FORTRAN implies that the data
                         * is indexed using the 1 based scheme, that array[1] is the first element
                         * of array.  
                         */
                        
                        int indexing() const;

                        /* constant Static integer values which enumberate what the values of the
                         * order function mean.  These values are publically available so the
                         * application programmer can test for their value.
                         */
  
                        static const int C;
                        static const int FORTRAN;
                        static const int ARB;


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
                        int _order; // The ordering of the data: C, FORTRAN, or ARB (default).
                        int _index; // The indexing of the data:  C (default) or FORTRAN

                        /* _apply is the generic apply utility funtion called by apply() and
                         * applyTranspose().  Walks through the non-zero elements of the Matrix and
                         * performs the proper calculation using the axpyin method defined by the
                         * field element above.
                         */

                        void _apply(Vector &, const Vector &, Index*, Index*) const;

                        /* _fyapply is an apply utility function with a few optimizations that take
                         * advantage of special ordering of the data in NAGSparse format.  It
                         * is called in a situation where the ordering is done over the index array
                         * used by the y vector in the apply calculation.
                         */

                        void _fyapply(Vector &, const Vector &, Index*,Index*) const;

                        /* _fxapply is an apply utility function with a few optimizations that take
                         * advantage of the special ordering of the data in NAGSparse format.  Same
                         * as _fyapply, but is used when ordering is done ov erthe index array used
                         * by the x vector in the apply calculation.
                         */

                        void _fxapply(Vector &, const Vector &, Index*, Index*) const;

        };

	/* Default constructor.  Not really useful, just kinda there.
	 */
	template<class Field, class Vector>
        NAGSparse<Field,Vector>::NAGSparse() : C(1), FORTRAN(2), ARB(3) {}


        /*  Constructor for the NAGSparse class.  This is the constructor that is
         * expected to be used.  To use it, you must pass in a field element that
         * will work over the data (F), pointers to the 3 arrays used by the NAGSparse
         * format (values, rowP, colP), the number of rows and columns (rows and
         * cols), the number of non-zero elements (NNz) and the ordering, which
         * defaults to 0 (no ordering implied).
         */

        template<class Field, class Vector>
        NAGSparse<Field, Vector>::NAGSparse(Field F, Element* values, Index* rowP, Index* colP, Index rows, Index cols, Index NNz, int order = 0, int index = 1):
        C(1), FORTRAN(2), ARB(0), _F(F), _values(values), _rowP(rowP), _colP(colP), _rows(rows), _cols(cols),_nnz(NNz), _order(order), _index(index) 
        { 
                if (_order != C && _order != FORTRAN) _order = ARB; 
                if (_index != C && _index != FORTRAN) _index = C; 
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
        BlackBoxArchetype<Vector>* NAGSparse<Field,Vector>::clone() const 
        {
                NAGSparse<Field,Vector>* p = new NAGSparse<Field,Vector>(_F,_values,_rowP,_colP,_rows,_cols,_nnz,_order,_index);
                return p;
        }

        /* BlackBoxArchetype rowdim function.  Not much to say here. Returns the number
         * of rows of the Matrix were it in dense format.
         */

        template<class Field, class Vector>
        size_t NAGSparse<Field,Vector>::rowdim() const
        {
                return _rows;
        }

        /* BlackBoxArchetype coldim function.  Not much to say here either. Returns the
         * number of columns of the Matrix were it in dense format.
         */

        template<class Field, class Vector>
        size_t NAGSparse<Field,Vector>::coldim() const
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
        Vector & NAGSparse<Field,Vector>::apply(Vector & y, const Vector & x) const 
        {
                switch(_order) {

                        case 1: // C Order - 
                                _fyapply(y,x,_rowP,_colP);
                                break;

                        case 2: // Fortran Order  
                                _fxapply(y, x, _rowP, _colP);
                                break;

                        case 0: // No valid ordering at input, use regular apply
                                _apply(y,x,_rowP,_colP);
                                break;
                }
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
        Vector & NAGSparse<Field,Vector>::applyTranspose(Vector & y, const Vector & x) const
        {
                switch(_order) {
                        
                        case 1: // C order
                                _fxapply(y, x, _colP, _rowP);
                                break;

                        case 2: // Fortran order
                                _fyapply(y,x,_colP,_rowP);
                                break;

                        case 0:  // No valid ordering input, use regular apply
                                _apply(y,x,_colP,_rowP);
                                break;
                }
                return y;
        }

        // Simple generic apply algorithm.  Works for completely unsorted NAG-sparse
        // arrays.  Doesn't take advantage of speed gains from use of pointers,
        // sorting of rows( C_order) or columns (Fortran_order).  Written to produce
        // something that would work.  Used by both apply and apply transpose

        template<class Field, class Vector>
        void NAGSparse<Field,Vector>::_apply(Vector & y, const Vector & x, Index* _i, Index* _j) const
        {
                typename Vector::iterator yp;
                typename Vector::const_iterator xp;
                Index* ip, *jp;
                Element* vp;
  
                // 0 out y.  Note, this implementation assumes a dense vector.
                for(yp = y.begin(); yp != y.end(); ++yp) 
                         _F.init(*yp , 0);

                switch(_index) {
                        
                        case 2:
                                yp = y.begin() - 1;
                                xp = x.begin() - 1;
                                break;
                        case 1:
                                yp = y.begin();
                                xp = x.begin();
                                break;
                }

                for(ip = _i, jp = _j, vp = _values; ip < _i + _nnz; ++ip, ++jp, ++vp) 
                        _F.axpyin( *(yp + *ip), *vp, *(xp + *jp)); 
        }

        // This function is a version of the _apply utility function that takes
        // Advantage of C or Fortran style ordering.  There is a version that
        // Uses the ordering for the y variable, and a version that uses it for
        // the x variable (depending if you are calling from apply() or
        // applyTranspose().  This is the y version.  Takes two vector, y & x,
        // returns y as y <- Ax (or ATx)

        template<class Field, class Vector>
        void NAGSparse<Field,Vector>::_fyapply(Vector & y, const Vector & x, Index* _i, Index* _j) const 
        {
                Index* ip = _i, *jp = _j, i_val;
                typename Vector::iterator yp, ybegin;
                typename Vector::const_iterator xp;
                Element* vp = _values;

                for(yp = y.begin(); yp != y.end(); ++yp) 
                        _F.init(*yp, 0);


                /** Switches on the indexing.  A workaround for the problem of attempting
                 * to work with row and column index arrays that use fortran style indexing
                 * Everything is 1 off.  This sets the starting point of the vector pointers
                 * 1 back, putting everything into place.
                 */

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

        // This function is a version of the _apply utility function that takes
        // Advantage of C or Fortran style ordering.  There is a version that
        // Uses the ordering for the x variable, and a version that uses it for
        // the y variable (depending if you are calling from apply() or
        // applyTranspose().  This is the x version.  Takes two vector, y & x,
        // returns y as y <- Ax (or ATx)       
        template<class Field, class Vector>
        void NAGSparse<Field,Vector>::_fxapply(Vector & y, const Vector & x, Index* _i, Index* _j) const 
        {
                Index* ip = _i, *jp = _j, j_val;
                typename Vector::iterator yp;
                typename Vector::const_iterator xp, xbegin;
                Element* vp = _values;

                // 0 out y
                for(yp = y.begin(); yp != y.end(); yp++)
                        _F.init(*yp, 0);

                /** Switches on the indexing.  A workaround for the problem of attempting
                 * to work with row and column index arrays that use fortran style indexing
                 * Everything is 1 off.  This sets the starting point of the vector pointers
                 * 1 back, evening everything else out.
                 */
                switch(_index) {

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

        template<class Field, class Vector>
        int NAGSparse<Field,Vector>::order() const 
        {
                return _order;
        }

        template<class Field, class Vector>
        int NAGSparse<Field,Vector>::indexing() const
        {
                return _index;
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
                cout << "NAGSparse blackbox suite" << endl;
                report << "NAGSparse blackbox suite" << endl;
    
                cout << "Form Test. . .";
                report << "Form Test. . .";
                report << "A:: Values: 1 3 5  7 9 11 13 15" << endl;
                report << "   columns: 0 0 0  1 2  3  3  3" << endl;
                report << "and   rows: 0 1 4  1 2  0  3  4" << endl;
                report << "A:: mxn = 4x5.  NNZ = 8.  " << endl;
                report << "Now creating Matrix." << endl;
        
                Element val[8];
                for(int k = 0; k < 8; k++) 
                        val[k] = F.init(blank, 2 * k - 1);
                
                Index cols[8] = {0,0,0,1,2,3,3,3};
                Index rows[8] = {0,1,4,1,2,0,3,4};
                NAGSparse<Field, Vector> testNAG(F, val, cols, rows, 4, 5, 8, 2, 1);

                int y_check[4][5];
                Vector x[4], y[5];
                report<< "Now checking each row:" << endl;
                report << "Create x as a m element vector w/ a 1 to look at each row" << endl;
            
                for(int i = 0; i < 4; ++i) 
                        for(int j = 0; j < 5; ++j) 
                                y_check[i][j] = 0;


                y_check[0][0] = 1;  y_check[0][1] = 3; y_check[0][4] = 4; y_check[1][1] = 7; y_check[2][2] = 9; y_check[3][0] = 11; y_check[3][3] = 13; y_check[3][4] = 15;
            
                for(int k = 0; k < testNAG.coldim(); ++k) {
                        report << "Checking row " << k + 1 << endl;
                        if(k > 0) x[k-1] = 0;
                        x[k] = 1;
                        testNAG.apply(y,x);
        
                        // Now checks y against y_check
                        for(int i = 0; i < 5; ++i) 
                                if(y[i] != y_check[k][i]) {
                                        res = false;
                                        report << "ERROR:  expected " << y_check[k][i] << " at position (" << k << "," << i << "), got " y[i] << ".  Not cool." << endl;
                                }
                }
        
                if(res == false) {
                        cout << "FAILED." << endl;
                        report << "Form Test. . .FAILED." << endl;
                } else {
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
                } else {
                        cout << "PASSED." << endl;
                        report << "FUNCTIONAL TEST. . .PASSED." << endl;
                }
                return res1 && res2;
        }
#endif // ifdef check
     
} // namespace LinBox

#endif // __NAG_SPARSE_H
