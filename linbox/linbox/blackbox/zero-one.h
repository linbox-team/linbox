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

#ifndef __ZERO_ONE_H
#define __ZERO_ONE_H

#include "linbox/blackbox/archetype.h"
#include "linbox/vector/vector-traits.h"
#include "linbox/util/debug.h"
#include "linbox/field/modular.h"

// For STL pair in RawIndexIterator
#include <utility>
#include <vector> // For vectors in _col2row and _row2col
#include <cstdlib> // For randomness in randomized quicksort
#include <ctime>

using std::vector;

namespace LinBox
{

        /** BlackBox representation of the 0-1's matrix
	 *  
         * A 0-1's matrix is a matrix with all 0's and 1's as entries.  In
	 * a nag-spasre format, applies could be performed with lightening speed
	 * When initalizing this class, you only need to build 2 arrays of equal length:
	 * an array of the row indices for the non-zero (1's) entries, and an array of the column
	 * indices for the non-zero (1's) entries.
         */

        template<class Field, class Vector>
        class ZeroOne: public BlackboxArchetype<Vector> {

                typedef typename Field::Element Element;
                typedef size_t Index;
	        typedef LinBox::uint32 uint32;
	        typedef LinBox::uint64 uint64;
                public:

                        // Default constructor, do nothing.
                        ZeroOne();
                        // The real constructor
                        ZeroOne(Field F, Index* rowP, Index* colP, Index rows, Index cols, Index NNz, bool rowSort = false, bool colSort = false);
                         // Destructor, once again do nothing
                         ~ZeroOne() {};


                        /** BlackBoxArchetype clone function.
                         * Creates a copy of the ZeroOne Matrix and passes a pointer to it.
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

                         /* Non blackbox function.  Tells the number of nonzero entries
                          */
                         size_t nnz() const;

                 	  void RowSort() const;
	                  void ColSort() const;


                         /* RawIterator class.  Iterates straight through the values of the matrix
                          */
	                       class RawIterator {
	                         public:
	                         	typedef Element value_type;

	                         	RawIterator(size_t pos, Element elem) :
	                         		_pos(pos), _elem(elem) {}

	                         	RawIterator(const RawIterator &In) :
	                         		_pos(In._pos), _elem(In._elem) {}

	                      		const RawIterator& operator=(const RawIterator& rhs) {
					        _pos = rhs._pos;
						_elem = rhs._elem;
		                      		return *this;
	                      		}


	                 		bool operator==(const RawIterator &rhs) {
		                 		return ( _pos == rhs._pos && _elem == rhs._elem);
	                 		}

	               			bool operator!=(const RawIterator &rhs) {
						return ( _pos != rhs._pos || _elem != rhs._elem );
					}

					RawIterator & operator++() {
						++_pos;
						return *this;
					}

					RawIterator operator++(int) {
						
					        RawIterator tmp = *this;
					        _pos++;
						return tmp;
					}

					value_type operator*() {
						return _elem;
					}

					const value_type operator*() const {
						return _elem;
					}

				private:
					value_type _elem;
			                size_t _pos;
			  };  

			 /* STL standard Begin and End functions.  Used to get
			  * the beginning and end of the data.  So that RawIterator
			  * can be used in algorithms like a normal STL iterator.
			  */

      	  	 	RawIterator rawBegin() { return RawIterator( 0, _F.init(_tmp, 1) ); }
			 RawIterator rawEnd() { return RawIterator( _nnz, _F.init(_tmp, 1) ); }
			 const RawIterator rawBegin() const { return RawIterator(0, _F.init(_tmp, 1) ); }
			 const RawIterator rawEnd() const { return RawIterator(_nnz, _F.init(_tmp, 1) ); } 

			/* RawIndexIterator - Iterates through the i and j of the current element
			 * and when accessed returns an STL pair containing the coordinates
			 */
			 class RawIndexIterator {
			    public:
				 typedef std::pair<size_t, size_t> value_type;

 			         RawIndexIterator() {}

				 RawIndexIterator(size_t* row, size_t* col):
				 	_row(row), _col(col) {}

				 RawIndexIterator(const RawIndexIterator &In):
				 	_row(In._row), _col(In._col) {}

				 const RawIndexIterator &operator=(const RawIndexIterator &rhs) {
					_row = rhs._row;
					_col = rhs._col;
					return *this;				
				}

				bool operator==(const RawIndexIterator &rhs) {
					_row == rhs._row && _col == rhs._col;	

				}

				bool operator!=(const RawIndexIterator &rhs) {
					return _row != rhs._row || _col != rhs._col;
				}

				const RawIndexIterator& operator++() {
					++_row; ++_col;	
					return *this;
				}

				const RawIndexIterator operator++(int) {
					RawIndexIterator tmp = *this;
					++_row; ++_col;
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
			 };

	                 RawIndexIterator indexBegin() {
			   return RawIndexIterator(_rowP, _colP);
			 }

                 	 const RawIndexIterator indexBegin() const {
			   return RawIndexIterator(_rowP, _colP);
			 }

	                 RawIndexIterator indexEnd() {
			   return RawIndexIterator(_rowP + _nnz, _colP + _nnz);
			 }

	                 const RawIndexIterator indexEnd() const {
                 	   return RawIndexIterator(_rowP + _nnz, _colP + _nnz);
			   } 

                    private:
                        Field _F; // The field used by this class
	  /* A temporary element used for initalization for the rawBegin() and
	   * rawEnd() methods of the ZeroOne class.  Is used to initalize a 1
	   * so that the RawIterator returned stores a 1 
	   */

	                Element _tmp; 


                        /* _rowP is a pointer to an array of row indexes.  _colP is a pointer
                         * to an array of column indexes. These two are the other arrays of a
                         * NAGSparse format Matrix.  _rows and _cols are the number of rows and
                         * columns of the Matrix if it were in dense format.  _nnz is the Number of
                         * Non-Zero elements in the Matrix.  It also happens to be the length of
                         * the three NAGSparse arrays.
                          */

                        Index _rows, _cols, _nnz;
	                mutable Index* _rowP, *_colP;
	                mutable bool _rowSort, _colSort; // status flags for sorting state

  
                        /* _apply has two versions:  general, and template specialized for Modular<uint32>
			 * part of a crazy scheme motivated by BlackBox's virtual nature
                         */
	                template<class theField>
                        Vector & _apply(theField &, Vector &, const Vector &, Index*, Index*) const;

	                Vector & _apply (LinBox::Modular<uint32> &, Vector &, const Vector &, Index*, Index*) const;

	                void _qsort(size_t start, size_t endp1, int &mode) const; // QuickSort function for when there is no sorting
	                size_t _part( size_t start, size_t endp1, int &mode) const; // Partition for quicksort
	                void _row2col() const; // Quick sorting from row sorted state to column sorted state
	                void _col2row() const; // Quick sorting from row sorted state to column sorted state

        };

	/* Default constructor.  Not really useful, just kinda there.
	 */

	template<class Field, class Vector>
        ZeroOne<Field,Vector>::ZeroOne() { srand( time(NULL) ); }


        /*  Constructor for the ZeroOne class.  This is the constructor that is
         * expected to be used.  To use it, you must pass in a field element that
         * will work over the data (F), pointers to the 3 arrays used by the NAGSparse
         * format (values, rowP, colP), the number of rows and columns (rows and
         * cols), the number of non-zero elements (NNz) and the ordering, which
         * defaults to 0 (no ordering implied).
         */

        template<class Field, class Vector>
        ZeroOne<Field, Vector>::ZeroOne(Field F, Index* rowP, Index* colP, Index rows, Index cols, Index NNz, bool rowSort, bool colSort):
        _F(F), _rows(rows), _cols(cols), _nnz(NNz), _rowP(rowP), _colP(colP), _rowSort(rowSort), _colSort(colSort) { srand(time(NULL)); }

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
        BlackboxArchetype<Vector>* ZeroOne<Field,Vector>::clone() const
        {
                ZeroOne<Field,Vector>* p = new ZeroOne<Field,Vector>(_F,_rowP,_colP,_rows,_cols,_nnz, _rowSort, _colSort);
                return p;
        }

        /* BlackBoxArchetype rowdim function.  Not much to say here. Returns the number
         * of rows of the Matrix were it in dense format.
         */

        template<class Field, class Vector>
        size_t ZeroOne<Field,Vector>::rowdim() const
        {
                return _rows;
        }

        /* BlackBoxArchetype coldim function.  Not much to say here either. Returns the
         * number of columns of the Matrix were it in dense format.
         */

        template<class Field, class Vector>
        size_t ZeroOne<Field,Vector>::coldim() const
        {
                return _cols;
        }

        /* BlackBoxArchetype apply function.  Performs the y = Ax calculation, where
         * y and x are vectors passed in apply(y,x), and A is the present matrix.
	 *
         */

        template<class Field, class Vector>
        Vector & ZeroOne<Field,Vector>::apply(Vector & y, const Vector & x) const
        {
	  return _apply(_F, y,x,_rowP,_colP);
        }

        /* BlackBoxArchetype applyTranspose function.  Performs the y = ATx, where
         * y and x are vectors passed in applyTranspose(y,x), and A is the present
         * Matrix.  Returns a reference to y.
         */

        template<class Field, class Vector>
        Vector & ZeroOne<Field,Vector>::applyTranspose(Vector & y, const Vector & x) const
        {
	  _apply(_F, y,x,_colP,_rowP);
	  return y;
        }

  template<class Field, class Vector>
  void ZeroOne<Field, Vector>::RowSort() const
  {
    int mode = 0;
    if( _rowSort) return;  // Already sorted, we're done
    if(_colSort) _col2row();
    else _qsort( (size_t) 0, _nnz, mode);
    _rowSort = true; _colSort = false;
    return;
  }

  template<class Field, class Vector>
  void ZeroOne<Field, Vector>::ColSort() const
  {
    int mode = 1;
    if( _colSort) return; // Already sorted, good to go
    if( _rowSort) _row2col(); // Just need to run the switch sorting algorithm
    else _qsort( (size_t) 0, _nnz, mode);
    _colSort = true; _rowSort = false;
    return;
  }

  template<class Field, class Vector>
  void ZeroOne<Field, Vector>::_row2col() const
  {
    size_t* rowp, *colp, j;
    vector<size_t>::iterator i;
    vector<vector<size_t> > R(_cols);
    for(rowp = _rowP, colp = _colP; rowp != _rowP + _nnz; ++rowp, ++colp) 
      R[*colp].push_back(*rowp);
    rowp = _rowP; colp = _colP;
    for(j = 0; j < _cols; ++j)
      for(i = R[j].begin(); i < R[j].end(); ++i, ++rowp, ++colp)
      {
	*rowp = *i;
	*colp = j;
      }
    return;
  }

  template<class Field, class Vector>
  void ZeroOne<Field, Vector>::_col2row() const
  {
    size_t* rowp, *colp, j;
    vector<size_t>::iterator i;
    vector<vector<size_t> > R(_rows);
    for(rowp = _rowP, colp = _colP; rowp != _rowP + _nnz; ++rowp, ++colp)
      R[*rowp].push_back(*colp);
    rowp = _rowP; colp = _colP;
    for(j = 0; j < _rows; ++j)
      for(i = R[j].begin(); i < R[j].end(); ++i, ++rowp, ++colp)
      {
	*colp = *i;
	*rowp = j;
      }

    return;
  }

  template<class Field, class Vector>
  void ZeroOne<Field, Vector>::_qsort(size_t p, size_t e, int &mode) const
  {
    int i;
    if( (e - p) <= 1) ;
    else {
      i = 1 + _part(p, e, mode);
      _qsort(p, i, mode);
      _qsort(i, e, mode);
    }
  }

  template<class Field, class Vector>
  size_t ZeroOne<Field, Vector>::_part(size_t p, size_t e, int &mode) const
  {
    size_t rtemp, ctemp, rowval, colval;
    int i = p + rand() % (e - p), j = e;
    rtemp = _rowP[p];
    ctemp = _colP[p];
    _rowP[p] = _rowP[i];
    _colP[p] = _colP[i];
    _rowP[i] = rtemp;
    _colP[i] = ctemp;
    rowval = _rowP[p];
    colval = _colP[p];
    i = p - 1;

    if(mode == 0) { // Row mode, go by row order, then column
      while(true) {
	do j--; while( _rowP[j] > rowval || ( _rowP[j] == rowval && _colP[j] > colval ));
	do i++; while( _rowP[i] < rowval || ( _rowP[i] == rowval && _colP[i] < colval ));
	if( i < j) {
	  rtemp = _rowP[j];
	  ctemp = _colP[j];
	  _rowP[j] = _rowP[i];
	  _colP[j] = _colP[i];
	  _rowP[i] = rtemp;
	  _colP[i] = ctemp;
	}
	else return j;
      }
    }
    else { // Col mode, go by col order, then row
      while(true) {
	do j--; while( _colP[j] > colval || ( _colP[j] == colval && _rowP[j] > rowval ));
	do i++; while( _colP[i] < colval || ( _colP[i] == colval && _rowP[i] < rowval ));
	if( i < j) {
	  rtemp = _rowP[j];
	  ctemp = _colP[j];
	  _rowP[j] = _rowP[i];
	  _colP[j] = _colP[i];
	  _rowP[i] = rtemp;
	  _colP[i] = ctemp;
	}
	else return j;
      }
    }
  }




        /* _apply utility function.  applies vector x to this matrix.  However, this Matrix is all
	 * 1's and 0's, so all we have to do is for every non-zero entry, add the vector entry (as
	 * normally say y_entry += x_entry * A_entry, but since A_entry is always 0 or 1, this formula
	 * becomes y_entry += x_entry).
	 */
           template<class Field, class Vector>
	   Vector & ZeroOne<Field,Vector>::_apply(LinBox::Modular<uint32> &F, Vector & y, const Vector & x, Index* _i, Index* _j) const
	 {
	   typename Vector::iterator yp = y.begin();
	   typename Vector::const_iterator xp = x.begin();
	   size_t *ip = _i, *jp = _j, i;
	   uint64 accum;

	   if( _i == _rowP) {
	     RowSort();
	     for(i = 1; i <= _rows; i++) {
	       accum = 0;
	       while( *ip == i) {
		 accum += *(xp + *jp);
		 ip++; jp++;
	       }
	       *(yp + i) = accum % F._modulus;
	     }
	   }
	   else {
	     ColSort();
	     for(i = 1; i <= _cols; i++) {
	       accum = 0;
	       while( *ip == i) {
		 accum += *(xp + *jp);
		 ip++; jp++;
	       }
	       *(yp + i) = accum % F._modulus;
	     }
	   }

	   return y;
	   
	 }
 


        template<class Field, class Vector>
	template<class theField>
        Vector & ZeroOne<Field,Vector>::_apply(theField & F, Vector & y, const Vector & x, Index* _i, Index* _j) const
        {
                typename Vector::iterator yp;
                typename Vector::const_iterator xp;
                Index* ip, *jp;
 
                // 0 out y.  Note, this implementation assumes a dense vector.
                for(yp = y.begin(); yp != y.end(); ++yp)
                         F.init(*yp , 0);
		
                for( yp = y.begin(), xp = x.begin(), ip = _i, jp = _j; ip < _i + _nnz; ++ip, ++jp)
		      F.addin( *(yp + *ip), *(xp + *jp) );

		return y;

        }

        template<class Field, class Vector>
	size_t ZeroOne<Field, Vector>::nnz() const
	{
	        return _nnz;
	}

       




#ifdef check
#include <iostream>
#include <fstream>
#include "tests/test_common.h"
#include "linbox/randiter/archetype.h"

        /* Linbox testing routine.  Tests the class to ensure that it functions
         * properly.  Uses the standard diagnostic model used by all Linbox BlackBox's
         * Runs 3 tests, a form test to ensure that results produced are proper, and a
         * function test to ensure that all BlackBoxArchetype API works properly.
         */

  template<class Vector,class Field>
  bool testZeroOneBlackBox(ostream &report, Field &F, const Vector &x, size_t n, size_t iter)
{
  typedef typename Field::Element Element;
  Element blank;
  bool res = true;
  int i, j, k, h;
  size_t *rows, *cols, n = 100, iter = 1;
  Vector y1, y2, x;
  typename Vector::iterator y1i, y2i;
  LinBox::RandIterArchetype gen(F);

  x.resize(n);
  for(xp = x.begin(); xp != x.end(); xp++) {
    gen.random(*xp);
  }


  report << "ZeroOne blackbox suite." << endl;
  
  report << "Form Test. . .";
  report << "Generating a " << n << "x" << n << " zero-one Matrix with 1s on the main diagonal, and" << endl;
  report << "the reverse diagonal." << endl;
  report << "Now creating Matrix." << endl;
  
  rows = new size_t[2 * n];
  cols = new size_t[2 * n];

  for(i = 0; i < n; i++) {
    rows[2 * i] = cols[2 * i] = i; // The diagonal
    
    rows[2 * i + 1] = i;         // The reverse diagonal
    cols[2 * i + 1] = 99 - i;
    
  }
  
  LinBox::ZeroOne<Field, Vector> testMatrix(F, rows, cols, n, n, 2 * n);

  for(h = 0; h < iter; h++) {
    report << "Now running test iteration: " << h + 1 << endl;

    report << "Testing rowdim: " << endl;
    if( testMatrix.rowdim() != n) {
      report << "rowdim FAIL!" << endl;
      report << "Was expecting " << n << ", recieved " << testMatrix.rowdim() << "." << endl; 
      res = false;
    }
  
    report << "Testing coldim: " << endl;
    if( testMatrix.coldim() != n) {
      report << "coldim FAIL!" <<  endl;
      report << "Was expecting " << n << ", recieved " << testMatrix.coldim() << "." << endl;
      res = false;
    }


    report << "Testing apply:" << endl;
    F.init(blank);
    y1.assign(100, blank); // ensure proper size of y values
    y2.assign(100,blank); 
  
    report << "Performing a manual application of Ax." << endl;
    for(i = 0; i < 100; i++) {
      for(j = 0; j < 100; j++) {
	if( i == j) {
	  F.addin( y1[i], x[j] );
	}
	else if(j == 99 - i) {
	  F.addin( y1[i], x[j]  );
	}
      }
    }
    // Now perform the apply
    report << "Performing application Ax using blackbox." << endl;
    testMatrix.apply(y2, x);
  
    // Now compare the two
    k = 0;
    for( i = 0; i < 100; i++ ) {
      if( y1[i] != y2[i]) {
	report << "Apply test FAIL!" << endl;
	report << "Was expecting " << y1[i] << ", but recieved " << y2[i] << "." << endl;
	report << "At location: " << i << endl;
	res = false;
	k++;
      }
    }
    if(k > 0) report << "Total number of erroneous vector entries: " << k << endl;

    report << "Testing applyTranspose" << endl;
    // 0's out y for the applyTranspose test
    y1.assign(100,blank);
    for(i = 0; i < 100; i++) {
      for(j = 0; j < 100; j++) {
	
	if(i == j) {
	  F.addin(y1[i], x[j]);
	}
      
	else if(j == 99 - i) {
	F.addin(y1[i], x[j]);
	}
      }
    }
  
    // Now perform the apply transpose
    testMatrix.applyTranspose(y2, x);
  
    k = 0;
    // Now go through and compare the two
    for( i = 0; i < 100; i++) {
      if( y1[i] != y2[i] ) {
	report << "applyTranspose test FAIL!" << endl;
	report << "Was expecting " << y1[i] << ", but recived " << y2[i] << " instead." << endl;
      res = false;
      report << "At location: " << i << endl; 
      k++;
      }
    }
    
    if( k > 0) report << "Total numberof erroneous vector entries: " << k << endl;

    report << "Checking Raw and Index Iterators." << endl;
    report << "Running pre-increment test." << endl;
    LinBox::ZeroOne<Field, Vector>::RawIterator ra_i = testMatrix.rawBegin();
    F.init(blank, 1); 
    for(i = 0; ra_i != testMatrix.rawEnd(); ++i, ++ra_i) {
      if(blank != *ra_i) {
	report << "Error, was expecting a 1, but got " << *ra_i << " instead." << endl;
      }
    }
    if( i != 2 * n) {
      report << "ERROR, should have processed " << 2 * n << " elements, processed " << i << " instead." << endl;
    }
  
    report << "Now running post-increment test." << endl;
    ra_i = testMatrix.rawBegin();
    for(i = 0; ra_i != testMatrix.rawEnd(); ++i, ra_i++) {
      if( blank != *ra_i) {
	report << "Error, was expecting a 1, but got " << *ra_i << " instead." << endl;
      }
    }
    if( i != 2 * n) {
      cout << "Error, should have processed " << 2 * n << " elements, processed " << i << " instead." << endl;
    }   

    // Now tests the RawIndexIterators using the same style of
    // linear tests as above.  This time checks the index pair
    // returned by the RawIndexIterator against the two row &
    // column initialization arrays used above
    report << "RawIterator test finished.  Now testing RawIndexIterator." << endl;
    LinBox::ZeroOne<Field, Vector>::RawIndexIterator ri_i = testMatrix.indexBegin();
    pair<size_t, size_t> testPair;
    report << "First running test using pre-increment." << endl;
    for(k = 0; ri_i != testMatrix.indexEnd(); ++k, ++ri_i) {
      testPair = *ri_i;
      if( testPair.first != rows[k] ) {
	res = false;
        report << "ERROR: expected row " << k << " of RawIndexIterator to be " << rows[k] << ", recieved " << testPair.first << "." << endl;
      }
      if( testPair.second != cols[k] ) {
	res = false;
	report << "ERROR: expect cols " << k << " of RawIndexIterator to be " << cols[k] << ", recieved " << testPair.second << "." << endl;
      }
    }
  
    if( k != 2 * n) {
      res = false;
      report << "ERROR: Supposed to process " << 2 * n << " indices.  Only processed " << k + 1 << "." << endl;
    }
    
    // Now repeats the test using, you guessed it, postincrement
    report << "Now running Index test using postincrement." <<endl;
    ri_i = testMatrix.indexBegin();
    for(k = 0; ri_i != testMatrix.indexEnd(); k++, ri_i++) {
      testPair = *ri_i;
      if( testPair.first != rows[k] ) {
	res = false;
	report << "RawIndexIterator test FAIL." << endl;
	report << "Expected row " << k << " of RawIndexIterator to be " << rows[k] << ", recieved " << testPair.first << "." << endl;
      }
      if( testPair.second != cols[k] ) {
	res = false;
	report << "RawIndexIterator test FAIL." << endl;
	report << "Expected col " << k << " of RawIndexIterator to be " << cols[k] << ", recieved " << testPair.second << "." << endl;
      }
    }
  
    if( k != 2 * n) {
      res = false;
      report << "ERROR: Supposed to process " << 2 * n << " indice pairs.  Only processed " << k << "." << endl;
    }
  

  }
  

  if(res == false) {
    report << "Form Test. . .FAILED." << endl;
  }
  else {
    report << "Form Test. . .PASSED." << endl;
  }
  delete [] rows;
  delete [] cols;

  return res;
}
#endif // ifdef check

} // namespace LinBox

#endif // __ZERO_ONE_H
