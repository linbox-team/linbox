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
                public:

                        // Default constructor, do nothing.
                        ZeroOne();
                        // The real constructor
                        ZeroOne(Field F, Index* rowP, Index* colP, Index rows, Index cols, Index NNz);
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
			     size_t* _row, _col;
			 };

	                 RawIndexIterator &indexBegin() {
			   return RawIndexIterator(_rowP, _colP);
			 }

                 	 const RawIndexIterator &indexBegin() const {
			   return RawIndexIterator(_rowP, _colP);
			 }

	                 RawIndexIterator &indexEnd() {
			   return RawIndexIterator(_rowP + _nnz, _colP + _nnz);
			 }

	                 const RawIndexIterator &indexEnd() const {
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

                        Index *_rowP, *_colP, _rows, _cols, _nnz;
  
                        /* _apply is the generic apply utility funtion called by apply() and
                         * applyTranspose().  Walks through the non-zero elements of the Matrix and
                         * performs the proper calculation using the axpyin method defined by the
                         * field element above.
                         */

                        void _apply(Vector &, const Vector &, Index*, Index*) const;

        };

	/* Default constructor.  Not really useful, just kinda there.
	 */

	template<class Field, class Vector>
        ZeroOne<Field,Vector>::ZeroOne() {}


        /*  Constructor for the ZeroOne class.  This is the constructor that is
         * expected to be used.  To use it, you must pass in a field element that
         * will work over the data (F), pointers to the 3 arrays used by the NAGSparse
         * format (values, rowP, colP), the number of rows and columns (rows and
         * cols), the number of non-zero elements (NNz) and the ordering, which
         * defaults to 0 (no ordering implied).
         */

        template<class Field, class Vector>
        ZeroOne<Field, Vector>::ZeroOne(Field F, Index* rowP, Index* colP, Index rows, Index cols, Index NNz):
        _F(F), _rowP(rowP), _colP(colP), _rows(rows), _cols(cols),_nnz(NNz) {}

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
                ZeroOne<Field,Vector>* p = new ZeroOne<Field,Vector>(_F,_rowP,_colP,_rows,_cols,_nnz);
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
	  _apply(y,x,_rowP,_colP);
	  return y;
        }

        /* BlackBoxArchetype applyTranspose function.  Performs the y = ATx, where
         * y and x are vectors passed in applyTranspose(y,x), and A is the present
         * Matrix.  Returns a reference to y.
         */

        template<class Field, class Vector>
        Vector & ZeroOne<Field,Vector>::applyTranspose(Vector & y, const Vector & x) const
        {
	  _apply(y,x,_colP,_rowP);
	  return y;
        }

        /* _apply utility function.  applies vector x to this matrix.  However, this Matrix is all
	 * 1's and 0's, so all we have to do is for every non-zero entry, add the vector entry (as
	 * normally say y_entry += x_entry * A_entry, but since A_entry is always 0 or 1, this formula
	 * becomes y_entry += x_entry).
	 */
 

        template<class Field, class Vector>
        void ZeroOne<Field,Vector>::_apply(Vector & y, const Vector & x, Index* _i, Index* _j) const
        {
                typename Vector::iterator yp;
                typename Vector::const_iterator xp;
                Index* ip, *jp;
 
                // 0 out y.  Note, this implementation assumes a dense vector.
                for(yp = y.begin(); yp != y.end(); ++yp)
                         _F.init(*yp , 0);
		
                for( yp = y.begin(), xp = x.begin(), ip = _i, jp = _j; ip < _i + _nnz; ++ip, ++jp)
		      _F.addin( *(yp + *ip), *(xp + *jp) );


        }

        template<class Field, class Vector>
	size_t ZeroOne<Field, Vector>::nnz() const
	{
	        return _nnz;
	}

        /** BlackBox representation of the 0-1's matrix
	 *  
         * A 0-1's matrix is a matrix with all 0's and 1's as entries.  In
	 * a nag-spasre format, applies could be performed with lightening speed
	 * When initalizing this class, you only need to build 2 arrays of equal length:
	 * an array of the row indices for the non-zero (1's) entries, and an array of the column
	 * indices for the non-zero (1's) entries.
         */
/*
        template<class Vector>
        class ZeroOne<LinBox::Modular<uint32>,Vector>: public BlackboxArchetype<Vector> {

                typedef LinBox::Modular<uint32>::Element Element;
                typedef size_t Index;
                public:

                        // Default constructor, do nothing.
                        ZeroOne() {}
                        // The real constructor
                        ZeroOne(LinBox::Modular<uint32> F, Index* rowP, Index* colP, Index rows, Index cols, Index NNz) : 
			  _F(F), _rowP(rowP), _colP(colP), _rows(rows), _cols(cols), _nnz(NNz) {}
                         // Destructor, once again do nothing
                         ~ZeroOne() {}
*/

                        /** BlackBoxArchetype clone function.
                         * Creates a copy of the ZeroOne Matrix and passes a pointer to it.
                         * In this case it isn't too helpful
                         * as this clone will of course suffer from the "siamese twin" problem.
                         * The clonse created will point to the same data as the parent.
                         * @return pointer to a new NAG format blackbox
                         */

 //                       BlackboxArchetype<Vector>* clone() const;

                        /** BlackBoxArchetype apply function.
                         * Take constant vector x and
                         * vector y, and perform the calculation y = Ax.  Uses one of the three
                         * private utility functions. It calls the generalized utility function
                         * _apply if there is no special ordering, _fyapply if there is C_ordering
                         * or _fxapply if there is fortran_ordering
                         */

  //                      Vector & apply(Vector &, const Vector &) const; // y = Ax;

                        /* BlackBoxArchetype applyTranspose function. Take constant vector x and
                         * vector y, and perform the calculation y = ATx.  Uses one of the three
                         * private utility functions, in the manner described above.  Worthy of
                         * note is the fact that applyTranspose works by passing the column
                         * positions to the _apply functions as if they were rows, and row positions
                         * as if they were columns, as if the matrix had been transposed.
                         */

   //                     Vector & applyTranspose(Vector &, const Vector &) const; // y = ATx

                        /* BlackBoxArchetype rowdim function.  Passes back the number of rows of
                         * the matrix represented.  Note that's the number of rows of the matrix
                         * as if it were in dense format, not in it's actual representation here.
                        */

    //                    size_t rowdim() const;

                        /* BlackBoxArchetype coldim function.  Passes back the number of columns
                         * of the matrix represented.  Not much more to say about this.
                         */

     //                   size_t coldim() const;

                         /* Non blackbox function.  Tells the number of nonzero entries
                          */
      //                   size_t nnz() const;

                         /* RawIterator class.  Iterates straight through the values of the matrix
                          */
/*	                       class RawIterator {
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
			  };*/  

			 /* STL standard Begin and End functions.  Used to get
			  * the beginning and end of the data.  So that RawIterator
			  * can be used in algorithms like a normal STL iterator.
			  */
/*
	  	 	RawIterator rawBegin() { return RawIterator( 0, _F.init(_tmp, 1) ); }
			 RawIterator rawEnd() { return RawIterator( _nnz, _F.init(_tmp, 1) ); }
			 const RawIterator rawBegin() const { return RawIterator(0, _F.init(_tmp, 1) ); }
			 const RawIterator rawEnd() const { return RawIterator(_nnz, _F.init(_tmp, 1) ); } 
*/
			/* RawIndexIterator - Iterates through the i and j of the current element
			 * and when accessed returns an STL pair containing the coordinates
			 */
/*			 class RawIndexIterator {
			    public:
				 typedef std::pair<size_t, size_t> value_type;

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
			     size_t* _row, _col;
			 };

	                 RawIndexIterator &indexBegin() {
			   return RawIndexIterator(_rowP, _colP);
			 }

                 	 const RawIndexIterator &indexBegin() const {
			   return RawIndexIterator(_rowP, _colP);
			 }

	                 RawIndexIterator &indexEnd() {
			   return RawIndexIterator(_rowP + _nnz, _colP + _nnz);
			 }

	                 const RawIndexIterator &indexEnd() const {
                 	   return RawIndexIterator(_rowP + _nnz, _colP + _nnz);
			   } 

                    private:
                        LinBox::Modular<uint32> _F; */ // The field used by this class
	  /* A temporary element used for initalization for the rawBegin() and
	   * rawEnd() methods of the ZeroOne class.  Is used to initalize a 1
	   * so that the RawIterator returned stores a 1 
	   */

	                //Element _tmp; 


                        /* _rowP is a pointer to an array of row indexes.  _colP is a pointer
                         * to an array of column indexes. These two are the other arrays of a
                         * NAGSparse format Matrix.  _rows and _cols are the number of rows and
                         * columns of the Matrix if it were in dense format.  _nnz is the Number of
                         * Non-Zero elements in the Matrix.  It also happens to be the length of
                         * the three NAGSparse arrays.
                          */

//                        Index *_rowP, *_colP, _rows, _cols, _nnz;
  
                        /* _apply is the generic apply utility funtion called by apply() and
                         * applyTranspose().  Walks through the non-zero elements of the Matrix and
                         * performs the proper calculation using the axpyin method defined by the
                         * field element above.
                         */

//                        void _apply(Vector &, const Vector &, Index*, Index*) const;

 //       };


        /* BlackBoxArchetype clone function.  Creates a another NAGSparse Matrix
         * and returns a pointer to it.  Very simple in construction, just uses the
         * new operator.  Of course needs to be deleted to prevent a memory leak.
         * Note, the BlackBox created by this clone function is not an independant
         * entity.  A NAGSparse is little more than a wrapper over a pre-created
         * NAG Sparse Matrix that allows the linbox algorithms to work on that matrix.
         * Thus, this constructor creates another wrapper for the same data pointed to
         * by this object.
         */

/*        template<class Field, class Vector>
        BlackboxArchetype<Vector>* ZeroOne<Field,Vector>::clone() const
        {
                ZeroOne<Field,Vector>* p = new ZeroOne<Field,Vector>(_F,_rowP,_colP,_rows,_cols,_nnz);
                return p;
        }
*/
        /* BlackBoxArchetype rowdim function.  Not much to say here. Returns the number
         * of rows of the Matrix were it in dense format.
         */

 /*       template<class Field, class Vector>
        size_t ZeroOne<Field,Vector>::rowdim() const
        {
                return _rows;
        }*/

        /* BlackBoxArchetype coldim function.  Not much to say here either. Returns the
         * number of columns of the Matrix were it in dense format.
         */

/*        template<class Field, class Vector>
        size_t ZeroOne<Field,Vector>::coldim() const
        {
                return _cols;
        }*/

        /* BlackBoxArchetype apply function.  Performs the y = Ax calculation, where
         * y and x are vectors passed in apply(y,x), and A is the present matrix.
	 *
         */

/*        template<class Field, class Vector>
        Vector & ZeroOne<Field,Vector>::apply(Vector & y, const Vector & x) const
        {
	  _apply(y,x,_rowP,_colP);
	  return y;
        }
*/
        /* BlackBoxArchetype applyTranspose function.  Performs the y = ATx, where
         * y and x are vectors passed in applyTranspose(y,x), and A is the present
         * Matrix.  Returns a reference to y.
         */

 /*       template<class Field, class Vector>
        Vector & ZeroOne<Field,Vector>::applyTranspose(Vector & y, const Vector & x) const
        {
	  _apply(y,x,_colP,_rowP);
	  return y;
        }
*/
        /* _apply utility function.  applies vector x to this matrix.  However, this Matrix is all
	 * 1's and 0's, so all we have to do is for every non-zero entry, add the vector entry (as
	 * normally say y_entry += x_entry * A_entry, but since A_entry is always 0 or 1, this formula
	 * becomes y_entry += x_entry).
	 */
 

/*        template<class Field, class Vector>
        void ZeroOne<Field,Vector>::_apply(Vector & y, const Vector & x, Index* _i, Index* _j) const
        {
                typename Vector::iterator yp;
                typename Vector::const_iterator xp;
                Index* ip, *jp;
 
                // 0 out y.  Note, this implementation assumes a dense vector.
                for(yp = y.begin(); yp != y.end(); ++yp)
                         _F.init(*yp , 0);
		
                for( yp = y.begin(), xp = x.begin(), ip = _i, jp = _j; ip < _i + _nnz; ++ip, ++jp)
		      _F.addin( *(yp + *ip), *(xp + *jp) );


        }

        template<class Field, class Vector>
	size_t ZeroOne<Field, Vector>::nnz() const
	{
	        return _nnz;
	}
*/		
#ifdef check
#include <iostream>
#include <fstream>
#include <test_common.h>
#include <ctime>
#include "linbox/randiter/archetype.h"

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
                bool res = true;
		int i, j, k;
		size_t rows[400], cols[400];
		Vector x, y1, y2;
		Vector::iterator y1i, y2i;
		LinBox::RandIterArchetype gen(F, 0, time(0) ); // initalize the random number generator

                cout << "ZeroOne blackbox suite." << endl;
                report << "ZeroOne blackbox suite." << endl;

                cout << "Form Test. . .";
                report << "Form Test. . .";
		report << "Generating a 100x100 zero-one Matrix with 1s on the main diagonal," << endl;
		report << "the reverse diagonal (100,1) up to (1,100), down a column (i,12)," << endl;
		report << "and across the center row (50,i)" << endl;
                report << "Now creating Matrix." << endl;

		for(i = 0; i < 100; i++) {
		  rows[4 * i] = cols[4 * i] = i; // The diagonal

		  rows[4 * i + 1] = i;         // The reverse diagonal
		  cols[4 * i + 1] = 99 - i;

		  rows[4 * i + 2] = 49;        // The center row
		  cols[4 * i + 2] = i;

		  rows[4 * i + 3] = i;        // The center column
		  cols[4 * i + 3] = 12;
		}
		  
		LinBox::ZeroOne<Field, Vector> testMatrix(F, rows, cols, 100, 100, 400);


                report << "Testing rowdim: " << endl;
		if( testMatrix.rowdim() != 100) {
		  report << "rowdim FAIL!" << endl;
		  report << "Was expecting 100, recieved " << testMatrix.rowdim() << "." << endl; 
		  res = false;
		}

		report << "Testing coldim: " << endl;
		if( testMatrix.coldim() != 100) {
		  report << "coldim FAIL!" <<  endl;
		  report << "Was expecting 100, recieved " << testMatrix.coldim() << "." << endl;
		  res = false;
		}


                report << "Testing apply:" << endl;
		report << "Creating random vector of elements." << endl;
		for(i = 0; i < 100; ++i) {
		  // Create a new random entry in from within the field
		  x.push_back( gen.random(blank) );
		}

		F.init(blank);
		y1.assign(100, blank); // ensure proper size of y values
		y2.assign(100,blank); 

		report << "Performing a manual application of Ax." << endl;
		for(i = 0; i < 100; i++) {
		  for(j = 0; j < 100; j++) {
		    if( i == j) {
		      F.addin(y1[i], x[j]);
		    }
		    else if(j == 99 - i) {
		      F.addin(y1[i], x[j]);
		    }
		    else if( i == 49 ) {
		      F.addin(y1[i], x[j]);
		    }
		    else if( j == 12) {
		      F.addin(y1[i], x[j]);
		    }
		  }
		}

		// Now perform the apply
		report << "Performing application Ax using blackbox." << endl;
		testMatrix.apply(y2, x);

		// Now compare the two
		for( y1i = y1.begin(), y2i = y2.begin(); y1i != y1.end(); ++y1i) {
		  if( *y1i != *y2i) {
		    report << "Apply test FAIL!" << endl;
		    report << "Was expecting " << *y1i << ", but recieved " << *y2i << "." << endl;
		    res = false;
		  }
		}

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
		    
		    else if(i == 12) {
		      F.addin(y1[i], x[j]);
		    }
		    
		    else if(j == 49) {
		      F.addin(y1[i], x[j]);
		    }
		  }
		}

		// Now perform the apply transpose
		testMatrix.applyTranspose(y2, x);

		// Now go through and compare the two
		for(Vector::iterator y1i = y1.begin(), Vector::iterator y2i = y2.begin(); y1i != y1.end(); ++y1i, ++y2i) {
		  if( *y1i != *y2i) {
		    report << "applyTranspose test FAIL!" << endl;
		    report << "Was expecting " << *y1i << ", but recived " << *y2i << " instead." << endl;
		    res = false;
		  }
		}

		/* Now tests Raw and Index Iterators 
		report << "Checking Raw and Index Iterators." << endl;
		report << "Running pre-increment test." << endl;

		... code goes here ...

		*/

		// Now tests the RawIndexIterators using the same style of
		// linear tests as above.  This time checks the index pair
		// returned by the RawIndexIterator against the two row &
		// column initialization arrays used above
		report << "RawIterator test finished.  Now testing RawIndexIterator." << endl;
		LinBox::ZeroOne<Field, Vector> ri_i = testMatrix.indexBegin();
		std::pair<size_t, size_t> testPair;
		report << "First running test using pre-increment." << endl;
		for(k = 0; ri_i != testNAG.indexEnd(); ++k, ++ri_i) {
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
		
		if( k != 7) {
		        res = false;
			report << "ERROR: Supposed to process 8 indices.  Only processed " << k + 1 << "." << endl;
		}

		// Now repeats the test using, you guessed it, postincrement
		report << "Now running Index test using postincrement." <<endl;
		for(k = 0; ri_i != testNAG.indexEnd(); k++, ri_i++) {
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
		
		if( k != 400) {
		        res = false;
			report << "ERROR: Supposed to process 400 indice pairs.  Only processed " << k + 1 << "." << endl;
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

                if( BBGeneralTest(testNAG, F, report) == false) {
		        res = false;
                        cout << "FAILED " << endl;
                        report << "FUNCTIONAL TEST. . .FAILED" << endl;
                } else {
                        cout << "PASSED." << endl;
                        report << "FUNCTIONAL TEST. . .PASSED." << endl;
                }
                return res;
        }
#endif // ifdef check

} // namespace LinBox

#endif // __ZERO_ONE_H
