/** -*- mode:C++ -*- */
/** File: random-matrix.h
 *  Author: Zhendong Wan
 */

/** randomMatrix generates random matrcies used in EGV and EGV+ algorithm
 */

#ifndef __LINBOX_RANDOM_MATRIX_H__
#define __LINBOX_RANDOM_MATRIX_H__

#include <linbox/blackbox/dense.h>
#include <linbox/blackbox/sparse.h>

namespace LinBox {
	
	class RandomMatrix {

	public:

		/** @memo General case.
		 */
		template <class Blackbox, class Field>
			static Blackbox*& randomMatrix (Blackbox* &, const Field& f, 
							int rowdim, int coldim);

		
		// constructor a random dense matrix, whose entries are random
		template<class Field>
			static DenseMatrix<Field>*& randomMatrix( DenseMatrix<Field>*& Ap,
								  const Field& f,
								  int rowdim, int coldim ) {
			
			Ap = new DenseMatrix<Field>(f, rowdim, coldim);
			
			typename DenseMatrix<Field>::RawIterator Ap_p;
			
			for (Ap_p = Ap -> rawBegin(); Ap_p != Ap -> rawEnd(); ++ Ap_p)
				f.init(*Ap_p, rand());
			

			return Ap;
		}

		// constructor a very special random sparse matrix
		template<class Field>
			static SparseMatrix<Field>*& randomMatrix( SparseMatrix<Field>*& Ap, 
								   const Field& f, 
								   int rowdim, int coldim) {
			
			Ap = new SparseMatrix<Field>(f, rowdim, coldim);
			
			const int m = rowdim < coldim ? rowdim : coldim;

			int i, j, k;

			typename Field::Element elt;

			f. init (elt, 1);

			for ( i = 0; i < m; ++ i) 
				
				Ap -> setEntry (i, i, elt);


			if ( m < rowdim ) {

				const int repeat = (rowdim - m) < 10 ? rowdim - m : 10;
				
				for ( i = m; i < rowdim; ++ i) {

					for ( k = 0; k < repeat; ++ k) {
						
						j = random() % coldim;
						
						f.init(elt, rand() % 3 + 1);
						
						Ap -> setEntry (i, j, elt);
					}
				}
			}

			else if ( m < coldim ) {

				int offset = coldim - m;

				const int repeat = offset < 10 ? offset : 10;
				
				for ( i = 0; i < rowdim; ++ i){
			       
					for ( k = 0; k < repeat; ++ k) {
						
						j = random() % offset + m;
							
						f.init(elt, random() % 3 + 1);
						
						Ap -> setEntry (i, j, elt);
					}
				}
			}

			else {}
			
			return Ap;

		}
	};
}

		
#endif
