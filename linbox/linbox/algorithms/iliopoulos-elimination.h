/** -*- mode:C++ -*- */

/** File iliopoulos-elimination.h
 *  Author: Zhendong Wan
 */
/** Compute Smith Form by elimination modulus S(n), the last invariant factor
 *  The elimination method is origanl described in Worst Case Complexity Bounds on Algorithms for computing the Canonical 
 *  Structure of Finit Abelian Groups and The Hermite and Smith Normal Fors of An Integer Matrix, by Costas Iliopoulos
 */

#ifndef __ILIOPOULOS_ELILIMINATION_H__
#define __ILIOPOULOS_ELILIMINATION_H__

#include <linbox/util/debug.h>
#include <linbox/vector/vector-domain.h>
#include <linbox/blackbox/dense.h>
#include <linbox/blackbox/submatrix-traits.h>

namespace LinBox {
	
	/** @memo Iliopoulos' algorithm
	 * Diagonalize
	 */

	class IliopoulosElimination {
		

	protected:
		/** @memo eliminationRow will make the first row (*, 0, ..., 0)
		 *  by col operations.
		 *  It is the implementation of Iliopoulos algorithm
		 */
		template<class Matrix, class Ring>
			static Matrix& eliminationRow (Matrix& A, const Ring& r,
						       const typename Ring::Element& d) {
#ifdef WANDEBUG
			std::cout << "In eliminationRow, Matrix:\n";
			A. write (std::cout);
#endif
			
			if (A. coldim() <= 1 || A. coldim()  == 0) return A;
			

			//typedef typename Matrix::Field Field;
			
			typedef typename Ring::Element Element;

			//Field r(A.field());

			VectorDomain<Ring> vd (r);

			// some tempory variable
			typename Matrix::RowIterator cur_r, tmp_r;
			
			typename Matrix::ColIterator cur_c, tmp_c;
			typename Matrix::Row::iterator row_p1, row_p2;
			//typename Matrix::Col::iterator col_p1, col_p2;
			
			cur_c = A. colBegin();

			cur_r = A. rowBegin();

			row_p1 = cur_r -> begin();

			Element g, s, t;
				
			r. xgcd (g, s, t, *row_p1, d);
			
			// if A[0][0] is coprime to d
			if (r. isUnit(g)) {

#ifdef WANDEBUG
				std::cout << "A[0][0] unit:\n";
#endif

				if (! r. isOne (s)) {
					
					vd. mulin(*cur_c, s);

					RemIn (*cur_c, r, d);
				}
				
			}
			
			// A[0][0] is not a unit
			else {
					
				// make A[0][0] = 0
				if ( !r.isZero(*row_p1)) {
					
					row_p2 = row_p1 + 1;
					
					Element y1, y2;

					Element g, s, t;
					
					r. xgcd (g, s, t, *row_p1, *row_p2);
					
					r. div (y1, *row_p2, g);
					
					r. div (y2, *row_p1, g);

					r.negin (y1);
										
#ifdef WANDEBUG
					std::cout <<"row_p2, y1 " << *row_p2 << " " << y1 <<"\n";
					
					std::cout <<"row_p1, y2 " << *row_p1 << " " << y2 <<"\n";
					
					std::cout <<" s row_p1 + t row_p2 should be " << g << " :" 
						  << (s * (*row_p1) + t * (*row_p2)) <<"\n";
					
					std::cout << "g, s, t, y1, y2, s y2 - t y1 (1): " <<
						g << " " << s << " " << t << " " <<
						y1 << " " << y2 << " " <<
						(s * y2 - t * y1) << "\n";

#endif
										
					tmp_c = cur_c + 1;
					
					std::vector<Element> tmp1 (A.rowdim()), tmp2 (A.rowdim());

					vd. mul (tmp1, *cur_c, y1);

					vd. axpyin (tmp1, y2, *tmp_c);

					vd. mul (tmp2, *cur_c, s);

					vd. axpyin (tmp2, t, *tmp_c);

					vd. copy (*cur_c,  tmp1);

					vd. copy (*tmp_c, tmp2);

					RemIn (*cur_c, r, d);

					RemIn (*tmp_c, r, d);
					
				}
#ifdef WANDEBUG
				std::cout << "A[0][0] should be 0:\n";
				A. write (std::cout);
#endif
				
				// matrix index is 0-based
				std::vector<Element> tmp_v(A.coldim());
				
				typename std::vector<Element>::iterator p1, p2;

				r. init(tmp_v[0], 1);

				p1 = tmp_v.begin() + 1;
				p2 = tmp_v.begin() + 1;
				
				row_p2 = row_p1 + 1; 
				
				Element g, s;
				
				r.assign(g, *row_p2); ++ row_p2;
				
				r.init(*p1,1); ++ p1;
				
				for (; row_p2 != cur_r -> end(); ++ row_p2, ++ p1) {
					
					r.xgcd(g,s,*p1,g,*row_p2);
					
					if (!r.isOne(s))
						for (p2 = tmp_v.begin() + 1; p2 != p1; ++ p2) {
							
							r. mulin (*p2, s);
							
							r. remin (*p2, d);
						}
					
				}
				
				// no pivot found
				if (r.isZero(g)) return A;
				
				Element tmp;

#ifdef WANDEBUG
				std::cout <<"g = " << g <<"\n";
				std::cout << "tmp_v:\n";
				
				for (p1 = tmp_v.begin(); p1 != tmp_v.end(); ++p1)
					std::cout << *p1 << " ";
				
				std::cout << '\n';
#endif
				
				for (tmp_r = cur_r; tmp_r != A.rowEnd(); ++ tmp_r) {
					
					vd. dot (tmp, *tmp_r, tmp_v);
					
					r. rem (*(tmp_r -> begin()), tmp, d);
				}

				
			}
				
			
			// after finding the pivot
			// column operation to make A[p][j] = 0, where k < j

#ifdef WANDEBUG

			std::cout << "After found the pivot:\n";
			A. write (std::cout);
#endif

			Element tmp;
			
			r. assign (g, *(cur_c -> begin()));

			for (tmp_c = cur_c + 1; tmp_c != A.colEnd(); ++ tmp_c) {
							       
				// test if needing to update 
				if (!r. isZero (*(tmp_c -> begin()))) {
					
					r.div (tmp, *(tmp_c -> begin()), g);
					
					r.negin(tmp);
					
					vd. axpyin (*tmp_c, tmp, *cur_c);

					RemIn (*tmp_c, r, d);
				}
			}

#ifdef WANDEBUG
			std::cout << "First row should be become (*, 0, ..., 0):\n";
			A.write(std::cout);
#endif
			return A;
		}
		

		
		/** @memo eliminationCol will make the first col (*, 0, ..., 0)
		 *  by elementary row operation.
		 *  It is the implementation of Iliopoulos algorithm
		 */
		template<class Matrix, class Ring>
			static Matrix& eliminationCol (Matrix& A, const Ring& r, 
						       const typename Ring::Element& d) {

#ifdef WANDEBUG
			std::cout << "In eliminationCol, Matrix:\n";
			A. write (std::cout);
#endif
			if((A.rowdim() <= 1) || (A.rowdim() == 0)) return A;
			
			//typedef typename Matrix::Field Field;
			typedef typename Ring::Element Element;
			
			//Field r (A.field());

			VectorDomain<Ring> vd(r);
			
			typename Matrix::ColIterator cur_c, tmp_c;
			typename Matrix::RowIterator cur_r, tmp_r;
			
			//typename Matrix::Row::iterator row_p1, row_p2;
			typename Matrix::Col::iterator col_p1, col_p2;
						
			
			cur_c = A.colBegin();
			cur_r = A.rowBegin();

			col_p1 = cur_c -> begin();

			Element g, s, t;

			r. xgcd (g, s, t, *col_p1, d);

			// If A[0][0] is coprime to d.
			if (r.isUnit (g) ) {
#ifdef WANDEBUG
				std::cout << "A[0][0] unit:\n";
#endif
				if (! r. isOne (s)) {
				
					
					vd. mulin (*cur_r, s);

					RemIn (*cur_r, r, d);

				}
				
			}
			else {
				// Make A[0][0] = 0;
				if (!r.isZero(*col_p1)) {
					
					Element g, s, t, y1, y2;
					
					std::vector<Element> tmp1(A.coldim()), tmp2(A.coldim());
					
					col_p2 = col_p1 + 1;
					
					r.xgcd(g, s, t, *col_p1, *col_p2);
					
					r.div (y1, *col_p2, g);
					
					r.div (y2, *col_p1, g);
					
					r.negin (y1);
#ifdef WANDEBUG
					std::cout <<"col_p2, y1 " << *col_p2 << " " << y1 <<"\n";
					
					std::cout <<"col_p1, y2 " << *col_p1 << " " << y2 <<"\n";
					
					std::cout <<" s row_p1 + t row_p2 should be " << g << " :" 
						  << (s * (*col_p1) + t * (*col_p2)) <<"\n";
					
					std::cout << "g, s, t, y1, y2, s y2 - t y1(1): " <<
						g << " " << s << " " << t << " " <<
						y1 << " " << y2 << " " 
						  << s * y2 - t * y1 << "\n";

#endif

					tmp_r = cur_r + 1;
					
					vd. mul (tmp1, *cur_r, y1);
					
					vd. axpyin (tmp1, y2, *tmp_r);

					vd. mul (tmp2, *cur_r, s);

					vd. axpyin (tmp2, t, *tmp_r);

					vd. copy (*cur_r, tmp1);

					vd. copy (*tmp_r, tmp2);

					RemIn (*cur_r, r, d);

					RemIn (*tmp_r, r, d);
				}

#ifdef WANDEBUG
				std::cout << "A[0][0] should be 0:\n";
				A. write (std::cout);
#endif				
				// matrix index is 0-based
				std::vector<Element> tmp_v (A.rowdim());
					
				typename std::vector<Element>::iterator p1, p2;
				
				Element g, s;
				
				col_p2 = col_p1 + 1;
				
				r.assign (g, *col_p2); ++ col_p2;
				
				r. init (tmp_v[0], 1);
				
				p1 = tmp_v.begin() + 1;
				
				r.init(*p1,1); ++ p1;

				for(; col_p2 != cur_c -> end(); ++ col_p2, ++ p1) {
					
					r.xgcd (g, s, *p1, g, *col_p2);
					
					if (! r.isOne(s)) 
						for (p2 = tmp_v.begin() + 1; p2 != p1; ++ p2) {
							
							r. mulin (*p2, s);

							r. remin (*p2, d);
						}
					
				}
				
				if (r.isZero(g))  return A;

				// no pivot found
#ifdef WANDEBUG
				std::cout << "g = " << g <<"\n";
				std::cout << "tmp_v:\n";
				
				for (p1 = tmp_v.begin(); p1 != tmp_v.end(); ++p1)
					std::cout << *p1 << " ";
				
				std::cout << '\n';
#endif


				
				Element tmp;
				
				for (tmp_c = cur_c; tmp_c != A.colEnd(); ++ tmp_c) {
					
					vd. dot(tmp, *tmp_c, tmp_v);
					
					r. rem (*(tmp_c -> begin()), tmp, d);

				}
			}			

#ifdef WANDEBUG

			std::cout << "After found the pivot:\n";
			A. write (std::cout);
#endif

			// A pivot is found

			Element tmp;

			r. assign (g, *( cur_r -> begin()));
			
			for (tmp_r = cur_r + 1; tmp_r != A.rowEnd(); ++ tmp_r) {
					
				if (! r.isZero(*(tmp_r -> begin() ) ) ) {
					
					r.div (tmp, *(tmp_r -> begin()), g);
					
					r.negin (tmp);
					
					vd. axpyin (*tmp_r, tmp, *cur_r);

					RemIn (*tmp_r, r, d);
				}
			}

#ifdef WANDEBUG
			std::cout << "First Col should become (*, 0, ... ,0):\n";
			A.write(std::cout);
#endif

			return A;

		}

		template<class Matrix, class Ring>				
			static bool check(const Matrix& A, const Ring& r, 
					  const typename Ring::Element& d) {
			
			//typedef typename Matrix::Ring Field;
			typedef typename Ring::Element Element;
			
			typename Matrix::ConstRowIterator cur_r;
			typename Matrix::ConstRow::const_iterator row_p;
			
			Element tmp, rem, g;
			
			cur_r = A.rowBegin();
			row_p = cur_r -> begin();
				
			tmp = * (A.rowBegin() -> begin());
			
			if (r.isZero(tmp)) return true;

			r. gcd (g, tmp, d);

			for (++ row_p; row_p != cur_r -> end(); ++ row_p ) {
				
				if (!r. isDivisor (g, *row_p))

					return false;
			}
			
			return true;
		}
				
		/** @memo Diagonalize the matrix A.
		 */
		template<class Matrix, class Ring>
			static Matrix& diagonalizationIn(Matrix& A, const Ring& r,
							 const typename Ring::Element& d) { 
			
			if (A.rowdim() == 0 || A.coldim() == 0) return A;

			
			do {
				
				eliminationRow (A, r, d);

				eliminationCol (A, r, d);
			}

			while (!check(A, r, d));

#ifdef WANDEBUG
			std::cout << "Matrix should be diagonal block:\n";
			A.write(std::cout);
#endif
			
			typename SubMatrixTraits<Matrix>::value_type 
				sub(A, (unsigned int)1, (unsigned int)1, 
				    A.rowdim() - 1, A.coldim() - 1);

			diagonalizationIn(sub, r, d);

			return A;
		}


	public:

		template<class Matrix>
			static  Matrix& smithIn(Matrix& A,
						const typename Matrix::Field::Element& d) {
			
			typedef typename Matrix::Field Ring;
			typedef typename Ring::Element Element;
			
			Ring r (A.field());

			typename Matrix::RowIterator row_p;
			
			for (row_p = A. rowBegin(); row_p != A. rowEnd();
			     ++ row_p) 
				
				RemIn (*row_p, r, d);
			
			Element tmp, zero, one;

			r. init (zero, 0);

			r. init (one, 1);
			
			diagonalizationIn(A, r, d);

			int min = A.rowdim() <= A.coldim() ? A.rowdim() : A.coldim();			

#ifdef WANDEBUG
			std::cout << "After Diagonalization:\n";
			A. write (std::cout);
#endif

			int i, j;
			
			Element g;
			
			for (i = 0; i < min; ++ i) {
				
				for ( j = i + 1; j < min; ++ j) {
					
					if (r. isUnit(A[i][i]))  break;
						
					else if (r. isZero (A[j][j])) continue;
					
					else if (r. isZero (A[i][i])) {
						std::swap (A[i][i], A[j][j]);
					}
					
					else {
						r. gcd (g, A[j][j], A[i][i]);
						
						r. divin (A[j][j], g);
						
						r. mulin (A[j][j], A[i][i]);

						r. assign (A[i][i], g);
					}
				}
				r. gcd (A[i][i], A[i][i], d);
			}
			
			return A;
			
		}
		
	private:

		/** @memo RemIn(v, r, d)
		 *  r. remIn (v, d);
		 */
		template<class Vector, class Ring>
			inline static Vector& RemIn (Vector& v, const Ring& r,
						     const typename Ring::Element& d) {
			
			typename Vector::iterator p;

			for (p = v. begin(); p != v. end(); ++ p)

				r. remin (*p, d);


			return v;
		}
	};
	
	
}
				
#endif	
