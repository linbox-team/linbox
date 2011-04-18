/* linbox-sage.C
 * Copyright (C) 2007 Martin Albrecht
 *               2008 Clement Pernet
 *
 * Written by Martin Albrecht
 *            Clement Pernet <clement.pernet@gmail.com>
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#include <iostream>
#include <gmp.h>

#include <cstdlib>
#include <vector>
#include <list>

#include "linbox-sage.h"

#include <linbox/util/commentator.h>

#include <linbox/blackbox/sparse.h>

//#include "linbox/element/givaro-polynomial.h"

#include <linbox/matrix/blas-matrix.h>
#include <linbox/matrix/sparse.h>
#include <linbox/vector/sparse.h>

#include <linbox/algorithms/blas-domain.h>
#include <linbox/algorithms/echelon-form.h>
#include "linbox/algorithms/gauss.h"
#include "linbox/algorithms/smith-form-adaptive.h"
#include "linbox/ffpack/ffpack.h"
#include <linbox/solutions/rank.h>          
#include <linbox/solutions/det.h>
#include <linbox/solutions/solve.h>
#include "linbox/solutions/methods.h"
#include <linbox/solutions/minpoly.h>
#include <linbox/solutions/charpoly.h>
#include "linbox/algorithms/double-det.h"
#include <linbox/integer.h>
#include <linbox/field/gmp-rational.h>
#include <linbox/ring/givaro-polynomial.h>
#include <linbox/field/modular.h>

using namespace LinBox;

/*************************************************************************
   dense modulo Z/nZ
*************************************************************************/

//we are using Modular<double> here as it seems to be best supported

/* NOTE: There are many echelon form functions, possible base rings, etc.  Strangely,
   most don't build.  This combination below does though.
*/

void linbox_modn_dense_delete_array_double(double * f) {delete[] f;}
void linbox_modn_dense_delete_array_float(float * f) {delete[] f;}

template <class Element>
unsigned long int linbox_modn_dense_echelonize (Element modulus, Element* matrix,
						size_t nrows, size_t ncols) {
	Modular<Element> F(modulus);

	size_t * P=new size_t[ncols];
	size_t * Q=new size_t[nrows];
	size_t r = FFPACK::ReducedRowEchelonForm(F, nrows, ncols, matrix, ncols, P,Q);

	for (size_t i=0; i<nrows;++i){
		for (size_t j=0; j<r; ++j){
			*(matrix+i*ncols+j) = 0;
		}
		if (i<r)
			*(matrix + i*(ncols+1)) = 1;
	}
	FFPACK::applyP (F, FFLAS::FflasRight, FFLAS::FflasNoTrans, nrows, 0,r, matrix, ncols, Q);
	delete[] P;
	delete[] Q;
	return r;
}

EXTERN unsigned long int linbox_modn_dense_echelonize_double (double modulus, double*matrix, size_t nrows, size_t ncols)
{return linbox_modn_dense_echelonize<double>(modulus, matrix, nrows, ncols);}

EXTERN unsigned long int linbox_modn_dense_echelonize_float (float modulus, float * matrix, size_t nrows, size_t ncols)
{return linbox_modn_dense_echelonize<float>(modulus,matrix,nrows,ncols);}

template<class Element>
unsigned long int linbox_modn_dense_rank (Element modulus,  Element* matrix, size_t nrows, size_t ncols) {

	Modular<Element> F (modulus);
	Element * Ad = new Element [nrows*ncols];
	for (size_t i=0; i< nrows; ++i)
		for (size_t j = 0; j < ncols; ++j)
			*(Ad + i * ncols + j) = *(matrix + i*ncols + j);
	
	size_t r = FFPACK::Rank (F, nrows, ncols, Ad, ncols);
	delete[] Ad;
	return r;
}

EXTERN unsigned long int linbox_modn_dense_rank_double (double modulus, double* matrix, size_t nrows, size_t ncols)
{return linbox_modn_dense_rank<double>(modulus,matrix,nrows,ncols);}
EXTERN unsigned long int linbox_modn_dense_rank_float (float modulus, float* matrix, size_t nrows, size_t ncols)
{return linbox_modn_dense_rank<float>(modulus,matrix,nrows,ncols);}

template<class Element>
Element linbox_modn_dense_det(Element modulus,Element*matrix,size_t nrows,size_t ncols){
	
	Modular<Element> F (modulus);
	Element * Ad = new Element [nrows*ncols];
	for (size_t i=0; i< nrows; ++i)
		for (size_t j = 0; j < ncols; ++j)
			*(Ad + i * ncols + j) = *(matrix+i*ncols+j);
	Element dd = FFPACK::Det (F, nrows, ncols, Ad, ncols);
	delete[] Ad;
	return  dd;
}

EXTERN double linbox_modn_dense_det_double (double modulus, double* matrix, size_t nrows, size_t ncols)
{return linbox_modn_dense_det<double>(modulus,matrix,nrows,ncols);}
EXTERN float linbox_modn_dense_det_float (float modulus, float* matrix, size_t nrows, size_t ncols)
{return linbox_modn_dense_det<float>(modulus,matrix,nrows,ncols);}

template<class Element>
Element* linbox_modn_dense_minpoly (Element modulus, Element ** mp, size_t* degree, 
				    size_t n, Element *matrix) {

	Modular<Element> F(modulus);
	    // Warning: super sketchy memory alloc here!!!!
	std::vector<Element> *minP=new std::vector<Element>(n);
	Element * X = new Element[n*(n+1)];
	size_t * P = new size_t[n];

	// FIXME: check the memory management: better to allocate mp in sage
	FFPACK::MinPoly (F, *minP, n, matrix, n, X, n, P);
	*degree=minP->size()-1;
	
	*mp = &(*minP)[0];
	delete[] P;
	delete[] X;
	return *mp;
}

EXTERN double* linbox_modn_dense_minpoly_double (double modulus, double ** mp, size_t* degree, size_t n, double*matrix)
{return linbox_modn_dense_minpoly<double> (modulus, mp, degree, n, matrix);}
EXTERN float* linbox_modn_dense_minpoly_float (float modulus, float ** mp, size_t* degree, size_t n, float *matrix)
{return linbox_modn_dense_minpoly<float> (modulus, mp, degree, n, matrix);}


template<class Polynomial, class Field>
Polynomial & mulpoly(const Field& F, Polynomial &res, const Polynomial & P1, const Polynomial & P2){
			size_t i,j;
			    // Warning: assumes that res is allocated to the size of the product
			    //res.resize(P1.size()+P2.size()-1);
			for (i=0;i<res.size();i++)
				F.assign(res[i], 0.0);
			for ( i=0;i<P1.size();i++)
				for ( j=0;j<P2.size();j++)
					F.axpyin(res[i+j],P1[i],P2[j]);
			return res;
			
		}

template<class Element>
Element* linbox_modn_dense_charpoly (Element modulus, Element *& cp, size_t n, Element *matrix) {

	Modular<Element> F(modulus);
	
	// FIXME: check the memory management: better to allocate mp in sage
	std::list<std::vector<Element> > P_list;

	
	FFPACK::CharPoly (F, P_list, n, matrix, n);

	std::vector<Element>* tmp = new std::vector<Element> (n+1);
	std::vector<Element> P;


	typename std::list<std::vector<Element> >::const_iterator it;
	it = P_list.begin();
	P = *(it++);
	while( it!=P_list.end() ){
		mulpoly (F,*tmp, P, *it);
		P = *tmp;
		//	delete &(*it);
		++it;
	}
	return cp = &(*tmp)[0];
	
}

EXTERN double* linbox_modn_dense_charpoly_double (double modulus, double ** cp, size_t n, double * matrix)
{return linbox_modn_dense_charpoly<double>(modulus,*cp,n,matrix);}
EXTERN float* linbox_modn_dense_charpoly_float (float modulus, float ** cp, size_t n, float * matrix)
{return linbox_modn_dense_charpoly<float>(modulus,*cp,n,matrix);}

template<class Element>
unsigned long linbox_modn_dense_col_rankprofile_submatrix_indices (Element modulus,
								   Element* matrix,
								   size_t ** row_idx,
								   size_t ** col_idx,
								   size_t* rank,
								   size_t nrows,
								   size_t ncols){
	Modular<Element> F (modulus);
	Element * Ad = new Element [nrows*ncols];
	for (size_t i=0; i< nrows; ++i)
		for (size_t j = 0; j < ncols; ++j)
			*(Ad + i * ncols + j) = *(matrix +i*ncols + j);
  
	size_t R;
	Element * X;
	FFPACK::ColRankProfileSubmatrixIndices (F, nrows, ncols, Ad, ncols,
						*row_idx, *col_idx, *rank);
	delete[] Ad;
	return *rank;
}
EXTERN unsigned long linbox_modn_dense_col_rankprofile_submatrix_indices_double (double modulus, double* matrix, size_t ** row_idx, size_t ** col_idx, size_t * rank, size_t nrows, size_t ncols)
{return linbox_modn_dense_col_rankprofile_submatrix_indices<double>(modulus,matrix,row_idx,col_idx,rank,nrows,ncols);}
EXTERN unsigned long linbox_modn_dense_col_rankprofile_submatrix_indices_float (float modulus, float* matrix, size_t ** row_idx, size_t ** col_idx, size_t * rank, size_t nrows,size_t ncols)
{return linbox_modn_dense_col_rankprofile_submatrix_indices<float>(modulus,matrix,row_idx,col_idx,rank,nrows,ncols);}

template<class Element>
unsigned long linbox_modn_dense_col_rankprofile_submatrix (Element modulus,
						 Element* matrix,
						 Element* ans,
						 size_t& rank,
						 size_t nrows, size_t ncols){
	Modular<Element> F (modulus);
	//FIXME: check the memory managmenent
	Element * Ad = new Element [nrows*ncols];
	for (size_t i=0; i< nrows; ++i)
		for (size_t j = 0; j < ncols; ++j)
			*(Ad + i * ncols + j) = *(matrix + i*ncols + j);
  
	size_t R = FFPACK::ColRankProfileSubmatrix (F, nrows, ncols, ans, ncols, ans, rank);
	delete[] Ad;
	return R;
}
EXTERN unsigned long linbox_modn_dense_col_rankprofile_submatrix_double (double modulus, double* matrix, double* out, size_t* rank, size_t nrows, size_t ncols)
{return linbox_modn_dense_col_rankprofile_submatrix<double>(modulus,matrix,out,*rank, nrows,ncols);}
EXTERN unsigned long linbox_modn_dense_col_rankprofile_submatrix_float (float modulus, float* matrix, float* out, size_t* rank, size_t nrows, size_t ncols)
{return linbox_modn_dense_col_rankprofile_submatrix<float>(modulus,matrix,out,*rank,nrows,ncols);}

template<class Element>
Element* linbox_modn_dense_matrix_matrix_multiply (Element modulus, Element *ans, 
					       Element *A, Element *B,
					       size_t m, size_t n, size_t k) {

	Modular<Element> F(modulus);

	FFLAS::fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, m,n,k, 1.0, 
		      A, k, B, n, 0.0, ans, n);
}

EXTERN double* linbox_modn_dense_matrix_matrix_multiply_double(double modulus, double * ans, double *A, double *B, size_t m,size_t n,size_t k)
{linbox_modn_dense_matrix_matrix_multiply<double>(modulus,ans,A,B,m,n,k);}
EXTERN float* linbox_modn_dense_matrix_matrix_multiply_float(float modulus, float * ans, float *A, float *B, size_t m,size_t n,size_t k)
{linbox_modn_dense_matrix_matrix_multiply<float>(modulus,ans,A,B,m,n,k);}

template<class Element>
Element *  linbox_modn_dense_matrix_matrix_general_multiply(Element modulus,
						      Element *ans,
						      Element alpha, Element beta,
						      Element *A, Element *B,
						      size_t m, size_t n, size_t k) {
	Modular<Element> F(modulus);
	FFLAS::fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, m,n,k, alpha, 
		      A, k, B, n, beta,ans, n);
	return ans;
	
}
EXTERN double* linbox_modn_dense_matrix_matrix_general_multiply_double (double modulus, double * ans, double alpha, double beta, double *A, double *B, size_t m, size_t n, size_t k)
{return linbox_modn_dense_matrix_matrix_general_multiply<double>(modulus,ans,alpha,beta,A,B,m,n,k);}
EXTERN float* linbox_modn_dense_matrix_matrix_general_multiply_float (float modulus, float * ans, float alpha, float beta, float *A, float *B, size_t m, size_t n, size_t k)
{return linbox_modn_dense_matrix_matrix_general_multiply<float>(modulus,ans,alpha,beta,A,B,m,n,k);}

/*************************************************************************
    dense over ZZ
*************************************************************************/

typedef PID_integer Integers;

template <class Field, class Polynomial>
void printPolynomial (const Field &F, const Polynomial &v) 
{
        for (int i = v.size () - 1; i >= 0; i--) {
                F.write (std::cout, v[i]);
                if (i > 0)
                        std::cout << " x^" << i << " + ";
        }
	std::cout << std::endl;
}

Integers ZZ;
SpyInteger spy;
typedef GivPolynomialRing<Integers,Dense> IntPolRing;

DenseMatrix<Integers> new_matrix(mpz_t** matrix, size_t nrows, size_t ncols) {
  DenseMatrix<Integers> A ( ZZ, nrows, ncols);

     size_t i, j, k;
     for (i=0; i < nrows; i++) {
	 for (j=0; j < ncols; j++) {
	     Integers::Element t;
	     mpz_set(spy.get_mpz(t), matrix[i][j]);
	     A.setEntry(i, j, t);
	 }
     }
     return A;
}

DenseMatrix<Integers> new_matrix_integers(mpz_t** matrix, size_t nrows, size_t ncols) {
     Integers ZZ;
     DenseMatrix<Integers> A ( ZZ,nrows, ncols);

     size_t i, j, k;
     for (i=0; i < nrows; i++) {
	 for (j=0; j < ncols; j++) {
	     Integers::Element t;
	     mpz_set(spy.get_mpz(t), matrix[i][j]);
	     A.setEntry(i, j, t);
	 }
     }
     return A;
}

template<class Element>
void set_matrix(mpz_t** matrix, BlasMatrix<Element>& A, size_t nrows, size_t ncols) {
     size_t i, j, k;
     for (i=0; i < nrows; i++) {
	 for (j=0; j < ncols; j++) {
	     mpz_set(matrix[i][j], spy.get_mpz(A.getEntry(i,j)));
	 }
     }
}

void linbox_integer_dense_minpoly_hacked(mpz_t* *mp, size_t* degree, size_t n, mpz_t** matrix, int do_minpoly) {
     /* We program around a bizarre bug in linbox, where minpoly doesn't work
	on matrices that are n x n with n divisible by 4!
     */
     size_t m;
     if (n % 4 == 0 || !do_minpoly) {
	 m = n + 1;
     } else {
	 m = n;
     }

     DenseMatrix<Integers> A( ZZ, m, m);

     size_t i, j;
     Integers::Element t;
     for (i=0; i < n; i++) {
	 for (j=0; j < n; j++) {
	     mpz_set(spy.get_mpz(t), matrix[i][j]);
	     A.setEntry(i, j, t);
	 }
     }

 //    vector<Integers::Element> m_A;
     IntPolRing::Element m_A;

     if (do_minpoly) 
	 minpoly(m_A, A);
     else
	 charpoly(m_A, A);

     if (n%4 == 0 || !do_minpoly) {
	 /* Program around the bug. 
	    It is OK that this code is crappy and redundant, since it will get replaced
	    when linbox gets fixed. */
	 int divide_by_x;

	 if (!do_minpoly)
	     divide_by_x = 1;
	 else {
	     long unsigned int r;
	     rank(r, A);
	     divide_by_x = (r==n);
	 }
	 if (divide_by_x) {
	     /* x was not a factor of the charpoly after all. */
	     (*mp) = new mpz_t[m_A.size()-1];
	     *degree = m_A.size() - 2;
	     for (size_t i=0; i <= *degree; i++) {
		 mpz_init((*mp)[i]);
		 mpz_set((*mp)[i], spy.get_mpz(m_A[i+1])); 
	     }
	     return;
	 }
     }

     (*mp) = new mpz_t[m_A.size()];
     *degree = m_A.size() - 1;
     for (size_t i=0; i <= *degree; i++) {
	 mpz_init((*mp)[i]);
	 mpz_set((*mp)[i], spy.get_mpz(m_A[i])); 
     }

}

 void linbox_integer_dense_charpoly(mpz_t* *mp, size_t* degree, size_t n, mpz_t** matrix) {
     /* THIS IS Broken when n % 4 == 0!!!!  Use above function instead. */
   /*    linbox_integer_dense_minpoly(mp, degree, n, matrix, 0); */

     DenseMatrix<Integers> A(new_matrix(matrix, n, n));
     IntPolRing::Element m_A;
     charpoly(m_A, A);

     (*mp) = new mpz_t[m_A.size()];
     *degree = m_A.size() - 1;
     for (size_t i=0; i <= *degree; i++) {
	 mpz_init((*mp)[i]);
	 mpz_set((*mp)[i], spy.get_mpz(m_A[i])); 
     }

 }

 void linbox_integer_dense_minpoly(mpz_t* *mp, size_t* degree, size_t n, mpz_t** matrix) {
     /* THIS IS Broken when n % 4 == 0!!!!  Use above function instead. */
   /*    linbox_integer_dense_minpoly(mp, degree, n, matrix, 0); */

   DenseMatrix<Integers> A(new_matrix(matrix, n, n));
     IntPolRing::Element m_A;
     minpoly(m_A, A);

     (*mp) = new mpz_t[m_A.size()];
     *degree = m_A.size() - 1;
     for (size_t i=0; i <= *degree; i++) {
	 mpz_init((*mp)[i]);
	 mpz_set((*mp)[i], spy.get_mpz(m_A[i])); 
     }

 }

void linbox_integer_dense_delete_array(mpz_t* f) {
	delete[] f;
 }

 int linbox_integer_dense_matrix_matrix_multiply(mpz_t** ans, mpz_t **A, mpz_t **B,
				   size_t A_nr, size_t A_nc, size_t B_nr, size_t B_nc) 
 {
   typedef PID_integer Integers;
   Integers ZZ;

     BlasMatrix<Integers::Element> AA(new_matrix_integers(A, A_nr, A_nc));
     BlasMatrix<Integers::Element> BB(new_matrix_integers(B, B_nr, B_nc));
     if (A_nc != B_nr)
	 return -1;   // error
     BlasMatrix<Integers::Element> CC( A_nr, B_nc);

     MatrixDomain<Integers> MD(ZZ);

     MD.mul(CC, AA, BB);

     set_matrix(ans, CC, A_nr, B_nc);

     return 0;
 }

 unsigned long linbox_integer_dense_rank(mpz_t** matrix, size_t nrows,
					 size_t ncols) {
     DenseMatrix<Integers> A(new_matrix(matrix, nrows, ncols));
     unsigned long r;
     rank(r, A);
     return r;
 }

 void linbox_integer_dense_det(mpz_t ans, mpz_t** matrix, size_t nrows,
			       size_t ncols) {
   commentator.setMaxDetailLevel(0);
   commentator.setMaxDepth (0);

   DenseMatrix<Integers> A(new_matrix_integers(matrix, nrows, ncols));
   Integers::Element d;
   det(d, A);
   mpz_set(ans, spy.get_mpz(d));
}


DenseMatrix<NTL_ZZ> new_matrix_integer_dense_ntl(mpz_t** matrix, size_t nrows, size_t ncols) {
     NTL_ZZ ZZ;
     DenseMatrix<NTL_ZZ> A (ZZ,nrows, ncols);
     size_t i, j, k;
     for (i=0; i < nrows; i++) {
	 for (j=0; j < ncols; j++) {
	   NTL_ZZ::Element t;
	   PID_integer::Element s;
	   mpz_set(spy.get_mpz(s), matrix[i][j]);
	   ZZ.init(t, s);
	   A.setEntry(i, j, t);
	 }
     }
     return A;
}

void linbox_integer_dense_double_det (mpz_t  ans1, mpz_t ans2, mpz_t **a, mpz_t ** b, mpz_t **c,
				      size_t n, int proof) {

  PID_integer ZZ;
  BlasBlackbox <PID_integer> A (ZZ, n+1, n);
  size_t i, j, k;
  for (i=0; i < n-1; i++) {
    for (j=0; j < n; j++) {
      PID_integer::Element t;
      mpz_set(spy.get_mpz(t), a[i][j]);
      A.setEntry(i, j, t);
    }
  }
  for (j=0; j < n; j++) {
    PID_integer::Element t;
    mpz_set(spy.get_mpz(t), b[0][j]);
    A.setEntry (n-1, j, t);
    mpz_set (spy.get_mpz(t), c[0][j]);
    A.setEntry (n, j, t);
  }
  PID_integer::Element d1,d2;
  doubleDet (d1, d2, A, proof);
  mpz_set(ans1, spy.get_mpz(d1));
  mpz_set(ans2, spy.get_mpz(d2));
}

/*
This won't build on OS X PPC.

void linbox_integer_dense_smithform(mpz_t **v, 
                                    mpz_t **matrix, size_t nrows, size_t ncols) {
  typedef NTL_ZZ Ints;
  Ints Z;
  DenseMatrix<Ints> M(new_matrix_integer_dense_ntl(matrix, nrows, ncols));
  vector<integer> w(ncols);
  SmithFormAdaptive::smithForm(w, M);

  (*v) = new mpz_t[ncols];
  for (size_t i=0; i < ncols; i++) {
    mpz_init((*v)[i]);
    mpz_set((*v)[i], spy.get_mpz(w[i]));
  }
}
*/

/*************************************************************************
   sparse modulo Z/nZ
*************************************************************************/

struct c_vector_modint_linbox {
  // copy of the declaration in vector_modn_sparse.pxi
  int *entries;
  int p;
  size_t *positions;
  size_t degree;
  size_t num_nonzero;
};
typedef unsigned int mod_int;
typedef Modular<unsigned int> GFp;
typedef GFp::Element  Element;
typedef std::vector <pair <size_t, Element> > SparseSeqVectorGFp;
typedef SparseMatrix<GFp, SparseSeqVectorGFp> SparseMatrixGFp;

static SparseMatrixGFp linbox_new_modn_sparse_matrix(mod_int modulus, size_t numrows, size_t numcols, void *rows) {
  GFp F(modulus);
  SparseMatrixGFp M(F, numrows, numcols);
  
  struct c_vector_modint_linbox *A = static_cast<struct c_vector_modint_linbox *>(rows);
  
  for(int i = 0; i < numrows; i++) {
    for(int j = 0; j < A[i].num_nonzero; j++) {
      M.setEntry(i, A[i].positions[j], A[i].entries[j]);
    } 
  }
  return M;
}

static vector<Element> linbox_new_modn_sparse_vector(mod_int modulus, size_t len, void *_vec) {
  GFp F(modulus);

  vector<GFp::Element> A(len);
  
  if (_vec==NULL) {
    return A;
  }

  struct c_vector_modint_linbox *vec = static_cast<struct c_vector_modint_linbox*>(_vec);
  for(int i = 0; i < vec->num_nonzero; i++) {
    A[vec->positions[i]] = vec->entries[i];
  }
  return A;
}

unsigned long linbox_modn_sparse_matrix_rank(mod_int modulus, 
					     size_t numrows, size_t numcols,
					     void *rows, int gauss) {
  GFp F(modulus);
  unsigned long M_rank;
  Element M_det;
  GaussDomain<GFp> dom(F);

  SparseMatrixGFp M( linbox_new_modn_sparse_matrix(modulus, numrows, numcols, rows) );

  if(!gauss) {
    dom.InPlaceLinearPivoting(M_rank, M_det, M, numrows, numcols);
  } else {
    dom.NoReordering(M_rank, M_det, M, numrows, numcols);
  }

  //*pivots = (int*)calloc(sizeof(int), dom.pivots.size());

//   int j=0;
//   for(vector<int>::const_iterator i= dom.pivots.begin(); i!= dom.pivots.end(); i++, j++){
//     (*pivots)[j] = *i;
//   }

  return M_rank;
}

vector<mod_int> linbox_modn_sparse_matrix_solve(mod_int p, size_t numrows, size_t numcols,
						void *_a, void *b, int method) {
  // solve ax = b, for x, a matrix, b vector, x vector
  GFp F(p);

  vector<mod_int> X( numrows);
  vector<mod_int> B( linbox_new_modn_sparse_vector(p, numcols, b));

  SparseMatrixGFp A(linbox_new_modn_sparse_matrix(p, numrows, numcols, _a));

  switch(method) {
  case 1:
    solve(X, A, B, Method::BlasElimination());
    break;

  case 2:
    solve(X, A, B, Method::Blackbox());
    break;

  case 3:
    solve(X, A, B, Method::Wiedemann());
    break;

  default:
    solve(X, A, B);
  }
  return X;
}
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
