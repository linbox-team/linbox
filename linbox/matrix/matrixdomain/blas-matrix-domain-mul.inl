/* linbox/matrix/blas-matrix-domain.inl
 * Copyright (C) 2004 Pascal Giorgi, Clément Pernet
 *               2020 Pascal Giorgi
 * Written by :
 *               Pascal Giorgi  pascal.giorgi@lirmm.fr
 *               Clément Pernet clement.pernet@imag.fr
 *               Brice Boyer (briceboyer) <boyer.brice@gmail.com>
 *
 * originally linbox/algorithms/blas-domain.inl
 *
 * ========LICENCE========
 * This file is part of the library LinBox.
 *
 * LinBox is free software: you can redistribute it and/or modify
 * it under the terms of the  GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 */

#ifndef __LINBOX_matrix_matrixdomain_blas_matrix_domain_mul_INL
#define __LINBOX_matrix_matrixdomain_blas_matrix_domain_mul_INL

#include "linbox/linbox-config.h"
#include "fflas-ffpack/ffpack/ffpack.h"
#include "fflas-ffpack/fflas/fflas.h"

namespace LinBox {
  
    
 

    
    /**********************************************
     **********************************************
     *****       MulAdd                       *****
     **********************************************
     ****         D= beta.C + alpha.A*B      *****
     **********************************************
     *********************************************/


    

    // Specialization for matrix container (i.e. BlasMatrix,  BlasSubMAtrix and BlasTransposedMatrix)
    template <typename Matrix1, typename Matrix2, typename Matrix3,typename Matrix4>
    struct BlasMatrixDomainMulAdd_specialized<Matrix1, Matrix2, Matrix3, Matrix4, ContainerCategories::Matrix, ContainerCategories::Matrix>
    {
        typedef typename Matrix1::Field::Element Element;
        
        static Matrix1&  muladdspe( Matrix1 &D,
                                    const Element &beta,   const Matrix2 &C,
                                    const Element & alpha, const Matrix3 &A, const Matrix4 &B) 
		{
			linbox_check( D.rowdim() == C.rowdim());
			linbox_check( D.coldim() == C.coldim());            
            D.copy(C);
            return muladdspe(beta, D, alpha, A, B);
		}
        
        static Matrix1&  muladdspe(const Element &beta,   Matrix1 &C,
                                   const Element & alpha, const Matrix3 &A, const Matrix4 &B) 
        {
            linbox_check( A.coldim() == B.rowdim()); linbox_check( C.rowdim() == A.rowdim());
			linbox_check( C.coldim() == B.coldim());            

            FFLAS::fgemm( C.field(), isTransposed<Matrix3>::value, isTransposed<Matrix4>::value,
                          C.rowdim(), C.coldim(), A.coldim(),
                          alpha, A.getPointer(), A.getStride() , B.getPointer(), B.getStride(),
                          beta,  C.getPointer(), C.getStride());
			return C;
        }
    };

    // Specialization for vector container (i.e. BlasVector,  BlasSubvector)
    template <typename Vector1,typename Vector2,typename Vector3,typename Matrix>
    struct BlasMatrixDomainMulAdd_specialized<Vector1,Vector2,Vector3,Matrix, ContainerCategories::Vector, ContainerCategories::Matrix>
    {
        typedef typename Matrix::Field::Element Element;        
        static Vector1&  muladdspe( Vector1 &d,
                                   const Element &beta,   const Vector2 &c,
                                   const Element & alpha, const Vector3 &a, const Matrix &B) 
		{
			linbox_check( d.size()   == c.size());
			d.copy(c);            
			return muladdspe(beta,d,alpha,a,B);
		}
        static Vector1& muladdspe(const Element &beta,   Vector1 &c,
                                  const Element & alpha, const Vector3 &a, const Matrix &B) 
        {
           	linbox_check( B.rowdim() == a.size());
			linbox_check( B.coldim() == c.size());
			FFLAS::fgemv( B.field(), isTransposed<TransposedBlasMatrix<Matrix>>::value,
                          B.rowdim(), B.coldim(),
                          alpha, B.getPointer(), B.getStride(),a.getPointer(), a.getInc(),
                          beta,  c.getPointer(), c.getInc());
			return c;
        }        
    };
    template <typename Vector1,typename Vector2,typename Matrix,typename Vector3>
    struct BlasMatrixDomainMulAdd_specialized<Vector1,Vector2,Matrix,Vector3,  ContainerCategories::Matrix, ContainerCategories::Vector>
    {
        typedef typename Matrix::Field::Element Element;        
        static Vector1&  muladdspe( Vector1 &d,
                                   const Element &beta,   const Vector2 &c,
                                   const Element & alpha, const Matrix &A, const Vector2 &b) 
		{
			linbox_check( d.size()   == c.size());
			d.copy(c);            
			return muladdspe(beta,d,alpha,A,b);
		}
        
        static Vector1&  muladdspe (const Element &beta,   Vector1 &c,
                                    const Element & alpha, const Matrix &A, const Vector3 &b) 
        {
          	linbox_check( A.coldim() == b.size());
			linbox_check( A.rowdim() == c.size()); 
			FFLAS::fgemv( A.field(), isTransposed<Matrix>::value,
                          A.rowdim(), A.coldim(),
                          alpha, A.getPointer(), A.getStride(),b.getPointer(), b.getInc(),
                          beta,  c.getPointer(), c.getInc());
			return c;
        }        
    };

    


    /**********************
     **********************
     *****      Mul   *****
     **********************
     **********************/
    
	/* does not fall back to MulAdd with BlasPermutation */
    
	// Matrix permutation product C = A*B (it should encompass both BlasMatrix and BlasVector)
	template<class Matrix1, class Matrix2>
	class BlasMatrixDomainMul<Matrix1, Matrix2, BlasPermutation<size_t> > {
	public:
		Matrix1 & operator()(Matrix1& C, const Matrix2& A, const BlasPermutation<size_t>& B) const
		{
			C.copy(A);
			return BlasMatrixDomainMulin<Matrix1, BlasPermutation<size_t> >()(C, B);
		}
	};
    template<class Matrix1, class Matrix2>
	class BlasMatrixDomainMul<Matrix1, BlasPermutation<size_t>, Matrix2> {
	public:
        Matrix1 & operator()(Matrix1& C, const BlasPermutation<size_t>& A, const Matrix2& B) const
		{
			C.copy(B);
			return BlasMatrixDomainMulin<Matrix1, BlasPermutation<size_t> >()(A, C);
		}
	};    
    template<class Matrix1, class Matrix2>
	class BlasMatrixDomainMul<Matrix1, Matrix2, TransposedBlasMatrix<BlasPermutation<size_t>> > {
	public:
		Matrix1 & operator()(Matrix1& C, const Matrix2& A, const TransposedBlasMatrix<BlasPermutation<size_t>>& B) const
		{
			C.copy(A);
			return BlasMatrixDomainMulin<Matrix1, TransposedBlasMatrix<BlasPermutation<size_t>> >()(C, B);
		}
	};
    template<class Matrix1, class Matrix2>
	class BlasMatrixDomainMul<Matrix1, TransposedBlasMatrix<BlasPermutation<size_t>>, Matrix2> {
	public:
        Matrix1 & operator()(Matrix1& C, const TransposedBlasMatrix<BlasPermutation<size_t>>& A, const Matrix2& B) const
		{
			C.copy(B);
			return BlasMatrixDomainMulin<Matrix1, TransposedBlasMatrix<BlasPermutation<size_t>> >()(A, C);
		}
	};

    

    /* does not fall back to MulAdd with TriangularBlasMatrix */
	template<class Matrix1,class Matrix2, class Matrix3>
	class BlasMatrixDomainMul<Matrix1, Matrix2, TriangularBlasMatrix<Matrix3> > {
	public:
		Matrix1& operator()(Matrix1& C, const Matrix2 & A, const TriangularBlasMatrix<Matrix3>& B) const
		{
			C.copy(A);
			return BlasMatrixDomainMulin<Matrix1, TriangularBlasMatrix<Matrix3> >()(C, B);
		}
	};
    template<class Matrix1,class Matrix2, class Matrix3>
	class BlasMatrixDomainMul<Matrix1, TriangularBlasMatrix<Matrix3>, Matrix2 > {
	public:
		Matrix1& operator()(Matrix1& C, const TriangularBlasMatrix<Matrix3>& B, const Matrix2& A) const
		{
          
			C.copy(A);
			return BlasMatrixDomainMulin<Matrix1, TriangularBlasMatrix<Matrix3> >()(B, C);
		}
	};

    /* does not fall back to MulAdd with TriangularBlasMatrix and BlasPermutation*/
	template<class Matrix1,class Matrix2>
	class BlasMatrixDomainMul<Matrix1, TriangularBlasMatrix<Matrix2>, BlasPermutation<size_t> > {
	public:
        Matrix1& operator()(Matrix1& C, const TriangularBlasMatrix<Matrix2>& A, const BlasPermutation<size_t>& B) const
		{
			C.copy(A);
			return BlasMatrixDomainMulin<Matrix1, BlasPermutation<size_t> >()(C, B);
		}
	};
	template<class Matrix1,class Matrix2>
	class BlasMatrixDomainMul<Matrix1, BlasPermutation<size_t>, TriangularBlasMatrix<Matrix2> > {
	public:
		Matrix1& operator()(Matrix1& C, const BlasPermutation<size_t>& B, const TriangularBlasMatrix<Matrix2>& A) const
		{
			C.copy(A);
			return BlasMatrixDomainMulin<Matrix1, BlasPermutation<size_t> >()(B, C);
		}
	};

    // special case for TriangularBlasMatrix by TriangularBlasMatrix multiplication (revert to trtrm case)
    template<class Matrix1, class Matrix2, class Matrix3>
	class BlasMatrixDomainMul<Matrix1,  TriangularBlasMatrix<Matrix2> , TriangularBlasMatrix<Matrix3> > {
	public:
        Matrix1 & operator()(Matrix1& C, const  TriangularBlasMatrix<Matrix2> & A, const TriangularBlasMatrix<Matrix3>& B) const
        {
            typename TriangularBlasMatrix<Matrix3>::constSubMatrixType Bplain(B);                        
            return BlasMatrixDomainMul<Matrix1, TriangularBlasMatrix<Matrix1> , typename TriangularBlasMatrix<Matrix3>::constSubMatrixType >()(C,A,Bplain);
		}
	};
    
    /**********************
     **********************
     *****    MulIn   *****
     **********************
     **********************/
    
    
    // In-place permutation of a matrix
    
	template<class Matrix>
	class BlasMatrixDomainMulin<Matrix, BlasPermutation<size_t> > {
	public:
		Matrix& operator()(Matrix& A, const BlasPermutation<size_t>& B) const
		{
			if (B.isIdentity()) return A ;
            if (isTransposed<Matrix>::value == FFLAS::FflasNoTrans){
                linbox_check(A.coldim() >= B.getSize() );
                FFPACK::applyP( A.field(), FFLAS::FflasRight, FFLAS::FflasNoTrans,
                                A.rowdim(), 0,(int) B.getOrder(),
                                A.getPointer(), A.getStride(), B.getPointer() );
            }
            else {
                linbox_check(A.rowdim() >= B.getSize() );
                FFPACK::applyP( A.field(), FFLAS::FflasLeft, FFLAS::FflasNoTrans,
                                A.coldim(), 0,(int) B.getOrder(),
                                A.getPointer(), A.getStride(), B.getPointer() );

            }
			return A;
		}
        Matrix& operator()(const BlasPermutation<size_t>& B, Matrix& A) const
		{
			if (B.isIdentity()) return A ;
            if (isTransposed<Matrix>::value == FFLAS::FflasNoTrans){
                linbox_check( A.rowdim() >= B.getSize() );
                FFPACK::applyP (A.field(), FFLAS::FflasLeft, FFLAS::FflasNoTrans,
                                A.coldim(), 0,(int) B.getOrder(), A.getPointer(), A.getStride(), B.getPointer() );
            }
            else {
                linbox_check( A.coldim() >= B.getSize() );
                FFPACK::applyP (A.field(), FFLAS::FflasRight, FFLAS::FflasNoTrans,
                                A.rowdim(), 0,(int) B.getOrder(), A.getPointer(), A.getStride(), B.getPointer() );
            
            }
			return A;
		}
	};
	template<class Matrix>
	class BlasMatrixDomainMulin<Matrix, TransposedBlasMatrix<BlasPermutation<size_t>>> {
	public:
		Matrix& operator()(Matrix& A, const TransposedBlasMatrix<BlasPermutation<size_t> >& B) const
		{
            BlasPermutation<size_t>& BB=B.getMatrix();
            if (BB.isIdentity()) return A ;
            if (isTransposed<Matrix>::value == FFLAS::FflasNoTrans){
                linbox_check(A.coldim() >= BB.getSize() );
                FFPACK::applyP( A.field(), FFLAS::FflasRight, FFLAS::FflasTrans,
                                A.rowdim(), 0,(int) BB.getOrder(),
                                A.getPointer(), A.getStride(), BB.getPointer() );
            }
            else {
                linbox_check(A.rowdim() >= BB.getSize() );
                FFPACK::applyP( A.field(), FFLAS::FflasLeft, FFLAS::FflasTrans,
                                A.coldim(), 0,(int) BB.getOrder(),
                                A.getPointer(), A.getStride(), BB.getPointer() );

            }
			return A;
		}
        Matrix& operator()(const TransposedBlasMatrix<BlasPermutation<size_t> >& B, Matrix& A) const
		{
            BlasPermutation<size_t>& BB=B.getMatrix();
            if (BB.isIdentity()) return A ;
            if (isTransposed<Matrix>::value == FFLAS::FflasNoTrans){
                linbox_check( A.rowdim() >= BB.getSize() );
                FFPACK::applyP (A.field(), FFLAS::FflasLeft, FFLAS::FflasTrans,
                                A.coldim(), 0,(int) BB.getOrder(), A.getPointer(), A.getStride(), BB.getPointer() );
            }
            else {
                linbox_check( A.coldim() >= BB.getSize() );
                FFPACK::applyP (A.field(), FFLAS::FflasRight, FFLAS::FflasTrans,
                                A.rowdim(), 0,(int) BB.getOrder(), A.getPointer(), A.getStride(), BB.getPointer() );
            
            }
			return A;
		}
	};

    // In-place permutation of a vector

	template<class Field>
	class BlasMatrixDomainMulin<BlasVector<Field>, BlasPermutation<size_t> > {
	public:
		BlasVector<Field>& operator()(BlasVector<Field>& A, const BlasPermutation<size_t>& B) const
		{           
			if (B.isIdentity()) return A ;
			linbox_check( A.size() >= B.getSize() );
			FFPACK::applyP( A.field(), FFLAS::FflasRight, FFLAS::FflasNoTrans,
                            1, 0,(int) B.getOrder(), A.getPointer(), 1, B.getPointer() );
			return A;
		}
		BlasVector<Field>& operator()(const BlasPermutation<size_t>& B, BlasVector<Field>& A)
		{
			if (B.isIdentity()) return A ;
			linbox_check( A.size() >= B.getSize() );
			FFPACK::applyP( A.field(), FFLAS::FflasLeft, FFLAS::FflasNoTrans,
                            1, 0,(int) B.getOrder(), A.getPointer(), 1, B.getPointer() );
			return A;
		}
	};

    template<class Field>
	class BlasMatrixDomainMulin<BlasVector<Field>, TransposedBlasMatrix<BlasPermutation<size_t>>> {
	public:
		BlasVector<Field>& operator()(BlasVector<Field>& A, const TransposedBlasMatrix<BlasPermutation<size_t>>& B) const
		{
             BlasPermutation<size_t>& BB=B.getMatrix();
			if (BB.isIdentity()) return A ;
			linbox_check( A.size() >= BB.getSize() );
			FFPACK::applyP( A.field(), FFLAS::FflasRight, FFLAS::FflasTrans,
                            1, 0,(int) BB.getOrder(), A.getPointer(), 1, BB.getPointer() );
			return A;
		}
		BlasVector<Field>& operator()(const TransposedBlasMatrix<BlasPermutation<size_t>>& B, BlasVector<Field>& A)
		{
            BlasPermutation<size_t>& BB=B.getMatrix();
			if (BB.isIdentity()) return A ;
			linbox_check( A.size() >= BB.getSize() );
			FFPACK::applyP( A.field(), FFLAS::FflasLeft, FFLAS::FflasTrans,
                            1, 0,(int) BB.getOrder(), A.getPointer(), 1, BB.getPointer() );
			return A;
		}
	};

	// In-place matrix*triangular matrix product
	template<class Matrix1, class Matrix2>
	class BlasMatrixDomainMulin<Matrix1, TriangularBlasMatrix<Matrix2>>{
    public:
        Matrix1& operator()(Matrix1& A, const TriangularBlasMatrix<Matrix2>& B) const
        {
            linbox_check( A.coldim() == B.rowdim() );

            FFLAS::ftrmm( A.field(), FFLAS::FflasRight, (FFLAS::FFLAS_UPLO) (B.getUpLo()),
                          FFLAS::FflasNoTrans,(FFLAS::FFLAS_DIAG) (B.getDiag()),
                          A.rowdim(), A.coldim(), A.field().one,
                          B.getPointer(), B.getStride(), A.getPointer(), A.getStride() );
            return A;
        }

        Matrix1& operator()(const TriangularBlasMatrix<Matrix2>& B, Matrix1& A) const
        {
            linbox_check( B.coldim() == A.rowdim() );
            FFLAS::ftrmm( A.field(), FFLAS::FflasLeft, (FFLAS::FFLAS_UPLO)(B.getUpLo()),
                          FFLAS::FflasNoTrans, (FFLAS::FFLAS_DIAG) (B.getDiag()),
                          A.rowdim(), A.coldim(), A.field().one,
                          B.getPointer(), B.getStride(),
                          A.getPointer(), A.getStride() );
            return A;
        }
    };


	/*! @internal In-place matrix*triangular matrix product with transpose.
     */
    template<class Matrix1, class Matrix2>
	class BlasMatrixDomainMulin<Matrix1, TransposedBlasMatrix<TriangularBlasMatrix<Matrix2>>>{
    public:
        Matrix1& operator()(Matrix1& A, const TransposedBlasMatrix<TriangularBlasMatrix<Matrix2>>& B) const
        {
            linbox_check( B.getMatrix().coldim() == A.coldim() );

            FFLAS::ftrmm( A.field(), FFLAS::FflasRight,
                          (FFLAS::FFLAS_UPLO)(B.getMatrix().getUpLo()),
                          FFLAS::FflasTrans,
                          (FFLAS::FFLAS_DIAG) (B.getMatrix().getDiag()),
                          A.rowdim(), A.coldim(),
                          A.field().one,
                          B.getMatrix().getPointer(), B.getMatrix().getStride(),
                          A.getPointer(), A.getStride() );
            return A;
        }

        Matrix1& operator()( const TransposedBlasMatrix<TriangularBlasMatrix<Matrix2>>& B, Matrix1& A) const
        {
            linbox_check( B.getMatrix().coldim() == A.rowdim() );
            FFLAS::ftrmm( A.field(), FFLAS::FflasLeft,
                          (FFLAS::FFLAS_UPLO) (B.getMatrix().getUpLo()),
                          FFLAS::FflasTrans,
                          (FFLAS::FFLAS_DIAG) (B.getMatrix().getDiag()),
                          A.rowdim(), A.coldim(), A.field().one,
                          B.getMatrix().getPointer(), B.getMatrix().getStride(),
                          A.getPointer(), A.getStride() );
            return A;
        }
    };



	/*
	 * specialization for Operand1 of type TriangularBlasMatrix<Field,_Rep> and Operand2 of type BlasPermutation
	 */

	// Matrix permutation product C = A*B


  



    


} // end of namespace LinBox

#endif //__LINBOX_matrix_matrixdomain_blas_matrix_domain_mul_INL


// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
