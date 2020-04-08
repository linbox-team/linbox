/* linbox/matrix/blas-matrix-domain.inl
 * Copyright (C) 2004 Pascal Giorgi, Clément Pernet
 *
 * Written by :
 *               Pascal Giorgi  pascal.giorgi@ens-lyon.fr
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



#ifndef __LINBOX_matrix_matrixdomain_blas_matrix_domain_INL
#define __LINBOX_matrix_matrixdomain_blas_matrix_domain_INL

#include "linbox/matrix/matrixdomain/blas-matrix-domain-mul.inl"


namespace LinBox
{ /* Det */

	template<class Matrix>
    typename Matrix::Field::Element
	BlasMatrixDomainDet<Matrix>::operator() (const Matrix& A) const
	{
        if (A.rowdim() != A.coldim())
            return A.field().zero;
        typename Matrix::matrixType Acopy(A);
        return 	BlasMatrixDomainDet<typename Matrix::matrixType>()(Acopy);
    }

	template<class Matrix>
	typename Matrix::Field::Element
	BlasMatrixDomainDet<Matrix>::operator() (Matrix& A) const
	{
        if (A.rowdim() != A.coldim())
            return A.field().zero;
        typename Matrix::Field::Element det; A.field().init(det);
        return FFPACK::Det(A.field(), det, A.coldim(), A.getPointer(), A.getStride());
	}

	template< class Matrix>
    class BlasMatrixDomainDet<TriangularBlasMatrix<Matrix> >{
    public:
        typename Matrix::Field::Element operator() (const TriangularBlasMatrix<Matrix> & A) const
        {
            if (A.rowdim() != A.coldim())
                return A.field().zero;
            if (A.getDiag() == Tag::Diag::Unit)
                return A.field().one;
            typename Matrix::Field::Element d ;
            A.field().init(d);A.field().assign(d,A.field().one);                
            for (size_t i=0;i<A.rowdim();i++)
                A.field().mulin(d,A.getEntry(i,i));
            return d;        
        }
        
        typename Matrix::Field::Element operator() (TriangularBlasMatrix<Matrix> & A) const
        {
            return (*this)(const_cast<TriangularBlasMatrix<Matrix> &>(A));
        }
    };
} // LinBox

namespace LinBox
{ /* Rank */

	template<class Matrix>
    size_t
	BlasMatrixDomainRank<Matrix>::operator() (const  Matrix  &A) const
	{
        typename Matrix::matrixType Acopy(A);
        return 	BlasMatrixDomainRank<typename Matrix::matrixType>()(Acopy);
	}

	template<class Matrix>
    size_t
	BlasMatrixDomainRank<Matrix>::operator() (Matrix        &A) const
	{
        return FFPACK::Rank(A.field(),A.rowdim(), A.coldim(), A.getPointer(), A.getStride());
	}

} // LinBox

namespace LinBox
{ /* Inverse */

	template<class Matrix1, class Matrix2>
	int BlasMatrixDomainInv<Matrix1, Matrix2>::operator() (Matrix1 &Ainv, const Matrix2 &A) const
	{
        typename Matrix2::matrixType Acopy(A);
        return 	BlasMatrixDomainInv<Matrix1,typename Matrix2::matrixType>()(Ainv, Acopy);
	}

	template<class Matrix1, class Matrix2>
	int BlasMatrixDomainInv<Matrix1, Matrix2>::operator() (Matrix1 &Ainv, Matrix2 &A) const
	{
        linbox_check( A.rowdim() == A.coldim());
        linbox_check( A.rowdim() == Ainv.rowdim());
        linbox_check( A.coldim() == Ainv.coldim());
        int nullity;
        FFPACK::Invert2 (A.field(), A.rowdim(), A.getPointer(), A.getStride(), Ainv.getPointer(), Ainv.getStride(),nullity);
        return nullity;
	}
} // LinBox

namespace LinBox
{ /* Add Sub */
	// Add
	template<class Matrix1, class Matrix2, class Matrix3>
	Matrix1&
	BlasMatrixDomainAdd<Matrix1, Matrix2, Matrix3 >::operator()(Matrix1& C, const Matrix2& A, const Matrix3& B) const
	{
		linbox_check( A.rowdim() == B.rowdim());
		linbox_check( C.rowdim() == A.rowdim());
		linbox_check( A.coldim() == B.coldim());
		linbox_check( C.coldim() == A.coldim());
        FFLAS::fadd (C.field(), C.rowdim(), C.coldim(),A.getPointer(), A.getStride(),B.getPointer(), B.getStride(),C.getPointer(), C.getStride());
		return C;
	}
	// AddIn
	template<class Matrix1, class Matrix3>
	Matrix1&
	BlasMatrixDomainAddin<Matrix1, Matrix3 >::operator()(Matrix1& C, const Matrix3& B) const
	{
		linbox_check( C.rowdim() == B.rowdim());
		linbox_check( C.coldim() == B.coldim());        
        FFLAS::faddin (C.field(), C.rowdim(), C.coldim(), B.getPointer(), B.getStride(), C.getPointer(), C.getStride());	
        return C;
	}
	// Sub
	template<class Matrix1, class Matrix2, class Matrix3>
	Matrix1&
	BlasMatrixDomainSub<Matrix1, Matrix2, Matrix3 >::operator()(Matrix1& C, const Matrix2& A, const Matrix3& B) const
	{
		linbox_check( A.rowdim() == B.rowdim());
		linbox_check( C.rowdim() == A.rowdim());
		linbox_check( A.coldim() == B.coldim());
		linbox_check( C.coldim() == A.coldim());
		FFLAS::fsub (A.field(), C.rowdim(), C.coldim(),A.getPointer(), A.getStride(),B.getPointer(), B.getStride(),C.getPointer(), C.getStride());
		return C;
	}
	// SubIn
	template<class Matrix1, class Matrix3>
	Matrix1&
	BlasMatrixDomainSubin<Matrix1, Matrix3 >::operator()(Matrix1& C, const Matrix3& B) const
	{
		linbox_check( C.rowdim() == B.rowdim());
		linbox_check( C.coldim() == B.coldim());
		FFLAS::fsubin (C.field(), C.rowdim(), C.coldim(),B.getPointer(), B.getStride(),C.getPointer(), C.getStride());
		return C;
	}
} // LinBox

namespace LinBox
{ /* Copy */
	//Copy
	template<class Matrix1, class Matrix2>
	Matrix1&
	BlasMatrixDomainCopy<Matrix1, Matrix2 >::operator()(Matrix1& B, const Matrix2& A) const
	{
		linbox_check( A.rowdim() == B.rowdim());
		linbox_check( A.coldim() == B.coldim());
        FFLAS::fassign (A.field(), A.rowdim(), A.coldim(), A.getPointer(), A.getStride(),B.getPointer(), B.getStride());
		return B;
	}
} // LinBox

namespace LinBox
{ /* Solve */

	/*
	 * Specialization for Operand of type BlasMatrix, BlasSubmatrix, BlasVector and BlasSubvector
	 */
	template <class Operand1, class Matrix, class Operand2>
	Operand1&
	BlasMatrixDomainLeftSolve<Operand1, Matrix, Operand2 >::operator() (Operand1& X, const Matrix& A, const Operand2& B) const
	{
		PLUQMatrix<typename Matrix::Field> PLUQ(A);
		PLUQ.left_solve(X, B);
		return X;
	}
	template <class Operand1, class Matrix, class Operand2>
	Operand1&
	BlasMatrixDomainLeftSolve<Operand1, Matrix, Operand2 >::operator() (const Matrix& A, Operand1& B) const
	{
		PLUQMatrix<typename Matrix::Field> PLUQ(A);
		PLUQ.left_solve(B);
		return B;
	}
	template <class Operand1, class Matrix, class Operand2>
	Operand1&
	BlasMatrixDomainRightSolve<Operand1, Matrix, Operand2 >::operator() (Operand1& X, const Matrix& A, const Operand2& B) const
	{
		PLUQMatrix<typename Matrix::Field> PLUQ(A);
		PLUQ.right_solve(X, B);
		return X;
	}
	template <class Operand1, class Matrix, class Operand2>
	Operand1&
	BlasMatrixDomainRightSolve<Operand1, Matrix, Operand2 >::operator() ( const Matrix& A, Operand1& B) const
	{
		PLUQMatrix<typename Matrix::Field> PLUQ(A);
		PLUQ.right_solve(B);
		return B;
	}


	/*
	 * ********************************************************
	 * *** Specialization for TriangularBlasMatrix<Field,_Rep> ***
	 * ********************************************************
	 */
  
    
	template <class Matrix1, class Matrix2, class Matrix3>
    class BlasMatrixDomainLeftSolve<Matrix1, TriangularBlasMatrix<Matrix2>, Matrix3> {
    public:
        Matrix1& operator() (Matrix1& X, const TriangularBlasMatrix<Matrix2>& A, const Matrix3& B) const
        {
            linbox_check( X.rowdim() == B.rowdim());
            linbox_check( X.coldim() == B.coldim());
            X.copy(B);        
            return (*this)(A, X);
        }        
        Matrix3& operator() (const TriangularBlasMatrix<Matrix2>& A, Matrix3& B) const
        {
            linbox_check( A.rowdim() == A.coldim());
            linbox_check( A.coldim() == B.rowdim());
            FFLAS::ftrsm( A.field(), FFLAS::FflasLeft, (FFLAS::FFLAS_UPLO) A.getUpLo(),isTransposed<Matrix3>::value,(FFLAS::FFLAS_DIAG) A.getDiag(),
                          A.rowdim(), B.coldim(), A.field().one, A.getPointer(), A.getStride(), B.getPointer(), B.getStride());
            return B;
        }
    };
    

    template <class Matrix1, class Matrix2, class Matrix3>
	class BlasMatrixDomainRightSolve<Matrix1, TriangularBlasMatrix<Matrix3>, Matrix2> {
    public:
        Matrix1&operator() (Matrix1& X, const TriangularBlasMatrix<Matrix3>& A, const Matrix2& B) const
        {
            linbox_check( X.rowdim() == B.rowdim());
            linbox_check( X.coldim() == B.coldim());
            X.copy(B);        
            return (*this)(A, X);
        }        
        Matrix2& operator() (const TriangularBlasMatrix<Matrix3>& A, Matrix2& B) const
        {
            linbox_check( A.rowdim() == A.coldim());
            linbox_check( A.coldim() == B.rowdim());
            FFLAS::ftrsm( A.field(), FFLAS::FflasRight, (FFLAS::FFLAS_UPLO) A.getUpLo(),isTransposed<Matrix2>::value,(FFLAS::FFLAS_DIAG) A.getDiag(),
                          A.rowdim(), B.coldim(), A.field().one, A.getPointer(), A.getStride(), B.getPointer(), B.getStride());
            return B;
        }
    };

    template <class Matrix>    
    class BlasMatrixDomainLeftSolve<BlasVector<typename Matrix::Field>, TriangularBlasMatrix<Matrix> > {
    public:
        typedef typename Matrix::Field Field;
        BlasVector<Field>& operator() (BlasVector<Field>& x,const TriangularBlasMatrix<Matrix>& A,  const BlasVector<Field>& b) const {
            linbox_check( A.rowdim() == b.size());
            linbox_check( A.coldim() == x.size());
            linbox_check( x.size() == b.size());
            x.copy(b);        
            return (*this)(A, x);
        }      
        BlasVector<Field>& operator() (const TriangularBlasMatrix<Matrix>& A,  BlasVector<Field>& b) const
        {
            linbox_check( A.rowdim() == A.coldim());
            linbox_check( A.rowdim() == b.size());
            FFLAS::ftrsv( A.field(),(FFLAS::FFLAS_UPLO)A.getUpLo(), FFLAS::FflasNoTrans, (FFLAS::FFLAS_DIAG)A.getDiag(),b.size(), A.getPointer(), A.getStride(),b.getPointer(),1);
            return b;
        }
    };
    

    template <class Matrix>
    class BlasMatrixDomainRightSolve<BlasVector<typename Matrix::Field>, TriangularBlasMatrix<Matrix> >{
    public:
        typedef typename Matrix::Field Field;
        BlasVector<Field>& operator() (BlasVector<Field>& x, const TriangularBlasMatrix<Matrix>& A, const BlasVector<Field>& b) const {
            linbox_check( A.rowdim() == b.size());
            linbox_check( A.coldim() == x.size());
            linbox_check( x.size() == b.size());
            x.copy(b);        
            return (*this)(A, x);
        }
        BlasVector<Field>& operator() (const TriangularBlasMatrix<Matrix>& A, BlasVector<Field>& b) const
        {
            linbox_check( A.rowdim() == A.coldim());
            linbox_check( A.rowdim() == b.size());
            FFLAS::ftrsv( A.field(),(FFLAS::FFLAS_UPLO)A.getUpLo(), FFLAS::FflasTrans, (FFLAS::FFLAS_DIAG)A.getDiag(),b.size(), A.getPointer(), A.getStride(),b.getPointer(),1);
            return b;
        }
    };

    
    template <class Matrix, class Vect>
    class  BlasMatrixDomainLeftSolve<BlasSubvector<Vect>, TriangularBlasMatrix<Matrix> >{
    public:
        BlasSubvector<Vect>& operator() (BlasSubvector<Vect>& x,const TriangularBlasMatrix<Matrix>& A,  const BlasSubvector<Vect>& b) const{
            linbox_check( A.rowdim() == b.size());
            linbox_check( A.coldim() == x.size());
            linbox_check( x.size() == b.size());
            x.copy(b);        
            return (*this)(A, x);
        }            
        BlasSubvector<Vect>& operator() (const TriangularBlasMatrix<Matrix>& A,  BlasSubvector<Vect>& b) const
        {
			linbox_check( A.rowdim() == A.coldim());
			linbox_check( A.rowdim() == b.size());
            FFLAS::ftrsv( A.field(),(FFLAS::FFLAS_UPLO)A.getUpLo(), FFLAS::FflasNoTrans, (FFLAS::FFLAS_DIAG)A.getDiag(),b.size(), A.getPointer(), A.getStride(),b.getPointer(),1);
			return b;
        }
    };
    
    
	template <class Matrix, class Vect>
    class BlasMatrixDomainRightSolve<BlasSubvector<Vect>, TriangularBlasMatrix<Matrix> >{
    public:
        BlasSubvector<Vect>& operator() (BlasSubvector<Vect>& x,const TriangularBlasMatrix<Matrix>& A,  const BlasSubvector<Vect>& b) const{
            linbox_check( A.rowdim() == b.size());
            linbox_check( A.coldim() == x.size());
            linbox_check( x.size() == b.size());
            x.copy(b);        
            return (*this)(A, x);
        }            
        BlasSubvector<Vect>& operator() (const TriangularBlasMatrix<Matrix>& A, BlasSubvector<Vect>& b) const
        {
			linbox_check( A.rowdim() == A.coldim());
			linbox_check( A.rowdim() == b.size());
            FFLAS::ftrsv( A.field(),A.getUpLo(), FFLAS::FflasTrans, A.getDiag(),b.size(), A.getPointer(), A.getStride(),b.getPointer(),1);
			return b;
        }
    };


} // LinBox

namespace LinBox
{ /* Minpoly Charpoly */

	template<class Polynomial, class Matrix>
	Polynomial&
	BlasMatrixDomainMinpoly<Polynomial,Matrix>::operator() (Polynomial& P, const Matrix& A) const
	{
		commentator().start ("Givaro::Modular Dense Minpoly ", "MDMinpoly");
		size_t n = A.coldim();
		linbox_check( n == A.rowdim());
		FFPACK::MinPoly( A.field(), P, n, A.getPointer(), A.getStride());
		commentator().report(Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION) << "minpoly with " << P.size() << " coefficients" << std::endl;
		commentator().stop ("done", NULL, "MDMinpoly");
		return P;
	}

	template<class Polynomial, class Matrix>
	Polynomial &
	BlasMatrixDomainCharpoly<Polynomial,Matrix>::operator() (Polynomial    &P, Matrix   &A) const
	{
		size_t n = A.coldim();
		P.clear();
		linbox_check( n == A.rowdim());
        typename Matrix::Field::RandIter G(A.field());
        typename Givaro::Poly1Dom< typename Matrix::Field> PolDom(A.field());
        FFPACK::CharPoly (PolDom, P, n, A.getPointer(), A.getStride(), G);
		return P;
	}

} //end of namespace LinBox

#endif // __LINBOX_matrix_matrixdomain_blas_matrix_domain_INL

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
