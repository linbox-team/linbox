/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
/* linbox/algorithms/blas-domain.inl
 * Copyright (C) 2004 Pascal Giorgi, Clément Pernet
 *
 * Written by :
 *               Pascal Giorgi  pascal.giorgi@ens-lyon.fr
 *               Clément Pernet clement.pernet@imag.fr
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



#ifndef __BLAS_MATRIX_DOMAIN_INL
#define __BLAS_MATRIX_DOMAIN_INL

namespace LinBox {


	/*
	 * **********************************************
	 * *** Specialization for BlasMatrix<Element> ***
	 * **********************************************
	*/	

	// Inversion
	template <class Field>
	template <>
	inline BlasMatrix<typename Field::Element>& 
	BlasMatrixDomain<Field>::inv<BlasMatrix<typename Field::Element> > (const BlasMatrix<typename Field::Element>& A,
									    BlasMatrix<typename Field::Element>& Ainv) const{}

	// Rank
	template <class Field>
	template <>	
	inline unsigned int 
	BlasMatrixDomain<Field>::rank<BlasMatrix<typename Field::Element> > (const BlasMatrix<typename Field::Element>& A) const{
		BlasMatrix<typename Field::Element> tmp(A);
		return rankin(tmp);
	}

	// in-place Rank (the matrix is modified)
	template <class Field>
	template <>	
	inline unsigned int 
	BlasMatrixDomain<Field>::rankin<BlasMatrix<typename Field::Element> > (BlasMatrix<typename Field::Element>& A) const{
		return FFLAPACK::rank(_F, A.rowdim(), A.coldim(),A.getPointer(), A.getStride());
	}

	// determinant
	template <class Field>
	template <>	
	inline Element& 
	BlasMatrixDomain<Field>::det<BlasMatrix<typename Field::Element> >(const BlasMatrix<typename Field::Element>& A) const{
		BlasMatrix<typename Field::Element> tmp(A);
		return detin(tmp);
	}

	//in-place Determinant (the matrix is modified
	template <class Field>
	template <>	
	inline Element& 
	BlasMatrixDomain<Field>::detin<BlasMatrix<typename Field::Element> > (BlasMatrix<typename Field::Element>& A) const{
		return FFLAPACK::det(_F, A.rowdim(), A.coldim(),A.getPointer(), A.getStride());
	}
		
	/*
	 * Solvers with Matrix right or left hand side
	 */ 

	// non-singular linear solve with matrix right hand side 
	// AX=B
	template <class Field>
	template <class Operand, class Matrix>
	inline Operand& BlasMatrixDomain<Field>::left_solve (Operand& X, const Matrix& A, const Operand& B) const {}
	
	// non-singular linear solve with matrix right hand side, the result is stored in-place in B
	// AX=B , (B<-X)
	template <class Field>
	template <class Operand,class Matrix>
	inline Operand& BlasMatrixDomain<Field>::left_solve (const Matrix& A, Operand& B) const{}
	
	// non-singular linear solve with matrix right hand side 
	// XA=B
	template <class Field>
	template <class Operand, class Matrix>
	inline Operand& BlasMatrixDomain<Field>::right_solve (Operand& X, const Matrix& A, const Operand& B) const{}
	
	// non-singular linear solve with matrix right hand side, the result is stored in-place in B
	// XA=B , (B<-X)
	template <class Field>
	template <class Operand, class Matrix>
	inline Operand& BlasMatrixDomain<Field>::right_solve (const Matrix& A, Operand& B) const{}
		

	/*
	 * ********************************************************
	 * *** Specialization for TriangularBlasMatrix<Element> ***
	 * ********************************************************
	 */ 


	/*
	 * specialization for Operand of type BlasMatrix<Element>
	 */
		
	template <class Field>
	template <>	
	inline  BlasMatrix<typename Field::Element>& 
	BlasMatrixDomain<Field>::left_solve<BlasMatrix<typename Field::Element>,
					    TriangularBlasMatrix<typename Field::Element> >	(BlasMatrix<typename Field::Element>& X,
												 const TriangularBlasMatrix<typename Field::Element>& A,
												 const BlasMatrix<typename Field::Element>& B) const{
		
		linbox_check( X.rowdim() == B.rowdim());
		linbox_check( X.coldim() == B.coldim());
		
		typename BlasMatrix<typename Field::Element>::ConstRawIterator  Biter =   B.rawBegin();
		typename BlasMatrix<typename Field::Element>::RawIterator       Xiter =   X.rawBegin();

		for (; Biter != B.rawEnd(); ++Biter,++Xiter)
			F.assign(*Xiter,*Biter);
		
		left_solve(A,X);
		return X;
	}
		
	template <class Field>
	template <>	
	inline  BlasMatrix<typename Field::Element>& 
	BlasMatrixDomain<Field>::left_solve<BlasMatrix<typename Field::Element>,
					    TriangularBlasMatrix<typename Field::Element> >  (const TriangularBlasMatrix<typename Field::Element>& A, 
											      BlasMatrix<typename Field::Element>& B) const{
			
		linbox_check( A.rowdim() == A.coldim());
		linbox_check( A.coldim() == B.rowdim());		

		switch (A.getUpLo()) {
		case up:
			switch(A.getDiag()) {
			case unit:
				FFLAS::ftrsm( _F, 
					      FFLAS::FflasLeft,FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasUnit,
					      A.rowdim(), B.coldim(),_One,A.getPointer(),A.getStride(),B.getPointer(),B.getStride());
			case nonunit:
				FFLAS::ftrsm( _F, 
					      FFLAS::FflasLeft,FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit,
					      A.rowdim(), B.coldim(),_One,A.getPointer(),A.getStride(),B.getPointer(),B.getStride());
			default:
				throw LinboxError ("Error in BlasMatrixDomain (triangular matrix not well defined)");
			}
		case low:
			switch(A.getDiag()) {
			case unit:
				FFLAS::ftrsm( _F, 
					      FFLAS::FflasLeft,FFLAS::FflasLower,FFLAS::FflasNoTrans,FFLAS::FflasUnit,
					      A.rowdim(), B.coldim(),_One,A.getPointer(),A.getStride(),B.getPointer(),B.getStride());
			case nonunit:
				FFLAS::ftrsm( _F, 
					      FFLAS::FflasLeft,FFLAS::FflasLower,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit,
					      A.rowdim(), B.coldim(),_One,A.getPointer(),A.getStride(),B.getPointer(),B.getStride());
			default:
				throw LinboxError ("Error in BlasMatrixDomain (triangular matrix not well defined)");
			}
		default:
			throw LinboxError ("Error in BlasMatrixDomain (triangular matrix not well defined)");
				
		}
		return B;
	}
		
	template <class Field>
	template <>	
	inline BlasMatrix<typename Field::Element>&
	BlasMatrixDomain<Field>::right_solve<BlasMatrix<typename Field::Element>,
					    TriangularBlasMatrix<typename Field::Element> > (BlasMatrix<typename Field::Element>& X,
											     const TriangularBlasMatrix<typename Field::Element>& A,
											     const BlasMatrix<typename Field::Element>& B) const{
			
		linbox_check( X.rowdim() == B.rowdim());
		linbox_check( X.coldim() == B.coldim());

		typename BlasMatrix<typename Field::Element>::ConstRawIterator  Biter =   B.rawBegin();
		typename BlasMatrix<typename Field::Element>::RawIterator       Xiter =   X.rawBegin();

		for (; Biter != B.rawEnd(); ++Biter,++Xiter)
			F.assign(*Xiter,*Biter);
		
		right_solve(A,X);
		return X;
	}
		
	template <class Field>
	template <>	
	inline BlasMatrix<typename Field::Element>&
	BlasMatrixDomain<Field>::right_solve<BlasMatrix<typename Field::Element>,
					    TriangularBlasMatrix<typename Field::Element> > (const TriangularBlasMatrix<typename Field::Element>& A, 
											     const BlasMatrix<typename Field::Element>& B) const{

		linbox_check( A.rowdim() == A.coldim());
		linbox_check( B.coldim() == A.rowdim());		
		
		switch (A.getUpLo()) {
		case up:
			switch(A.getDiag()) {
			case unit:
				FFLAS::ftrsm( _F, 
					      FFLAS::FflasRight,FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasUnit,
					      B.rowdim(), A.coldim(),_One,A.getPointer(),A.getStride(),B.getPointer(),B.getStride());
			case nonunit:
				FFLAS::ftrsm( _F, 
					      FFLAS::FflasRight,FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit,
					      B.rowdim(), A.coldim(),_One,A.getPointer(),A.getStride(),B.getPointer(),B.getStride());
			default:
				throw LinboxError ("Error in BlasMatrixDomain (triangular matrix not well defined)");
			}
		case low:
			switch(A.getDiag()) {
			case unit:
				FFLAS::ftrsm( _F, 
					      FFLAS::FflasRight,FFLAS::FflasLower,FFLAS::FflasNoTrans,FFLAS::FflasUnit,
					      B.rowdim(), A.coldim(),_One,A.getPointer(),A.getStride(),B.getPointer(),B.getStride());
			case nonunit:
				FFLAS::ftrsm( _F, 
					      FFLAS::FflasRight,FFLAS::FflasLower,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit,
					      B.rowdim(), A.coldim(),_One,A.getPointer(),A.getStride(),B.getPointer(),B.getStride());
			default:
				throw LinboxError ("Error in BlasMatrixDomain (triangular matrix not well defined)");
			}
		default:
			throw LinboxError ("Error in BlasMatrixDomain (triangular matrix not well defined)");
				
		}
		return B;
	}



	/*
	 * specialization for Operand of type std::vector<Element>
	 */

	template <class Field>
	template <>	
	inline std::vector<typename Field::Element>&
	BlasMatrixDomain<Field>::left_solve<TriangularBlasMatrix<typename Field::Element>, 
					    std::vector<typename Field::Element>  (std::vector<typename Field::Element>& x,
										   const TriangularBlasMatrix<typename Field::Element>& A, 
										   const std::vector<typename Field::Element>& b) const{

		linbox_check (X.size() == B.size());
		std::vector<typename Field::Element>::const_iterator biter = b.begin();
		std::vector<typename Field::Element>::iterator       xiter = x/begin();   
		for (;biter!=b.end();++biter,++xiter)
			_F.assign(*xiter,*biter);
		left_solve(A,x);		
	}
		
	
	template <class Field>
	template <>	
	inline std::vector<typename Field::Element>&
	BlasMatrixDomain<Field>::left_solve<TriangularBlasMatrix<typename Field::Element>, 
					    std::vector<typename Field::Element> (const TriangularBlasMatrix<typename Field::Element>& A,
										  std::vector<typename Field::Element>& b) const{

		linbox_check( A.rowdim() == A.coldim());
		linbox_check( A.rowdim() == b.size());


		switch (A.getUpLo()) {
		case up:
			switch(A.getDiag()) {
			case unit:
				FFLAS::ftrsv( _F, 
					      FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasUnit,
					      b.size(),A.getPointer(),A.getStride(),&b[0],1);
			case nonunit:
				FFLAS::ftrsv( _F, 
					      FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit,
					      b.size(),A.getPointer(),A.getStride(),&b[0],1);
			default:
				throw LinboxError ("Error in BlasMatrixDomain (triangular matrix not well defined)");
			}
		case low:
			switch(A.getDiag()) {
			case unit:
				FFLAS::ftrsv( _F, 
					      FFLAS::FflasLower,FFLAS::FflasNoTrans,FFLAS::FflasUnit,
					      b.size(),A.getPointer(),A.getStride(),&b[0],1);
			case nonunit:
				FFLAS::ftrsv( _F, 
					      FFLAS::FflasLower,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit,
					      b.size(),A.getPointer(),A.getStride(),&b[0],1);
			default:
				throw LinboxError ("Error in BlasMatrixDomain (triangular matrix not well defined)");
			}
		default:
			throw LinboxError ("Error in BlasMatrixDomain (triangular matrix not well defined)");
				
		}
		
	}
		
	
	template <class Field>
	template <>	
	inline std::vector<typename Field::Element>&
	BlasMatrixDomain<Field>::right_solve< std::vector<typename Field::Element>,
					      TriangularBlasMatrix<typename Field::Element> > (std::vector<typename Field::Element>& x,
											       const TriangularBlasMatrix<typename Field::Element>& A,
											       const std::vector<typename Field::Element>& b) const{

		linbox_check (X.size() == B.size());
		std::vector<typename Field::Element>::const_iterator biter = b.begin();
		std::vector<typename Field::Element>::iterator       xiter = x/begin();   
		for (;biter!=b.end();++biter,++xiter)
			_F.assign(*xiter,*biter);
		right_solve(A,x);				
	}
		
	// non-singular linear solve with matrix right hand side, the result is stored in-place in b
	template <class Field>
	template <>	
	inline std::vector<typename Field::Element>&
	BlasMatrixDomain<Field>::right_solve< std::vector<typename Field::Element>,
					      TriangularBlasMatrix<typename Field::Element> > (const TriangularBlasMatrix<typename Field::Element>& A, 
											       std::vector<typename Field::Element>& b) const{

		linbox_check( A.rowdim() == A.coldim());
		linbox_check( A.coldim() == b.size());


		switch (A.getUpLo()) {
		case up:
			switch(A.getDiag()) {
			case unit:
				FFLAS::ftrsv( _F, 
					      FFLAS::FflasUpper,FFLAS::FflasTrans,FFLAS::FflasUnit,
					      b.size(),A.getPointer(),A.getStride(),&b[0],1);
			case nonunit:
				FFLAS::ftrsv( _F, 
					      FFLAS::FflasUpper,FFLAS::FflasTrans,FFLAS::FflasNonUnit,
					      b.size(),A.getPointer(),A.getStride(),&b[0],1);
			default:
				throw LinboxError ("Error in BlasMatrixDomain (triangular matrix not well defined)");
			}
		case low:
			switch(A.getDiag()) {
			case unit:
				FFLAS::ftrsv( _F, 
					      FFLAS::FflasLower,FFLAS::FflasTrans,FFLAS::FflasUnit,
					      b.size(),A.getPointer(),A.getStride(),&b[0],1);
			case nonunit:
				FFLAS::ftrsv( _F, 
					      FFLAS::FflasLower,FFLAS::FflasTrans,FFLAS::FflasNonUnit,
					      b.size(),A.getPointer(),A.getStride(),&b[0],1);
			default:
				throw LinboxError ("Error in BlasMatrixDomain (triangular matrix not well defined)");
			}
		default:
			throw LinboxError ("Error in BlasMatrixDomain (triangular matrix not well defined)");
				
		}
	}
	       
	/*
	 *  Method to apply Permutation
	 */
	// Apply a BlasPermutation matrix P to a dense matrix A: 
	// B = A.P 
	template <class Field>
	template <class Operand>	
	inline Operand& BlasMatrixDomain<Field>::applyRight(  Operand& B, const Operand& A, const BlasPermutation& P){}

	// B = A.P^t
	template <class Field>
	template <class Operand>	
	inline Operand& BlasMatrixDomain<Field>::applyRightTranspose(  Operand& B, const Operand& A, const BlasPermutation& P){}

	// B = P.A 
	template <class Field>
	template <class Operand>	
	inline Operand& BlasMatrixDomain<Field>::applyLeft(  Operand& B, const Operand& A, const BlasPermutation& P){}
		
	// B = A.P^t
	template <class Field>
	template <class Operand>	
	inline Operand& BlasMatrixDomain<Field>::applyLeftTranspose(  Operand& B, const Operand& A, const BlasPermutation& P){}
		
	// In place apply.
	// A = A.P 
	template <class Field>
	template <class Operand>	
	inline Operand& BlasMatrixDomain<Field>::applyinRight( Operand& A, const BlasPermutation& P){}
		
	// A = A.P^t
	template <class Field>
	template <class Operand>	
	inline Operand& BlasMatrixDomain<Field>::applyinRightTranspose( Operand& A, const BlasPermutation& P){}       

	// A = P.A 
	template <class Field>
	template <class Operand>	
	inline Operand& BlasMatrixDomain<Field>::applyinLeft( Operand& A, const BlasPermutation& P){}
		
	// A = A.P^t
	template <class Field>
	template <class Operand>	
	inline Operand& BlasMatrixDomain<Field>::applyinLeftTranspose( Operand& A, const BlasPermutation& P){}

	// Conversion from BlasPermutation to BlackBoxPermutation 
	template <class Field>
	inline Permutation& BlasMatrixDomain<Field>::convert ( Permutation& P, const BlasPermutation& BP ){}


} //end of namespace LinBox
	
#endif
