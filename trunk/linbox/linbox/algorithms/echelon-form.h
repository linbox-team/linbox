/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/algorithms/sigma-basis.h
 * Copyright (C) 2006 Pascal Giorgi
 *
 * Written by Pascal Giorgi pascal.giorgi@univ-perp.fr
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



#ifndef __ECHELON_FORM_H
#define __ECHELON_FORM_H

#include <linbox/matrix/blas-matrix.h>
#include <linbox/algorithms/blas-domain.h>

#include <linbox/matrix/matrix-domain.h>
#include <linbox/matrix/factorized-matrix.h>

 

namespace LinBox {



	template<class Field>
	class EchelonFormDomain{
	
	private:
	
		Field                      _F;
		BlasMatrixDomain<Field>  _BMD;
		MatrixDomain<Field>       _MD;
	
	public:
		typedef typename Field::Element Element;

		// constructor
		EchelonFormDomain(const Field &F) : _F(F), _BMD(F), _MD(F) {}
		

		//  row echelon form 
		// E is supposed to be the zero matrix
		template<class Matrix>
		int rowEchelon(Matrix &E, const Matrix& A){

			size_t m,n, rank;
			m = A.rowdim();
			n = A.coldim();
			
			// get the transposed of A
			BlasMatrix<Element> At(n, m);			
			for (size_t i=0;i<m;++i)
				for (size_t j=0;j<n;++j)
					At.setEntry(j,i,A.getEntry(i,j));

			// compute the LQUP of At
			LQUPMatrix<Field> LQUP(_F, At);

			// get L
			BlasMatrix<Element> L(n,n);
			TriangularBlasMatrix<Element> LL(L, BlasTag::low, BlasTag::unit );
			LQUP.getL(LL);

			
			// permute L with Q
			BlasPermutation Qt = LQUP.getQ();		
			
			TransposedBlasMatrix<BlasPermutation> Q(Qt);
			_BMD.mulin_left(L, Q);		

			// get the rank
			rank = LQUP.getrank();
			//std::cout<<"get transpose :"<<m<<" "<<n<<" "<<rank<<"\n";;
			
			// read the transpose of the echelon form from the rank 1st column of L
			for (size_t i=0; i<rank;++i)
				for (size_t j=0;j<n;++j){					
					E.setEntry(i,j, L.getEntry(j,i));
				}
			return rank;
		}

		// row reduced echelon form (using copy)
		template<class Matrix>
		int rowReducedEchelon(Matrix &E, const Matrix& A){

			size_t m,n, rank;
			m = A.rowdim();
			n = A.coldim();
			Element zero, one;
			_F.init(zero,0);
			_F.init(one,1);
			
			// get the transposed of A
			BlasMatrix<Element> At(n, m);			
			for (size_t i=0;i<m;++i)
				for (size_t j=0;j<n;++j)
					At.setEntry(j,i,A.getEntry(i,j));

			// compute the LQUP of At
			LQUPMatrix<Field> LQUP(_F, At);

			// get L		
			BlasMatrix<Element> L(n,n);	
			TriangularBlasMatrix<Element> LL(L, BlasTag::low, BlasTag::unit );
			LQUP.getL(LL);

			// permute L such that L <- Q^t.L.Q
			// -> gather all the pivots in the top left corner
			BlasPermutation Qt = LQUP.getQ();
			TransposedBlasMatrix<BlasPermutation> Q(Qt);
			_BMD.mulin_left(L, Q);
			_BMD.mulin_right(Qt,L);


			// get the rank
			rank = LQUP.getrank();
			
			// Update the first r columns of L by Lrr^(-1)
			BlasMatrix<Element> Lr(L,0,0,rank,rank);
			TriangularBlasMatrix<Element> Lrr(Lr, BlasTag::low, BlasTag::unit);
			BlasMatrix<Element> Ln(L,rank,0,n-rank,rank);

			_BMD.right_solve(Lrr, Ln);
			for (size_t i=0;i<rank;++i){
				for (size_t j=0;j<i;++j)
					L.setEntry(i,j,zero);
			}
			
			// permute L such that L<-Q.L
			_BMD.mulin_right(Q,L);
			

			// read the transpose of the echelon form from the rank 1st column of L
			for (size_t i=0; i<rank;++i)
				for (size_t j=0;j<n;++j)
					E.setEntry(i,j, L.getEntry(j,i));
			return rank;
		}
			

		//  column echelon form (using copy)
		template<class Matrix>
		int columnEchelon(Matrix &E, const Matrix& A){

			size_t m,n, rank;
			m = A.rowdim();
			n = A.coldim();
			Element zero, one;
			_F.init(zero,0);
			_F.init(one,1);

			// copy  A in E
			for (size_t i=0;i<m;++i)
				for (size_t j=0;j<n;++j)
					E.setEntry(i,j,A.getEntry(i,j));

			// compute the LQUP of E
			LQUPMatrix<Field> LQUP(_F, E);

			// get the rank
			rank = LQUP.getrank();
					
			// permute L with Q
			BlasPermutation Qt = LQUP.getQ();				

			// Zero out upper triangular part of E
			for (size_t i=0;i<m;++i)
				for (size_t j=i;j<n;++j)
					E.setEntry(i,j,zero);

			// put one inplace of pivot
			for (size_t i=0;i<rank;++i)
				E.setEntry(i,*(Qt.getPointer()+i),one);
			
			return rank;
		}

		// column reduced echelon form (using copy)
		template<class Matrix>
		int columnReducedEchelon(Matrix &E, const Matrix& A){
			size_t m,n, rank;
		
			m = A.rowdim();
			n = A.coldim();
			Element zero, one;
			_F.init(zero,0);
			_F.init(one,1);

			// copy  A in E
			for (size_t i=0;i<m;++i)
				for (size_t j=0;j<n;++j)
					E.setEntry(i,j,A.getEntry(i,j));

			// compute the LQUP of E
			LQUPMatrix<Field> LQUP(_F, E);

			// get the rank
			rank = LQUP.getrank();
					
			// permute L with Q
			BlasPermutation Qt = LQUP.getQ();		
			TransposedBlasMatrix<BlasPermutation> Q(Qt);

			// Zero out upper triangular part of E
			for (size_t i=0;i<m;++i)
				for (size_t j=i;j<n;++j)
					E.setEntry(i,j,zero);

			// put one inplace of pivot
			for (size_t i=0;i<rank;++i)
				E.setEntry(i,*(Qt.getPointer()+i),one);
			
			
			// Update the first r columns of E by Err^(-1)
			BlasMatrix<Element> Er(E,0,0,rank,rank);
			TriangularBlasMatrix<Element> Err(Er, BlasTag::low, BlasTag::unit);
			BlasMatrix<Element> En(E,rank,0,m-rank,rank);

			_BMD.right_solve(Err, En);

			for (size_t i=0;i<rank;++i){
				for (size_t j=0;j<i;++j)
					E.setEntry(i,j,zero);
			}
			
			// permute L such that L<-Q.E
			_BMD.mulin_right(Q,E);


			return rank;

		}



		template<class Matrix>
		void write_maple(const char* name, const Matrix& A) {
			
			size_t m,n;
			m = A.rowdim();
			n = A.coldim();
			std::cout<<name<<":= Matrix([";

			for (size_t i=0;i<m-1;++i){
					std::cout<<"[";
					for (size_t j=0;j<n-1;++j)
						_F.write(std::cout,A.getEntry(i,j))<<",";
					_F.write(std::cout, A.getEntry(i,n-1))<<"] , ";
				}
				std::cout<<"[";
				for (size_t j=0;j<n-1;++j)
					_F.write(std::cout,A.getEntry(m-1,j))<<",";				
				_F.write(std::cout, A.getEntry(m-1,n-1))<<"]]);\n ";
		}

	};
 
}

#endif
