/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/algorithms/lsp.h
 * Copyright (C) 2003 Pascal Giorgi
 *
 * Written by Pascal Giorgi pascal.giorgi@ens-lyon.fr
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

#ifndef __LSP_H
#define __LSP_H

#include <iostream>
#include <vector>
#include "linbox/algorithms/lsp-tools.h"

#ifdef __CHECK_LSP
#include "linbox/matrix/matrix-domain.h"
#endif



namespace LinBox {

/** This class provide the decomposition M = LSP of a dense matrix which is stored contiguously by row
 * and where row size is lower or equal than the column size. 
 *
 * The class is templatized by a Matrix type which must have the function FullIterator(). This function returns a pointer
 * to the 1st element of the matrix.
 */

	template <class Field,class Matrix>
	class lsp {

		typedef typename Field::Element Element;
		typedef std::vector<int>  Perm;

	private:
    
		Field  _F;
		int _m;
		int _n;
		Element Zero;
		Element One;
		int _rank;

		Matrix _ML;
		Matrix _MS;
		Perm   _MP;
#ifdef __CHECK_LSP
		Matrix _MM;
#endif	       

	public:

		// Constructor of LSP.
		// Initialize the matrix L,S and the permutation P and some constant
		lsp (const Field& F, const Matrix& M) : 
			_F(F),
			_m(M.rowdim()),
			_n(M.coldim()),
			_MS(M) 
#ifdef __CHECK_LSP
			, _MM(M)
#endif
		{				
			_MP.reserve(_n);
			_MP.resize(_n);
			for (int i=0;i<_n;i++)
				_MP[i]=i;
			_F.init(Zero,0UL);
			_F.init(One,1UL);
		
			Matrix Identity(M.rowdim(),M.rowdim());
			for (unsigned int i=0;i<M.rowdim();++i)
				Identity.setEntry(i,i,One);
			_ML=Identity;			
			//_rank = LSPCompute (_m,_n,_ML.FullIterator(),_m, _MS.FullIterator(),_n, _MP);
		}
	    
		// Copy constructor
		lsp (const lsp& fact) :
			_F(fact._F),
			_ML(fact._ML),
			_MS(fact._MS),
			_MP(fact._MP),
			Zero(fact.Zero),
			One(fact.One),
			_rank(fact._rank)
		{}

		// Function to get the matrix L.
		const Matrix& get_L () const {
			return _ML;
		}
		// Function to print out the matrix L.
		std::ostream& write_L (std::ostream& os, const Field& F) {
			return _ML.write(os,F);
		}
		// Function to get the matrix S.
		const Matrix& get_S () const {
			return _MS;
		}
		// Function to print out the matrix S.
		std::ostream& write_S (std::ostream& os, const Field& F) {
			return _MS.write(os,F);
		}
		// function to get the permutation P.
		const Perm& get_P () const {
			return _MP;
		}

		// function to print out the matrix of permutation P.
		std::ostream& write_P (std::ostream& os, const Field& F) {
			Matrix tmp(_n,_n);
			for (int i=0;i<_n;i++)
				tmp.setEntry(i,_MP[i],One);			
			return tmp.write(os,F);
		}
		
		
		// function to get the rank of the matrix passed as argument of LSP.
		const int& rank () const {
			return _rank;
		}


		// function to get the vector of non-zero-row's index  in S.
		std::vector<int> get_PivotIndex() {
			std::vector<int> P(_rank); 
			int i=0;
			int j=0;
			while (i< _rank) {
				if (_F.isZero(*(_MS.FullIterator()+i+j*_n)))
					j++;
				else {
					P[i]=j;
					i++;
					j++;
				}
			}
			

			return P;
		}

		// function to get the vector of zero-rows'index in S.
		std::vector<int> get_ZeroIndex() {
			std::vector<int> P(_m-_rank);			
			int k=0;
			for (int i=0;i<_m;i++) {
				bool zero=true;
				for (int j=0;j<_n;j++)
					if (!_F.isZero(*(_MS.FullIterator()+j+i*_n)))
						zero=false;
				if (zero) {
					P[k]=i;
					k++;
				}
			}
			
			return P;
		}

		
		// launcher of the computation of the LSP.
		void compute() 
		{  			
			_rank = (_n > _m)?  
				LSPCompute_moreCols (_m,_n,_ML.FullIterator(),_m, _MS.FullIterator(),_n, _MP):
				LSPCompute_moreRows (_m,_n,_ML.FullIterator(),_m, _MS.FullIterator(),_n, _MP);
#ifdef __CHECK_LSP
			Element one;
			_F.init(one,1UL);
			Matrix PP(_n,_n);
			for (int i=0;i<_n;i++)
				PP.setEntry(i,_MP[i],one);
			
			MatrixDomain<Field> MD(_F);
			Matrix MM(_m,_n);			
			MD.mul(MM,_ML,_MS);						
			MD.mulin(MM,PP);			
			if (! MD.areEqual(_MM,MM))
				cerr<<"LSP COMPUTED IS WRONG !!! \n";
			else{
				cerr<<"LSP IS CORRECT !!! \n";
				
			}
#endif
		}
			

	protected:

		unsigned int LSPCompute_moreCols (size_t m, size_t n,
						  Element* L, int ldl,
						  Element* S, int lds,
						  Perm& P) {
			

			unsigned int rank=0;
			unsigned int rank_high;
			
			if ( m > 1) { 
			
				int m_up   = m >>1;
				int m_down = m - m_up;
				
				Element * L1 = L;
				Element * L2 = L + m_up*ldl;
				Element * S1 = S;
				Element * S2 = S + m_up*lds;
				// splitting and call of lower recursion level

				Perm P1(n);
				for (unsigned int i=0;i<n;i++) P1[i]=i;
			
				// Computing the LSP decomposition with the first half of rows from the entry matrix M.
				rank_high= LSPCompute_moreCols (m_up,n,L1,ldl,S1,lds,P1);

				
				rank+=rank_high;
			
				if (rank_high != 0) {									

					// Application of transposed permutation of P1 with size(n*n) on A2 
					ApplyColPermTrans (_F,S2,m_down,n,lds,P1);					
					ComputeG (_F, S1,m_up,rank_high,lds,S2,m_down,lds,L2,ldl);

					// updating of S2 with the first part of L2, computed above, and with S1. S2= S2 - L2*S1     			      

					for (int i=0;i<m_down;i++)
						for (unsigned int j=0;j<rank_high;j++)
							_F.assign(*(S2+j+i*lds), Zero);
					
					Field_dgemm (_F,m_down,n-rank_high,m_up,-1,L2,ldl,S1+rank_high,lds,1,S2+rank_high,lds);
									
				}
				
				Perm P2(n-rank_high);
				for (unsigned int i=0;i<n-rank_high;i++) 
					P2[i]=i;					
				// Computing the LSP decomposition with the second half of row from 
				// the entry matrix M which have been modified above , according to 

				// the LSP decomposition of the first half of rows from the entry matrix M.					
				rank+= LSPCompute_moreCols (m_down,n-rank_high,L2+m_up,ldl,S2+rank_high,lds,P2);					

				// Application of transposed permutation of P2  on S1 with size(m_up,n) on columns from rank to the last. 
				ApplyColPermTrans (_F,S1+rank_high,m_up,n-rank_high,lds,P2);							        									       							
				unsigned int i=0;
				for (;i<rank_high;i++){
					P[i]=P1[i];}
				for (;i<n;i++){
					P[i]=P1[P2[i-rank_high]+rank_high];}
				
							
				return rank;
			}
			else { //last recusion level				
				unsigned int idx=0;
				while ( (idx < n) && (_F.isZero(*(S+idx))))
					idx++;
	
				// the row is zero
				if ( idx != n) {				
					// notify the permutation for the pivot according to the entire matrix M given at the first recursion level	
					P[0]= idx;
					P[idx]=0;
					if (idx != 0)
						{
							// swap the first element of S (wich is equal to zero) and the idx'th element of S. 
							Element tmp;
							_F.assign(tmp,*S);
							_F.assign(*S,*(S+idx));
						 	_F.assign(*(S+idx),tmp);
						}
					
					return 1;
				}
				else 
					return 0;				
			}
		}	// end of the function LSPCompute_moreCols

		
		unsigned int LSPCompute_moreRows (size_t m, size_t n,
						  Element* L, int ldl,
						  Element* S, int lds,
						  Perm& P) {
		
		unsigned int rank=0;
		unsigned int rank_high;
		
		if ( m > 1) { 
			
			int m_up   = m >>1;
			int m_down = m - m_up;
			
			Element * L1 = L;
			Element * L2 = L + m_up*ldl;
			Element * S1 = S;
			Element * S2 = S + m_up*lds;
			// splitting and call of lower recursion level
			
			Perm P1(n);
			for (unsigned int i=0;i<n;i++) P1[i]=i;
			
			// Computing the LSP decomposition with the first half of rows from the entry matrix M.
			rank_high= LSPCompute_moreRows (m_up,n,L1,ldl,S1,lds,P1);
				
			rank+=rank_high;
			
			if (rank_high != 0) {									
				// Application of transposed permutation of P1 with size(n*n) on A2 
				ApplyColPermTrans (_F,S2,m_down,n,lds,P1);					
				ComputeG (_F, S1,m_up,rank_high,lds,S2,m_down,lds,L2,ldl);
				// updating of S2 with the first part of L2, computed above, and with S1. S2= S2 - L2*S1     			      
				for (int i=0;i<m_down;i++)
					for (unsigned int j=0;j<rank_high;j++)
						_F.assign(*(S2+j+i*lds), Zero);
				if (n - rank_high > 0) {
					if (n - rank_high <= 0) cout<<"error: "<<_n<<" "<<rank<<" "<<n<<endl;;
					Field_dgemm (_F,m_down,n-rank_high,m_up,-1,L2,ldl,S1+rank_high,lds,1,S2+rank_high,lds);
				}				
			}
			if (n - rank_high > 0) {
				Perm P2(n-rank_high);
				for (unsigned int i=0;i<n-rank_high;i++) 
					P2[i]=i;					
				// Computing the LSP decomposition with the second half of row from 
				// the entry matrix M which have been modified above , according to 
				// the LSP decomposition of the first half of rows from the entry matrix M.					
				rank+= LSPCompute_moreRows (m_down,n-rank_high,L2+m_up,ldl,S2+rank_high,lds,P2);					

				// Application of transposed permutation of P2  on S1 with size(m_up,n) on columns from rank to the last. 
				ApplyColPermTrans (_F,S1+rank_high,m_up,n-rank_high,lds,P2);							        									       							
				unsigned int i=0;
				for (;i<rank_high;i++){
					P[i]=P1[i];}
				for (;i<n;i++){
					P[i]=P1[P2[i-rank_high]+rank_high];}
			}
			else {
				unsigned int i=0;
				for (;i<rank_high;i++){
					P[i]=P1[i];}	
			}				
			return rank;
		}
		else { //last recusion level				
			unsigned int idx=0;
			while ( (idx < n) && (_F.isZero(*(S+idx))))
				idx++;
			
			// the row is zero
			if ( idx != n) {				
				// notify the permutation for the pivot according to the entire matrix M given at the first recursion level	
				P[0]= idx;
				P[idx]=0;
				if (idx != 0)
					{
						// swap the first element of S (wich is equal to zero) and the idx'th element of S. 
						Element tmp;
						_F.assign(tmp,*S);
						_F.assign(*S,*(S+idx));
						_F.assign(*(S+idx),tmp);
					}
				
				return 1;
			}
			else 
				return 0;				
		}
	}	// end of the function LSPCompute_moreRows
	}; /* end of class lsp */

} /* end of namespace LinBox */

#endif /* __LSP_H */
