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

#include <iostream.h>
#include <vector>
#include "linbox/algorithms/lsp-tools.h"



/* This class provide the decomposition LSP of a dense matrix which is stored contiguously by row
 * and where is row size is lower or equal than the column size. 
 *
 * The class is templatized by a Matrix type which must have the function FullIterator(). This function returns a pointer
 * to the 1st element of the matrix.
 */

namespace LinBox {

	template <class Field,class Matrix>
	class lsp {

		typedef typename Field::Element Element;
		typedef std::vector<int>  Perm;

	private:
    
		Matrix _L;
		Matrix _S;
		Perm   _P;
		Field  _F;
		int _m;
		int _n;
		Element Zero;
		Element One;
		int _rank;

	public:

		// Constructor of LSP.
		// Initialize the matrix L,S and the permutation P and some constant
		lsp (const Field& F, const Matrix& M) : 
			_F(F),
			_m(M.rowdim()),
			_n(M.coldim()),
			_S(M) 
		{				
			_P.reserve(_n);
			_P.resize(_n);
			for (int i=0;i<_n;i++)
				_P[i]=i;
			_F.init(Zero,0UL);
			_F.init(One,1UL);
		
			Matrix Identity(M.rowdim(),M.rowdim());
			for (int i=0;i<M.rowdim();i++)
				Identity.setEntry(i,i,One);
			_L=Identity;			
			//_rank = LSPCompute (_m,_n,_L.FullIterator(),_m, _S.FullIterator(),_n, _P);
		}
	    
		// Copy constructor
		lsp (const lsp& fact) :
			_F(fact._F),
			_L(fact._L),
			_S(fact._S),
			_P(fact._P),
			Zero(fact.Zero),
			One(fact.One),
			_rank(fact._rank)
		{}

		// Function to get the matrix L.
		const Matrix& get_L () const {
			return _L;
		}
		// Function to print out the matrix L.
		ostream& write_L (ostream& os, const Field& F) {
			return _L.write(os,F);
		}
		// Function to get the matrix S.
		const Matrix& get_S () const {
			return _S;
		}
		// Function to print out the matrix S.
		ostream& write_S (ostream& os, const Field& F) {
			return _S.write(os,F);
		}
		// function to get the permutation P.
		const Perm& get_P () const {
			return _P;
		}

		// function to print out the matrix of permutation P.
		ostream& write_P (ostream& os, const Field& F) {
			Matrix tmp(_n,_n);
			for (int i=0;i<_n;i++)
				tmp.setEntry(i,_P[i],One);			
			return tmp.write(cout,F);
		}
		
		
		// function to get the rank of the matrix passed as argument of LSP.
		const int& rank () const {
			return _rank;
		}
		
		// launcher of the computation of the LSP.
		void compute() 
		{  			
			_rank = LSPCompute (_m,_n,_L.FullIterator(),_m, _S.FullIterator(),_n, _P);			
		}
			

	protected:

		int LSPCompute (size_t m, size_t n,
				Element* L, int ldl,
				Element* S, int lds,
				Perm& P) {
      
			int rank=0;
			int rank_high;
			
			if ( m > 1) { 
			
				int m_up   = m >>1;
				int m_down = m - m_up;
				
				Element * L1 = L;
				Element * L2 = L + m_up*ldl;
				Element * S1 = S;
				Element * S2 = S + m_up*lds;
				// splitting and call of lower recursion level

				Perm P1(n);
				for (int i=0;i<n;i++) P1[i]=i;
			
				// Computing the LSP decomposition with the first half of rows from the entry matrix M.
				rank_high= LSPCompute (m_up,n,
						   L1,ldl,
						   S1,lds,
						   P1);
				
				rank+=rank_high;

				if (rank_high != 0) {
					// Application of transposed permutation of P1 with size(n*n) on A2 
					ApplyColPermTrans (_F,S2,m_down,n,lds,P1);
					
					ComputeG (_F, S1,m_up,rank_high,lds,S2,m_down,lds,L2,ldl);
										
					// updating of S2 with the first part of L2, computed above, and with S1. S2= S2 - L2*S1
					
					for (int i=0;i<m_down;i++)
						for (int j=0;j<rank_high;j++)
							_F.assign(*(S2+j+i*lds), Zero);
					
					Field_dgemm (_F,m_down,n-rank_high,m_up,-1,L2,ldl,S1+rank_high,lds,1,S2+rank_high,lds);
				
				}
				
				Perm P2(n-rank_high);
				for (int i=0;i<n-rank_high;i++) 
					P2[i]=i;
												
				// Computing the LSP decomposition with the second half of row from 
				// the entry matrix M which have been modified above , according to 
				// the LSP decomposition of the first half of rows from the entry matrix M.
				rank+= LSPCompute (m_down,n-rank_high,
						   L2+m_up,ldl,
						   S2+rank_high,lds,
						   P2);
				
				// Application of transposed permutation of P2 with size(n*n) on S1 with size(m_up,n) on columns from rank to the last. 
				ApplyColPermTrans (_F,S1+rank_high,m_up,n-rank_high,lds,P2);							        									       		
				
				int i=0;
				for (;i<rank_high;i++){
					P[i]=P1[i];}
				for (;i<n;i++){
					P[i]=P1[P2[i-rank_high]+rank_high];}
				
				
				return rank;
			}
			else { //last recusion level
				
				int idx=0;
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
		}	// end of the function LSPCompute


	}; //end of class lsp

} // end of namespace LinBox

#endif
