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
#include "linbox/fflapack/fflapack.h"
#include "linbox/blackbox/permutation.h"
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

	template <class Field>
	class lsp {

		typedef typename Field::Element Element;
		typedef std::vector<int>  Perm;

	private:
    
		Field  _F;

	public:

		// Constructor of LSP.
		// Initialize the matrix L,S and the permutation P and some constant
		lsp (const Field& F ) :
			_F(F)
			{}
	    
		// Copy constructor
		lsp (const lsp& fact) :
			_F(fact._F)
			{}

	        // function to get the vector of non-zero-row's index  in S.
		template<class Matrix>
		std::vector<int>& get_PivotIndex( const Matrix& S, std::vector<int>& P) {
			
			int m = S.rowdim();
			int n = S.coldim();
			P.resize( (m<n)?m:n ); 
			int i=0, ind =0;
			int j=0;
			Element* ms_it = _MS.FullIterator();
			while (i < m) {
				if (_F.isZero(*ms_it)){
					ms_it += n;
					i++;
				}
				else {
					P[ind++]=i++;
					ms_it += n+1;
				}
			}
			P.resize(ind);
			return P;
		}

		// function to get the vector of zero-rows'index in S.
		template<class Matrix>
		std::vector<int>& get_ZeroIndex(const Matrix& S, std::vector<int>& P ) {
			int m = S.rowdim();
			int n = S.coldim();
			P.resize( (m<n)?m:n ); 
			int i=0, ind =0;
			int j=0;
			Element* ms_it = _MS.FullIterator();
			while (i < m) {
				if (_F.isZero(*ms_it)){
					ms_it += n;
					P[ind++] = i++;
				}
				else {
					i++;
					ms_it += n+1;
				}
			}
			P.resize(ind);
			return P;
		}
		
		// launcher of the computation of the LSP. 
		// Computation is inplace in M
		// P is a Lapack-style permutation vector
		template<class Matrix>
		size_t computein ( Matrix& M, Perm& P, size_t stride = 0  ) 
		{
			size_t n = M.coldim();
			size_t m = M.rowdim();
			size_t r;
			P.resize(n);
			size_t Permut[n];
			size_t rowP[m];
			r = FFLAPACK::LUdivine( _F, FFLAS::FflasNonUnit, m, n, M.FullIterator(),
						(stride?stride:n), 
                                                Permut, FFLAPACK::FflapackLQUP, rowP );

			Perm::iterator it = P.begin();
			size_t* Pi=Permut;
			for (;it!=P.end();it++, ++Pi)
				*it = *Pi;
                        
			return r;
		}

		// launcher of the computation of the LSP. 
		// Computation is inplace in M
		// P is a BlackBox permutation matrix
		template<class Matrix, class Vector>
		size_t computein ( Matrix& M, Permutation<Vector>& P, size_t stride = 0 ) 
		{
			size_t n = M.coldim();
			size_t m = M.rowdim();
			size_t r;
			size_t Plapack[n];
			size_t rowP[m];
			r = FFLAPACK::LUdivine( _F, FFLAS::FflasUnit, m, n, M.FullIterator(),
						(stride?stride:n), 
                                                Plapack, FFLAPACK::FflapackLQUP, rowP );
			
			Perm Pbb(n);
			size_t tmp;
			for (size_t i=0;i<n;i++)
				Pbb[i] = i;
			for (size_t i=0;i<n;i++)
				if ( Plapack[i] != i ) {
					tmp = Pbb[i];
					Pbb[i] = Pbb[Plapack[i]];
					Pbb[Plapack[i]] =  tmp;
				}
			
			P = Permutation<Vector>(Pbb);
			return r;
		}


		
		// launcher of the computation of the LSP of M in the matrices L and S
		// P is a Lapack-style permutation vector
		template<class Matrix>
		size_t compute ( const Matrix& M, Matrix& L, Matrix& S, Perm& P, size_t stride = 0) 
		{  			
			size_t r;
			size_t m = M.rowdim();
			size_t n = M.coldim(); 
			S = M;
			r = computein( S, P, stride );
			
			Element * ms = S.FullIterator();
			Element * ml = L.FullIterator();
			Element one, zero;
			_F.init( one, 1UL );
			_F.init( zero, 0UL );
			size_t pivot = 0;
			size_t i,j,jl=0;
			for ( j=0;j<r;++j){
				while ( _F.isZero ( ms[j+n*pivot]) && jl<m){
					// inserting 0 column in ml
					for ( i = 0; i<m; ++i )
						_F.assign( ml[i*m+jl], zero );
					pivot++;
					jl++;
				}
				for ( i=0; i<pivot;++i)
					_F.assign(ml[i*m+jl], zero);
				_F.assign( ml[pivot*m+jl], one );
				for ( i=pivot+1;i<m;++i){
					_F.assign( ml[i*m+jl], ms[i*n+j] );
				}
				// search for the next non zero pivot
				pivot++;
				jl++;
				
			}
			
			pivot = 0;
			for ( i=0;i<n;i++){
				if (!_F.isZero(ms[i*n+pivot])){
					for (j=0;j<pivot;++j)
						ms[i*n+j]=zero;
					pivot++;
				}
				else{
					for (j=0;j<=pivot;++j)
						ms[i*n+j]=zero;
				}
			}
						
		return r;
		}
		
		// launcher of the computation of the LSP of M in the matrices L and S
		// P is a BlackBox permutation matrix
		template<class Matrix, class Vector>
		size_t compute ( const Matrix& M, Matrix& L, Matrix& S, Permutation<Vector>& P) 
		{  			
			// To be corrected ( construction of L is not correct ) 
			size_t r;
			size_t m = M.rowdim();
			size_t n = M.coldim(); 
			S = M;
			
			r = computein( S, P );
			
			Element * ms_it = S.FullIterator();
			Element * ml_it = L.FullIterator();
			Element one, zero;
			_F.init( one, 1UL );
			_F.init( zero, 0UL );
			for ( size_t i=0; i<m; ++i){
				for ( size_t j=0; j<i; ++j){
					_F.assign( *ml_it, *ms_it );
					_F.assign( *ms_it, zero );
					ml_it++;
					ms_it++;
				}
				_F.assign( *ml_it, *ms_it );
				_F.assign( *ms_it, one );
				ml_it++;
				ms_it += n-i;
				for ( size_t j=i+1; j<n; ++j){
					_F.assign( *ml_it, zero );
					ml_it++;
				}
			}
			
			return r;
		}
	}; /* end of class lsp */

} /* end of namespace LinBox */

#endif /* __LSP_H */
