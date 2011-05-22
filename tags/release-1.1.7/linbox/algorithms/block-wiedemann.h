/* linbox/algorithms/block-wiedemann.h
 * Copyright (C) 2004 Pascal Giorgi
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


#ifndef __LINBOX_block_wiedemann_H
#define __LINBOX_block_wiedemann_H

#include <vector>

#include <linbox/integer.h>
#include <linbox/matrix/blas-matrix.h>
#include <linbox/algorithms/blas-domain.h>
#include <linbox/algorithms/blackbox-block-container.h>
#include <linbox/algorithms/block-massey-domain.h>
#include <linbox/vector/vector-domain.h>
#include <linbox/blackbox/transpose.h>

#include <linbox/util/error.h>
#include <linbox/util/debug.h>

namespace LinBox 
{

	template <class _Field>
	class BlockWiedemannSolver{

	public:
		typedef _Field                          Field;
		typedef typename Field::Element       Element;
		typedef typename Field::RandIter     RandIter;
		typedef std::vector<Element>         Vector;
		typedef BlasMatrix<Element>           Block;

	protected:
		Field                         _F;
		BlasMatrixDomain<Field>     _BMD;
		VectorDomain<Field>         _VDF;
		RandIter                   _rand;

	public:
		BlockWiedemannSolver (const Field &F) : _F(F), _BMD(F), _VDF(F), _rand(F) {}

		BlockWiedemannSolver (const Field &F, const RandIter &rand) : _F(F), _BMD(F), _VDF(F), _rand(rand) {}


		template <class Blackbox>
		Vector &solveNonSingular (Vector &x, const Blackbox &B, const Vector &y) const {

			
			Transpose<Blackbox> A(B);

			size_t m,n;
			m = A.rowdim();
			n = A.coldim();

			size_t p,q;
			integer tmp;
			tmp = m;
			//p = tmp.bitsize()-1;
			p=sqrt(tmp);
			tmp = n;
			//q = tmp.bitsize()-1;
			q=sqrt(tmp);
			cout<<"row block: "<<p<<endl;
			cout<<"col block: "<<q<<endl;


			Block U(p,m), UA(p-1,m), V(n,q);
			
			for (size_t i=0;i<n;++i)
				for (size_t j=0;j<q;++j)
					_rand.random(V.refEntry(i,j));
			
			for (size_t i=0;i<p-1;++i)
				for (size_t j=0;j<m;++j)
					_rand.random(UA.refEntry(i,j));
			
			Block::RowIterator        iter_U  = U.rowBegin();
			Block::ConstRowIterator   iter_UA = UA.rowBegin();
			++iter_U;
			for (; iter_UA != UA.rowEnd(); ++iter_UA, ++iter_U) 
				A.applyTranspose( *iter_U , *iter_UA );
				
			for (size_t i=0;i<m;++i)
				U.setEntry(0,i,y[i]);

			BlackboxBlockContainer<Field,Transpose<Blackbox> > Sequence (&A,_F,U,V);			
			BlockMasseyDomain <Field,BlackboxBlockContainer<Field,Transpose<Blackbox> > > MBD(&Sequence);
			
			std::vector<Block> minpoly;
			std::vector<size_t> degree;
			MBD.left_minpoly(minpoly,degree);
			

			size_t idx=0;
			if ( _F.isZero(minpoly[0].getEntry(0,0))) {

				size_t i=1;
				while ( _F.isZero(minpoly[0].getEntry(i,0)))
					++i;
				if (i == m)
					throw LinboxError(" block minpoly: matrix seems to be singular - abort");
				else 
					idx=i	;			
			}
			

			bool classic = false;
			if ( classic) {
				/*
				 * Compute the solution according to the polynomial combination
				 * given by each column of the idx-th row of MinPoly such that the constant term of
				 * the first element in this row is non zero.
				 * we use y and UA as projection (UA= U.A)
				 */
				size_t deg = degree[idx];		
				std::vector<Vector> combi(p,Vector(deg+1));
				for (size_t i=0;i<p;++i) 
					for (size_t k=0;k<deg+1;++k)
						combi[i][k]=minpoly[k].getEntry(idx,i);
					
				Vector lhs(n);			
				A.applyTranspose(lhs,y);
				_VDF.mulin(lhs,combi[0][deg]);
				Vector lhsbis(lhs);
				for (int i = deg-1 ; i > 0;--i) {
					_VDF.axpy (lhs, combi[0][i], y, lhsbis);
					A.applyTranspose (lhsbis, lhs);			
				}   
		
				Vector accu (lhs);
				for (size_t k=1;k<p;++k){			
					Vector row(m);
					for (size_t j=0;j<m;++j)
						row[j]=UA.getEntry(k-1,j);				
					A.applyTranspose(lhs,row);			
					_VDF.mulin(lhs,combi[k][deg]);
					Vector lhsbis(lhs);
					for (int i = deg-1 ; i >= 0;--i) {												
						_VDF.axpy (lhs, combi[k][i], row, lhsbis);
						A.applyTranspose (lhsbis, lhs);
					}   
				
				
					_VDF.addin(accu,lhs);
				}
			
				Element scaling;
				_F.init(scaling);
				_F.neg(scaling,combi[0][0]);
				_F.invin(scaling);
				_VDF.mul(x,accu,scaling);
				
			}
			else {
				/*
				 * Compute the solution according to the polynomial combination
				 * given by the product of the idx-th row of MinPoly and UA.
				 * this should decrease the number of sparse apply but increase memory requirement.
				 */
				size_t deg = degree[idx];		
				BlasMatrix<Element> idx_poly(deg+1,p-1);
				for (size_t i=0;i<deg+1;++i) 
					for (size_t j=0;j<p-1;++j)
						idx_poly.setEntry(i,j,minpoly[i].getEntry(idx,j+1));

				BlasMatrix<Element> Combi(deg+1,m);
				_BMD.mul(Combi,idx_poly,UA);
					

				Vector lhs(n),row(m);	
				for (size_t i=0;i<m;++i)
					row[i]= Combi.getEntry(deg,i);
					
				A.applyTranspose(lhs,row);					
				Vector lhsbis(lhs);
				for (int i = deg-1 ; i >= 0;--i) {
					for (size_t j=0;j<m;++j)
						row[j]= Combi.getEntry(i,j);
					_VDF.add (lhs,row,lhsbis);
					A.applyTranspose (lhsbis, lhs);			
				}   
					
				Vector accu (lhs);
					
					
				A.applyTranspose(lhs,y);
				_VDF.mulin(lhs,minpoly[deg].getEntry(idx,0));
				lhsbis=lhs;
				for (size_t i = deg-1 ; i > 0;--i) {
					_VDF.axpy (lhs,minpoly[i].getEntry(idx,0) , y, lhsbis);
					A.applyTranspose (lhsbis, lhs);			
				}  
					
				_VDF.addin(accu,lhs);
				Element scaling;
				_F.init(scaling);
				_F.neg(scaling,minpoly[0].getEntry(idx,0));
				_F.invin(scaling);
				_VDF.mul(x,accu,scaling);
			}
	
			return x;	
		}



	}; // end of class BlockWiedemannSolver
    

    
}// end of namespace LinBox

#endif //__LINBOX_block_wiedemann_H
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
