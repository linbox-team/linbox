/* linbox/algorithms/coppersmith.h
 * evolved from block-wiedemann.h by George Yuhasz
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


#ifndef __LINBOX_coppersmith_H
#define __LINBOX_coppersmith_H

#include <vector>
#include <iostream>
using namespace std;


#include "linbox/integer.h"
#include "linbox/algorithms/matrix-domain.h"
#include "linbox/algorithms/blackbox-block-container.h"
#include "linbox/algorithms/bm-seq.h"
#include "linbox/vector/vector-domain.h"

#include "linbox/util/error.h"
#include "linbox/util/debug.h"

namespace LinBox
{

	template <class _Field>
	class CoppersmithSolver{

	public:
		typedef _Field                          Field;
		typedef typename Field::Element       Element;
		typedef typename Field::RandIter     RandIter;
		typedef std::vector<Element>           Vector;
		typedef MatrixDomain<Field>::Matrix 	Block;

	protected:
		Field                         _field;
		MatrixDomain<Field>     _MD;
		VectorDomain<Field>         _VD;
		RandIter                   _rand;

	public:
		CoppersmithSolver(const Field &F) :
			_field(F), _MD(F), _VD(F), _rand(F)
		{}

		CoppersmithSolver (const Field &F, const RandIter &rand) :
			_field(F), _MD(F), _VD(F), _rand(rand)
		{}

		template <class Blackbox>
		Vector &solveNonSingular (Vector &x, const Blackbox &B, const Vector &y) const
		{ 
			// stub to be replaced by real code
			return y = x;
		}

#if 0 //  Pascal's left side code to be used as a template
		template <class Blackbox>
		Vector &solveNonSingular (Vector &x, const Blackbox &B, const Vector &y) const
		{
			Transpose<Blackbox> A(B);

			size_t m,n;
			m = A.rowdim();
			n = A.coldim();

			size_t p,q;
			integer tmp;
			tmp = m;
			p = tmp.bitsize()-1;
			//p=sqrt(tmp);
			tmp = n;
			q = tmp.bitsize()-1;
			//q=sqrt(tmp);
			//std::cout<<"row block: "<<p<<std::endl;
			//std::cout<<"col block: "<<q<<std::endl;


			Block U(_field,p,m), UA(_field,p-1,m), V(_field,n,q);

			for (size_t i=0;i<n;++i)
				for (size_t j=0;j<q;++j)
					_rand.random(V.refEntry(i,j));

			for (size_t i=0;i<p-1;++i)
				for (size_t j=0;j<m;++j)
					_rand.random(UA.refEntry(i,j));

			typename Block::RowIterator        iter_U  = U.rowBegin();
			typename Block::ConstRowIterator   iter_UA = UA.rowBegin();
			++iter_U;
			for (; iter_UA != UA.rowEnd(); ++iter_UA, ++iter_U)
				A.applyTranspose( *iter_U , *iter_UA );

			for (size_t i=0;i<m;++i)
				U.setEntry(0,i,y[i]);

			BlackboxBlockContainer<Field,Transpose<Blackbox> > Sequence (&A,_field,U,V);
			BlockMasseyDomain <Field,BlackboxBlockContainer<Field,Transpose<Blackbox> > > MBD(&Sequence);

			std::vector<Block> minpoly;
			std::vector<size_t> degree;
			MBD.left_minpoly(minpoly,degree); 
			//MBD.printTimer();

                        //cout<<"minpoly is: \n";
                        //write_maple(_field,minpoly);
                        //cout<<endl;

			size_t idx=0;
			if ( _field.isZero(minpoly[0].getEntry(0,0))) {

				size_t i=1;
				while (i<p && _field.isZero(minpoly[0].getEntry(i,0)))
					++i;
				if (i == p)
					throw LinboxError(" block minpoly: matrix seems to be singular - abort");
				else
					idx=i	;
			}


			bool classic = true;
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
				_field.init(scaling);
				_field.neg(scaling,combi[0][0]);
				_field.invin(scaling);
				_VDF.mul(x,accu,scaling);

			}
			else {
				/*
				 * Compute the solution according to the polynomial combination
				 * given by the product of the idx-th row of MinPoly and UA.
				 * this should decrease the number of sparse apply but increase memory requirement.
				 */
				size_t deg = degree[idx];
				Block idx_poly(_field,deg+1,p-1);
				for (size_t i=0;i<deg+1;++i)
					for (size_t j=0;j<p-1;++j)
						idx_poly.setEntry(i,j,minpoly[i].getEntry(idx,j+1));

				Block Combi(_field,deg+1,m);
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
				_field.init(scaling);
				_field.neg(scaling,minpoly[0].getEntry(idx,0));
				_field.invin(scaling);
				_VDF.mul(x,accu,scaling);
			}

			return x;
		}
#endif



	}; // end of class CoppersmithSolver



}// end of namespace LinBox

#endif //__LINBOX_coppersmith_H

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

