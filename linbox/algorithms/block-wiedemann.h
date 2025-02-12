/* linbox/algorithms/block-wiedemann.h
 * Copyright (C) 2004 Pascal Giorgi
 *
 * Written by Pascal Giorgi pascal.giorgi@ens-lyon.fr
 * modified by Pascal Giorgi (pascal.giorgi@lirmm.fr) 2011
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


#ifndef __LINBOX_block_wiedemann_H
#define __LINBOX_block_wiedemann_H

#include <vector>
#include <iostream>


#include "linbox/integer.h"
#include "linbox/matrix/dense-matrix.h"
#include "linbox/matrix/matrix-domain.h"
#include "linbox/algorithms/blackbox-block-container.h"
#include "linbox/algorithms/block-massey-domain.h"
#include "linbox/vector/vector-domain.h"
#include "linbox/blackbox/transpose.h"
#include "linbox/util/commentator.h"

#include "linbox/util/error.h"
//#include "linbox/util/debug.h"

namespace LinBox
{

	template <class Context_>
	class BlockWiedemannSolver{

	public:
		typedef typename Context_::Field                 Field;
		typedef typename Field::Element       Element;
		typedef typename Field::RandIter     RandIter;
		typedef BlasVector<Field>           Vector;
		typedef BlasMatrix<Field>               Block;

	protected:
		Context_                    _BMD;
		VectorDomain<Field>         _VDF;
        RandIter                   _rand;
        size_t                 _left_blockdim;
        size_t                 _right_blockdim; 


#define BW_BLOCK_DEFAULT 8UL
        //#define BW_MAX_TRY 100
#define BW_MAX_TRY 10
	public:
		const Field & field() const { return _BMD.field(); }

		BlockWiedemannSolver (const Context_ &C, size_t lblock=BW_BLOCK_DEFAULT, size_t rblock=BW_BLOCK_DEFAULT+1) :
			_BMD(C.field()), _VDF(C.field()), _rand(const_cast<Field&>(C.field())), _left_blockdim(lblock), _right_blockdim(rblock)
		{
            if (_left_blockdim ==0) _left_blockdim=BW_BLOCK_DEFAULT;
            if (_right_blockdim ==0) _right_blockdim=BW_BLOCK_DEFAULT;
        }

		BlockWiedemannSolver (const Field &F, RandIter &rand, size_t lblock=BW_BLOCK_DEFAULT, size_t rblock=BW_BLOCK_DEFAULT+1) :
			_BMD(F), _VDF(F), _rand(rand) , _left_blockdim(lblock), _right_blockdim(rblock)
		{
            if (_left_blockdim ==0) _left_blockdim=BW_BLOCK_DEFAULT;
            if (_right_blockdim ==0) _right_blockdim=BW_BLOCK_DEFAULT;
        }

		template <class Blackbox>
		Vector &solve (Vector &x, const Blackbox &B, const Vector &y) const {
            try {
                solveNonSingular(x,B,y);
            }
            catch (LinboxError& e) {
                std::cerr<<e<<std::endl;
            }
            return x;
        }

		template <class Blackbox>
		Vector &solveNonSingular (Vector &x, const Blackbox &B, const Vector &y) const
		{

			Transpose<Blackbox> A(B);

			size_t m,n,bw_try=0;
			m = A.rowdim();
			n = A.coldim();
            Vector z(field(),y.size());
                        
            if (_left_blockdim >m/2 || _right_blockdim >n/2) 
                std::cerr<<"BlockWiedemannSolver (Warning) : block size too large, number of tries might be large"<<std::endl;

            //std::cout<<"row block: "<<_left_blockdim<<std::endl;
			//std::cout<<"col block: "<<_right_blockdim<<std::endl;

			Block U(field(),_left_blockdim,m), UA(field(),_left_blockdim-1,m), V(field(),n,_right_blockdim);

            // LAS VEGAS VERSION
            do
                {
                                
                    // set V at random
                    for (size_t i=0;i<n;++i)
                        for (size_t j=0;j<_right_blockdim;++j)
                            _rand.random(V.refEntry(i,j));

                    for (size_t i=0;i<_left_blockdim-1;++i)
                        for (size_t j=0;j<m;++j)
                            _rand.random(UA.refEntry(i,j));

                    typename Block::RowIterator        iter_U  = U.rowBegin();
                    typename Block::ConstRowIterator   iter_UA = UA.rowBegin();
                    ++iter_U;
                    for (; iter_UA != UA.rowEnd(); ++iter_UA, ++iter_U)
                        A.applyTranspose( *iter_U , *iter_UA );

                    for (size_t i=0;i<m;++i)
                        U.setEntry(0,i,y[(size_t)i]);

                    BlackboxBlockContainer<Field,Transpose<Blackbox> > Sequence (&A,field(),U,V);
                    BlockMasseyDomain <Field,BlackboxBlockContainer<Field,Transpose<Blackbox> > > MBD(&Sequence);
                                        
                    std::vector<Block> minpoly;
                    std::vector<size_t> degree;
                    MBD.left_minpoly_rec(minpoly,degree);
                    //MBD.printTimer();

                    // std::cout<<"U:=";
                    // U.write(std::cout,Tag::FileFormat::Maple)<<";\n";
                    // std::cout<<"V:=";
                    // V.write(std::cout,Tag::FileFormat::Maple)<<";\n";

                    
                    // std::cout<<"minpoly is: \n";
                    // write_maple(field(),minpoly);
                    // std::cout<<std::endl;
                                        
                    size_t idx=0;
                    if ( field().isZero(minpoly[0].getEntry(0,0))) {

                        size_t i=1;
                        while (i<_left_blockdim && field().isZero(minpoly[0].getEntry(i,0)))
                            ++i;
                        if (i == _left_blockdim){
                            std::cerr<<"BW: matrix is singular \n";
                            throw LinboxError(" block minpoly: matrix seems to be singular - abort");
                        }
                        else
                            idx=i	;
                    }


                    typename Blackbox::Field F = A.field();

                    bool classic = true;
                    if ( classic) {
                        /*
                         * Compute the solution according to the polynomial combination
                         * given by each column of the idx-th row of MinPoly such that the constant term of
                         * the first element in this row is non zero.
                         * we use y and UA as projection (UA= U.A)
                         */
                        size_t deg = degree[(size_t)idx];
                        std::vector<Vector> combi(_left_blockdim,Vector(F,deg+1));
                        for (size_t i=0;i<_left_blockdim;++i)
                            for (size_t k=0;k<deg+1;++k)
                                combi[(size_t)i][k]=minpoly[k].getEntry(idx,i);
                        Vector lhs(F,n);
                        A.applyTranspose(lhs,y);
                        _VDF.mulin(lhs,combi[0][deg]);
                        Vector lhsbis(lhs);
                        for (int i = (int)deg-1 ; i > 0;--i) {
                            _VDF.axpy (lhs, combi[0][(size_t)i], y, lhsbis);
                            A.applyTranspose (lhsbis, lhs);
                        }
                        Vector accu (lhs);
                        for (size_t k=1;k<_left_blockdim;++k){
                            Vector row(F,m);
                            for (size_t j=0;j<m;++j)
                                row[j]=UA.getEntry(k-1,j);
                            A.applyTranspose(lhs,row);
                            _VDF.mulin(lhs,combi[k][deg]);
                            Vector lhsbis_loc(lhs);
                            for (int i = (int)deg-1 ; i >= 0;--i) {
                                _VDF.axpy (lhs, combi[k][(size_t)i], row, lhsbis_loc);
                                A.applyTranspose (lhsbis_loc, lhs);
                            }
                            _VDF.addin(accu,lhs);
                        }
                        Element scaling;
                        field().init(scaling);
                        field().neg(scaling,combi[0][0]);
                        field().invin(scaling);
                        _VDF.mul(x,accu,scaling);                                
                    }
                    else {
                        /*
                         * Compute the solution according to the polynomial combination
                         * given by the product of the idx-th row of MinPoly and UA.
                         * this should decrease the number of sparse apply but increase memory requirement.
                         */
                        size_t deg = degree[(size_t)idx];
                        Block idx_poly(field(),deg+1,_left_blockdim-1);
                        for (size_t i=0;i<deg+1;++i)
                            for (size_t j=0;j<_left_blockdim-1;++j)
                                idx_poly.setEntry(i,j,minpoly[(size_t)i].getEntry(idx,j+1));

                        Block Combi(field(),deg+1,m);
                        _BMD.mul(Combi,idx_poly,UA);

                        Vector lhs(F,n),row(F,m);
                        for (size_t i=0;i<m;++i)
                            row[(size_t)i]= Combi.getEntry(deg,i);

                        A.applyTranspose(lhs,row);
                        Vector lhsbis(lhs);
                        for (int i = (int)deg-1 ; i >= 0;--i) {
                            for (size_t j=0;j<m;++j)
                                row[j]= Combi.getEntry((size_t)i,j);
                            _VDF.add (lhs,row,lhsbis);
                            A.applyTranspose (lhsbis, lhs);
                        }

                        Vector accu (lhs);

                        A.applyTranspose(lhs,y);
                        _VDF.mulin(lhs,minpoly[deg].getEntry(idx,0));
                        lhsbis=lhs;
                        for (size_t i = deg-1 ; i > 0;--i) {
                            _VDF.axpy (lhs,minpoly[(size_t)i].getEntry(idx,0) , y, lhsbis);
                            A.applyTranspose (lhsbis, lhs);
                        }

                        _VDF.addin(accu,lhs);
                        Element scaling;
                        field().init(scaling);
                        field().neg(scaling,minpoly[0].getEntry(idx,0));
                        field().invin(scaling);
                        _VDF.mul(x,accu,scaling);
                    }

                    B.apply(z,x); // checking result
                    if ( bw_try>  BW_MAX_TRY ) throw LinboxError("BlockWiedemann solve: LasVegas maximum tries reached");
                    bw_try++;

                } while (!_VDF.areEqual(z,y));
#ifdef _BW_LASVEGAS_COUNT
            std::cerr<<"BlockWiedemannSolver: nbr of try: "<<bw_try<<std::endl;
#endif 
                        commentator().report()<<"BlockWiedemannSolver: nbr of tries: "<<bw_try<<std::endl;
			return x;
		}



	}; // end of class BlockWiedemannSolver



}// end of namespace LinBox

#endif //__LINBOX_block_wiedemann_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
