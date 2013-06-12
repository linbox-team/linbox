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
#include "linbox/util/commentator.h"
#include "linbox/matrix/matrix-category.h"
#include "linbox/matrix/matrix-domain.h"
#include "linbox/algorithms/blackbox-block-container.h"
#include "linbox/algorithms/block-coppersmith-domain.h"
//#include "linbox/algorithms/bm-seq.h"
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
		typedef typename MatrixDomain<_Field>::Matrix 	Block;

		inline const Field & field() const { return *_field; }
	protected:
		const Field                         *_field;
		MatrixDomain<Field>     _MD;
		VectorDomain<Field>         _VD;
		RandIter                   _rand;

	public:
		CoppersmithSolver(const Field &F) :
			_field(&F), _MD(F), _VD(F), _rand(F)
		{}

		CoppersmithSolver (const Field &F, const RandIter &rand) :
			_field(&F), _MD(F), _VD(F), _rand(rand)
		{}

		template <class Vector, class Blackbox>
		Vector &solveNonSingular (Vector &x, const Blackbox &B, const Vector &y) const
		{
			commentator().start ("Coppersmith solveNonSingular", "solveNonSingular");
			std::ostream& report = commentator().report(Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);

			//Set up the projection matrices and their dimensions
			size_t d = B.coldim();
			size_t r,c;
			integer tmp = d;

			//Set the blocking size, Using Pascal Giorgi's convention
			r=tmp.bitsize()-1;
			c=tmp.bitsize()-1;

			//Create the blocks
			Block U(field(),r,d);
			Block W(field(),d,c-1);
			Block V(field(),d,c);

			//Pick random entries for U and W. W will become the last c-1 columns of V
			for(size_t i=0; i<r;i++)
				for(size_t j=0; j<d; j++)
					_rand.random(U.refEntry(i,j));
			for(size_t i=0; i<d;i++)
				for(size_t j=0; j<c-1; j++){
					_rand.random(W.refEntry(i,j));
			}

			//Multiply W by B on the left and place it in the last c-1 columns of V
			typename MatrixDomain<Field>::Submatrix V2(V,0,1,d,c-1);
			_MD.mul(V2,B,W);

			//Make the first column of V a copy of the right side of the system, y
			for(size_t i=0; i<d; i++)
				V.setEntry(i,0,y[i]);

			//Create the sequence container and its iterator that will compute the projection
			BlackboxBlockContainer<Field, Blackbox > blockseq(&B,field(),U,V);

			//Get the generator of the projection using the Coppersmith algorithm (slightly modified by Yuhasz)
			BlockCoppersmithDomain<Field, BlackboxBlockContainer<Field, Blackbox> > BCD(&blockseq,d);
			std::vector<Block> gen;
			std::vector<size_t> deg;
			deg = BCD.right_minpoly(gen);

			//Reconstruct the solution
			//Pick a column of the generator with a nonzero element in the first row of the constant coefficient
			size_t idx = 0;
			if(field().isZero(gen[0].getEntry(0,0))){
				size_t i = 1;
				while(i<c && field().isZero(gen[0].getEntry(0,i)))
					i++;
				if(i==c)
					throw LinboxError(" block minpoly: matrix seems to be singular - abort");
				else
					idx=i;
			}

			//from 1 to the degree of the index column, multiply A^(i-1)V times the idx column of the generator coefficient x^i
			//Accumulate these results in xm
			size_t mu = deg[idx];
			Block BVo(V);
			Block BVe(field(),d,c);
			Block xm(field(),d,1);
			bool odd = true;
			for(size_t i = 1; i < mu+1; i++){
				typename MatrixDomain<Field>::Submatrix gencol(gen[i],0,idx,c,1); // BB changed d,1 to c,1
				Block BVgencol(field(),d,1);
				if(odd){
					_MD.mul(BVgencol,BVo,gencol);
					_MD.addin(xm, BVgencol);
					_MD.mul(BVe,B,BVo);
					odd=false;
				}
				else{
					_MD.mul(BVgencol,BVe,gencol);
					_MD.addin(xm, BVgencol);
					_MD.mul(BVo,B,BVe);
					odd=true;
				}

			}

			//For the constant coefficient, loop over the elements in the idx column except the first row
			//Multiply the corresponding column of W (the last c-1 columns of V before application of B) by the generator element
			//Accumulate the results in xm
			for(size_t i = 1; i < c; i++){
				typename MatrixDomain<Field>::Submatrix Wcol(W,0,i-1,d,1);
				Block Wcolgen0(field(),d,1);
				_MD.mul(Wcolgen0, Wcol, gen[0].getEntry(i,idx));
				_MD.addin(xm,Wcolgen0);
			}

			//Multiply xm by -1(move to the correct side of the equation) and divide the the 0,idx entry of the generator constant
			_MD.negin(xm);
			typename Field::Element gen0inv;
			_MD.mulin(xm, field().inv(gen0inv, gen[0].getEntry(0,idx)));

			//Test to see if the answer works with U
			Block Bxm(field(),d,1), UBxm(field(),r,1), Uycol(field(), r,1);
			typename MatrixDomain<Field>::Submatrix ycol(V,0,0,d,1);
			_MD.mul(Uycol, U, ycol);
			_MD.mul(Bxm, B, xm);
			_MD.mul(UBxm, U, Bxm);

			if(_MD.areEqual(UBxm, Uycol))
				report << "The solution matches when projected by U" << endl;
			else
				report << "The solution does not match when projected by U" << endl;

			//Copy xm into x (Change type from 1 column matrix to Vector)
			for(size_t i =0; i<d; i++)
				x[i]=xm.getEntry(i,0);

			commentator().stop ("done", NULL, "solveNonSingular");
			return x;
		}




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

