/* linbox/blackbox/block-hankel-inverse.h
 * Copyright (C) 2005 Pascal Giorgi
 *
 * Written by Pascal Giorgi pgiorgi@uwaterlo.ca
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


#ifndef __LINBOX_block_hankel_inverse_H
#define __LINBOX_block_hankel_inverse_H

//#define  _BM_TIMING

//#define __CHECK_SIGMA_BASIS


#include "linbox/matrix/polynomial-matrix.h"
#include "linbox/algorithms/polynomial-matrix/order-basis.h"
#include "linbox/blackbox/block-hankel.h"
#include "linbox/blackbox/compose.h"
#include "linbox/vector/vector-domain.h"

//#define PADEMATRIX



namespace LinBox
{


	template<class _Field>
	class BlockHankelInverse {
	public:
		typedef _Field                    Field;
		typedef typename Field::Element Element;
		typedef BlasMatrix<Field> Coefficient;
        typedef PolynomialMatrix<Field,PMType::matfirst> PMatrix;

	private:
		const Field           *_field;
		VectorDomain<Field>    _VD;
		BlockHankel<Field>    *_H1;
		BlockHankel<Field>    *_T1;
		BlockHankel<Field>    *_H2;
		BlockHankel<Field>    *_T2;
		size_t _row, _col;
		BlasMatrixDomain<Field> _BMD;
		size_t   _numblock;
		size_t      _block;
	public:

		// Constructor from a stl vector of BlasMatrix representing
		// all different elements in the Hankel representation
		// vector is of odd size and represent the 1st column and last row append together
		BlockHankelInverse(const Field &F, const std::vector<BlasMatrix<Field> > &P) :
			_field(&F), _VD(F), _BMD(F)
		{
			//write_maple("UAV",P);

			size_t block    = P[0].rowdim();
			size_t deg      = P.size();
			size_t rowblock = (deg+1)>>1;
			size_t colblock = rowblock;
			_row = _col = colblock*block;
			Element one=field().one;
			_numblock=rowblock;
			_block= block;


			// compute right and left matrix Pade Approximant of P of order deg+2 and deg


			// construct the matrix power series
                        // LeftPowerSerie = [ P(x)^T  I ]^T
                        // RighPowerSerie = [ P(x)   I ]^T --> this should be transposed but we use left order basis instead of right one

                        PMatrix LeftPowerSerie (field(),2*block,block,deg+2);
                        PMatrix RightPowerSerie(field(),block,2*block,deg+2);
			for (size_t i=0;i< deg+2; ++i){
				if (i <deg) {
					for (size_t j=0;j<block;++j)
						for (size_t k=0;k<block;++k){
				                        field().assign(LeftPowerSerie.ref(j,k,i),  P[i].getEntry(j,k));
							field().assign(RightPowerSerie.ref(k,j,i), P[i].getEntry(j,k));
                                                }
				}
			}
			for (size_t j=0;j<block;++j){
                                field().assign(LeftPowerSerie.ref (block+j, j ,0), one);
                                field().assign(RightPowerSerie.ref(block+j, j, 0), one);
			}

			OrderBasis<Field> OB(F);

                        PMatrix LP1(field(),2*block,2*block,deg), RP1(field(),2*block,2*block,deg);
                        PMatrix LP2(field(),2*block,2*block,deg+2), RP2(field(),2*block,2*block,deg+2);

			size_t two_n= block<<1;
			std::vector<size_t> dlp1(two_n,0), drp1(two_n,0), dlp2(two_n,0), drp2(two_n,0);

			for (size_t i=block;i<two_n;++i){
				dlp1[i]=1;
				drp1[i]=1;
				dlp2[i]=1;
				drp2[i]=1;
			}

			// Compute the sigma basis
			//Timer chrono;
			//chrono.start();
                        OB.PM_Basis(LP1, LeftPowerSerie, deg-1, dlp1);
                        OB.PM_Basis(LP2, LeftPowerSerie, deg+1, dlp2); // MUST BE OPTIMIZED -> modify power serie thanks to LP1 and compute small basis
                        OB.PM_Basis(RP1, RightPowerSerie, deg-1, dlp1);
                        OB.PM_Basis(RP2, RightPowerSerie, deg+1, dlp2); // MUST BE OPTIMIZED -> modify power serie thanks to RP1 and compute small basis


			std::vector<BlasMatrix<Field> > SLP1, SLP2, SRP1, SRP2;
			extractLeftSigma  (SLP1, LP1, dlp1, block);
			extractLeftSigma  (SLP2, LP2, dlp2, block);
			extractTransposedRightSigma (SRP1, RP1, drp1, block);
			extractTransposedRightSigma (SRP2, RP2, drp2, block);


			BlasMatrix<Field> Res(field(),block,block), Inv(field(),block, block);

			int singular;
			// Normalization of SLP2 (V*)
			_BMD.inv(Inv, SLP2[0],singular);
			if (singular)
				throw LinboxError("BLock Hankel Inversion failed\n;");
			for (size_t i=0;i<SLP2.size();++i){
				Res = SLP2[i];
				_BMD.mul(SLP2[i], Inv, Res);
			}

			// Normalization of SRP2 (V)
			_BMD.inv(Inv, SRP2[0],singular);
			if (singular)
				throw LinboxError("BLock Hankel Inversion failed\n;");
			for (size_t i=0;i<SRP2.size();++i){
				Res = SRP2[i];
				_BMD.mul(SRP2[i], Res, Inv);
			}


			// Normalization of SLP1 (Q*)
			//_BMD.mul(Res, SLP1[0], LeftPowerSerie[deg-1]);
			_BMD.mul(Res, SLP1[0], P[deg-1]);
			for (size_t i=1;i<SLP1.size(); ++i)
				//_BMD.axpyin(Res, SLP1[i], LeftPowerSerie[deg-1-i]);
				_BMD.axpyin(Res, SLP1[i], P[deg-1-i]);

			_BMD.inv(Inv, Res,singular);
			if (singular)
				throw LinboxError("BLock Hankel Inversion failed\n;");
			for (size_t i=0;i<SLP1.size();++i){
				Res=SLP1[i];
				_BMD.mul(SLP1[i], Inv, Res);
			}


			// Normalization of SRP1 (Q)
			//_BMD.mul(Res, RightPowerSerie[deg-1], SRP1[0]);
			_BMD.mul(Res, P[deg-1], SRP1[0]);
			for (size_t i=1;i<SRP1.size(); ++i)
				//_BMD.axpyin(Res, RightPowerSerie[deg-1-i], SRP1[i]);
				_BMD.axpyin(Res, P[deg-1-i], SRP1[i]);

			_BMD.inv(Inv, Res, singular);
			if (singular)
				throw LinboxError("BLock Hankel Inversion failed\n;");
			for (size_t i=0;i<SRP1.size();++i){
				Res=SRP1[i];
				_BMD.mul(SRP1[i], Res, Inv);
			}



			/*
			   write_maple("QStar",SLP1);
			   write_maple("Q",SRP1);
			   write_maple("VStar",SLP2);
			   write_maple("V",SRP2);
			   */

			Coefficient Zero2(field());
			std::vector<Coefficient> rev_poly(SRP2.size(),Zero2);
			for (size_t i=0;i<SRP2.size();++i)
				rev_poly[i]= SRP2[SRP2.size()-1-i];

			rev_poly.erase(rev_poly.begin());


			_H1 = new BlockHankel<Field>  (field(), rev_poly, BlockHankelTag::up); // V

			rev_poly.resize(SRP2.size());
			const BlasMatrix<Field> Zero(field(),block,block);
			for (size_t i=0;i<SRP1.size();++i)
				rev_poly[i]= SRP1[SRP1.size()-1-i];
			rev_poly[SRP1.size()]= Zero;
			rev_poly.erase(rev_poly.begin());



			_H2 = new BlockHankel<Field>  (field(), rev_poly, BlockHankelTag::up); // Q

			_T1 = new BlockHankel<Field>  (field(), SLP1, BlockHankelTag::up); // Qstar

			SLP2.erase(SLP2.begin());
			_T2 = new BlockHankel<Field>  (field(), SLP2, BlockHankelTag::up); // Vstar

		}

		//!@bug copy constructor ??

		~BlockHankelInverse()
		{
			delete _H1;
			delete _T1;
			delete _H2;
			delete _T2;

		}


		//template<class Vector1, class Vector2>
		std::vector<Element>& apply (std::vector<Element> &x, const std::vector<Element> &y) const
		{

			std::vector<Element> z1(_row), z2(_row);
			std::vector<Element> rev_y(_row);
			// reverse y according to block structure
			for (size_t i=0; i< _numblock; ++i)
				for (size_t j=0;j<_block;++j){
					field().assign(rev_y[(_numblock-i-1)*_block+j], y[i*_block+j]);
				}

			_T1->apply(z1, rev_y);
			_H1->apply(z2, z1);
			_T2->apply(z1, rev_y);
			_H2->apply(rev_y, z1) ;
			_VD.sub(x, z2, rev_y);

			return x;
		}

		template<class Vector1, class Vector2>
		Vector1& applyTranspose (Vector1 &x, const Vector2 &y) const
		{
			// not handled yet
			return x;
		}

		size_t rowdim() const { return _row;}

		size_t coldim() const { return _col;}

		const Field& field() const { return *_field;}

	protected:

		void extractLeftSigma(std::vector<BlasMatrix<Field> >                &S,
				      PMatrix                                &SigmaBase,
				      std::vector<size_t>                       &defect,
				      size_t                                      block) const
		{

			// take the block rows which have lowest defect
			// compute permutation such that first block rows have lowest defect
			std::vector<size_t> Perm(2*block);
			for (size_t i=0;i<2*block;++i)
				Perm[i]=i;
			for (size_t i=0;i<2*block;++i) {
				size_t idx_min=i;
				for (size_t j=i+1;j<2*block;++j)
					if (defect[j]< defect[idx_min])
						idx_min=j;
				std::swap(defect[i],defect[idx_min]);
				Perm[i]=idx_min;
			}
			BlasPermutation<size_t>  BPerm(Perm);

			// Apply BPerm to the Sigma Base
			for (size_t i=0;i<SigmaBase.size();++i)
				_BMD.mulin_right(BPerm,SigmaBase[i]);

			size_t max=defect[0];
			for (size_t i=0;i<block;++i)
				if (defect[i] > max)
					//max=defect[notnull[i]];
					max=defect[i];

			// prepare S to receive the sigma base
			const BlasMatrix<Field> Zero(field(),block,block);
			S.resize(max+1, Zero);

			// extract the sigma base
			for (size_t k=0;k<S.size();++k){
				for(size_t i=0;i<block;++i)
					for (size_t j=0;j<block;++j)
						S[k].setEntry(i,j, SigmaBase[k].getEntry(i,j));
			}
		}


		void extractTransposedRightSigma(std::vector<BlasMatrix<Field> >            &S,
                                                 std::vector<BlasMatrix<Field> >    &SigmaBase,
                                                 std::vector<size_t>                     &defect,
                                                 size_t                                    block) const
		{

			// take the m rows which have lowest defect
			// compute permutation such that first block rows have lowest defect
			std::vector<size_t> Perm(2*block);
			for (size_t i=0;i<2*block;++i)
				Perm[i]=i;
			for (size_t i=0;i<2*block;++i) {
				size_t idx_min=i;
				for (size_t j=i+1;j<2*block;++j)
					if (defect[j]< defect[idx_min])
						idx_min=j;
				std::swap(defect[i],defect[idx_min]);
				Perm[i]=idx_min;
			}
			BlasPermutation<size_t>  BPerm(Perm);

			// Apply BPerm to the Sigma Base
			for (size_t i=0;i<SigmaBase.size();++i)
				_BMD.mulin_right(BPerm,SigmaBase[i]);

			size_t max=defect[0];
			for (size_t i=0;i<block;++i)
				if (defect[i] > max)
					max=defect[i];

			// prepare S to receive the sigma base
			const BlasMatrix<Field> Zero(field(),block,block);
			S.resize(max+1, Zero);

			// extract the sigma base
			for (size_t k=0;k<S.size();++k){
				for(size_t i=0;i<block;++i)
					for (size_t j=0;j<block;++j)
						//S[k].setEntry(i,j, SigmaBase[k].getEntry(notnull[i],j));
						S[k].setEntry(i,j, SigmaBase[k].getEntry(i,j));
			}

			// transpose the sigma base since it comes from a left sigma basis
			for (size_t k=0;k<S.size();++k){
				for (size_t i=0;i<block;++i)
					for (size_t j=i;j<block;++j)
						std::swap(S[k].refEntry(i,j), S[k].refEntry(j,i));
			}
		}


	protected:

		void write_maple(const char* name, const std::vector<Coefficient> & P)
		{
			size_t m,n;
			m = P[0].rowdim();
			n = P[0].coldim();
			std::cout<<name<<":=[";
			for (size_t k=0;k<P.size()-1;++k){
				std::cout<<"Matrix([";
				for (size_t i=0;i<m-1;++i){
					std::cout<<"[";
					for (size_t j=0;j<n-1;++j)
						field().write(std::cout,P[k].getEntry(i,j))<<",";
					field().write(std::cout, P[k].getEntry(i,n-1))<<"] , ";
				}
				std::cout<<"[";
				for (size_t j=0;j<n-1;++j)
					field().write(std::cout,P[k].getEntry(m-1,j))<<",";
				field().write(std::cout, P[k].getEntry(m-1,n-1))<<"]]) , ";
			}

			std::cout<<"Matrix([";
			for (size_t i=0;i<m-1;++i){
				std::cout<<"[";
				for (size_t j=0;j<n-1;++j)
					field().write(std::cout,P[P.size()-1].getEntry(i,j))<<",";
				field().write(std::cout, P[P.size()-1].getEntry(i,n-1))<<"] , ";
			}
			std::cout<<"[";
			for (size_t j=0;j<n-1;++j)
				field().write(std::cout,P[P.size()-1].getEntry(m-1,j))<<",";
			field().write(std::cout, P[P.size()-1].getEntry(m-1,n-1))<<"]])]; \n";
		}


	}; //end of class BlockHankelInverse



}// end of namespace LinBox


#endif //__LINBOX_block_hankel_inverse_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
