/* linbox/algorithms/lattice.inl
 * Copyright (C) 2011 The LinBox group
 * Written by Brice Boyer (briceboyer) <boyer.brice@gmail.com>
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
#ifndef __LINBOX_algorithms_lattice_INL
#define __LINBOX_algorithms_lattice_INL


/*!@internal
 * @file algorithms/lattice.inl
 * @brief  LLL reduction
 * @ingroup algorithms
 * @ingroup lattice
 *
 * Implements the various wrappers
 */

/*  Interface to NTL LLL */
namespace LinBox
{
#ifdef __LINBOX_HAVE_NTL

	//! @todo we should use mat_ZZ here instead of BlasMatrix<Ring>
	template<class Ring, bool withU>
	void
	lllReduceInBase(BlasMatrix<Ring>                      & H,
			BlasMatrix<Ring>                      & UU,
			const latticeMethod::latticeNTL_LLL   & meth)
	{
		// convert to mat_ZZ
		NTL::mat_ZZ B ;
		NTL::mat_ZZ U ;
		B.SetDims(H.rowdim(),H.coldim());
		NTL::clear(B);
		// XXX Maybe Bij, maybe Bji...
		NTL_ZZ Ints ;
		NTL_ZZ::Element Bij ;
		for (size_t i = 0 ; i < H.rowdim(); ++i) {
			for (size_t j = 0  ; j < H.coldim() ; ++j) {
				Ints.init(Bij,H.getEntry(i,j));
				B[(long)i][(long)j] = Bij ;
			}
		}
		if (withU) {
			U.SetDims(UU.rowdim(),UU.coldim());
			NTL::clear(U);
		}
		// do NTL's LLL
		if (withU) {
			switch (meth.getMeth()) {
			case (latticeMethod::latticeNTL_LLL::FP) :
				{
					if (meth.givens()) {
						G_LLL_FP(B,U,meth.getDelta(),meth.getDepth());
					}
					else {
						LLL_FP(B,U,meth.getDelta(),meth.getDepth());
					}
				}
				break;
			case (latticeMethod::latticeNTL_LLL::RR) :
				{
					if (meth.givens()) {
						G_LLL_RR(B,U,meth.getDelta(),meth.getDepth());
					}
					else {
						LLL_RR(B,U,meth.getDelta(),meth.getDepth());
					}
				}
				break;
			case (latticeMethod::latticeNTL_LLL::XD) :
				{
					if (meth.givens()) {
						G_LLL_XD(B,U,meth.getDelta(),meth.getDepth());
					}
					else {
						LLL_XD(B,U,meth.getDelta(),meth.getDepth());
					}
				}
				break;
			case (latticeMethod::latticeNTL_LLL::QP) :
				{
					if (meth.givens()) {
						G_LLL_QP(B,U,meth.getDelta(),meth.getDepth());
					}
					else {
						LLL_QP(B,U,meth.getDelta(),meth.getDepth());
					}
				}
				break;

			}
		}
		else {
			switch (meth.getMeth()) {
			case (latticeMethod::latticeNTL_LLL::FP) :
				{
					if (meth.givens()) {
						G_LLL_FP(B,meth.getDelta(),meth.getDepth());
					}
					else {
						LLL_FP(B,meth.getDelta(),meth.getDepth());
					}
				}
				break;
			case (latticeMethod::latticeNTL_LLL::RR) :
				{
					if (meth.givens()) {
						G_LLL_RR(B,meth.getDelta(),meth.getDepth());
					}
					else {
						LLL_RR(B,meth.getDelta(),meth.getDepth());
					}
				}
				break;
			case (latticeMethod::latticeNTL_LLL::XD) :
				{
					if (meth.givens()) {
						G_LLL_XD(B,meth.getDelta(),meth.getDepth());
					}
					else {
						LLL_XD(B,meth.getDelta(),meth.getDepth());
					}
				}
				break;
			case (latticeMethod::latticeNTL_LLL::QP) :
				{
					if (meth.givens()) {
						G_LLL_QP(B,meth.getDelta(),meth.getDepth());
					}
					else {
						LLL_QP(B,meth.getDelta(),meth.getDepth());
					}
				}
				break;

			}
		}
		// convert from mat_ZZ
		Integer Hij;
		for (size_t i = 0 ; i < H.rowdim(); ++i) {
			for (size_t j = 0  ; j < H.coldim() ; ++j) {
				Ints.convert(Hij,B[(long)i][(long)j]);
				H.setEntry( i,j,Hij );
			}
		}
		if (withU) {
			Integer Uij;
			for (size_t i = 0 ; i < H.rowdim(); ++i) {
				for (size_t j = 0  ; j < H.coldim() ; ++j) {
					Ints.convert(Uij,U[(long)i][(long)j]);
					UU.setEntry(i,j,Uij );
				}
			}
		}

	}
#endif // __LINBOX_HAVE_NTL
}

/*  Interface to NTL BKZ */
namespace LinBox
{
#ifdef __LINBOX_HAVE_NTL

	//! @todo we should use mat_ZZ here instead of BlasMatrix<Ring>
	template<class Ring, bool withU>
	void
	lllReduceInBase(BlasMatrix<Ring>                      & H,
			BlasMatrix<Ring>                      & UU,
			const latticeMethod::latticeNTL_BKZ   & meth)
	{
		// convert to mat_ZZ
		NTL::mat_ZZ B ;
		NTL::mat_ZZ U ;
		B.SetDims(H.rowdim(),H.coldim());
		NTL::clear(B);
		// XXX Maybe Bij, maybe Bji...
		NTL_ZZ Ints ;
		NTL_ZZ::Element Bij ;
		for (size_t i = 0 ; i < H.rowdim(); ++i) {
			for (size_t j = 0  ; j < H.coldim() ; ++j) {
				Ints.init(Bij, H.getEntry(i,j));
				B[(long)i][(long)j] = Bij ;
			}
		}
		if (withU) {
			U.SetDims(UU.rowdim(),UU.coldim());
			NTL::clear(U);
		}
		// do NTL's BKZ
		if (withU) {
			switch (meth.getMeth()) {
			case (latticeMethod::latticeNTL_BKZ::FP) :
				{
					if (meth.givens()) {
						G_BKZ_FP(B,U,meth.getBlockSize(),meth.getPrune());
					}
					else {
						BKZ_FP(B,U,meth.getBlockSize(),meth.getPrune());
					}
				}
				break;
			case (latticeMethod::latticeNTL_BKZ::RR) :
				{
					if (meth.givens()) {
						G_BKZ_RR(B,U,meth.getBlockSize(),meth.getPrune());
					}
					else {
						BKZ_RR(B,U,meth.getBlockSize(),meth.getPrune());
					}
				}
				break;
			case (latticeMethod::latticeNTL_BKZ::XD) :
				{
					if (meth.givens()) {
						G_BKZ_XD(B,U,meth.getBlockSize(),meth.getPrune());
					}
					else {
						BKZ_XD(B,U,meth.getBlockSize(),meth.getPrune());
					}
				}
				break;
			case (latticeMethod::latticeNTL_BKZ::QP) :
				{
					if (meth.givens()) {
						G_BKZ_QP(B,U,meth.getBlockSize(),meth.getPrune());
					}
					else {
						BKZ_QP(B,U,meth.getBlockSize(),meth.getPrune());
					}
				}
				break;

			}
		}
		else {
			switch (meth.getMeth()) {
			case (latticeMethod::latticeNTL_BKZ::FP) :
				{
					if (meth.givens()) {
						G_BKZ_FP(B,meth.getBlockSize(),meth.getPrune());
					}
					else {
						BKZ_FP(B,meth.getBlockSize(),meth.getPrune());
					}
				}
				break;
			case (latticeMethod::latticeNTL_BKZ::RR) :
				{
					if (meth.givens()) {
						G_BKZ_RR(B,meth.getBlockSize(),meth.getPrune());
					}
					else {
						BKZ_RR(B,meth.getBlockSize(),meth.getPrune());
					}
				}
				break;
			case (latticeMethod::latticeNTL_BKZ::XD) :
				{
					if (meth.givens()) {
						G_BKZ_XD(B,meth.getBlockSize(),meth.getPrune());
					}
					else {
						BKZ_XD(B,meth.getBlockSize(),meth.getPrune());
					}
				}
				break;
			case (latticeMethod::latticeNTL_BKZ::QP) :
				{
					if (meth.givens()) {
						G_BKZ_QP(B,meth.getBlockSize(),meth.getPrune());
					}
					else {
						BKZ_QP(B,meth.getBlockSize(),meth.getPrune());
					}
				}
				break;

			}
		}
		// convert from mat_ZZ
		Integer Hij;
		for (size_t i = 0 ; i < H.rowdim(); ++i) {
			for (size_t j = 0  ; j < H.coldim() ; ++j) {
				Ints.convert(Hij,B[(long)i][(long)j]);
				H.setEntry( i,j,Hij );
			}
		}
		if (withU) {
			Integer Uij;
			for (size_t i = 0 ; i < H.rowdim(); ++i) {
				for (size_t j = 0  ; j < H.coldim() ; ++j) {
					Ints.convert(Uij,U[(long)i][(long)j]);
					UU.setEntry(i,j,Uij );
				}
			}
		}

	}
#endif // __LINBOX_HAVE_NTL

}

/* Interface to FPLLL */
namespace LinBox
{
#ifdef __LINBOX_HAVE_FPLLL

#if FPLLL_MAJOR_VERSION < 5
#define lll_reduction lllReduction
#endif

	//! @bug we suppose Ring and mpz_t understand eachother...
	template<class Ring, bool withU>
	void
	lllReduceInBase(BlasMatrix<Ring>                   & H,
			BlasMatrix<Ring>                   & UU,
			const latticeMethod::latticeFPLLL  & meth)
	{
		typedef mpz_t ZT ;
		if (withU)
			throw NotImplementedYet("not U");
		// Convert H
		fplll::ZZ_mat<ZT> B(H.rowdim(),H.coldim()) ;
		for (size_t i = 0 ; i < H.rowdim() ; ++i) {
			for (size_t j = 0 ; j < H.coldim() ; ++j) {
                                fplll::Z_NR<ZT> x = H.getEntry(i,j);
                                B(i,j) = x;
			}
		}
		// LLL()
		switch (meth.getMeth()) {
		case (latticeMethod::latticeFPLLL::P) :
			{
				lll_reduction(B, meth.getDelta(), meth.getEta(), fplll::LM_PROVED, fplll::FT_DOUBLE, meth.getPrecision());
			}
			break;
		case (latticeMethod::latticeFPLLL::W) :
			{
                                // precision is *ignored*
				lll_reduction(B, meth.getDelta(), meth.getEta(), fplll::LM_WRAPPER);
			}
			break;
		case (latticeMethod::latticeFPLLL::H) :
			{
                                int flags = meth.getSiegel() ? fplll::LLL_SIEGEL : fplll::LLL_DEFAULT;
				lll_reduction(B, meth.getDelta(), meth.getEta(), fplll::LM_HEURISTIC, fplll::FT_DOUBLE, meth.getPrecision(), flags);
			}
			break;
		}



		// Convert back H
		for (size_t i = 0 ; i < H.rowdim() ; ++i) {
			for (size_t j = 0 ; j < H.coldim() ; ++j) {
				H.setEntry(i,j, B(i,j) );
			}
		}

	}
#endif // __LINBOX_HAVE_FPLLL


}

#endif // __LINBOX_algorithms_lattice_INL

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
