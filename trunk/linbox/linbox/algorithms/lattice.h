/* linbox/algorithms/lattice.h
 * Copyright (C) 2011 The LinBox group
 * Written by Brice Boyer <bboyer@imag.fr>
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
#ifndef __LINBOX_algorithms_lattice_H
#define __LINBOX_algorithms_lattice_H


/*! @file algorithms/lattice.h
 * @brief  LLL reduction
 * @ingroup algorithms
 * @ingroup lattice
 *
 * This is an interface to NTL/FPLLL.
 *
 * @todo Create a BlasMatrix<NTL_ZZ> that is just like a mat_ZZ !
 * @todo Create a BlasMatrix<FPLLL_ZZ> that is just like a IntMatrix !
 * @todo This will avoid copy back/forth a BlasMatrix<PID_integer>
 */

#if !defined(__LINBOX_HAVE_FPLLL) && !defined(__LINBOX_HAVE_NTL)
#error "you need either FPLLL or NTL here"
#endif


#ifdef __LINBOX_HAVE_NTL
#include <NTL/LLL.h>
#endif



#ifdef __LINBOX_HAVE_FPLLL
// this is a damn FPLLL bug !!!
#define round
#define trunc
#include <fplll/fplll.h>
#include <fplll/heuristic.h>
#include <fplll/proved.h>
#include <fplll/wrapper.h>
#undef round
#undef trunc

#endif


namespace LinBox
{ /*  Methods */

	/*! NTL methods.
	 *
	 * This lists the methods implemented.
	 */
	class latticeMethod {
	public:
		struct genericMethod {};
#ifdef __LINBOX_HAVE_NTL
		/*! NTL_LLL.
		 * This is NTL's LLL
		 * The Defaults are NTL's
		 */
		struct latticeNTL_LLL : public virtual genericMethod {
		public :
			enum  localMeth { FP , XD , QP, RR } ;
		private :
			double _delta ;
			long   _deep ;
			// LLLCheckFct
			enum  localMeth  _met ;
			bool   _givens ; // true for Givens
		public :
			latticeNTL_LLL() :
				_delta(0.99),
				_deep(0),
				_met(XD),
				_givens(false)
			{}
			void setDelta(const double & delta)
			{
				_delta = delta ;
			}
			void setDepth(const long & deep)
			{
				_deep = deep ;
			}
			void setMet( enum localMeth met )
			{
				_met = met ;
			}
			void setGivens (bool giv )
			{
				_givens = giv ;
			}
			double getDelta() const
			{
				return _delta ;
			}
			long getDepth() const
			{
				return _deep ;
			}
			enum localMeth getMeth() const
			{
				return _met ;
			}
			bool givens() const
			{
				return _givens ;
			}

		};

		/*! NTL_BKZ.
		 * This is NTL's BKZ.
		 * The Defaults are NTL's
		 */
		struct latticeNTL_BKZ : public virtual genericMethod {
			enum localMeth { FP , XD , QP, RR } ;
		private :
			double  _delta ;
			long    _bsize ;
			long    _prune ;
			// LLLCheckFct
			enum  localMeth _met ;
			bool   _givens ;
		public :
			latticeNTL_BKZ() :
				_delta(0.99),
				_bsize(10),
				_prune(0),
				_met(XD),
				_givens(false)
			{}
			void setDelta(const double & delta)
			{
				_delta = delta ;
			}
			void setBlockSize(const long & bsize)
			{
				_bsize = bsize ;
			}
			void setPrune( const long & prune)
			{
				_prune = prune ;
			}
			void setMet( enum localMeth met )
			{
				_met = met ;
			}
			void setGivens (bool giv )
			{
				_givens = giv ;
			}
			double getDelta() const
			{
				return _delta ;
			}
			long getBlockSize() const
			{
				return _bsize ;
			}
			long getPrune() const
			{
				return _prune ;
			}
			enum localMeth getMeth() const
			{
				return _met ;
			}
			bool givens() const
			{
				return _givens ;
			}

		};
#endif // __LINBOX_HAVE_NTL

#ifdef __LINBOX_HAVE_FPLLL

		/*! FPLLL LLL.
		 * Wrapper to fplll
		 * Babai,GetBase not implemented,
		 * so class is not stored there yet.
		 */
		// template<class ZT, class FT>
		struct latticeFPLLL : public virtual genericMethod {
			enum localMeth {
				P    //!< proved
				, H  //!< heuristic
				, W  //!< wrapper
			} ;
		private :
			int     _pres ;
			double   _eta ;
			double _delta ;
			int   _siegel ;
			enum localMeth _met;
			// proved<ZT,FT>    * _meth1 ;
			// heuristic<ZT,FT> * _meth2 ;
			// wrapper          * _meth3 ;
		public :
			latticeFPLLL() :
				_pres (0),
				_eta(0.51),
				_delta(0.99),
				_siegel(0),
				_met(P)
			{}
			void setPrecision(int pres)
			{
				_pres = pres ;
			}
			void setEta(double eta)
			{
				_eta = eta ;
			}
			void setDelta(double delta)
			{
				_delta = delta ;
			}
			void setSiegel(int siegel)
			{
				_siegel = siegel ;
			}
			void setMeth(enum localMeth met)
			{
				_met = met ;
			}
			int getPrecision() const
			{
				return _pres;
			}
			double getEta() const
			{
				return _eta;
			}
			double getDelta() const
			{
				return _delta;
			}
			int getSiegel() const
			{
				return _siegel;
			}
			enum localMeth getMeth() const
			{
				return _met;
			}
		};
#endif // __LINBOX_HAVE_FPLLL


	};

} // LinBox

#include "linbox/algorithms/lattice.inl"


#ifdef __LINBOX_HAVE_FPLLL
#define defaultLllMeth latticeMethod::latticeFPLLL
#else
#define defaultLllMeth latticeMethod::latticeNTL_LLL
#endif

namespace LinBox
{
	template<class Ring, class myMethod>
	void lllReduceIn(BlasMatrix<Ring> & H,
			 const myMethod   & meth = defaultLllMeth())
	{
		Ring Z ;
		BlasMatrix<Ring> U(Z,0,0);
		const bool withU = false;
		lllReduceInBase<Ring,withU>(H,U,meth);

	}

	template<class Ring, class myMethod>
	void lllReduceIn(BlasMatrix<Ring> & H,
			 BlasMatrix<Ring> & U,
			 const myMethod   & meth = defaultLllMeth())
	{
		const bool withU = true;
		lllReduceInBase<Ring,withU>(H,U,meth);
	}

	template<class Ring, class myMethod>
	void lllReduce(BlasMatrix<Ring>       & H,
		       const BlasMatrix<Ring> & A,
		       const myMethod         & meth = defaultLllMeth())
	{
		H = A ;
		lllReduceIn<Ring>(H,meth);
	}

	template<class Ring, class myMethod>
	void lllReduce(BlasMatrix<Ring>       & H,
			 BlasMatrix<Ring>     & U,
		       const BlasMatrix<Ring> & A,
			 const myMethod       & meth = defaultLllMeth())
	{
		H = A ;
		lllReduceIn<Ring>(H,U,meth);
	}

}

#undef defaultLllMeth

#endif // __LINBOX_algorithms_lattice_H

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

