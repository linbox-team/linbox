/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* linbox/blackbox/nag-sparse.h
 * Copyright (C) 2002 Rich Seagraves
 *
 * Written by Rich Seagraves <seagrave@cis.udel.edu>
 * Modified by Zhendong Wan, -bds
 * ------------------------------------
 *
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 *.
 */

#ifndef __LINBOX_bb_zoi_INL
#define __LINBOX_bb_zoi_INL

#if 0
#ifndef _ZOI_OUTERLOOP_CHUNK_SIZE
#include "test-zo-openmp-chunksize.h" // BB: where is this file ?
#endif
#endif

namespace LinBox
{

	/*
	// for GF(310), apply to packvec will do adds with periodic normalizations
	template<>
	template<Field>
	Packvec<typename Field::Accumulator>& ZeroOne<Field>::apply<Packvec<typename Field::Accumulator>, Packvec<typename Field::Accumulator> >
	(Pacvec<typename Field::Accumulator>& y,
	Pacvec<typename Field::Accumulator>& x)
	{	return y; }

	// more stuff for paley graphs
	template<class OutVector, class InVector>
	template<>
	OutVector & ZeroOne<Special3_10Field>::apply(OutVector & y, const InVector & x) const
	{
	typedef Special3_10Field Field;
	typedef typename Packvec<Field::Accumulator> Vector;
	Vector xx, yy;
	field().convert(xx, x);
	apply(yy, xx);
	field().convert(y, yy);
	}
	*/

	template<class Field>
	template<class OutVector, class InVector>
	OutVector & ZeroOne<Field>::apply(OutVector & y, const InVector & x) const
	//OutVector & ZeroOne<Field>::apply(OutVector & y, const InVector & x)
	{

		//std::cout << endl;
		//std::cout << " new apply pass " << endl;
		//std::cout << this->rowdim() << " " << this->coldim() << endl;
		//std::cout << " inside apply " << endl;
		//std::cout << x.size() << " " << y.size() << endl;

		linbox_check((y.size()==rowdim())&&(x.size()==coldim()));

		typename OutVector::iterator yp;
		typename InVector::const_iterator xp;
		//PointerVector::const_iterator ip;
		int i; // for index version of loop suitable for openMP
		int ilimit; // for index version of loop suitable for openMP
		IndexVector::const_iterator jp;


		//std::cout << " before sorting " << endl;

		if( !sorted )
		{
			//std::cout << " -- apply: switch sort (" << _rowdim << ", " << _coldim << ", " << nnz() << ", " << sorted << " " << &(this->sorted) << ") -- " << std::endl;
			switch_sort();
			//std::cout << " -- apply: switch sort (" << _rowdim << ", " << _coldim << ", " << nnz() << ", " << sorted << ") -- " << std::endl;
		}
		//std::cout << " -- zo apply: " << this << "(" << _rowdim << ", " << _coldim << ", " << nnz() << ")" << std::endl;

		xp=x.begin();
		yp=y.begin();
		//    accum.reset();


		// outer loop to parallelize

		//  Timer timer;
		//  timer.start();
		ilimit = _indexP.size() -1;
		//  std::cout << std::endl << "Using ilimit: " << ilimit << std::endl;


#ifdef _RUNOPENMP
		//     std:: cout<< "_ZOI_OUTERLOOP_CHUNK_SIZE is:" << _ZOI_OUTERLOOP_CHUNK_SIZE << std:: endl;
		//  std::cout<<"number of threads is " << omp_get_num_threads() << std ::endl;
#pragma omp parallel
		{
#pragma omp for firstprivate(xp, ilimit) private(jp, i) schedule(static, _ZOI_OUTERLOOP_CHUNK_SIZE)
#endif

			for(i = 0; i < ilimit; ++i)

				//for(i = 0; i < _indexP.size()-1; ++i)
				//for(ip = _indexP.begin(); ip !=_indexP.end()-1; ++ip, ++yp)
			{

				FieldAXPY<Field> accum (_field);

#ifdef _RUNOPENMP
				//#pragma omp critical
				//	{
				//	  if (i%5000==0)
				//	    std::cout<< ":i/thread no.:" << i << ":" << omp_get_thread_num() << std::endl;
				//	}
#endif

				for(jp = _indexP[i]; jp !=_indexP[i+1]; ++jp)
				{
					//std::cout << jp - _index.begin() << endl;

					accum.accumulate_special( *(xp + *jp) );

				}

				//std::cout << " accumulate is done for one iteration " << endl;
				//std::cout << " before accum.get " << yp - y.begin() << endl;
				accum.get(y[i]);
				//std::cout << " before accum.reset " << endl;
				//	        accum.reset();
			}
#ifdef _RUNOPENMP
		}
#endif

		// timer.stop();

		// std::cout << "time:" << timer << std::endl;
		// std::cout << "user/sys/real/threads/chunksize:" << timer.usertime() << "/" <<
		// timer.systime() << "/" << timer.realtime() <<
		// "/" << omp_get_max_threads() << "/" << _ZOI_OUTERLOOP_CHUNK_SIZE << std::endl;

		return y;
	}

	/* if you want to keep two copies for the matrix, one of which is sorted by row,
	 * the other by column, then you want to use this applyTranspose function. In
	 * this case, un-comment this one and comment out the applyTranspose further down
	 */
	/*
	   template<class Field>
	   template<class OutVector, class InVector>
	   OutVector & ZeroOne<Field>::applyTranspose(OutVector & y, const InVector & x) const
	   {
	   linbox_check((y.size()==coldim())&&(x.size()==rowdim()));

	   FieldAXPY<Field> accum (_field);

	   typename OutVector::iterator yp;
	   typename InVector::const_iterator xp;
	   PointerVector::const_iterator ip;
	   IndexVector::const_iterator jp;

	   xp=x.begin();
	   yp=y.begin();
	   accum.reset();

	   for(ip = _colP.begin(); ip < _colP.end()-1; ++ip, ++yp)
	   {
	   for(jp = *ip; jp <*(ip + 1); ++jp)
	   accum.accumulate_special( *(xp + *jp) );
	   accum.get(*yp);
	   accum.reset();
	   }

	   return y;

	   }
	   */

	template<class Field>
	template<class OutVector, class InVector>
	OutVector & ZeroOne<Field>::applyTranspose(OutVector & y, const InVector & x) const
	{
		linbox_check((y.size()==coldim())&&(x.size()==rowdim()));

		FieldAXPY<Field> accum (_field);

		typename OutVector::iterator yp;
		typename InVector::const_iterator xp;
		PointerVector::const_iterator ip;
		IndexVector::const_iterator jp;

		xp=x.begin();
		yp=y.begin();
		accum.reset();

		if( sorted )
			//{
			//std::cout << " -- in apply transpose, before switch sort -- " << std::endl;
			switch_sort();
		//std::cout << " -- in apply transpose, after switch sort -- " << std::endl;
		//}

		for(ip = _indexP.begin(); ip < _indexP.end()-1; ++ip, ++yp)
		{
			for(jp = *ip; jp <*(ip + 1); ++jp)
				accum.accumulate_special( *(xp + *jp) );
			accum.get(*yp);
			accum.reset();
		}
		return y;
	}

}//End of LinBox

//#undef _ZOI_OUTERLOOP_CHUNK_SIZE

#endif //__LINBOX_bb_zoi_INL

