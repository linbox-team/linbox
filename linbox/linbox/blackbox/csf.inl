/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* -*- mode: C++; style: linux -*- */
/*
 * Copyright (c) LinBox
 * Written by Bryan Youse
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

namespace LinBox
{
	template<class Field>
	template<class OutVector, class InVector>
	OutVector & CSF<Field>::apply(OutVector & y, const InVector & x) const
	{
		linbox_check((y.size()==rowdim())&&(x.size()==coldim()));

		FieldAXPY<Field> accum (_field);

		typename OutVector::iterator yp;
		typename InVector::const_iterator xp;
		PtrVector::const_iterator ip;

		/*
		if( !sorted )
			switch_sort();
			*/

		xp=x.begin();
		yp=y.begin();
		accum.reset();

		for(Index i = _ptrs[0]; (size_t)i < _ptrs.size()-1; ++i, ++yp) {
			for(Index j = _ptrs[i]; j < _ptrs[i+1]; ++j)
			{
				// y = a*x
				accum.mulacc(_vals[j], *(xp + _inds[j]) );
			}
			accum.get(*yp);
			accum.reset();
		}

		return y;
	}

	/*
	template<class Field>
	template<class OutVector, class InVector>
	OutVector & CSF<Field>::apply(OutVector & y, const InVector & x) const
	//OutVector & CSF<Field>::apply(OutVector & y, const InVector & x)
	{
		//std::cout << endl;
		//std::cout << " new apply pass " << endl;
		//std::cout << this->rowdim() << " " << this->coldim() << endl;
		//std::cout << " inside apply " << endl;
		//std::cout << x.size() << " " << y.size() << endl;

		linbox_check((y.size()==rowdim())&&(x.size()==coldim()));

		FieldAXPY<Field> accum (_field);

		typename OutVector::iterator yp;
		typename InVector::const_iterator xp;
		PtrVector::const_iterator ip;
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
		accum.reset();

		//std::cout << " before for loop " << endl;
		//std::cout << _ptrs.end() - _ptrs.begin() << " ";
		//std::cout << _inds.end() - _inds.begin() << endl;

		//std::cout << " stupid test here " << (*_ptrs.begin()) - _inds.begin() << endl;
		//cout << &_inds << "  of size ";
		//cout << sizeof(_inds) << " ";
		//cout  << &(*_inds.begin()) << endl;

		for(ip = _ptrs.begin(); ip !=_ptrs.end()-1; ++ip, ++yp)
			//for(ip = _ptrs.begin(); ip < _ptrs.end()-2; ++ip, ++yp)  // zigzag way
		{
			//std::cout << " inside the outer for loop " << ip - _ptrs.begin() << endl;
			for(jp = *ip; jp !=*(ip + 1); ++jp)
			{
				//std::cout << jp - _inds.begin() << endl;
				accum.accumulate_special( *(xp + *jp) );
			}
			//std::cout << " accumulate is done for one iteration " << endl;
			//std::cout << " before accum.get " << yp - y.begin() << endl;
			accum.get(*yp);
			//std::cout << " before accum.reset " << endl;
			accum.reset();

			// /// zigzag way
			   ++ip;++yp;

			   for(jp = *(ip + 1); jp >*ip; --jp)
			   accum.accumulate_special( *(xp + *(jp-1)) );
			   accum.get(*yp);
			   accum.reset();


		}

		return y;
	}
*/


	/* if you want to keep two copies for the matrix, one of which is sorted by row,
	 * the other by column, then you want to use this applyTranspose function. In
	 * this case, un-comment this one and comment out the applyTranspose further down
	 */
	/*
	   template<class Field>
	   template<class OutVector, class InVector>
	   OutVector & CSF<Field>::applyTranspose(OutVector & y, const InVector & x) const
	   {
	   linbox_check((y.size()==coldim())&&(x.size()==rowdim()));

	   FieldAXPY<Field> accum (_field);

	   typename OutVector::iterator yp;
	   typename InVector::const_iterator xp;
	   PtrVector::const_iterator ip;
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

#if 0
	template<class Field>
	template<class OutVector, class InVector>
	OutVector & CSF<Field>::applyTranspose(OutVector & y, const InVector & x) const
	{
		linbox_check((y.size()==coldim())&&(x.size()==rowdim()));

		FieldAXPY<Field> accum (_field);

		typename OutVector::iterator yp;
		typename InVector::const_iterator xp;
		PtrVector::const_iterator ip;
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

		for(ip = _ptrs.begin(); ip < _ptrs.end()-1; ++ip, ++yp)
		{
			for(jp = *ip; jp <*(ip + 1); ++jp)
				accum.accumulate_special( *(xp + *jp) );
			accum.get(*yp);
			accum.reset();
		}
		return y;
	}
#endif
}//End of LinBox
