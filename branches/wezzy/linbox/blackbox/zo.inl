/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* -*- mode: C++; style: linux -*- */

/* linbox/blackbox/nag-sparse.h
 * Copyright (C) 2002 Rich Seagraves
 *
 * Written by Rich Seagraves <seagrave@cis.udel.edu>
 * Modified by Zhendong Wan, -bds
 * ------------------------------------
 *
 * See COPYING for license information.
 */

namespace LinBox
{

	/*
	// for GF(3^10), apply to packvec will do adds with periodic normalizations
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

		FieldAXPY<Field> accum (_field);

		typename OutVector::iterator yp;
		typename InVector::const_iterator xp;
		PointerVector::const_iterator ip;
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
		//std::cout << _indexP.end() - _indexP.begin() << " ";
		//std::cout << _index.end() - _index.begin() << endl;

		//std::cout << " stupid test here " << (*_indexP.begin()) - _index.begin() << endl;
		//cout << &_index << "  of size ";
		//cout << sizeof(_index) << " ";
		//cout  << &(*_index.begin()) << endl;

		for(ip = _indexP.begin(); ip !=_indexP.end()-1; ++ip, ++yp)
			//for(ip = _indexP.begin(); ip < _indexP.end()-2; ++ip, ++yp)  // zigzag way
		{
			//std::cout << " inside the outer for loop " << ip - _indexP.begin() << endl;
			for(jp = *ip; jp !=*(ip + 1); ++jp)
			{
				//std::cout << jp - _index.begin() << endl;
				accum.accumulate_special( *(xp + *jp) );
			}
			//std::cout << " accumulate is done for one iteration " << endl;
			//std::cout << " before accum.get " << yp - y.begin() << endl;
			accum.get(*yp);
			//std::cout << " before accum.reset " << endl;
			accum.reset();

			/* // zigzag way
			   ++ip;++yp;

			   for(jp = *(ip + 1); jp >*ip; --jp)
			   accum.accumulate_special( *(xp + *(jp-1)) );
			   accum.get(*yp);
			   accum.reset();
			   */

		}

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
