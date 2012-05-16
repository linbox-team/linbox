/* linbox/field/givaro-zpz.inl
 * Written by Pascal Giorgi <pascal.giorgi@ens-lyon.fr>
 *
 * ------------------------------------
 *
 * Copyright (c) LinBox
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
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 *.
 */

/*!@internal
 * @file linbox/field/Givaro/givaro-zpz.inl
 * @brief implementation of givaro-zpz.h.
 */


#ifndef __LINBOX_givaro_zpz_INL
#define __LINBOX_givaro_zpz_INL


#include <iostream>
#include "linbox/integer.h"

namespace LinBox
{



	template <class Vector1, class Vector2>
	inline GivaroZpz< Givaro::Std32>::Element &DotProductDomain<GivaroZpz< Givaro::Std32> >::dotSpecializedDD
	(GivaroZpz< Givaro::Std32>::Element &res, const Vector1 &v1, const Vector2 &v2) const
	{
		typedef typename GivaroZpz<Givaro::Std32>::Element Elem;
		uint64_t inter,best ;
		inter=best=0;
		if (v1.size()==0)
			return res=Elem(0);
		else {

			uint64_t size      = v1.size();
			uint64_t min       = Max / (Corr+ (uint64_t)(_field.characteristic()-1)*(uint64_t)(_field.characteristic()-1));
			uint64_t min_size  =  (size < min ? size : min);
			uint64_t good1     = (size > min_size ?  size - min_size: 0);
			uint64_t good2     = (long)(size / min_size)* min_size ;
			uint64_t good_size = (good1 > good2 ? good1 : good2 );

			typename Vector1::const_iterator i=v1.begin();
			typename Vector2::const_iterator j=v2.begin();

			unsigned long k=0;

			for (;k<min_size;i++,j++,k++)
				best+=(uint64_t)*i * (uint64_t)*j;

			for (inter=best;k<good_size;inter=best) {
				for (unsigned long l=0;l<min_size;i++,j++,k++,l++)
					best+= (uint64_t)*i * (uint64_t)*j;
				if (inter > best) best+=Corr;
			}

			for (;k<size;i++,j++,k++)
				best+= (uint64_t)*i * (uint64_t)*j;
			if (inter > best) best+=Corr;


			return res =  Elem(best % (uint64_t)_field.characteristic());
		}
	}

	template <class Vector1, class Vector2>
	inline GivaroZpz< Givaro::Std32>::Element &DotProductDomain<GivaroZpz< Givaro::Std32> >::dotSpecializedDSP
	(GivaroZpz< Givaro::Std32>::Element &res, const Vector1 &v1, const Vector2 &v2) const
	{
		uint64_t inter,best ;
		inter=best=0;
		if ((v1.first).size()== 0)
			return res=GivaroZpz< Givaro::Std32>::Element(0);
		else {
			uint64_t size      = (v1.first).size();
			uint64_t min       = Max / (Corr+ (uint64_t)(_field.characteristic()-1)*(uint64_t)(_field.characteristic()-1));
			uint64_t min_size  =  (size < min ? size : min);
			uint64_t good1     = (size > min_size ?  size - min_size: 0);
			uint64_t good2     = (long)(size / min_size)* min_size ;
			uint64_t good_size = (good1 > good2 ? good1 : good2 );

			typename Vector1::first_type::const_iterator i_idx  =  v1.first.begin ();
			typename Vector1::second_type::const_iterator i_elt =  v1.second.begin ();

			unsigned long k=0;

			for (;k<min_size;i_idx++,i_elt++,k++)
				best+=(uint64_t)*i_elt * (uint64_t)v2[*i_idx];

			for (inter=best;k<good_size;inter=best) {
				for (unsigned long l=0;l<min_size;i_idx++,i_elt++,k++,l++)
					best+= (uint64_t)*i_elt * (uint64_t)v2[*i_idx];
				if (inter > best) best+=Corr;
			}

			for (;k<size;i_idx++,i_elt++,k++)
				best+= (uint64_t)*i_elt * (uint64_t)v2[*i_idx];
			if (inter > best) best+=Corr;

			return res =  best % _field.characteristic();
		}
	}



	template <class Vector1, class Vector2>
	inline GivaroZpz< Givaro::Std16>::Element &DotProductDomain<GivaroZpz< Givaro::Std16> >::dotSpecializedDD
	(GivaroZpz< Givaro::Std16>::Element &res, const Vector1 &v1, const Vector2 &v2) const
	{
		uint32_t inter,best ;
		inter=best=0;
		if (v1.size() == 0)
			return  res=GivaroZpz< Givaro::Std16>::Element(0);
		else {
			uint32_t size      = v1.size();
			uint32_t min       = Max / (Corr+ ((uint32_t)_field.characteristic()-1)*(uint32_t)(_field.characteristic()-1));
			uint32_t min_size  =  (size < min ? size : min);
			uint32_t good1     = (size > min_size ?  size - min_size: 0);
			uint32_t good2     = (long)(size / min_size)* min_size ;
			uint32_t good_size = (good1 > good2 ? good1 : good2 );


			typename Vector1::const_iterator i=v1.begin();
			typename Vector2::const_iterator j=v2.begin();

			uint32_t k=0;

			for (;k<min_size;i++,j++,k++)
				best+=(uint32_t)*i * (uint32_t)*j;

			for (inter=best;k<good_size;inter=best) {
				for (unsigned long l=0;l<min_size;i++,j++,k++,l++)
					best+= (uint32_t)*i * (uint32_t)*j;
				if (inter > best) best+=Corr;
			}

			for (;k<size;i++,j++,k++)
				best+= (uint32_t)*i * (uint32_t)*j;
			if (inter > best) best+=Corr;

			return res = best % (uint32_t)_field.characteristic();
		}
	}

	template <class Vector1, class Vector2>
	inline GivaroZpz< Givaro::Std16>::Element &DotProductDomain<GivaroZpz< Givaro::Std16> >::dotSpecializedDSP
	(GivaroZpz< Givaro::Std16>::Element &res, const Vector1 &v1, const Vector2 &v2) const
	{
		uint32_t inter,best ;
		inter=best=0;
		if ((v1.first).size()==0)
			return  res=GivaroZpz< Givaro::Std16>::Element(0);
		else {
			uint32_t size      = (v1.first).size();
			uint32_t min       = Max / (Corr+ (uint32_t)(_field.characteristic()-1)*(uint32_t)(_field.characteristic()-1));
			uint32_t min_size  =  (size < min ? size : min);
			uint32_t good1     = (size > min_size ?  size - min_size: 0);
			uint32_t good2     = (long)(size / min_size)* min_size ;
			uint32_t good_size = (good1 > good2 ? good1 : good2 );

			typename Vector1::first_type::const_iterator i_idx  =  v1.first.begin ();
			typename Vector1::second_type::const_iterator i_elt =  v1.second.begin ();
			uint32_t k=0;

			for (;k<min_size;i_idx++,i_elt++,k++)
				best+=(uint32_t)*i_elt * (uint32_t)v2[*i_idx];

			for (inter=best;k<good_size;inter=best) {
				for (unsigned long l=0;l<min_size;i_idx++,i_elt++,k++,l++)
					best+= (uint32_t)*i_elt * (uint32_t)v2[*i_idx];
				if (inter > best) best+=Corr;
			}

			for (;k<size;i_idx++,i_elt++,k++)
				best+= (uint32_t)*i_elt * (uint32_t)v2[*i_idx];
			if (inter > best) best+=Corr;

			return res =  best % (uint32_t)_field.characteristic();
		}
	}


#ifdef XMLENABLED

	template<>
	bool GivaroZpz< Givaro::Std16>::toTag(Writer &W) const
	{
		string s;
		int16_t m = ZpzDom< Givaro::Std16>::residu();

		W.setTagName("field");
		W.setAttribute("implDetail", "givaro-zpz-std16");
		W.setAttribute("cardinality", Writer::numToString(s, m));

		W.addTagChild();
		W.setTagName("finite");

		W.addTagChild();
		W.setTagName("characteristic");
		W.addNum(m);
		W.upToParent();
		W.upToParent();

		return true;
	}

	template <>
	bool GivaroZpz< Givaro::Std32>::toTag(Writer &W) const
	{
		string s;
		int32_t m = ZpzDom< Givaro::Std32>::residu();

		W.setTagName("field");
		W.setAttribute("implDetail", "givaro-zpz-std32");
		W.setAttribute("cardinality", Writer::numToString(s, m));

		W.addTagChild();
		W.setTagName("finite");

		W.addTagChild();
		W.setTagName("characteristic");
		W.addNum(m);
		W.upToParent();

		W.upToParent();

		return true;
	}

	template <>
	bool GivaroZpz<Log16>::toTag(Writer &W) const
	{
		string s;
		int16_t m = ZpzDom<Log16>::residu();

		W.setTagName("field");
		W.setAttribute("implDetail", "givaro-zpz-log16");
		W.setAttribute("cardinality", Writer::numToString(s, m));

		W.addTagChild();
		W.setTagName("finite");

		W.addTagChild();
		W.setTagName("characteristic");
		W.addNum(m);
		W.upToParent();

		W.upToParent();

		return true;
	}
#endif


}

#endif //__LINBOX_givaro_zpz_INL


// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

