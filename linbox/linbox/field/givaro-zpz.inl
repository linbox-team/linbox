/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/field/givaro-zpz.inl
 * Written by Pascal Giorgi <pascal.giorgi@ens-lyon.fr>
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */




#ifndef __GIVARO_ZPZ_INL 
#define __GIVARO_ZPZ_INL


#include <iostream>
#include "linbox/integer.h"

namespace LinBox {



template <class Vector1, class Vector2>
inline GivaroZpz<Std32>::Element &DotProductDomain<GivaroZpz<Std32> >::dotSpecializedDD
	(GivaroZpz<Std32>::Element &res, const Vector1 &v1, const Vector2 &v2) const
{
        uint64 inter,best ;
	inter=best=0;

	uint64 size      = v1.size();
	uint64 min       = Max / (Corr+ (uint64)(_F._p-1)*(uint64)(_F._p-1));
	uint64 min_size  =  (size < min ? size : min);
	uint64 good1     = (size > min_size ?  size - min_size: 0);
	uint64 good2     = (long)(size / min_size)* min_size ;
	uint64 good_size = (good1 > good2 ? good1 : good2 );
	
	typename Vector1::const_iterator i=v1.begin();
	typename Vector2::const_iterator j=v2.begin();
		
	unsigned long k=0;
	
	for (;k<min_size;i++,j++,k++) 
		best+=(uint64)*i * (uint64)*j;

	for (inter=best;k<good_size;inter=best) {
		for (unsigned long l=0;l<min_size;i++,j++,k++,l++)
			best+= (uint64)*i * (uint64)*j;
		if (inter > best) best+=Corr;	
	 }
	
	for (;k<size;i++,j++,k++)
		best+= (uint64)*i * (uint64)*j;
	if (inter > best) best+=Corr;

	
	return res =  best % (uint64)_F._p;
}

template <class Vector1, class Vector2>
inline GivaroZpz<Std32>::Element &DotProductDomain<GivaroZpz<Std32> >::dotSpecializedDSP
	(GivaroZpz<Std32>::Element &res, const Vector1 &v1, const Vector2 &v2) const
{
	uint64 inter,best ;
	inter=best=0;
	
	uint64 size      = (v1.first).size();
	uint64 min       = Max / (Corr+ (uint64)(_F._p-1)*(uint64)(_F._p-1));
	uint64 min_size  =  (size < min ? size : min);
	uint64 good1     = (size > min_size ?  size - min_size: 0);
	uint64 good2     = (long)(size / min_size)* min_size ;
	uint64 good_size = (good1 > good2 ? good1 : good2 );

	typename Vector1::first_type::const_iterator i_idx;
	typename Vector1::second_type::const_iterator i_elt;
  
	unsigned long k=0;
	
	for (;k<min_size;i_idx++,i_elt++,k++) 
		best+=(uint64)*i_elt * (uint64)v2[*i_idx];

	for (inter=best;k<good_size;inter=best) {
		for (unsigned long l=0;l<min_size;i_idx++,i_elt++,k++,l++)
			best+= (uint64)*i_elt * (uint64)v2[*i_idx];
		if (inter > best) best+=Corr;	
	 }
	
	for (;k<size;i_idx++,i_elt++,k++)
		best+= (uint64)*i_elt * (uint64)v2[*i_idx];
	if (inter > best) best+=Corr;
	
	return res =  best % _F._p;
}



template <class Vector1, class Vector2>
inline GivaroZpz<Std16>::Element &DotProductDomain<GivaroZpz<Std16> >::dotSpecializedDD
	(GivaroZpz<Std16>::Element &res, const Vector1 &v1, const Vector2 &v2) const
{
        uint32 inter,best ;
	inter=best=0;
	
	uint32 size      = v1.size();
	uint32 min       = Max / (Corr+ ((uint32)_F._p-1)*(uint32)(_F._p-1));
	uint32 min_size  =  (size < min ? size : min);
	uint32 good1     = (size > min_size ?  size - min_size: 0);
	uint32 long good2     = (long)(size / min_size)* min_size ;
	uint32 long good_size = (good1 > good2 ? good1 : good2 );
	

	typename Vector1::const_iterator i=v1.begin();
	typename Vector2::const_iterator j=v2.begin();
		
	unsigned long k=0;
	
	for (;k<min_size;i++,j++,k++) 
		best+=(uint32)*i * (uint32)*j;

	for (inter=best;k<good_size;inter=best) {
		for (unsigned long l=0;l<min_size;i++,j++,k++,l++)
			best+= (uint32)*i * (uint32)*j;
		if (inter > best) best+=Corr;	
	 }
	
	for (;k<size;i++,j++,k++)
		best+= (uint32)*i * (uint32)*j;
	if (inter > best) best+=Corr;
	
	return res = best % (uint32)_F._p;
}

template <class Vector1, class Vector2>
inline GivaroZpz<Std16>::Element &DotProductDomain<GivaroZpz<Std16> >::dotSpecializedDSP
	(GivaroZpz<Std16>::Element &res, const Vector1 &v1, const Vector2 &v2) const
{
	uint32 inter,best ;
	inter=best=0;
	
	uint32 size      = (v1.first).size();
	uint32 min       = Max / (Corr+ (uint32)(_F._p-1)*(uint32)(_F._p-1));
	uint32 min_size  =  (size < min ? size : min);
	uint32 good1     = (size > min_size ?  size - min_size: 0);
	uint32 good2     = (long)(size / min_size)* min_size ;
	uint32 good_size = (good1 > good2 ? good1 : good2 );

	typename Vector1::first_type::const_iterator i_idx;
	typename Vector1::second_type::const_iterator i_elt;
  
	unsigned long k=0;
	
	for (;k<min_size;i_idx++,i_elt++,k++) 
		best+=(uint32)*i_elt * (uint32)v2[*i_idx];

	for (inter=best;k<good_size;inter=best) {
		for (unsigned long l=0;l<min_size;i_idx++,i_elt++,k++,l++)
			best+= (uint32)*i_elt * (uint32)v2[*i_idx];
		if (inter > best) best+=Corr;	
	 }
	
	for (;k<size;i_idx++,i_elt++,k++)
		best+= (uint32)*i_elt * (uint32)v2[*i_idx];
	if (inter > best) best+=Corr;
	
	return res =  best % (uint32)_F._p;
}




}




#endif
