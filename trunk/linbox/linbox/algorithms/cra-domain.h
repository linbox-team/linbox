/* linbox/algorithms/cra-domain.h
 * Copyright (C) 1999-2010 The LinBox group
 *
 * Selector for ChineseRemainder
 * Parallel versions are transparent to the user
 * Time-stamp: <30 Mar 10 15:11:42 Jean-Guillaume.Dumas@imag.fr> 
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */
#ifndef __LINBOX_cra_domain_H
#define __LINBOX_cra_domain_H

//#if defined(OMP_H)
#ifdef _OPENMP

#include "linbox/algorithms/cra-domain-omp.h"
namespace LinBox 
{
    template<class CRABase>
    struct ChineseRemainder : public ChineseRemainderOMP<CRABase> {
        typedef typename CRABase::Domain	Domain;
        typedef typename CRABase::DomainElement	DomainElement;
        template<class Param>
        ChineseRemainder(const Param& b) : ChineseRemainderOMP<CRABase>(b) {}
        
	ChineseRemainder(const CRABase& b) : ChineseRemainderOMP<CRABase>(b) {}
    };
}

#else

#include "linbox/algorithms/cra-domain-seq.h"
namespace LinBox 
{
    template<class CRABase>
    struct ChineseRemainder : public ChineseRemainderSeq<CRABase> {
        typedef typename CRABase::Domain	Domain;
        typedef typename CRABase::DomainElement	DomainElement;
        template<class Param>
        ChineseRemainder(const Param& b) : ChineseRemainderSeq<CRABase>(b) {}
        
	ChineseRemainder(const CRABase& b) : ChineseRemainderSeq<CRABase>(b) {}
    };
}

#endif


#endif
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
