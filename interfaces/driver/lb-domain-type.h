/* lb-domain-type.h
 * Copyright (C) 2005 Pascal Giorgi
 *
 * Written by Pascal Giorgi <pgiorgi@uwaterloo.ca>
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

#ifndef __LINBOX_lb_domain_type_H
#define __LINBOX_lb_domain_type_H


#include "givaro/zring.h"
#include "linbox/field/gmp-rational.h"
#include "givaro/givrational.h"
#include "givaro/modular.h"

#ifdef __LINBOX_HAVE_NTL
//#include "linbox/ring/ntl.h"
#endif

/**************************************
 * Define the list of all Domain Type *
 **************************************/

#define __LINBOX_MINIMIZE_DOMAIN
#define __LINBOX_DOMAIN_ONLY

typedef LinBoxTypelist < Givaro::Modular<double>          , LinBoxDumbType> DL1;
typedef LinBoxTypelist < Givaro::ZRing<Givaro::Integer>   , DL1> DL2;
typedef LinBoxTypelist < Givaro::Modular<Givaro::Integer> , DL2> DL3;
typedef LinBoxTypelist < Givaro::QField<Givaro::Rational> , DL3> DL4;
typedef LinBoxTypelist < Givaro::Modular<int64_t>         , DL4> DL5;


#ifdef __LINBOX_MINIMIZE_DOMAIN
typedef DL4 linbox_domain_list;
#else
typedef DL6 linbox_domain_list;
#endif

//#ifdef __LINBOX_HAVE_NTL
//typedef LinBoxTypelist < LinBox::NTL_ZZ_p, LinBoxDumbType> DN1;
//typedef LinBoxTypelist < LinBox::NTL_zz_p, DN1> DN2;
//typedef LinBoxTypelist < LinBox::NTL_ZZ  , DN2> DN3;

// #ifdef __LINBOX_MINIMIZE_DOMAIN
// typedef DN1 ntl_domain_list;
// #else
// typedef DN3 ntl_domain_list;
// #endif
// #else
typedef LinBoxDumbType ntl_domain_list;
//#endif


// define DomainList to be the list of all domains
#ifdef __LINBOX_DOMAIN_ONLY
typedef linbox_domain_list DomainList;
#else
typedef LinBoxTL::Append< linbox_domain_list, ntl_domain_list>::Result DomainList;
#endif








#endif

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
