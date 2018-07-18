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

#include "linbox/ring/modular.h"
#include "givaro/zring.h"
#include "linbox/field/gmp-rational.h"
#include "linbox/field/givaro.h"

#ifdef __LINBOX_HAVE_NTL
#include "linbox/ring/ntl.h"
#endif

/**************************************
 * Define the list of all Domain Type *
 **************************************/

#define __LINBOX_MINIMIZE_DOMAIN
#define __LINBOX_DOMAIN_ONLY

typedef LinBoxTypelist < Givaro::Modular<double>          , LinBoxDumbType> DL1;
typedef LinBoxTypelist < LinBox::Givaro::ZRing<Integer>              , DL1> DL2;
typedef LinBoxTypelist < LinBox::GMPRationalField         , DL2> DL3;
typedef LinBoxTypelist < Givaro::Modular<int32_t>         , DL3> DL4;
typedef LinBoxTypelist < LinBox::Givaro::Modular<int64_t>         , DL4> DL5;
typedef LinBoxTypelist < Givaro::Modular<LinBox::integer> , DL5> DL6;

#ifdef __LINBOX_MINIMIZE_DOMAIN
typedef DL3 linbox_domain_list;
#else
typedef DL6 linbox_domain_list;
#endif

#ifdef __LINBOX_HAVE_NTL
typedef LinBoxTypelist < LinBox::NTL_ZZ_p, LinBoxDumbType> DN1;
typedef LinBoxTypelist < LinBox::NTL_zz_p, DN1> DN2;
typedef LinBoxTypelist < LinBox::NTL_ZZ  , DN2> DN3;

#ifdef __LINBOX_MINIMIZE_DOMAIN
typedef DN1 ntl_domain_list;
#else
typedef DN3 ntl_domain_list;
#endif
#else
typedef LinBoxDumbType ntl_domain_list;
#endif

typedef LinBoxTypelist < LinBox::Givaro::GFq, LinBoxDumbType> DG1;
typedef LinBoxTypelist < Givaro::Modular< int32_t>, DG1> DG2;
#ifdef __LINBOX_MINIMIZE_DOMAIN
typedef DG1 givaro_domain_list;
#else
typedef DG2 givaro_domain_list;
#endif

// define DomainList to be the list of all domains
#ifdef __LINBOX_DOMAIN_ONLY
typedef linbox_domain_list DomainList;
#else
typedef LinBoxTL::Append< linbox_domain_list, LinBoxTL::Append< ntl_domain_list, givaro_domain_list>::Result>::Result DomainList;
#endif





/*******************************************
 * Update the Factory with all domain type *
 *******************************************/
extern Domain_Factory linbox_domain;

void UpdateDomain(){
	linbox_domain.add("linbox_field_dbl"      , constructDomain<Givaro::Modular<double> >);
	linbox_domain.add("linbox_field_rational" , constructDomain<LinBox::GMPRationalField>);
	linbox_domain.add("linbox_ring_integer"   , constructDomain<LinBox::Givaro::ZRing<Integer> >);
#ifndef __LINBOX_MINIMIZE_DOMAIN
	linbox_domain.add("linbox_field_32"       , constructDomain<Givaro::Modular<int32_t> >);
	linbox_domain.add("linbox_field_64"       , constructDomain<LinBox::Givaro::Modular<int64_t> >);
	linbox_domain.add("linbox_field_mp"       , constructDomain<Givaro::Modular<LinBox::integer> >);
#endif
#ifdef __LINBOX_HAVE_NTL
	linbox_domain.add("ntl_field_ZZ_p"      , constructDomain<LinBox::NTL_ZZ_p>);
#ifndef __LINBOX_MINIMIZE_DOMAIN
	linbox_domain.add("ntl_field_zz_p"      , constructDomain<LinBox::NTL_zz_p >);
	linbox_domain.add("ntl_ring_integer"    , constructDomain<LinBox::NTL_ZZ>);
#endif
#endif

	linbox_domain.add("givaro_field_gfq"    , constructDomain<LinBox::Givaro::GFq>);
#ifndef  __LINBOX_MINIMIZE_DOMAIN
	linbox_domain.add("givaro_field_32"     , constructDomain<Givaro::Modular< int32_t> >);
#endif
}


/****************************
 * Default type for Domains *
 ****************************/

// definition of the default type for prime field
#define default_prime_field  "linbox_field_dbl"

// definition of the default type for rational field
#define default_rational_field "linbox_field_rational"

// definition of the default type for integer ring
#define default_integer_ring "linbox_ring_integer"



#endif

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
