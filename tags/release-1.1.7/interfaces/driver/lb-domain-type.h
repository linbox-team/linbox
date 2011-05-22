/* lb-domain-type.h
 * Copyright (C) 2005 Pascal Giorgi
 *
 * Written by Pascal Giorgi <pgiorgi@uwaterloo.ca>
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

#ifndef __LINBOX_lb_domain_type_H
#define __LINBOX_lb_domain_type_H

#include <linbox/field/modular.h>
#include <linbox/field/PID-integer.h>
#include <linbox/field/gmp-rational.h>

#ifdef __LINBOX_HAVE_NTL
#include <linbox/field/ntl-lzz_p.h>
#include <linbox/field/ntl-ZZ_p.h>
#include <linbox/field/ntl-ZZ.h>
#endif

#ifdef __LINBOX_HAVE_GIVARO
#include <linbox/field/givaro-zpz.h>
#include <linbox/field/givaro-gfq.h>
#endif

/**************************************
 * Define the list of all Domain Type *
 **************************************/

#define __LINBOX_MINIMIZE_DOMAIN
#define __LINBOX_DOMAIN_ONLY

typedef LinBoxTypelist < LinBox::Modular<double>          , LinBoxDumbType> DL1;
typedef LinBoxTypelist < LinBox::PID_integer              , DL1> DL2;
typedef LinBoxTypelist < LinBox::GMPRationalField         , DL2> DL3;
typedef LinBoxTypelist < LinBox::Modular<LinBox::int32>   , DL3> DL4; 
typedef LinBoxTypelist < LinBox::Modular<LinBox::int64>   , DL4> DL5;
typedef LinBoxTypelist < LinBox::Modular<LinBox::integer> , DL5> DL6;

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

#ifdef __LINBOX_HAVE_GIVARO
typedef LinBoxTypelist < LinBox::GivaroGfq, LinBoxDumbType> DG1;
typedef LinBoxTypelist < LinBox::GivaroZpz<Std32>, DG1> DG2;
#ifdef __LINBOX_MINIMIZE_DOMAIN
typedef DG1 givaro_domain_list;
#else
typedef DG2 givaro_domain_list;
#endif
#else
typedef LinBoxDumbType givaro_domain_list;
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
	linbox_domain.add("linbox_field_dbl"      , constructDomain<LinBox::Modular<double> >);
	linbox_domain.add("linbox_field_rational" , constructDomain<LinBox::GMPRationalField>);
	linbox_domain.add("linbox_ring_integer"   , constructDomain<LinBox::PID_integer>);
#ifndef __LINBOX_MINIMIZE_DOMAIN	
	linbox_domain.add("linbox_field_32"       , constructDomain<LinBox::Modular<LinBox::int32> >);
	linbox_domain.add("linbox_field_64"       , constructDomain<LinBox::Modular<LinBox::int64> >);
	linbox_domain.add("linbox_field_mp"       , constructDomain<LinBox::Modular<LinBox::integer> >);
#endif
#ifdef __LINBOX_HAVE_NTL
	linbox_domain.add("ntl_field_ZZ_p"      , constructDomain<LinBox::NTL_ZZ_p>);
#ifndef __LINBOX_MINIMIZE_DOMAIN	
	linbox_domain.add("ntl_field_zz_p"      , constructDomain<LinBox::NTL_zz_p >);
	linbox_domain.add("ntl_ring_integer"    , constructDomain<LinBox::NTL_ZZ>);
#endif
#endif
#ifdef __LINBOX_HAVE_GIVARO
	linbox_domain.add("givaro_field_gfq"    , constructDomain<LinBox::GivaroGfq>);
#ifndef  __LINBOX_MINIMIZE_DOMAIN
	linbox_domain.add("givaro_field_32"     , constructDomain<LinBox::GivaroZpz<Std32> >);
#endif
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
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
