/* lb-domain.C
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

#ifndef __LINBOX_lb_domain_C
#define __LINBOX_lb_domain_C

#include "lb-domain.h"
#include "lb-domain-data.h"
#include "lb-domain-function.h"
#include "lb-domain-type.h"


/****************************
 * Allocate global variable *
 ****************************/

// global hash table for allocated domains
DomainTable domain_hashtable;

// global variable for the factory
Domain_Factory linbox_domain;





/*******************************************
 * Update the Factory with all domain type *
 *******************************************/
//extern Domain_Factory linbox_domain;

void UpdateDomain(){
	linbox_domain.add("linbox_field_dbl"      , constructDomain<Givaro::Modular<double> >);
	//linbox_domain.add("linbox_field_rational" , constructDomain<LinBox::GMPRationalField>);
	linbox_domain.add("linbox_ring_integer"   , constructDomain<Givaro::ZRing<Givaro::Integer> >);
#ifndef __LINBOX_MINIMIZE_DOMAIN
	linbox_domain.add("linbox_field_32"       , constructDomain<Givaro::Modular<int32_t> >);
	linbox_domain.add("linbox_field_64"       , constructDomain<Givaro::Modular<int64_t> >);
	linbox_domain.add("linbox_field_mp"       , constructDomain<Givaro::Modular<Givaro::Integer> >);
#endif
#ifdef __LINBOX_HAVE_NTL
        //linbox_domain.add("ntl_field_ZZ_p"      , constructDomain<LinBox::NTL_ZZ_p>);
#ifndef __LINBOX_MINIMIZE_DOMAIN
	//linbox_domain.add("ntl_field_zz_p"      , constructDomain<LinBox::NTL_zz_p >);
	//linbox_domain.add("ntl_ring_integer"    , constructDomain<LinBox::NTL_ZZ>);
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

// global variable for current prime field type
const char* current_prime_field  = default_prime_field;

// global variable for current rational field type
const char* current_rational_field  = default_rational_field;

// global variable for current integer ring type
const char* current_integer_ring = default_integer_ring;




/****************************
 * API to contruct domains  *
 ****************************/

const DomainKey& createDomain( const LinBox::integer characteristic, const char *name){
	const char* type=name;
	if (name == NULL){
		if (characteristic == 0)
			type = current_integer_ring;
		else
			type = current_prime_field;
	}

	DomainKey key(characteristic, type);

	// check if the domain is already constructed in domain hashtable
	// if so return the pointer to the domain
	DomainTable::const_iterator it= domain_hashtable.find(key) ;
	if (it != domain_hashtable.end()){
		it->first.copy();
		return it->first;
	}

	// need to create a new field
	DomainAbstract *domain = linbox_domain.create(type, characteristic);


	std::pair<DomainTable::const_iterator, bool> status= domain_hashtable.insert(std::pair<DomainKey, DomainAbstract*> (key, domain));

	if (status.second)
		return status.first->first;
	else
			throw lb_runtime_error("LinBox ERROR: domain creation failed \n");// throw an exception
}


/************************
 * API to copy domains  *
 ************************/
const DomainKey copyDomain( const DomainKey &k){
	return k;
}



/**********************************************
 * API to modify the current prime field type *
 **********************************************/
void setPrimeField(const char* t){
	current_prime_field= t;
}


/*************************************************
 * API to modify the current rational field type *
 *************************************************/
void setRationalField(const char* t){
	current_rational_field= t;
}


/***********************************************
 * API to modify the current integer ring type *
 ***********************************************/
void setIntegerRing(const char* t){
	current_integer_ring= t;
}

/*********************************
 * API to write info on a domain *
 *********************************/
void writeDomainInfo(const DomainKey &key, std::ostream& os){
	os<<"[LinBox Domain (type = "<<key.Type()<<", charact = "<<key.Characteristic()<<")]\n";
}

#endif // end of file


// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

