/* lb-domain.C
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

#ifndef __LINBOX_lb_domain_C
#define __LINBOX_lb_domain_C

#include <lb-domain.h>
#include <lb-domain-function.h>

	
/****************************
 * Allocate global variable *
 ****************************/

// global hash table for allocated domains
DomainTable domain_hashtable;

// global variable for the factory
Domain_Factory linbox_domain;
	
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
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
