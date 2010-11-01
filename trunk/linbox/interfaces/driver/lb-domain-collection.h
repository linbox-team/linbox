/* lb-domain-collection.h
 * Copyright (C) 2005 Pascal Giorgi
 *
 * Written by Pascal Giorgi <pgiorgi@uwaterloo.ca>
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 du -h tes *
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

#ifndef __LINBOX_lb_domain_collection_H
#define __LINBOX_lb_domain_collection_H

#include <linbox/integer.h>

#include <map>
#include <lb-domain-abstract.h>


/*************************
 * Collection of Domains *
 *************************/

// definition of a unique key (using embedded  reference counting to manage domain)
class DomainKey {
private:
	LinBox::integer    *pfirst;
	const char       **psecond;
	size_t            *counter;
	mutable bool        autogc;


public:
	DomainKey(const LinBox::integer &p, const char* name) 
		: pfirst(new LinBox::integer(p)), psecond(new const char*(name)), counter(new size_t(0)), autogc(false) 
	{}

	DomainKey(const DomainKey &k, bool gc = false) 
		: pfirst(k.pfirst), psecond(k.psecond), counter(k.counter), autogc(gc) {(*counter)++;}

	~DomainKey() {
		if (autogc){
			extern void deleteDomain(const DomainKey &key);
			if ((*counter) == 0) 
				deleteDomain(*this);
			else 
				(*counter)--;	
			
		}
		else{ 		
			if ((*counter) == 0) { delete counter; delete pfirst; delete psecond;}
			else (*counter)--;	
		}
	}	
	
	const DomainKey& operator= (const DomainKey &k) {
		if (autogc){
			extern void deleteDomain(const DomainKey &key);
			if ((*counter) == 0) 
				deleteDomain(*this);
			else 
				(*counter)--;	
			
		}
		else{ 		
			if ((*counter) == 0) {  delete counter; delete pfirst; delete psecond;}
			else (*counter)--;	
		}
		 
		pfirst  = k.pfirst;
		psecond = k.psecond;
		counter = k.counter;
		autogc  = k.autogc;
		(*counter)++;
		return *this;
	}

	bool free() const {return ((*counter) == 0);}
	
	void dispose() const { (*counter)--; }

	void copy() const { (*counter)++; }

	void set_autogc() const {autogc=true;}

	bool lessThan(const DomainKey &k) const {
		return (strcmp(*psecond, *(k.psecond))<0) || ((strcmp(*psecond, *(k.psecond))==0) && ( *pfirst < *(k.pfirst )));
	}

	LinBox::integer Characteristic() const {return *pfirst;}

	const char* Type() const {return *psecond;}

	size_t getcounter() const  {return *counter;}
};


 
// comparison functor on key
struct DomainKeyLessThan{
	bool operator()(const DomainKey& k1, const DomainKey &k2) 
	{ return k1.lessThan(k2);}
};

// definition of a hash table type
typedef std::map<DomainKey, DomainAbstract*, DomainKeyLessThan>  DomainTable;



#endif
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
