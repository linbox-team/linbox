/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* Copyright (C) 1999 LinBox
 * Written by <Jean-Guillaume.Dumas@imag.fr>
 *
 *
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

// ================================================================
// Black Box iterator and container
// For symmetric matrix with same left and right vector
// the sequence is u^t v, u^t A v, ...,  u^t A^n v,
// ================================================================


#ifndef __LINBOX_bbcontainer_symmetric_H
#define __LINBOX_bbcontainer_symmetric_H

#include <LinBox/lin_rand.h>
#include <LinBox/lin_base_bbit.h>

template<class BlackBoxDomain, class Vecteur = typename BlackBoxDomain::PreferredInMatrix_t, class RandIter = Random>
class BB_Symmetric_Container : public Base_BB_Container< BlackBoxDomain, Vecteur > {
public:
	BB_Symmetric_Container() {}

	BB_Symmetric_Container(BlackBoxDomain_t * D, const Vecteur& u0) :
		Base_BB_Container< BlackBoxDomain, Vecteur>(D) { init(u0,u0); }

	BB_Symmetric_Container(BlackBoxDomain_t * D, RandIter& g ) :
		Base_BB_Container< BlackBoxDomain, Vecteur>(D) { init(g); }

protected:
	void _launch () {
		if (casenumber > 0) {
			if (casenumber == 1) {
				casenumber = 2;
				_BB_domain->Apply( v, u);  // v <- B(B^i u_0) = B^(i+1) u_0
				DOTPROD(_value,u,v);       // t <- u^t v = u_0^t B^(2i+1) u_0
			} else {
				casenumber = -1;
				DOTPROD(_value,v,v);       // t <- v^t v = u_0^t B^(2i+2) u_0
			}
		} else {
			if (casenumber == 0) {
				casenumber = 1;
				DOTPROD(_value,u,u);       // t <- u^t u = u_0^t B^(2i+4) u_0
			} else {
				casenumber = 0;
				_BB_domain->Apply( u, v);  // u <- B(B^(i+1) u_0) = B^(i+2) u_0
				DOTPROD(_value,v,u);       // t <- v^t u = u_0^t B^(2i+3) u_0
			}
		}
	}

	void _wait () {}
};


#endif // __LINBOX_bbcontainer_symmetric_H

