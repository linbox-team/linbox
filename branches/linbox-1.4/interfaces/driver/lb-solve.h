/* lb-solve.h
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

#ifndef __LINBOX_lb_solve_H
#define __LINBOX_lb_solve_H

#include <lb-blackbox-collection.h>
#include <lb-vector-collection.h>


/***********************************************************
 * API for solving linear systems                          *
 * vector solution  is returned through a given vector key *
 ***********************************************************/
void lb_solve(const VectorKey &res, const BlackboxKey &Bkey, const VectorKey &Vkey);


/*****************************************************
 * API for solving linear systems                    *
 * vector solution  is returned through a vector key *
 *****************************************************/
const VectorKey&  lb_solve(const BlackboxKey &Bkey, const VectorKey &Vkey);



#endif

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

