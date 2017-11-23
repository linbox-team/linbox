/* lb-mul.h
 * Copyright (C) 2017 Pascal Giorgi
 *
 * Written by Pascal Giorgi <pascal.giorgi@lirmm.fr
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

#ifndef __LINBOX_lb_mul_H
#define __LINBOX_lb_mul_H

#include "lb-blackbox-collection.h"


/*******************************************************
 * API for matrix multiplication                       *
 * matrix result is given through a blackbox key       *
 *******************************************************/
void lb_mul(const BlackboxKey &Ckey, const BlackboxKey &Akey, const BlackboxKey &Bkey);


/*******************************************************
 * API for matrix multiplication                       *
 * matrix result is given through a blackbox key       *
 *******************************************************/
const BlackboxKey & lb_mul(const BlackboxKey &Akey, const BlackboxKey &Bkey);


#endif
