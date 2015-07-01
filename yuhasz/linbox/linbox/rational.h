/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/rational.h
 * Copyright (C) 2002 Bradford Hovinen
 *
 * Written by Bradford Hovinen <bghovinen@math.uwaterloo.ca>
 *
 * -------------------------------------
 *
 * See COPYING for license information
 */

#ifndef __RATIONAL_H
#define __RATIONAL_H

#include "linbox-config.h"

#include "linbox/field/gmp-rational.h"

namespace LinBox
{
	typedef GMPRationalElement rational;
}

#endif // __RATIONAL_H
