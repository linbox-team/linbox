/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/*
 * examples/integer-mul.C
 *
 * Copyright (C) 2002, 2005, 2010 G Villard, D. Saunders
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 */

/** \file examples/integer-mul.C
 * @example examples/integer-mul.C
 * \author Gilles Villard
 * \brief The LinBox arbitrary precision integer type illustrated.
 * \ingroup examples
 *
 * The class `integer' is a wrapper of <a href=http://gmplib.org>GMP</a> integers.
 */

// ---------------------------------------------
#include <iostream>
#include <fstream>
// ---------------------------------------------

#include "linbox/linbox-config.h"

// Use of Gmp based LinBox integers
#include "linbox/integer.h"

using namespace LinBox;
using namespace std;

// ---------------------------------------------

/// no command line args.  Prompts for two integers.
int main()
{

	integer a,b;

	cout << "1st integer > ";
	cin >> a;
	cout << "2nd integer > ";
	cin >> b;

	cout << "The product " << a*b << "\n";

	return 0;
};
