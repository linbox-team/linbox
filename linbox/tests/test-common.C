/* linbox/tests/test-common.C
 * Copyright (C) 2001, 2002 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
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


/*! @file  tests/test-common.C
 * @ingroup tests
 * @brief  no doc
 */



#ifndef __LINBOX_test_common_C
#define __LINBOX_test_common_C


#include "linbox-config.h"

#include <iostream>
#include <fstream>

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <list>

#include "linbox/util/commentator.h"
#include "linbox/field/archetype.h"

#include "test-common.h"

bool isPower (LinBox::integer n, LinBox::integer m)
{
	return (n == 1) || (((n % m) == 0) && isPower (n/m, m));
}

inline double incompleteGamma (double a, double x, double tol)
{
	double xa_ex = pow (x, a) * exp (-x);
	double pi = 1.0;
	double xn = 1.0;
	double sigma = 0.0;
	double last_sigma;

	int n = 0;

	do {
		pi *= a + n;
		last_sigma = sigma;
		sigma += xn / pi;
		xn *= x;
		++n;
	} while (abs (sigma - last_sigma) >= tol) ;

	return sigma * xa_ex;
}

double chiSquaredCDF (double chi_sqr, double df)
{
	return incompleteGamma (df / 2.0, chi_sqr / 2.0, 1e-10) / exp (lgamma (df / 2.0));
}
#endif // __LINBOX_test_common_C


// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

