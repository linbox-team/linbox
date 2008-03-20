/*  Author Z. Wan
 *  Modified to fit in linbox
 *  Using GMP directly
 *  Implement thte idea in the manuscript, available at http://www.cis.udel.edu/~wan/jsc_wan.ps
 */
#ifndef __RATIONAL_RECONSTRUCTION2_H__
#define __RATIONAL_RECONSTRUCTION2_H__

#include <linbox/integer.h>

namespace LinBox{
static int rational_reconstruction (integer& a, integer& b, const integer& n0, const integer& d0, const integer& B) {

	integer p0, p1, p2; integer q0, q1, q2;
	integer q, r; integer n, d;

	p0 = 0; p1 = 1; p2 = 0; q0 = 1; q1 = 0; q2 = 0;
	n = abs (n0); d = abs (d0);

	int sgn = sign(n0) * sign(d0);

	if (d0 < B) {
		if (sgn > 0)
			 a = n0;
		else
			a = -n0;

		b= d0;
		return 0;
	}
	
	do {
		integer::divmod (q, r, n, d);
		p2 = p0 + q * p1;
		q2 = q0 + q * q1;
		p0 = p1; p1 = p2;
		q0 = q1; q1 = q2;
		n = d; d = r;

	} while ((q2 < B) && (sign(d) != 0));

	if (q2 > B) {
		if (sgn >= 0) 
			a = p0;
		else
			a = - p0;
		
		b = q0;
	}
	else {
		if (sgn >= 0)
			a = p1;
		else
			a = -p1;
		b = q1;
	}
	
	return 0;
}
}// LinBox

#endif
