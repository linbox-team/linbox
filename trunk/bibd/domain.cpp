/* domain.c
 * Written by: Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * Various utility functions for the base field we are using
 */

#include "domain.h"

/* ext_euclid
 *
 * Extended Euclidean algorithm implementation
 */

static void ext_euclid (unsigned int a, unsigned int b, unsigned int *d, 
			int *u, int *v) 
{
	long long r[2], s[2], t[2], q;
	int bit = 0, lastbit = 0;

	r[0] = MAX (a, b);
	r[1] = MIN (a, b);
	s[0] = 1; t[0] = 0;
	s[1] = 0; t[1] = 1;

	while (r[bit] != 0) {
		bit = lastbit;
		lastbit = (bit + 1) % 2;
		q = r[bit] / r[lastbit];
		r[bit] = r[bit] % r[lastbit];
		s[bit] = s[bit] - q * s[lastbit];
		t[bit] = t[bit] - q * t[lastbit];
	}

	*d = r[lastbit];
	*u = (a > b) ? s[lastbit] : t[lastbit];
	*v = (a > b) ? t[lastbit] : s[lastbit];
}

/* inv_mod
 *
 * Find the modular inverse of a number
 */

unsigned int inv_mod (unsigned int x, unsigned int n) 
{
	unsigned d;
	int u, v;

	ext_euclid (x, n, &d, &u, &v);
	if (d != 1) return 0;
	return ((u >= 0) ? u : u + n) % n;
}

GMP_Rational_Field field;
