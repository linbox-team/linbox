/* domain.h
 * Written by: Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * Declarations for the base field we are using
 */

#ifndef __DOMAIN_H
#define __DOMAIN_H

#include "LinBox/gmp-rational-field.h"
#include "LinBox/gmp-rational-number.h"
#include "LinBox/gmp-rational-random.h"

#ifndef MIN
#  define MIN(a, b) ((a < b) ? a : b)
#endif
#ifndef MAX
#  define MAX(a, b) ((a > b) ? a : b)
#endif

using namespace LinBox;

typedef GMP_Rational_Number Element;
typedef GMP_Rational_Random Random;

extern GMP_Rational_Field field;

#endif /* __DOMAIN_H */
