/* -bds 7/00 
The class integer is used throughout linbox where "casual" use of 
arbitrary precision integers is needed.  Examples are cardinalities
of sets that may be very large.  In such casual use, operator forms of 
arithmetic are used.  This is distinct from using a domain representation
of the integers and it's functions such as add, mul, axpy, etc.
*/
#ifndef __lin_integer__
#define __lin_integer__
#include <stdlib.h>
/* we use a Givaro based wrapper of gmp with the above location. 
 */
#include <GMP++/gmp++.C>

/* Now gmp 4.0 supports the C++ operator forms, 
 * so we may want to use it directly.
 * ...however, it is immature, so we may wait and try it later.
#include <gmpxx.h>
 */

namespace LinBox{

// these are unused defs - size_t is being used where int_32 might have been.
typedef int int_32; // should guarantee 32 bits though
typedef long long int_64; // should guarantee 64 bits though

/* the gmp++.C based def:
 */
typedef Integer integer;

/* the gmp 4.0 based def:
typedef mpz_class integer;
*/

} // namespace LinBox
#endif
