// this is a very temporary test to check if integers are working. -bds

#include <iostream.h>
#include "LinBox/integer.h"      // integer = long
/*
#include "LinBox/lin_integer.h"  // integer = Givaro Integer -> gmp
*/

using namespace LinBox;
main()
{
integer a, b;
a = 84; b = 60;
cout << a + b << " " << a * b << endl;
integer g, u, v;
g = gcd(a,b,u,v);
cout << g << " = " << a << " * " << u << " + " << b << " * " << v << endl;
}
