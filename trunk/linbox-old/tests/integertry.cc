# this is a very temporary test to check if integers are working.

#include <iostream.h>
/*
#include "LinBox/integer.h"      // integer = long
*/
#include "LinBox/lin_integer.h"  // integer = Givaro Integer -> gmp

using namespace LinBox;
main()
{
integer a, b;
a = 2; b = 3;
cout << a + b << " " << a * b << endl;
}
