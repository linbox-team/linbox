
#include <NTL/quad_float.h>

main()
{
   quad_float a, b, c, d;

   quad_float::SetOutputPrecision(25);

   if (PrecisionOK())
      cout << "Precision OK\n";
   else
      cout << "Precision not OK\n";


   cin >> a;
   cout << a << "\n";

   cin >> b;
   cout << b << "\n";

   c = a + b;
   d = a;
   d += b;
   cout << c << "\n";
   cout << d << "\n";

   c = a - b;
   d = a;
   d -= b;
   cout << c << "\n";
   cout << d << "\n";

   c = a * b;
   d = a;
   d *= b;
   cout << c << "\n";
   cout << d << "\n";

   c = a / b;
   d = a;
   d /= b;
   cout << c << "\n";
   cout << d << "\n";

   c = -a;
   cout << c << "\n";

   c = sqrt(a);
   cout << c << "\n";

   power(c, to_quad_float(10), 20);
   cout << c << "\n";

   {

   long n, n1;
   int shamt = min(NTL_BITS_PER_LONG,2*NTL_DOUBLE_PRECISION);

   n = (1L << (shamt-1)) - 1;
   c = to_quad_float(n);
   n1 = to_long(c);

   if (n1 == n)
      cout << "long conversion OK\n";
   else
      cout << "long conversion not OK\n";

   n = (1L << (shamt-1));
   c = to_quad_float(n);
   n1 = to_long(c);

   if (n1 == n)
      cout << "long conversion OK\n";
   else
      cout << "long conversion not OK\n";

   }

   {

   int n, n1;
   int shamt = min(NTL_BITS_PER_INT,2*NTL_DOUBLE_PRECISION);

   n = (1 << (shamt-1)) - 1;
   c = to_quad_float(n);
   n1 = to_long(c);

   if (n1 == n)
      cout << "int conversion OK\n";
   else
      cout << "int conversion not OK\n";

   n = (1 << (shamt-1));
   c = to_quad_float(n);
   n1 = to_long(c);

   if (n1 == n)
      cout << "int conversion OK\n";
   else
      cout << "int conversion not OK\n";

   }

}
