
#include <NTL/ZZ_pX.h>





main()
{


   cerr << "This is NTL version " << NTL_VERSION << "\n"; 
   cerr << "configuration flags: ";

#ifdef NTL_LONG_LONG
   cerr << "NTL_LONG_LONG ";
#endif

#ifdef NTL_LONG_LONG_TYPE
   cerr << "NTL_LONG_LONG_TYPE ";
#endif

#ifdef NTL_CPLUSPLUS_ONLY
   cerr << "NTL_CPLUSPLUS_ONLY ";
#endif


#ifdef NTL_X86_FIX
   cerr << "NTL_X86_FIX ";
#endif

#ifdef NTL_NO_X86_FIX
   cerr << "NTL_NO_X86_FIX ";
#endif

#ifdef NTL_AVOID_FLOAT
   cerr << "NTL_AVOID_FLOAT ";
#endif

#ifdef NTL_AVOID_BRANCHING
   cerr << "NTL_AVOID_BRANCHING ";
#endif

#ifdef NTL_FFT_PIPELINE
   cerr << "NTL_FFT_PIPELINE ";
#endif

#ifdef NTL_SINGLE_MUL
   cerr << "NTL_SINGLE_MUL ";
#endif

#ifdef NTL_FAST_INT_MUL
   cerr << "NTL_FAST_INT_MUL ";
#endif

#ifdef NTL_TBL_REM
   cerr << "NTL_TBL_REM ";
#endif

#ifdef NTL_NO_INIT_TRANS
   cerr << "NTL_NO_INIT_TRANS ";
#endif

#ifdef NTL_RANGE_CHECK
   cerr << "NTL_RANGE_CHECK ";
#endif


   cerr << "\n";


   long n, k;

   n = 200;
   k = 200;

   ZZ p;

   GenPrime(p, k);


   ZZ_p::init(p);         // initialization

   ZZ_pX f, g, h, r1, r2, r3;

   random(g, n);    // g = random polynomial of degree < n
   random(h, n);    // h =             "   "
   random(f, n);    // f =             "   "

   // SetCoeff(f, n);  // Sets coefficient of X^n to 1
   
   ZZ_p lc;

   do {
      random(lc);
   } while (IsZero(lc));

   SetCoeff(f, n, lc);


   // For doing arithmetic mod f quickly, one must pre-compute
   // some information.

   ZZ_pXModulus F;
   build(F, f);

   PlainMul(r1, g, h);  // this uses classical arithmetic
   PlainRem(r1, r1, f);

   MulMod(r2, g, h, F);  // this uses the FFT

   MulMod(r3, g, h, f);  // uses FFT, but slower

   // compare the results...

   if (r1 != r2) {
      cerr << "r1 != r2!!\n";
      return 1;
   }
   else if (r1 != r3) {
      cerr << "r1 != r3!!\n";
      return 1;
   }

   cerr << "test is OK\n";

   ZZ x1, x2, x3, x4;
   double t;
   long i;

   RandomLen(x1, 512);
   RandomBnd(x2, x1);
   RandomBnd(x3, x1);

   mul(x4, x2, x3);

   t = GetTime();
   for (i = 0; i < 100000; i++)
      mul(x4, x2, x3);
   t = GetTime()-t;

   cerr << "time for 512-bit mul: " << t*10 << "us\n";

   t = GetTime();
   for (i = 0; i < 100000; i++)
      rem(x2, x4, x1);
   t = GetTime()-t;

   cerr << "time for 1024/512-bit rem: " << t*10 << "us\n";

   GenPrime(p, 512);
   RandomBnd(x1, p);
   if (IsZero(x1)) set(x1);

   t = GetTime();
   for (i = 0; i < 10000; i++)
      InvMod(x2, x1, p);
   t = GetTime()-t;

   cerr << "time for 512-bit modular inverse: " << t*100 << "us\n";


   // test modulus switching
   
   n = 1024;
   k = 1024;
   RandomLen(p, k);

   ZZ_p::init(p);

   ZZ_pX j1, j2, j3;

   random(j1, n);
   random(j2, n);

   t = GetTime();

   mul(j3, j1, j2);

   t = GetTime()-t;

   cerr << "time to multiply degree 1023 polynomials mod a 1024-bit number: ";
   cerr << t << "s\n";

   return 0;
}
