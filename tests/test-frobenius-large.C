#include "linbox/linbox-config.h"

#include <algorithm>
#include <iostream>
#include <vector>

#include "linbox/ring/modular.h"
#include "linbox/util/commentator.h"
#include "linbox/ring/ntl.h"

#include "linbox/algorithms/frobenius-large.h"

#include "test-frobenius.h"

using namespace LinBox;

int main(int argc, char** argv) {
	uint64_t p = 10007;
	uint64_t e = 1;
	size_t k = 0;
	int seed = time(NULL);

	static Argument args[] = {
		{ 'k', "-k K", "Number of invariant factors to compute (default k=0 for all)", TYPE_INT, &k},
		{ 'p', "-p P", "Set characteristic of field GF(p^e)", TYPE_INT, &p},
		{ 'e', "-e E", "Set extension degree of field GF(p^e)", TYPE_INT, &e},
		{ 'r', "-r R", "Random seed", TYPE_INT, &seed},
		END_OF_ARGUMENTS
	};

	parseArguments(argc,argv,args);

	srand(seed);

   //typedef Givaro::Modular<double> Field;
   typedef NTL_zz_p Field;
   Field F(p,e);

   typedef NTL_zz_pX PolyRing;
	PolyRing R(F);
	
	FrobeniusLarge<PolyRing> FLD(R);
   bool pass = testFrobenius(FLD, F, R, k);
   return pass ? 0 : -1;
}
