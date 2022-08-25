#include "linbox/linbox-config.h"

#include <algorithm>
#include <iostream>
#include <vector>

#include "linbox/ring/modular.h"
#include "linbox/util/commentator.h"
#include "linbox/ring/ntl.h"

#include "linbox/algorithms/frobenius-small.h"

#include "test-frobenius.h"

using namespace LinBox;

int main(int argc, char** argv) {
	uint64_t p = 3;
	size_t k = 0;
	int seed = time(NULL);

	static Argument args[] = {
		{ 'k', "-k K", "Number of invariant factors to compute (default k=0 for all)", TYPE_INT, &k},
		{ 'p', "-p P", "Set the field GF(p)", TYPE_INT, &p},
		{ 'r', "-r R", "Random seed", TYPE_INT, &seed},
		END_OF_ARGUMENTS
	};

	parseArguments(argc,argv,args);

	srand(seed);

   typedef Givaro::Modular<double> Field;
   Field F(p);

   typedef NTL_zz_pX PolyRing;
	PolyRing R(p);
	
	FrobeniusSmall<Field, PolyRing> FSD(F, R);
   bool pass = testFrobenius(FSD, F, R, k);
   return pass ? 0 : -1;
}
