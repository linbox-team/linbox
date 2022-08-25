#include "linbox/linbox-config.h"

#include <algorithm>
#include <iostream>

#include "linbox/util/commentator.h"

#include "linbox/matrix/sparse-matrix.h"

#include "givaro/givtimer.h"

using namespace LinBox;

// This does a very basic 8x8 case.  A fuller suite of test matrices would be good.
template<class FrobeniusObject, class Field, class PolyRing>
bool testFrobenius(FrobeniusObject& FO, Field & F, PolyRing & R, size_t k=0) {
   typedef typename PolyRing::Element Polynomial;
   typedef typename PolyRing::Coeff Coeff;
   typedef SparseMatrix<Field, SparseMatrixFormat::CSR> SparseMat;

	
   std::ostream &report = commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);

   //typedef typename FrobeniusObject::PolynomialRing::CoeffField Field;
  //typedef typename FrobeniusObject::Field Field;

   SparseMat M(F);                                           
   M.resize(8,8);
   M.setEntry(1,0,F.one);
   M.setEntry(2,1,F.one);
   M.setEntry(3,3,F.one);
   M.setEntry(5,4,F.one);
   M.setEntry(7,6,F.one);
	M.finalize(); // invariants: (x^3)*(x-1), x^2, x^2.

   Polynomial g(3); // x^2
   R.setCoeff(g, 0, (Coeff)0);
   R.setCoeff(g, 1, (Coeff)0);
   R.setCoeff(g, 2, (Coeff)1);

   Polynomial f(5); // (x^3)(x-1)
   R.setCoeff(f, 0, (Coeff)0);
   R.setCoeff(f, 1, (Coeff)0);
   R.setCoeff(f, 2, (Coeff)0);
   R.setCoeff(f, 3, (Coeff)-1);
   R.setCoeff(f, 4, (Coeff)1);

	std::vector<Polynomial> fl(3);
   fl[0] = f;
   fl[1] = g;
   fl[2] = g;
	
	std::vector<Polynomial> fs;

   Givaro::Timer T;

	T.clear(); T.start();
	FO.frobeniusInvariants(fs, M, k);
	T.stop();

	report << "usertime: " << T.usertime() << std::endl;
	
   bool valid = true;
   if (fl.size() > fs.size())
      valid = false;
   else
      for (size_t i = 0; i < fl.size(); ++i) 
         if (fl[i] != fs[i]) valid = false;
	
   report << "computed\n";
	std::for_each(fs.begin(), fs.end(), [&](const Polynomial &v) {
		R.write(report, v) << std::endl;
	});
   if (not valid){
      report << "expected\n";
   std::for_each(fl.begin(), fl.end(), [&](const Polynomial &v) {
      R.write(report, v) << std::endl;
	   });
   }

	return valid;
}
