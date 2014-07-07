
#include "linbox/linbox-config.h"

#include <omp.h>

#define LINBOX_USES_OMP 1
#include "linbox/field/modular.h"
#include "linbox/matrix/sparse-matrix.h"
#include "linbox/matrix/dense-matrix.h"
#include "linbox/matrix/matrix-domain.h"

#include "linbox/algorithms/poly-interpolation.h"
#include "linbox/algorithms/poly-det.h"

#include <givaro/givpoly1.h>
#include <givaro/givzpz.h>

using namespace LinBox;

typedef Modular<double> Field;
typedef typename Field::Element FieldElt;
typedef Givaro::Poly1Dom<Field,Givaro::Dense> PolyDom;
typedef GivaroPoly<PolyDom> Ring;
typedef typename Ring::Element RingElt;
typedef MatrixDomain<PolyDom> PolyMatDom;
typedef PolyMatDom::OwnMatrix PolyMat;



int main(int argc, char** argv)
{
	static Argument args[] = {
		END_OF_ARGUMENTS
	};
	parseArguments (argc, argv, args);


	bool pass=true;

	std::vector<integer> x,y;
	x.push_back(1);
	x.push_back(2);
	x.push_back(3);
	x.push_back(4);
	x.push_back(5);
	y.push_back(0);
	y.push_back(2);
	
	int q=101;
	RingElt P1,P2,P4;
	Field F(q);
	PolyDom PD(F,"x");
	Ring R(PD);
	R.init(P1,x);
	R.init(P4,y);


	std::vector<FieldElt> vals,pts;
	pts.push_back(1);
	pts.push_back(2);
	pts.push_back(3);
	pts.push_back(4);
	pts.push_back(5);
	pts.push_back(6);
	pts.push_back(7);
	pts.push_back(8);


	FieldElt d;
	PolyInterpolation<Field,PolyDom> PO;
	PolyInterpolation<Field,PolyDom>::ProductTree mtree;
	PO.productTree(mtree,pts,PD);
	PO.evaluate(vals,P1,mtree,PD,F);
	for (int i=0;i<vals.size();++i) {
		PD.eval(d,P1,pts[i]);
		pass=pass&&F.areEqual(d,vals[i]);
	}
	PO.interpolate(P2,pts,vals,PD,F);

	int n=3,m=3;
	typename MatrixDomain<PolyDom>::OwnMatrix A(PD,m,n);
	for (int i=0;i<m;++i) {
		for (int j=0;j<n;++j) {
			if (i==j) {
				A.setEntry(i,j,P4);
			} else {
				A.setEntry(i,j,R.zero);
			}
		}
	}
	A.setEntry(0,1,P4);


	//A.setEntry(0,0,P1);
	for (int i=0;i<m;++i) {
		for (int j=0;j<n;++j) {
			PolyDom::Element tempP;
			A.getEntry(tempP,i,j);
			PD.write(std::cout,tempP);
			std::cout << " ";
		}
		std::cout << std::endl;
	}
	PolyDom::Element P3;
	computePolyDetExtension(P3,F,A);
	R.write(std::cout,P3);
	std::cout << std::endl;




	R.write(std::cout,P1);
	std::cout << std::endl;
	R.write(std::cout,P2);
	std::cout << std::endl;

	return pass?0:-1;
}

