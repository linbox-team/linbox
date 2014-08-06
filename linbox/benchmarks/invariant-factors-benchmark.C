#include "linbox/linbox-config.h"

#include <iostream>
#include <vector>
#include <omp.h>

#define __LINBOX_USE_OPENMP 1
#include "linbox/field/modular.h"
#include "linbox/matrix/sparse-matrix.h"
#include "linbox/matrix/dense-matrix.h"
#include "linbox/matrix/matrix-domain.h"

#include "linbox/algorithms/blackbox-block-container.h"
#include "linbox/algorithms/block-coppersmith-domain.h"

#include "linbox/solutions/det.h"

#include <givaro/givzpz.h>
#include <givaro/givpoly1.h>
#include <linbox/ring/givaro-poly.h>
#include <linbox/algorithms/smith-form-kannan-bachem.h>

// Computes the invariant factors of a sparse matrix (given in Matrix Market Format)
// Effectively times: TPL_omp, BlockCoppersmithDomain and KannanBachem

using namespace LinBox;

typedef Modular<double> Field;
typedef typename Field::Element Element;
typedef SparseMatrix<Field, SparseMatrixFormat::TPL_omp> SparseMat;
typedef MatrixDomain<Field> Domain;
typedef typename Domain::OwnMatrix Block;


typedef Givaro::ZpzDom<Givaro::Std32> BaseDom;
typedef Givaro::Poly1Dom<BaseDom,Givaro::Dense> PolyDom;
typedef GivaroPoly<PolyDom> Ring;
typedef MatrixDomain<Ring> PolyMatDom;
typedef BlasMatrix<Ring> PolyMat;


int randRange(int start, int end)
{
        double rval = rand();
        static const double NORMALIZING_CONSTANT = 1.0/(1.0+RAND_MAX);
        double normedRVal = rval*NORMALIZING_CONSTANT;
        double rangeSize = end-start;
        int offset = rangeSize*normedRVal;
        return start+offset;
}


void benchmarkBCD(Field& F,
                  Domain& MD,
                  SparseMat& M,
                  Block& U,
                  Block& V,
                  std::vector<Block>& gen,
                  std::vector<size_t>& deg,
                  int t,
                  int p)
{
	BlackboxBlockContainer<Field,SparseMat> blockseq(&M,F,U,V);
	BlockCoppersmithDomain<Domain,BlackboxBlockContainer<Field,SparseMat> >
		BCD(MD,&blockseq,t);

	double start=omp_get_wtime();
	deg=BCD.right_minpoly(gen);

	BaseDom BD(p);
	PolyDom PD(BD, "x");
	Ring R(PD);
	PolyMatDom PMD(R);
	int b=U.rowdim();
	int d=gen.size();
	PolyMat MM(R,b,b);
	Ring::Element temp;
	temp.resize(d);
	for (int i=0;i<b;++i) {
		for (int j=0;j<b;++j) {
			for (int k=0;k<d;++k) {
				temp[k]=gen[k].getEntry(i,j);
			}
			MM.setEntry(i,j,temp);
		}
	}

	SmithFormKannanBachemDomain<PolyMatDom> SFKB(PMD);
	BlasVector<Ring> diag(R,b,R.zero);
	SFKB.solve(diag,MM);
	PolyDom::Type_t lcoef;
	for (size_t i=0;i<diag.size();++i) {
		PD.leadcoef(lcoef, diag[i]);
		PD.divin(diag[i],lcoef);
	}
	double time=omp_get_wtime()-start;
	std::cout << time << std::endl;
	for (size_t i=0;i<diag.size();++i) {
		R.write(std::cout,diag[i]);
		std::cout << std::endl;
	}

	testCorrectness(F,p,M,diag[diag.size()-1],PD);
}

bool testCorrectness(Field& F,
                     int p,
                     SparseMat& M,
                     Ring& minPoly
                     PolyDom& PD)
{
	int n=M.rowdim();
	Element d,e;
	F.init(d,0);
	BlasVector<Element> v(F,n,d),sum(F,n,d),vPrime(F,n,d);
	for (int i=0;i<n;++i) {
		F.init(v[i],randRange(0,p));
	}
	for (int i=0;i<=PD.degree(minPoly);++i) {
		PD.getEntry(d,i,minPoly);
		for (int j=0;j<n;++j) {
			F.mul(e,d,v[j]);
			F.addin(sum[j],e);
		}
		M.apply(vPrime,v);
		for (int j=0;j<n;++j) {
			v[j]=vPrime[j]
		}
	}
	for (int i=0;i<n;++i) {
		if (!F.isZero(v[i])) {
			return false;
		}
	}
	return true;
}

int main(int argc, char** argv)
{
	int earlyTerm;
	int p;
	std::string uFname,vFname,mFname;

	static Argument args[] = {
		{ 'p', "-p P", "Set the field GF(p)", TYPE_INT, &p},
		{ 't', "-t T", "Early term threshold", TYPE_INT, &earlyTerm},
		{ 'm', "-m M", "Name of file for matrix M", TYPE_STR, &mFname},
		{ 'u', "-u U", "Name of file for matrix U", TYPE_STR, &uFname},
		{ 'v', "-v V", "Name of file for matrix V", TYPE_STR, &vFname},
		END_OF_ARGUMENTS
	};

	parseArguments(argc,argv,args);

	Field F(p);
	Domain MD(F);
	SparseMat M(F);
	Block U(F),V(F);

	{
		ifstream iF(mFname);
		M.read(iF);
		M.finalize();
		iF.close();
	}
	{
		ifstream iF(uFname);
		U.read(iF);
		iF.close();
	}
	{
		ifstream iF(vFname);
		V.read(iF);
		iF.close();
	}

	std::vector<Block> gen;
	std::vector<size_t> deg;

	benchmarkBCD(F,MD,M,U,V,gen,deg,earlyTerm,p);

	return 0;
}
