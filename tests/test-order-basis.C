#include <iostream>
#include "linbox/matrix/polynomial-matrix.h"
#include "linbox/randiter/random-fftprime.h"
#include "linbox/randiter/random-prime.h"
#include "linbox/ring/modular.h"

#include "linbox/algorithms/polynomial-matrix/polynomial-matrix-domain.h"
#include "linbox/algorithms/polynomial-matrix/order-basis.h"
#include "linbox/algorithms/block-coppersmith-domain.h"

using namespace LinBox;
using namespace std;


//ostream& report = commentator().report();
//ostream& report = std::cout;

template<typename Field, typename Mat>
bool check_sigma(const Field& F, const Mat& sigma,  Mat& serie, size_t ord, string& msg){
	ostream &report = commentator().report ();//Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	Mat T(F,sigma.rowdim(),serie.coldim(),sigma.size()+serie.size()-1);
	PolynomialMatrixMulDomain<Field> PMD(F);
	PMD.mul(T,sigma,serie);
	MatrixDomain<Field> MD(F);
	size_t i=0;
	msg= string(".....");
	bool nul_sigma=true;
	while(i<ord && MD.isZero(T[i])){
		if (!MD.isZero(sigma[i])) nul_sigma=false;		
		i++;
	}
	if (i<ord){
		report<<"error at degree="<<i<<endl;
		T[i].write(report, Tag::FileFormat::Plain);
		report<<"***"<<endl;
		report<<serie<<endl;
		report<<sigma<<endl;	
	}
	
	
	if (i==ord && !nul_sigma)
		{msg+="done"; return true;}
	else
		{msg+="error"; return false;}
}

template<typename MatPol>
bool operator==(const MatPol& A, const MatPol& B){
	ostream &report = commentator().report ();//Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	MatrixDomain<typename MatPol::Field> MD(A.field());
	if (A.real_degree()!=B.real_degree()|| A.rowdim()!= B.rowdim() || A.coldim()!=B.coldim()){
		report<<A.size()<<"("<<A.rowdim()<<"x"<<A.coldim()<<") <> "
		    <<B.size()<<"("<<B.rowdim()<<"x"<<B.coldim()<<") <> "<<endl;
		return false;
	}
	size_t i=0;
	while (i<=A.real_degree() && MD.areEqual(A[i],B[i]))
		i++;

	if (i<=A.real_degree() && A.rowdim()<10 && A.coldim()<10){
		report<<"first:"<<endl<<A<<endl;
		report<<"second:"<<endl<<B<<endl;
	}

	return i>A.real_degree();
}
 

template<typename Field, typename RandIter>
bool check_sigma(const Field& F, RandIter& Gen, size_t m, size_t n, size_t d) {
	ostream &report = commentator().report ();//Commentator::LEVEL_ALWAYS, INTERNAL_DESCRIPTION);
	//typedef typename Field::Element Element;
	typedef PolynomialMatrix<Field, PMType::matfirst> MatrixP;
	//typedef PolynomialMatrix<Field, PMType::polfirst> MatrixP;
	MatrixP Serie(F, m, n,  d);
	MatrixP Sigma1(F, m, m, d+1),Sigma2(F, m, m, d+1),Sigma3(F, m, m, d+1);

	// set the Serie at random
	for (size_t k=0;k<d;++k)
		for (size_t i=0;i<m;++i)
			for (size_t j=0;j<n;++j)
				Gen.random(Serie.ref(i,j,k));

	//report<<"Serie:="<<Serie<<std::endl;
	
	// define the shift
	vector<size_t> shift(m,0);
	vector<size_t> shift2(shift),shift3(shift);

	OrderBasis<Field> SB(F);
    bool passed(true); string msg;
    // MBasis check
	SB.M_Basis(Sigma3, Serie, d, shift3);
    passed&=check_sigma(F,Sigma3,Serie,d, msg);    
	report << "M-Basis       : " <<msg<<endl;
    // PMBasis check
	SB.PM_Basis(Sigma1,Serie, d, shift);
    passed&=check_sigma(F,Sigma1,Serie,d, msg);    
	report << "PM-Basis      : " <<msg<<endl;

    // PMBasis online check
	// SB.oPM_Basis(Sigma2, Serie, d, shift2);
    // passed&=check_sigma(F,Sigma2,Serie,d, msg);    
	// report << "PM-Basis iter : " <<msg<<endl;

	report<<endl;
    return passed;
}

bool runTest(uint64_t m,uint64_t n, uint64_t d, long seed){

    commentator().start ("Testing order basis computation", "testOrderBasis", 1);
    
	bool ok=true,passed;
    size_t bits= (53-integer(n).bitsize())/2;

    typedef Givaro::Modular<double>              SmallField;	
	typedef Givaro::Modular<Givaro::Integer>      LargeField1;
    typedef Givaro::Modular<RecInt::ruint128,RecInt::ruint256> LargeField2;

	ostream &report = commentator().report ();//Commentator::LEVEL_ALWAYS, INTERNAL_DESCRIPTION);
	report<<"###  matrix series is of size "<<m<<" x "<<n<<" of degree "<<d<<std::endl;
    
    
	// fourier prime < 2^(53--log(n))/2
	{
		integer p;
		RandomFFTPrime::seeding (seed);
		if (!RandomFFTPrime::randomPrime (p, 1<<bits, integer(d).bitsize()+1))
			throw LinboxError ("RandomFFTPrime::randomPrime failed");
		SmallField  F((int32_t)p);
        typename SmallField::RandIter G(F,seed);
        report<<"   - checking with small FFT prime p="<<p<<endl;
        ok&=passed=check_sigma (F,G,m,n,d);
        report<<"   ---> "<<(passed?"done":"error")<<std::endl<<std::endl;
		
	}
	// normal prime < 2^(53--log(n))/2
	{
		typedef Givaro::Modular<double> Field;
		PrimeIterator<IteratorCategories::HeuristicTag> Rd(FieldTraits<Field>::bestBitSize(n),seed);
		integer p;
		p=*Rd;
        SmallField  F((int32_t)p);
        typename SmallField::RandIter G(F,seed);
        report<<"   - checking with small generic prime p="<<p<<std::endl;
		ok&=passed=check_sigma (F,G,m,n,d);
        report<<"   ---> "<<(passed?"done":"error")<<std::endl<<std::endl;
	}

	// multi-precision prime
	 {
	 	size_t bits=114;
	 	PrimeIterator<IteratorCategories::HeuristicTag> Rd(bits,seed);
	 	integer p= *Rd;
        // Modular<integer>
        LargeField1  F1(p);
        typename LargeField1::RandIter G1(F1,seed);
        report<<"   - checking with multiprecision prime (Modular<integer>) p="<<p<<std::endl;
		ok&=passed=check_sigma (F1,G1,m,n,d);
        report<<"   ---> "<<(passed?"done":"error")<<std::endl<<std::endl;
        // Modular<recint<128 ,256 >>
        LargeField2  F2(p);
        typename LargeField2::RandIter G2(F2,seed);
        report<<"   - checking with multiprecision prime (Modular<recint>) p="<<p<<std::endl;
		ok&=passed=check_sigma (F2,G2,m,n,d);
        report<<"   ---> "<<(passed?"done":"error")<<std::endl<<std::endl;
	 }

     commentator().stop (MSG_STATUS (ok), (const char *) 0, "testOrderBasis"); 
     return ok;
}

int main(int argc, char** argv){
	static size_t  m = 8; // matrix dimension
	static size_t  n = 4; // matrix dimension
    //	static size_t  b = 20; // entries bitsize
	static size_t  d = 80;  // matrix degree
	static long    seed = time(NULL);

	static Argument args[] = {
		{ 'm', "-m M", "Set row dimension of matrix series to M.", TYPE_INT,     &m },
		{ 'n', "-n N", "Set column dimension of matrix series to N.", TYPE_INT,     &n },
		{ 'd', "-d D", "Set degree of  matrix series to D.", TYPE_INT,     &d },
        //		{ 'b', "-b B", "Set bitsize of the matrix entries", TYPE_INT, &b },
		{ 's', "-s s", "Set the random seed to a specific value", TYPE_INT, &seed},
		END_OF_ARGUMENTS
	};
	
	parseArguments (argc, argv, args);

	return (runTest(m,n,d,seed)?0:-1);
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
