//#define __FFLASFFPACK_SEQUENTIAL

#include <iostream>
#include <iomanip>
size_t getPeakRSS( );
size_t getCurrentRSS( );
//#define MEMINFO std::right<<std::setw(20)<<"                     ---->   Max Mem: "<<getPeakRSS()/1000000.<<"Mo"
#define MB(x) ((x)/(double)(1<<20))
//#define MB(x) ((x)>>20)
#define MEMINFO std::right<<" ---->   Mem: "<<MB(getCurrentRSS())<<" Mo  (Max: "<<MB(getPeakRSS())<<" Mo)"  
#include "linbox/matrix/polynomial-matrix.h"
#include "linbox/randiter/random-fftprime.h"
#include "linbox/randiter/random-prime.h"
#include "linbox/ring/modular.h"
#include "linbox/algorithms/polynomial-matrix/polynomial-matrix-domain.h"
#include "linbox/algorithms/polynomial-matrix/order-basis.h"
#include "linbox/algorithms/block-coppersmith-domain.h"

/* MEMORY INFO */
#if defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
#include <unistd.h>
#include <sys/types.h>
#include <sys/param.h>
#if defined(BSD)
#include <sys/sysctl.h>
#endif
#else
#error "Unable to define getMemorySize( ) for an unknown OS."
#endif

#if defined(_WIN32)
#include <windows.h>
#include <psapi.h>

#elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
#include <unistd.h>
#include <sys/resource.h>

#if defined(__APPLE__) && defined(__MACH__)
#include <mach/mach.h>

#elif (defined(_AIX) || defined(__TOS__AIX__)) || (defined(__sun__) || defined(__sun) || defined(sun) && (defined(__SVR4) || defined(__svr4__)))
#include <fcntl.h>
#include <procfs.h>

#elif defined(__linux__) || defined(__linux) || defined(linux) || defined(__gnu_linux__)
#include <stdio.h>

#endif

#else
#error "Cannot define getPeakRSS( ) or getCurrentRSS( ) for an unknown OS."
#endif
/* END MEMORY INFO */

size_t getPeakRSS( )
{
#if defined(_WIN32)
	/* Windows -------------------------------------------------- */
	PROCESS_MEMORY_COUNTERS info;
	GetProcessMemoryInfo( GetCurrentProcess( ), &info, sizeof(info) );
	return (size_t)info.PeakWorkingSetSize;

#elif (defined(_AIX) || defined(__TOS__AIX__)) || (defined(__sun__) || defined(__sun) || defined(sun) && (defined(__SVR4) || defined(__svr4__)))
	/* AIX and Solaris ------------------------------------------ */
	struct psinfo psinfo;
	int fd = -1;
	if ( (fd = open( "/proc/self/psinfo", O_RDONLY )) == -1 )
		return (size_t)0L;		/* Can't open? */
	if ( read( fd, &psinfo, sizeof(psinfo) ) != sizeof(psinfo) )
		{
			close( fd );
			return (size_t)0L;		/* Can't read? */
		}
	close( fd );
	return (size_t)(psinfo.pr_rssize * 1024L);

#elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
	/* BSD, Linux, and OSX -------------------------------------- */
	struct rusage rusage;
	getrusage( RUSAGE_SELF, &rusage );
#if defined(__APPLE__) && defined(__MACH__)
	return (size_t)rusage.ru_maxrss;
#else
	return (size_t)(rusage.ru_maxrss * 1024L);
#endif

#else
	/* Unknown OS ----------------------------------------------- */
	return (size_t)0L;			/* Unsupported. */
#endif
}
/**
 * Returns the current resident set size (physical memory use) measured
 * in bytes, or zero if the value cannot be determined on this OS.
 */
size_t getCurrentRSS( )
{
#if defined(_WIN32)
	/* Windows -------------------------------------------------- */
	PROCESS_MEMORY_COUNTERS info;
	GetProcessMemoryInfo( GetCurrentProcess( ), &info, sizeof(info) );
	return (size_t)info.WorkingSetSize;

#elif defined(__APPLE__) && defined(__MACH__)
	/* OSX ------------------------------------------------------ */
	struct mach_task_basic_info info;
	mach_msg_type_number_t infoCount = MACH_TASK_BASIC_INFO_COUNT;
	if ( task_info( mach_task_self( ), MACH_TASK_BASIC_INFO,
		(task_info_t)&info, &infoCount ) != KERN_SUCCESS )
		return (size_t)0L;		/* Can't access? */
	return (size_t)info.resident_size;

#elif defined(__linux__) || defined(__linux) || defined(linux) || defined(__gnu_linux__)
	/* Linux ---------------------------------------------------- */
	long rss = 0L;
	FILE* fp = NULL;
	if ( (fp = fopen( "/proc/self/statm", "r" )) == NULL )
		return (size_t)0L;		/* Can't open? */
	if ( fscanf( fp, "%*s%ld", &rss ) != 1 )
	{
		fclose( fp );
		return (size_t)0L;		/* Can't read? */
	}
	fclose( fp );
	return (size_t)rss * (size_t)sysconf( _SC_PAGESIZE);

#else
	/* AIX, BSD, Solaris, and Unknown OS ------------------------ */
	return (size_t)0L;			/* Unsupported. */
#endif
}

/**
 * Returns the size of physical memory (RAM) in bytes.
 */
size_t getMemorySize( )
{
#if defined(_WIN32) && (defined(__CYGWIN__) || defined(__CYGWIN32__))
	/* Cygwin under Windows. ------------------------------------ */
	/* New 64-bit MEMORYSTATUSEX isn't available.  Use old 32.bit */
	MEMORYSTATUS status;
	status.dwLength = sizeof(status);
	GlobalMemoryStatus( &status );
	return (size_t)status.dwTotalPhys;

#elif defined(_WIN32)
	/* Windows. ------------------------------------------------- */
	/* Use new 64-bit MEMORYSTATUSEX, not old 32-bit MEMORYSTATUS */
	MEMORYSTATUSEX status;
	status.dwLength = sizeof(status);
	GlobalMemoryStatusEx( &status );
	return (size_t)status.ullTotalPhys;

#elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
	/* UNIX variants. ------------------------------------------- */
	/* Prefer sysctl() over sysconf() except sysctl() HW_REALMEM and HW_PHYSMEM */

#if defined(CTL_HW) && (defined(HW_MEMSIZE) || defined(HW_PHYSMEM64))
	int mib[2];
	mib[0] = CTL_HW;
#if defined(HW_MEMSIZE)
	mib[1] = HW_MEMSIZE;            /* OSX. --------------------- */
#elif defined(HW_PHYSMEM64)
	mib[1] = HW_PHYSMEM64;          /* NetBSD, OpenBSD. --------- */
#endif
	int64_t size = 0;               /* 64-bit */
	size_t len = sizeof( size );
	if ( sysctl( mib, 2, &size, &len, NULL, 0 ) == 0 )
		return (size_t)size;
	return 0L;			/* Failed? */

#elif defined(_SC_AIX_REALMEM)
	/* AIX. ----------------------------------------------------- */
	return (size_t)sysconf( _SC_AIX_REALMEM ) * (size_t)1024L;

#elif defined(_SC_PHYS_PAGES) && defined(_SC_PAGESIZE)
	/* FreeBSD, Linux, OpenBSD, and Solaris. -------------------- */
	return (size_t)sysconf( _SC_PHYS_PAGES ) *
		(size_t)sysconf( _SC_PAGESIZE );

#elif defined(_SC_PHYS_PAGES) && defined(_SC_PAGE_SIZE)
	/* Legacy. -------------------------------------------------- */
	return (size_t)sysconf( _SC_PHYS_PAGES ) *
		(size_t)sysconf( _SC_PAGE_SIZE );

#elif defined(CTL_HW) && (defined(HW_PHYSMEM) || defined(HW_REALMEM))
	/* DragonFly BSD, FreeBSD, NetBSD, OpenBSD, and OSX. -------- */
	int mib[2];
	mib[0] = CTL_HW;
#if defined(HW_REALMEM)
	mib[1] = HW_REALMEM;		/* FreeBSD. ----------------- */
#elif defined(HW_PYSMEM)
	mib[1] = HW_PHYSMEM;		/* Others. ------------------ */
#endif
	unsigned int size = 0;		/* 32-bit */
	size_t len = sizeof( size );
	if ( sysctl( mib, 2, &size, &len, NULL, 0 ) == 0 )
		return (size_t)size;
	return 0L;			/* Failed? */
#endif /* sysctl and sysconf variants */

#else
	return 0L;			/* Unknown OS. */
#endif
}

using namespace LinBox;
using namespace std;

template<typename Field, typename Mat>
string check_sigma(const Field& F, const Mat& sigma,  Mat& serie, size_t ord){
	Mat T(F,sigma.rowdim(),serie.coldim(),sigma.size()+serie.size()-1);
	PolynomialMatrixMulDomain<Field> PMD(F);
	PMD.mul(T,sigma,serie);
	MatrixDomain<Field> MD(F);
	size_t i=0;
	string msg(".....");
	bool nul_sigma=true;
	while(i<ord && MD.isZero(T[i])){
		if (!MD.isZero(sigma[i])) nul_sigma=false;		
		i++;
	}
	if (i<ord){
		cout<<"error at degree="<<i<<endl;
		T[i].write(std::cout, Tag::FileFormat::Plain);
		cout<<"***"<<endl;
		cout<<serie<<endl;
		cout<<sigma<<endl;	
	}
	
	
	if (i==ord && !nul_sigma)
		msg+="done";
	else
		msg+="error";
	return msg;
}

template<typename MatPol>
bool operator==(const MatPol& A, const MatPol& B){
	MatrixDomain<typename MatPol::Field> MD(A.field());
	if (A.size()!=B.size()|| A.rowdim()!= B.rowdim() || A.coldim()!=B.coldim()){
		cout<<A.size()<<"("<<A.rowdim()<<"x"<<A.coldim()<<") <> "
		    <<B.size()<<"("<<B.rowdim()<<"x"<<B.coldim()<<") <> "<<endl;
		return false;
	}
	size_t i=0;
	while (i<A.size() && MD.areEqual(A[i],B[i]))
		i++;

	if (i<A.size() && A.rowdim()<10 && A.coldim()<10){
		cout<<"first:"<<endl<<A<<endl;
		cout<<"second:"<<endl<<B<<endl;
	}

	return i==A.size();
}
 

template<typename Field, typename RandIter>
void bench_sigma(const Field& F,  RandIter& Gen, size_t m, size_t n, size_t d, string target) {
	//typedef typename Field::Element Element;
	//typedef PolynomialMatrix<Field,PMType::matfirst> MatrixP;
	typedef PolynomialMatrix<Field,PMType::polfirst> MatrixP;
	std::cout<<"Order Basis computation over ";F.write(cout)<<endl;
	integer p;
	F.characteristic(p);
	size_t memp=length(p)+(p.bitsize()>=64?16:0);
	//size_t data_in=3*m*n*d*memp;
	//size_t data_out=2*m*m*(d+1)*memp;
	//size_t data_comp= 2*m*m*d*(length(uint64_t(m*d)*p*p)+(p.bitsize()>26?8:0));
	std::cout<<"**************************"<<std::endl;
	std::cout<<"mem(p)        : "<<memp<<std::endl;
	//std::cout<<"Projected Memory : "<< MB(data_in+data_out+data_comp)<<"Mo"<<std::endl;
	std::cout<<"Available memory : "<<MB(getMemorySize())<<std::endl;
	std::cout<<"**************************"<<std::endl;
	std::cout<<"**************************"<<std::endl<<std::endl<<std::endl;
	std::cout<<"[begin ] : "<<MEMINFO<<std::endl; 

	
	MatrixP *Serie = new MatrixP(F, m, n, d);	
	// set the Serie at random
	for (size_t k=0;k<d;++k)
		for (size_t i=0;i<m;++i)
			for (size_t j=0;j<n;++j)
				Gen.random(Serie->ref(i,j,k));
	std::cout<<"[initial sequence] : "<<MB(Serie->realmeminfo())<<"Mo"<<MEMINFO<<std::endl;
			
	// define the shift
	vector<size_t> shift(m,0);
	
	OrderBasis<Field> SB(F);
	Timer chrono;
#ifdef BENCH_MBASIS
	if (target=="ALL"){
		MatrixP Sigma1(F, m, m, d+1);
		vector<size_t> shift2(m,0);
		chrono.start();
		SB.M_Basis(Sigma1, *Serie, d, shift2);
		chrono.stop();
		std::cout << "M-Basis       : " <<chrono.usertime()<<" s"<<std::endl;
	}
#endif


#ifndef  LOW_MEMORY_PMBASIS
	MatrixP Sigma2(F, m, m, d+1);
	std::cout<<"[output sigma    ] : "<<MB(Sigma2.realmeminfo())<<"Mo"<<MEMINFO<<std::endl;	
	chrono.clear();		
	chrono.start();
	SB.PM_Basis(Sigma2, *Serie, d, shift);
	chrono.stop();
	std::cout << "PM-Basis      : " <<chrono.usertime()<<" s"<<std::endl;
	chrono.clear();
	delete Serie;
#else
	MatrixP* sigma_ptr;
	chrono.clear();		
	chrono.start();
	SB.PM_Basis_low(sigma_ptr, Serie, d, shift);
	// Serie is deleted within PM_Basis_low
	chrono.stop();
	std::cout << "PM-Basis      : " <<chrono.usertime()<<" s"<<std::endl;
	chrono.clear();
	delete sigma_ptr;
#endif

	
	// MatrixP Sigma3(F, m, m, d+1);
	//vector<size_t> shift3(m,0);
	// chrono.start();
	// SB.oPM_Basis(Sigma3, Serie, d, shift3);
	// chrono.stop();
	// std::cout << "PM-Basis iter : " <<chrono.usertime()<<" s"<<std::endl;
	std::cout<<endl;
	std::cout<<"[end] :  "<<MEMINFO<<std::endl;
}

int main(int argc, char** argv){
	
	static size_t  m = 64; // matrix dimension
	static size_t  n = 32; // matrix dimension
	static size_t  b = 20; // entries bitsize
	static size_t  d = 32;  // matrix degree
	static long    seed = time(NULL);
	static string target="BEST";

	static Argument args[] = {
		{ 'm', "-m M", "Set row dimension of matrix series to M.", TYPE_INT,     &m },
		{ 'n', "-n N", "Set column dimension of matrix series to N.", TYPE_INT,     &n },
		{ 'd', "-d D", "Set degree of  matrix series to D.", TYPE_INT,     &d },
		{ 'b', "-b B", "Set bitsize of the matrix entries", TYPE_INT, &b },
		{ 's', "-s s", "Set the random seed to a specific value", TYPE_INT, &seed},
		{ 't', "-t T", "Set the targeted benchmark {ALL, BEST}.",            TYPE_STR , &target },
		END_OF_ARGUMENTS
	};

	parseArguments (argc, argv, args);
	
	typedef Givaro::Modular<double>              SmallField;	
	//typedef Givaro::Modular<Givaro::Integer>      LargeField;
	typedef Givaro::Modular<RecInt::ruint128,RecInt::ruint256>  LargeField;

	size_t logd=integer((uint64_t)d).bitsize();

	
	std::cout<<"###  matrix series is of size "<<m<<" x "<<n<<" of degree "<<d<<std::endl;
	if (b < 26){
#ifdef FFT_PROFILER		
		FFT_PROF_LEVEL=1;
#endif

		if (logd>b-4){
			std::cout<<"degree is to large for field bitsize: "<<b<<std::endl;
			exit(0);
		}
        integer p;
		RandomFFTPrime::seeding (seed);
		if (!RandomFFTPrime::randomPrime (p, 1<<b, logd+1))
			throw LinboxError ("RandomFFTPrime::randomPrime failed");
		std::cout<<"# starting sigma basis computation over Fp[x] with p="<<p<<endl;;		
		SmallField F(p);
		typename SmallField::RandIter G(F,0,seed);
		bench_sigma(F,G,m,n,d,target);
	}
	else {
#ifdef FFT_PROFILER		
			FFT_PROF_LEVEL=2;
#endif

		PrimeIterator<IteratorCategories::HeuristicTag> Rd(b,seed);
		integer p = *Rd;
		std::cout<<"# starting sigma basis computation over Fp[x] with p="<<p<<endl;;		
		LargeField F(p);		
		typename LargeField::RandIter G(F,b,seed);

		
		bench_sigma(F,G,m,n,d,target);
	}
	
	
	return 0;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
