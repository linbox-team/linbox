AC_DEFUN([LB_OPT],
[
AC_MSG_CHECKING([whether to use run time optimization])

AC_ARG_ENABLE(optimization, [--enable-optimization  Enable run time optimization in LinBox code],
[
AC_MSG_RESULT(yes)

	
BACKUP_CXXFLAGS=${CXXFLAGS}
BACKUP_LIBS=${LIBS}
	
if test "x$HAVE_ATLAS" = "xyes" ;then
AC_MSG_CHECKING([best threshold for Strassen matrix multiplication])


CXXFLAGS="${BACKUP_CXXFLAGS} -I`pwd` ${ATLAS_CFLAGS} ${GMP_CFLAGS}  ${GIVARO_CFLAGS} " 
LIBS="${BACKUP_LIBS} ${ATLAS_LIBS} ${GMP_LIBS}" 

AC_TRY_RUN([	#define LinBoxSrcOnly
		#include <iostream>
		#include <linbox/fflas/fflas.h>
		#include <linbox/field/modular-double.h>
		#include <linbox/util/timer.h>

		using namespace LinBox;
		int main () {  
  
		  Modular<double> F(17);
		  size_t n=300, nmax=5000, prec=512, nbest=0, count=0;
		  LinBox::Timer chrono;
		  double basetime, time;
		  bool bound=false;
  
		  double *A, *B, *C;	 
		  A = new double[nmax*nmax];
		  B = new double[nmax*nmax];
		  C = new double[nmax*nmax];
		  for (size_t i=0; i<nmax*nmax;++i){
		    A[i]=2.;
		    B[i]=3.;
	  	  }
    
		  do {    
		    chrono.start();	
		    FFLAS::fgemm(F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
				 n, n, n, 1., A, n, B, n, 0., C, n, 0);
		    chrono.stop();
		    basetime= chrono.usertime();
		    chrono.clear();
		    chrono.start();	
		    FFLAS::fgemm(F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
				 n, n, n, 1., A, n, B, n, 0., C, n, 1);
		    chrono.stop();
		    time= chrono.usertime();
        
		    if (basetime < time ){ 
		      count=0;	    
		      if (bound)
			prec=prec>>1;
		      n+=prec;
		    }
		    else{
		      count++;
		      if (count > 2){	
	    		 nbest=n;
		         bound=true;
		         prec=prec>>1;
		         n-=prec;     
		      }
		    }
		  } while ((prec > 32 ) && (n < nmax));

		  ofstream out("WinoThreshold");
		  out<<nbest;
		  out.close();

		  delete[] A;
		  delete[] B;
		  delete[] C;  
  
		  return 0; 
		}
	],[	
	AC_MSG_RESULT(done)
	WT="`cat WinoThreshold`"
	if test "$WT" != "0"; then 
	 AC_DEFINE(STRASSEN_OPTIMIZATION,,[Define if optimized  threshold for Strassen matrix multiplication is available])
	 AC_DEFINE_UNQUOTED(WINOTHRESHOLD, $WT, [optimized threshold for switching to strassen matrix multiplication])
	fi
	],[
	AC_MSG_RESULT(problem)
	strassen_opti="no"
	break
	],[
	AC_MSG_RESULT(cross compilation)
	strassen_opti="no"
	break
	])

fi;
],[
AC_MSG_RESULT(problem)
])

])