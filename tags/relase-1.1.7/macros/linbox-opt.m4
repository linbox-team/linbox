AC_DEFUN([LB_OPT],
[
AC_MSG_CHECKING([whether to use run time optimization])

AC_ARG_ENABLE(optimization, 
[  --enable-optimization  Enable run time optimization in LinBox code],
[
AC_MSG_RESULT(yes)

	
BACKUP_CXXFLAGS=${CXXFLAGS}
BACKUP_LIBS=${LIBS}
	
if test "x$HAVE_BLAS" = "xyes" ;then
AC_MSG_CHECKING([best threshold for Strassen-Winograd matrix multiplication])


CXXFLAGS="${BACKUP_CXXFLAGS} -I`pwd` -I`pwd`/linbox ${BLAS_CFLAGS} ${GMP_CFLAGS}  ${GIVARO_CFLAGS} ${CBLAS_FLAG}" 
LIBS="${BACKUP_LIBS} ${BLAS_LIBS} ${GIVARO_LIBS} ${GMP_LIBS} " 


echo   " #define __LINBOX_INT8  $LINBOX_INT8  	 
	 #define __LINBOX_INT16 $LINBOX_INT16 	 
	 #define __LINBOX_INT32 $LINBOX_INT32 	 
	 #define __LINBOX_INT64 $LINBOX_INT64 	 
" > linbox/linbox-config.h 


AC_TRY_RUN([	#define LinBoxSrcOnly
		#include <iostream>
		#include <fstream>
		#define _LINBOX_LINBOX_CONFIG_H
		#define __LINBOX_CONFIGURATION
		#include <linbox/config-blas.h>
		#include <linbox/linbox-config.h>
		#include <linbox/field/modular-double.h>
		#include <linbox/fflas/fflas.h>
		#include <linbox/util/timer.h>

		using namespace LinBox;
		int main () {  
  
		  Modular<double> F(17);
		  size_t n=300, nmax=5000, prec=512, nbest=0, count=0;
		  LinBox::Timer chrono;
		  double basetime, time;
		  bool bound=false;
  
		  double *A, *C;	 
		  A = new double[nmax*nmax];
		  C = new double[nmax*nmax];
		  for (size_t i=0; i<nmax*nmax;++i){
		    A[i]=2.;		  
	  	  }
    
		  std::ofstream outlog;
		  outlog.open("config.log", std::ofstream::out | std::ofstream::app);
		  outlog << std::endl 
			 << "Threshold for finite field Strassen-Winograd matrix multiplication"
			 << std::endl;
		  do {    
		    chrono.start();	
		    FFLAS::fgemm(F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
				 n, n, n, 1., A, n, A, n, 0., C, n, 0);
		    chrono.stop();
		    std::cout << std::endl 
			      << "fgemm " << n << "x" << n << ": " 
			      << chrono.usertime() << " s, " 
                              << (2.0/chrono.usertime()*n/100.0*n/100.0*n/100.0) << " Mffops" 
			      << std::endl;
		    outlog << std::endl 
			      << "fgemm " << n << "x" << n << ": " 
			      << chrono.usertime() << " s, " 
                              << (2.0/chrono.usertime()*n/100.0*n/100.0*n/100.0) << " Mffops" 
			      << std::endl;
		    basetime= chrono.usertime();
		    chrono.clear();
		    chrono.start();	
		    FFLAS::fgemm(F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
				 n, n, n, 1., A, n, A, n, 0., C, n, 1);
		    chrono.stop();
		    std::cout << "1Wino " << n << "x" << n << ": " 
			      << chrono.usertime() << " s, " 
			      << (2.0/chrono.usertime()*n/100.0*n/100.0*n/100.0) << " Mffops" 
			      << std::endl;
		    outlog << "1Wino " << n << "x" << n << ": " 
			      << chrono.usertime() << " s, " 
			      << (2.0/chrono.usertime()*n/100.0*n/100.0*n/100.0) << " Mffops" 
			      << std::endl;
		    time= chrono.usertime();
        
		    if (basetime > time ){ 
		      count++;
		      if (count > 1){	
	    		 nbest=n;
		         bound=true;
		         prec=prec>>1;
		         n-=prec;     
		      }
		    }
		    else{
		      count=0;	    
		      if (bound)
			prec=prec>>1;
		      n+=prec;
		    }
		  } while ((prec > 64 ) && (n < nmax));

		  std::ofstream out("WinoThreshold");
		  out<<nbest;
		  out.close();

		  outlog << "defined __LINBOX_STRASSEN_OPTIMIZATION" << std::endl
			 << "defined __LINBOX_WINOTHRESHOLD to " << nbest << "" << std::endl;
	          outlog.close();

		  delete[] A;		 
		  delete[] C;  
  
		  return 0; 
		}
	],[	
	AC_MSG_RESULT(done)
	WT="`cat WinoThreshold`"
	if test "$WT" != "0"; then 
	 AC_DEFINE(STRASSEN_OPTIMIZATION,,[Define if optimized  threshold for Strassen-Winograd matrix multiplication is available])
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
	\rm -f linbox/linbox-config,h 2>&1 > /dev/null
	])

fi;
],[
AC_MSG_RESULT(no)
\rm -f linbox/linbox-config.h  2>&1 > /dev/null
])

])
