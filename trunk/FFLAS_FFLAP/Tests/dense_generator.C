#include <stdlib.h>
#include <stdio.h>
#include <iostream.h>



template<class T>
T& myrand (T& r, long size) {
	  if (size < 0)
		     return r = T( (lrand48() % (-size-size)) + size );
	    else
		       return r = T(  lrand48() % size ) ;
};


int main(int argc, char ** argv) {

	 long ni=10,nj=10,max=100;
	 int offset = 0;
		    
	 if (argc > ++offset)
	          ni = atoi( argv[offset] );
	 if (argc > ++offset)
	       nj = atoi( argv[offset] );
	 if (argc > ++offset)
	       max = atoi( argv[offset] );
				        
	 long tmp;
	 printf("%ld %ld M\n", ni, nj);
	 for (long i = 0; i < ni; ++i) 
	   for (long j = 0; j < nj; ++j)
	     printf("%ld %ld %ld\n", i+1, j+1, myrand(tmp, max));

	printf("0 0 0\n");

return 0;
}
