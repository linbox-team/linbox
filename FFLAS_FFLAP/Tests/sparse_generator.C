#include <stdlib.h>
#include <stdio.h>
#include <iostream.h>
#include <vector.h>
#include <set.h>
#include <algo.h>


size_t vssize(const vector< set<long> >& vs) {
    size_t t = 0;
    for (vector< set<long> >::const_iterator i = vs.begin(); i != vs.end(); ++i)
        t += (*i).size();
    return t;
}

template<class T>
T& myrand (T& r, long size) {
  if (size < 0)
   return r = T( (lrand48() % (-size-size)) + size );
  else
   return r = T(  lrand48() % size ) ;
};



int main(int argc, char ** argv) {
    long ni=10,nj=10,nz=100,max=100;
    int offset = 0;
    
    if (argc > ++offset)
        ni = atoi( argv[offset] );
    if (argc > ++offset)
        nj = atoi( argv[offset] );
    if (argc > ++offset)
        nz = atoi( argv[offset] );
    if (nz > (ni*nj))
        nz = ni*nj;
    if (argc > ++offset)
        max = atoi( argv[offset] );
    
    vector< set<long> > Matrix (ni);

    while (vssize(Matrix) < nz) {
        long placei = lrand48() % ni;
        long placej = lrand48() % nj;
        Matrix[placei].insert( placej );
    }
    
    printf("%ld %ld M\n", ni, nj);
    long tmp;

    for (long i = 0; i < ni; ++i) {
        vector<long> vv(Matrix[i].size());
        set<long>::const_iterator si = Matrix[i].begin();
        for(vector<long>::iterator vi = vv.begin(); vi != vv.end();)
            *vi++ = *si++;
        
        sort( vv.begin(), vv.end() );
        for(vector<long>::const_iterator j = vv.begin(); j != vv.end() ; ++j)
            printf("%ld %ld %ld\n", i+1, (*j)+1, myrand(tmp, max) );
    }
    
     printf("0 0 0\n");
   
    
}
