/*
(Message inbox:353)
Received: from chaplin.csd.uwo.ca by mail.eecis.udel.edu id aa21643;
          25 Jan 1999 14:15 EST
Message-Id: <199901251915.TAA28926@heffalump.csd.uwo.ca>
Date: Mon, 25 Jan 1999 19:15:08 GMT
From: Mark Giesbrecht <mwg@csd.uwo.ca>
MIME-Version: 1.0
Content-Type: text/plain; charset=us-ascii
Content-Transfer-Encoding: 7bit
To: Dave Saunders <saunders@mail.eecis.udel.edu>
Subject: Re: code 
In-Reply-To: <199901251346.aa19973@mail.eecis.udel.edu>
References: <199901251631.QAA28032@heffalump.csd.uwo.ca>
	<199901251346.aa19973@mail.eecis.udel.edu>
X-Mailer: VM 6.22 under 19.15 XEmacs Lucid


Hi Dave,

So I reversed the order that refs in the yy array were allocated and
the factor when up to 3x slower for the references.  Suggests caching
to me.  I'm sure more pathelogical screwing around (i.e., arranging it
so that every reference hits the same cache line) would probably bump
this further but is this "interesting".  Also, it is hard to imagine
seeing 20x slowdown without really trying for it (if then).  Something
else is going on.  Did Erich show you the code?
 
Mark

*/

#include <iostream>
#include <vector>
#include <time.h>

void nada(int&); // an external routine which does nothing
                 // i.e., in a separate file: void nada(int& z) {}
                 // You may or may not need this...

int main() {

   long starttime;   
   int size, iters;
   cin >> size >> iters;
   // int size=100000, iters=1000;

   vector<int> x((vector<int>::size_type)size,1), 
               y((vector<int>::size_type)size,1);
   int z;
   
   starttime=clock();
   for (int j=0; j<iters; j++) {
     z=0;
     for (int i=0; i<size; i++)
       z+= x[i]*y[i];
     nada(z);  // just so the optimizer doesn't eat my loop!
   }
   int inttime=clock()-starttime;



   cout << " Array of " << size << " ints ("
        << iters << " iterations) : Time = " << inttime/1000 << endl;


   vector<int*> xx(size), yy(size);
   for (int i=0; i<size; i++) {
     xx[i]=new int; *(xx[i])=1;
     yy[size-1-i]=new int; *(yy[size-1-i])=1;
   }

   starttime=clock();
   for (int j=0; j<iters; j++) {
      z=0;
      for (int i=0; i<size; i++)
        z+=  *(xx[i]) * *(yy[i]);
      nada(z);  // just so the optimizer doesn't eat my loop!
   }
   int reftime=clock()-starttime;

   cout << " Array of " << size << " refs to spread out ints ("
        << iters << " iterations) : Time = " << reftime/1000 << endl;

   cout << endl;
   cout << "Ratio (reftime/inttime): " << ((double)reftime)/(double)inttime 
        << endl;
   
} 
