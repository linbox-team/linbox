/*
(Message inbox:351)
Received: from chaplin.csd.uwo.ca by mail.eecis.udel.edu id aa08400;
          24 Jan 1999 21:39 EST
Message-Id: <199901250239.CAA04965@heffalump.csd.uwo.ca>
Date: Mon, 25 Jan 1999 02:39:13 GMT
From: Mark Giesbrecht <mwg@csd.uwo.ca>
MIME-Version: 1.0
Content-Type: text/plain; charset=us-ascii
Content-Transfer-Encoding: 7bit
To: Dave Saunders <saunders@mail.eecis.udel.edu>
Subject: back?
In-Reply-To: <199901241115.aa23010@mail.eecis.udel.edu>
References: <199901241115.aa23010@mail.eecis.udel.edu>
X-Mailer: VM 6.22 under 19.15 XEmacs Lucid


Hi Dave,

You back yet?  I got your results on field arithmetic.  I also did some
very naive inner product tests): inner product of two vector<int> and
of two vector<int*>.  The former was about 2.3 times as fast.  I'll
dig up the code and send it you you.

I arrived back and promptly got a rather nasty resperatory infection
(after working myself silly paying for the sin of going to France) so
I've spent the weekend lying around and taking antibiotics.  Have you
though about your travel plans?

Cheers, Mark

Hmm.  Here's the code with timing routines for solaris. Pretty rough,
but you get the idea...

*/
#include <iostream>
#include <vector>
#include <time.h>

int main() {

   long starttime;   
   int size, iters;
   int q = 65521;
   int stride;
   cout << "Vector size, iterations, prime, stride: ";
   cin >> size >> iters >> q >> stride;

   vector<int> x((size_t)size,1), y((size_t)size,1);
   int z;

   starttime=clock();
   for (int j=0; j<iters; j++) {
     z=0;
     for (int i=0; i<size; i++)
       z+= x[i]*y[i];
   }

   cout << z << endl;
   cout << " Array of " << size << " int 1's: Time = " << (clock() -starttime)/1000 << endl;

   for (int i=0; i<size; i++) {
     x[i]=i&q;
     y[i]=(i+1)&q;
   }

   starttime=clock();
   for (int j=0; j<iters; j++) {
     z=0;
     for (int i=0; i<size; i++)
       z+= x[i]+y[i];
   }

   cout << z << endl;
   cout << " Array of " << size << " ints with add not mul: Time = " << (clock() -starttime)/1000 << endl;

   
   starttime=clock();
   for (int j=0; j<iters; j++) {
     z=0;
     for (int i=0; i<size; i++)
       z+= x[i]*y[i];
   }

   cout << z << endl;
   cout << " Array of " << size << " ints: Time = " << (clock() -starttime)/1000 << endl;

   starttime=clock();
   for (int j=0; j<iters; j++) {
     z=0;
     for (int i=0; i<size; i++)
       z = (z + (x[i]*y[i])%q)%q;
   }

   cout << z << endl;
   cout << " Array of " << size << " ints mod "<<q<<": Time = " << (clock() -starttime)/1000 << endl;

// references to (presumably) adjacent values (adjacent in malloc time at least)
   vector<int*> xx(size), yy(size);
   for (int i=0; i<size; i++) {
     xx[i]=new int; *(xx[i])=i&q;
     yy[i]=new int; *(yy[i])=(i+1)&q;
   }

   starttime=clock();
   for (int j=0; j<iters; j++) {
      z=0;
      for (int i=0; i<size; i++)
        z+=  *(xx[i]) * *(yy[i]);
   }

   cout << z << endl;
   cout << " Array of " << size << " refs to ints: Time = " << (clock() -starttime)/1000 << endl;
   
   starttime=clock();
   for (int j=0; j<iters; j++) {
      z=0;
      for (int i=0; i<size; i++)
       z = (z + (*(xx[i])* *(yy[i]))%q)%q;
   }

   cout << z << endl;
   cout << " Array of " << size << " refs to ints mod "<<q<<": Time = " << (clock() -starttime)/1000 << endl;
   
// references to strided locations
   int ii = 0;
   for (int i=0; i<size; i++) {
     xx[i]=&x[ii];
     yy[i]=&y[ii];
     ii = (ii+stride)%size;
   }

   starttime=clock();
   for (int j=0; j<iters; j++) {
      z=0;
      for (int i=0; i<size; i++)
        z+=  *(xx[i]) * *(yy[i]);
   }

   cout << z << endl;
   cout << " Array of " << size << " refs to spaced ints: Time = " << (clock() -starttime)/1000 << endl;
   
   starttime=clock();
   for (int j=0; j<iters; j++) {
      z=0;
      for (int i=0; i<size; i++)
       z = (z + (*(xx[i])* *(yy[i]))%q)%q;
   }

   cout << z << endl;
   cout << " Array of " << size << " refs to spaced ints mod "<<q<<": Time = " << (clock() -starttime)/1000 << endl;
   
   starttime=clock();
   for (int j=0; j<iters; j++) {
      z=0;
      ii = 0;
      for (int i=0; i<size; i++){
       z = (z + (*(xx[ii])* *(yy[ii]))%q)%q;
	ii = (ii+stride)%size;
      }
   }

   cout << z << endl;
   cout << " Array of " << size << " spaced refs to spaced ints mod "<<q<<": Time = " << (clock() -starttime)/1000 << endl;

// references to spread locations, broken stride   
   int st = stride;
   ii = 0;
   for (int i=0; i<size; i++) {
     xx[i]=&x[ii];
     yy[i]=&y[ii];
     ii = (ii+st)%size;
     st = (st+1)%size;
   }

   starttime=clock();
   for (int j=0; j<iters; j++) {
      z=0;
      for (int i=0; i<size; i++)
        z+=  *(xx[i]) * *(yy[i]);
   }

   cout << z << endl;
   cout << " Array of " << size << " refs to variably spaced ints: Time = " << (clock() -starttime)/1000 << endl;
   
   starttime=clock();
   for (int j=0; j<iters; j++) {
      z=0;
      for (int i=0; i<size; i++)
       z = (z + (*(xx[i])* *(yy[i]))%q)%q;
   }

   cout << z << endl;
   cout << " Array of " << size << " refs to variably spaced ints mod "<<q<<": Time = " << (clock() -starttime)/1000 << endl;
   
     
   st = stride;
   starttime=clock();
   for (int j=0; j<iters; j++) {
      z=0;
      ii = 0;
      for (int i=0; i<size; i++){
       z = (z + (*(xx[ii])* *(yy[ii]))%q)%q;
	ii = (ii+st)%size;
	st = (1+st)%size;
      }
   }

   cout << z << endl;
   cout << " Array of " << size << " variably spaced refs to variably spaced ints mod "<<q<<": Time = " << (clock() -starttime)/1000 << endl;
   
     
   
   
} 

