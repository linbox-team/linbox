/* mapleLB.cc - July 11, 2002
 *  Maple - LinBox Interface
 *  Copyright(C) 2002 David Saunders
 *               2002 Rich Seagraves
 *
 */

/* Now then.  This file contains the various functions used to interface
 * between Maple and the current solutions availble in the linbox package.
 * At core there are 3 LinBox methods currently available:  rank, to cacluate
 * the rank of a Matrix, det, to calculate the determinant of a matrix, and
 * minpoly, to calculate the Minimal Polynomial.  The functions provided below
 * translate input from a number of Maple procedures and call the linbox
 * procedures on that input.
 *
 * Note: All functions that are wrapped in the "extern "C" { ... } "
 * declarations are done so that they recieve C name mangling as opposed to
 * C++ name mangling.  This is done for all functions called by Maple, as
 * Maple can perform C name mangling but not C++ name mangling
 */

#include "maplec.h"
#include "MapleBB.h"
#include "linbox/field/modular.h"
#include "linbox/integer.h"
#include "linbox/solutions/rank.h"
#include "linbox/solutions/det.h"
#include "linbox/solutions/minpoly.h"
#include <vector>

using LinBox::Modular;
using LinBox::MapleBB;
using LinBox::integer;
using std::vector;


extern "C" ALGEB rank1(MKernelVector,ALGEB*);
extern "C" ALGEB rank2(MKernelVector, ALGEB*);
extern "C" ALGEB rank3(MKernelVector,ALGEB*);

extern "C" ALGEB det1(MKernelVector,ALGEB*);
extern "C" ALGEB det2(MKernelVector,ALGEB*);
extern "C" ALGEB det3(MKernelVector,ALGEB*);

extern "C" ALGEB minpoly1(MKernelVector,ALGEB*);
extern "C" ALGEB minpoly2(MKernelVector,ALGEB*);
extern "C" ALGEB minpoly3(MKernelVector,ALGEB*);

extern "C" ALGEB apply(MKernelVector, ALGEB*);
extern "C" ALGEB appplyt(MKernelVector, ALGEB*);

void LItoM(integer &, int*, int length);
void MtoLI(integer &, int*, int length);

extern "C" { 
  ALGEB rank1(MKernelVector kv, ALGEB* args) 
  {
    unsigned long result;
    size_t* rowP, *colP, m, n, nonzeros;
    long* data, prime;

    prime = MapleToInteger32(kv,args[1]);
    nonzeros = RTableNumElements(kv,args[2]);
    m = RTableUpperBound(kv,args[2],1);
    n = RTableUpperBound(kv,args[2],2);
    rowP = (size_t*) RTableSparseIndexRow(kv,args[2],1);
    colP = (size_t*) RTableSparseIndexRow(kv,args[2],2);
    data = (long*) RTableDataBlock(kv,args[2]);

    Modular<long> field(prime);
    MapleBB< Modular<long>, vector<long> > BB( field, data, rowP, colP, m, n, nonzeros);
    LinBox::rank(result, BB, field);
    return ToMapleInteger(kv, result);
  }
} 

/* rank2 - called when the user calls the LBrank command with a sparse
 * matrix, specifically when that Matrix contains entries of a wordsize or
 * smaller.  Calls the same function as rankLB above, and returns the integer
 * result of the rank calculation
 */

extern "C" {
  ALGEB rank2(MKernelVector kv, ALGEB* args) 
  {
    long prime, *data;
    size_t nonzeros, m, n, *rowP, *colP;
    unsigned long result;

    prime = MapleToInteger32(kv,args[1]);
    rowP = (size_t*) RTableDataBlock(kv,args[2]);
    nonzeros  = RTableUpperBound(kv, args[2],1) - RTableLowerBound(kv,args[2],1) + 1;
    colP = (size_t*) RTableDataBlock(kv,args[3]);
    m = MapleToInteger32(kv,args[5]);
    n = MapleToInteger32(kv,args[6]);
    data = (long*) RTableDataBlock(kv,args[4]);

    Modular<long> field(prime);
    MapleBB< Modular<long>, vector<long> > BB( field, data, rowP, colP, m, n, nonzeros);

    LinBox::rank(result, BB, field);
    return ToMapleInteger(kv, result);
  }
}  

/* rank3 - Function called by the LBrank procedure in maple when Matrix
 * entries are greater than 1 wordsize long.  Translates between Maple large
 * integers and GMP integers.  Returns the rank of the Matrix
 */

extern "C" {
  ALGEB rank3(MKernelVector kv, ALGEB* args) 
  {
    unsigned long result;
    integer prime, *iArray;
    size_t nonzeros, m, n, *rowP, *colP;
    ALGEB chunks;

    MtoLI(prime, (int*) RTableDataBlock(kv, args[1]), RTableNumElements(kv,args[1]) );
    rowP = (size_t*) RTableDataBlock(kv,args[2]);
    nonzeros  = RTableUpperBound(kv, args[2],1) - RTableLowerBound(kv,args[2],1) + 1;
    colP = (size_t*) RTableDataBlock(kv,args[3]);
    m = (size_t) MapleToInteger32(kv,args[5]);
    n = (size_t) MapleToInteger32(kv,args[6]);

    iArray = new integer[nonzeros];
    for(int i = 0; i < nonzeros; i++) {
      chunks = MapleListSelect(kv,args[4],(M_INT) i + 1 );
      if( IsMapleInteger(kv,chunks)) 
	iArray[i] = integer(MapleToInteger32(kv,chunks));
      else 
	MtoLI(iArray[i], (int*) RTableDataBlock(kv,chunks), RTableNumElements(kv,chunks));
    }

    Modular<integer> field(prime);
    MapleBB< Modular<integer>, vector<integer> > BB( field, iArray, rowP, colP, m, n, nonzeros);

    LinBox::rank(result, BB, field);
    delete [] iArray;
    return ToMapleInteger(kv, result);
  }
}  

/* det1 - Called by the LBdet maple procedure when LBdet receives a NAG
 * sparse Matrix as input.  Calls the linboxdet function and returns the 
 * determinant.
 */

extern "C" {
  ALGEB det1(MKernelVector kv, ALGEB* args) 
  {

    long prime, *data, result;
    size_t* rowP, *colP, m, n, nonzeros;

    prime = MapleToInteger32(kv,args[1]);
    nonzeros = RTableNumElements(kv,args[2]);
    m = RTableUpperBound(kv,args[2],1);
    n = RTableUpperBound(kv,args[2],2);
    rowP = (size_t*) RTableSparseIndexRow(kv,args[2],1);
    colP = (size_t*) RTableSparseIndexRow(kv,args[2],2);
    data = (long*) RTableDataBlock(kv,args[2]);

    Modular<long> field(prime);
    MapleBB< Modular<long>, vector<long> > BB( field, data, rowP, colP, m, n, nonzeros);

    LinBox::det( result, BB, field);
    return ToMapleInteger(kv, result);
  }
}

/* det2 - Called by the LBdet maple procedure when LBdet recieves a sparse
 * Matrix with entries that are a wordsize long.  calls teh linboxdet
 * function and returns the determinant
 */

extern "C" {
  ALGEB det2(MKernelVector kv, ALGEB* args) 
  {

    long result, prime, *data;
    size_t nonzeros, m, n, *rowP, *colP;

    prime = MapleToInteger32(kv,args[1]);
    rowP = (size_t*) RTableDataBlock(kv,args[2]);
    nonzeros  = RTableUpperBound(kv, args[2],1) - RTableLowerBound(kv,args[2],1) + 1;
    colP = (size_t*) RTableDataBlock(kv,args[3]);
    m = (size_t) MapleToInteger32(kv,args[5]);
    n = (size_t) MapleToInteger32(kv,args[6]);
    data = (long*) RTableDataBlock(kv,args[4]);

    Modular<long> field(prime);
    MapleBB< Modular<long>, vector<long> > BB(field, data, rowP, colP, m, n, nonzeros);

    LinBox::det( result, BB, field);
    return ToMapleInteger(kv, result);
  }
}  

/* det3 - Called by the LBdet maple procedure when LBdet recieves as input
 * a sparse Matrix with entries larger than a wordsize integer.  Converts
 * all integers to GMP integers and calls the linbox command.  The answer
 * is converted into an array of wordsize bytes which are returned to Maple
 */

extern "C" {
  ALGEB det3(MKernelVector kv, ALGEB* args) 
  {
    // Numerous variable declarations
    integer prime, *iArray, detres;
    INTEGER32 *data;
    size_t nonzeros, m, n, *rowP, *colP;
    ALGEB chunks, rtable;
    RTableSettings s;
    M_INT bound[2];

    // Extract data from Maple call to external code
    MtoLI(prime, (int*) RTableDataBlock(kv, args[1]), RTableNumElements(kv,args[1]) );
    rowP = (size_t*) RTableDataBlock(kv,args[2]);
    nonzeros  = RTableUpperBound(kv, args[2],1) - RTableLowerBound(kv,args[2],1) + 1;
    colP = (size_t*) RTableDataBlock(kv,args[3]);
    m = (size_t) MapleToInteger32(kv,args[5]);
    n = (size_t) MapleToInteger32(kv,args[6]);

    // Create array of GMP integers
    iArray = new integer[nonzeros];
    for(int i = 0; i < nonzeros; i++) {
      chunks = MapleListSelect(kv,args[4],(M_INT) i + 1 );
      if( IsMapleInteger(kv,chunks)) 
	iArray[i] = integer(MapleToInteger32(kv,chunks));
      else 
        MtoLI(iArray[i], (int*) RTableDataBlock(kv,chunks), RTableNumElements(kv,chunks));
    }
    
    Modular< integer > field(prime);
    MapleBB< Modular<integer>, vector<integer> > BB(field, iArray, rowP, colP, m, n, nonzeros);
    LinBox::det(detres, BB, field);
    delete [] iArray;

    // Create Array  for result
    bound[0] = 1; bound[1] = detres.size();
    kv->rtableGetDefaults(&s);
    s.data_type = RTABLE_INTEGER32;
    s.subtype = RTABLE_ARRAY;
    s.storage = RTABLE_RECT;
    s.num_dimensions = 1;
    rtable = kv->rtableCreate(&s,NULL,bound);

    // Initialize the array
    data = (INTEGER32*) RTableDataBlock(kv, rtable);
    LItoM(detres, data, bound[1]);
    
    // Return the result
    return rtable;
  }
}  

/* minpoly1 - Called by the LBminpoly maple procedure when it recieves a
 * NAG sparse format matrix as input.  Calls the linbox minpoly function,
 * recieves a vector of coefficients (note this is a dense vector of
 * coefficients), and turns this vector into a Maple list that is returned.
 * This list is then turned into a polynomial by the maple procedure that
 * called minpolyLB.
 */

extern "C" { 
  ALGEB minpoly1(MKernelVector kv, ALGEB* args) 
  {

    long prime, *data;
    size_t *rowP, *colP, m, n, nonzeros;
    vector<long> result;
    vector<long>::iterator r_i;
    ALGEB retList;
    int i;

    prime = MapleToInteger32(kv,args[1]);
    nonzeros = RTableNumElements(kv,args[2]);
    m = RTableUpperBound(kv,args[2],1);
    n = RTableUpperBound(kv,args[2],2);
    rowP = (size_t*) RTableSparseIndexRow(kv,args[2],1);
    colP = (size_t*) RTableSparseIndexRow(kv,args[2],2);
    data = (long*) RTableDataBlock(kv,args[2]);
    
    Modular< long> field(prime);
    MapleBB< Modular<long>, vector<long> > BB( field, data, rowP, colP, m, n, nonzeros);

    LinBox::minpoly(result, BB, field);

    retList = MapleListAlloc(kv, (M_INT) result.size() );
    for(i = 1, r_i = result.begin(); r_i != result.end(); ++i, ++r_i) {
      MapleListAssign(kv,retList, i, ToMapleInteger(kv, *r_i));
    }
    return retList;
  }
} 

/* minpoly2 - Called by the maple procedure LBminpoly when LBminpoly recieves
 * a Matrix in sparse format, whose entries are all wordsize int's.  Operates
 * in the samne manner as minpolyLB above.
 */

extern "C" {
  ALGEB minpoly2(MKernelVector kv, ALGEB* args) 
  {

    long prime, *data;
    size_t nonzeros, m, n, *rowP, *colP;
    vector<long> result;
    vector<long>::iterator r_i;
    ALGEB retList;
    int i;
    
    prime = MapleToInteger32(kv,args[1]);
    rowP = (size_t*) RTableDataBlock(kv,args[2]);
    nonzeros  = RTableUpperBound(kv, args[2],1) - RTableLowerBound(kv,args[2],1) + 1;
    colP = (size_t*) RTableDataBlock(kv,args[3]);
    m = (size_t) MapleToInteger32(kv,args[5]);
    n = (size_t) MapleToInteger32(kv,args[6]);
    data = (long*) RTableDataBlock(kv,args[4]);

    Modular<long> field(prime);
    MapleBB< Modular<long>, vector<long> > BB( field, data, rowP, colP, m, n, nonzeros);

    LinBox::minpoly(result, BB, field);
    
    retList = MapleListAlloc(kv, result.size() );
    for( i = 1, r_i = result.begin(); r_i != result.end(); ++i, ++r_i) {
      MapleListAssign(kv, retList, i, ToMapleInteger(kv, *r_i));
    }
    
    return retList;
  }
}  

/* minpoly3 - Called by the LBminpoly maple procedure when LBminpoly rceives
 * a Matrix in sparse format with entries greater than 1 wordsize integer.
 * Performs the most data conversions of all, as it converts all maple
 * integers into GMP integers, then converts all coefficients back.
 */


extern "C" {
  ALGEB minpoly3(MKernelVector kv, ALGEB* args) 
  {
    // Numerous variable declarations
    integer prime, *iArray;
    size_t nonzeros, m, n, *rowP, *colP;
    INTEGER32 *data;
    ALGEB chunks, retList;
    RTableSettings s;
    M_INT i, bound[2];
    vector<integer> result;
    vector<integer>::iterator r_i;

    // Extract data from Maple call to external code
    MtoLI(prime, (int*) RTableDataBlock(kv, args[1]), RTableNumElements(kv,args[1]) );
    rowP = (size_t*) RTableDataBlock(kv,args[2]);
    nonzeros  = RTableUpperBound(kv, args[2],1) - RTableLowerBound(kv,args[2],1) + 1;
    colP = (size_t*) RTableDataBlock(kv,args[3]);
    m = (size_t) MapleToInteger32(kv,args[5]);
    n = (size_t) MapleToInteger32(kv,args[6]);

    // Create array of GMP integers
    iArray = new integer[nonzeros];
    for(i = 0; i < nonzeros; i++) {
      chunks = MapleListSelect(kv,args[4],(M_INT) i + 1 );
      if( IsMapleInteger(kv,chunks)) 
	iArray[i] = integer(MapleToInteger32(kv,chunks));
      else 
        MtoLI(iArray[i], (int*) RTableDataBlock(kv,chunks), RTableNumElements(kv,chunks));
    }
    
    
    Modular< integer> field(prime);
    MapleBB< Modular<integer>, vector<integer> > BB( field, iArray, rowP, colP, m, n, nonzeros);

    LinBox::minpoly(result, BB, field);
    delete [] iArray;

    // Create settings for GMP->Maple Arrays
    bound[0] = 1;
    kv->rtableGetDefaults(&s);
    s.data_type = RTABLE_INTEGER32;
    s.subtype = RTABLE_ARRAY;
    s.num_dimensions = 1;

    // Allocate list of returned values
    retList = MapleListAlloc(kv, result.size());
    // Loop through each value, create the array then add it to the list
    for(i = 1, r_i = result.begin(); r_i != result.end(); ++i, ++r_i) {
      bound[1] = r_i->size();
      chunks = kv->rtableCreate(&s, NULL, bound);
      data = (INTEGER32*) RTableDataBlock(kv, chunks);
      LItoM( *r_i, data, bound[1]);
      MapleListAssign(kv, retList, i, chunks);
    }
    // Return the result
    return retList;
  }
}  

extern "C" {
  ALGEB apply(MKernelVector kv, ALGEB* args) 
  {
    ALGEB rvec, indList, datList;
    int bounds[2], i, j, index, nnzV;
    long prime, *data;
    size_t *rowP, *colP, m, n, nonzeros;
    RTableSettings s;
    RTableData d;
    vector<long> x, y;

    prime = MapleToInteger32(kv,args[1]);
    data = (long*) RTableDataBlock(kv, args[2]);
    rowP = (size_t*) RTableDataBlock(kv, args[3]);
    colP = (size_t*) RTableDataBlock(kv, args[4]);
    nonzeros = RTableNumElements(kv, args[2]);

    m = MapleToInteger32(kv, args[5]);
    n = MapleToInteger32(kv, args[6]);

    Modular<long> field(prime);
    MapleBB< Modular<long>, vector<long> > BB(field, data, rowP, colP, m, n, nonzeros);

    x.reserve(m);
    indList = args[7];
    datList = args[8];
    nnzV = MapleToInteger32(kv, args[9]);

    for(i = 1, j = 1 ;i <= nnzV; ++i, ++j) {

      index = MapleToInteger32(kv, MapleListSelect(kv, indList, i) );
      while(j < index) {
	x.push_back(0L);
	++j;
      }

      x.push_back( (long) MapleToInteger32(kv, MapleListSelect(kv, datList, i) ) );
      
    }
    while(j <= m) {
      x.push_back(0L);
      ++j;
    }


    y.reserve(n);
    for(i = 0; i < n; i++)
      y.push_back(0L);

    BB.apply(y, x);
    
    bounds[0] = 1; bounds[1] = n;
    kv->rtableGetDefaults(&s);
    s.data_type = RTABLE_INTEGER32;
    s.subtype = RTABLE_COLUMN;
    s.storage = RTABLE_SPARSE;
    s.num_dimensions = 1;
    rvec = kv->rtableCreate(&s, NULL, bounds);

    vector<long>::iterator v_i;
    for(index = 1, v_i = y.begin(); v_i != y.end(); ++index, ++v_i) {
      if(*v_i != 0) {
	d.int32 = *v_i;
	RTableAssign(kv, rvec, &index, d);
      }
    }

    return rvec;
  }
}

extern "C" {
  ALGEB applyt(MKernelVector kv, ALGEB* args) 
  {
    ALGEB rvec, indList, datList;
    int bounds[2], i, j, index, nnzV;
    long prime, *data;
    size_t m, n, nonzeros, *rowP, *colP;
    RTableSettings s;
    RTableData d;
    vector<long> x, y;

    prime = MapleToInteger32(kv,args[1]);
    data = (long*) RTableDataBlock(kv, args[2]);
    rowP = (size_t*) RTableDataBlock(kv, args[3]);
    colP = (size_t*) RTableDataBlock(kv, args[4]);
    nonzeros = RTableNumElements(kv, args[2]);

    m = MapleToInteger32(kv, args[5]);
    n = MapleToInteger32(kv, args[6]);

    Modular<long> field(prime);
    MapleBB< Modular<long>, vector<long> > BB(field, data, rowP, colP, m, n, nonzeros);

    x.reserve(n);
    indList = args[7];
    datList = args[8];
    nnzV = MapleToInteger32(kv, args[9]);

    for(i = 1, j = 1; i <= nnzV; ++i, ++j) {

      index = MapleToInteger32(kv, MapleListSelect(kv, indList, i));
      while(j < index) {
	x.push_back(0L);
	++j;
      }

      x.push_back( (long) MapleToInteger32(kv, MapleListSelect(kv, datList, i) ) );
      
    }
    while(j <= n) {
      x.push_back(0L);
      ++j;
    }


    y.reserve(m);
    for(i = 0; i < m; i++)
      y.push_back(0L);

    BB.applyTranspose(y, x);
    
    bounds[0] = 1; bounds[1] = m;
    kv->rtableGetDefaults(&s);
    s.data_type = RTABLE_INTEGER32;
    s.subtype = RTABLE_ROW;
    s.storage = RTABLE_SPARSE;
    s.num_dimensions = 1;
    rvec = kv->rtableCreate(&s, NULL, bounds);

    vector<long>::iterator v_i;
    for(index = 1, v_i = y.begin(); v_i != y.end(); ++index, ++v_i) {
      if(*v_i != 0) {
	d.int32 = *v_i;
	RTableAssign(kv, rvec, &index, d);
      }
    }

    return rvec;
  }
}




/* MtoLI - Conversion from Maple Array to GMP integer.  This function takes
 * an array of mod 10000 integers (the basic Maple storage scheme for large
 * integers) and returns a GMP integer.  This function will probably change
 * in future versions in favor of a faster method.
 */

void MtoLI(integer &result, int* digits, int length) {

  int i;
  static integer base(10000);
  result = 0;
  for(i = 0; i < length; i++) {
    result += ( pow(base, i) * digits[i] );
  }
}

/* LItoM - Conversion from a GMP integer to a Maple Array.  Creates a Maple
 * array of wordsize integers, which are returned to Maple and reconstructed
 * Note in future versions this code will most likely change in favor of
 * something more efficient
 */

void LItoM(integer & In, int* data, int length) 
{
  size_t i;
  for( i = 0; i < length; i++) {
    data[i] = In[i];
  }

  return;
}

