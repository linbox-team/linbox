/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* Copyright (C) LinBox
 *
 *
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin Street - Fifth Floor,
 * Boston, MA 02110-1301, USA.
 */


#include "linbox/field/modular.h"
#include "linbox/integer.h" // <- Wrapper for gmp BIG int support
// #include "linbox/field/integer.h" <- When linbox supports computations
//                                         over the whole integers
#include "linbox/blackbox/triplesbb.h" // Special blackbox for this interface
#include "linbox/solutions/rank.h"
#include "linbox/solutions/det.h"
#include "linbox/solutions/minpoly.h"
// #include "linbox/solutions/ssolve.h" <- When linbox supports System solver

#include <vector>
#include <map>
#include <utility>
#include <cstring>
#include <cstdio>

#include "maplec.h"

using LinBox::integer;
typedef std::vector<long> Vectorl;
typedef std::vector<integer> VectorI;
typedef LinBox::TriplesBB<LinBox::Modular<long> > TriplesBBi;
typedef LinBox::TriplesBB<LinBox::Modular<integer> > TriplesBBI;

/* Type references: Used to identify the type of object pointed to by
 * elements in the hash table.
 * Types:  BlackBoxi - Black Box matrix, single-word size entries
 *           Class:  LinBox::TriplesBB<LinBox::Modular<long>, std::vector<long> > (TriplesBBi)
 *         BlackBoxI - Black Box matrix, multi-word size entries
 *           Class:  LinBox::TriplesBB<LinBox::Modular<integer>, std::vector<integer> > (TriplesBBI)
 *         SmallV - STL vector, single-word size entries
 *           Class:  std::vector<long> (Vectorl)
 *         LargeV - STL vector, multi-word size entries
 *           Class:  std::vector<integer> (VectorI)
 */

const int BlackBoxi = 1;
const int BlackBoxI = 2;
const int SmallV = 3;
const int LargeV = 4;

extern "C"
{
  ALGEB End(MKernelVector kv, ALGEB* args);
  ALGEB initBB(MKernelVector kv, ALGEB* args);
  ALGEB getMatrix(MKernelVector kv, ALGEB* args);
  ALGEB killMatrix(MKernelVector kv, ALGEB* args);
  ALGEB initV(MKernelVector kv, ALGEB* args);
  ALGEB getVector(MKernelVector kv, ALGEB* args);
  ALGEB killVector(MKernelVector kv, ALGEB* args);
  ALGEB apply(MKernelVector kv, ALGEB* args);
  ALGEB applyT(MKernelVector kv, ALGEB* args);
  ALGEB rank(MKernelVector kv, ALGEB* args);
  ALGEB det(MKernelVector kv, ALGEB* args);
  ALGEB minpoly(MKernelVector kv, ALGEB* args);
  //  ALGEB ssolve(MKernelVector kv, ALGEB* args);
}

ALGEB & LiToM(MKernelVector &, const integer &, ALGEB &);
integer & MtoLI(MKernelVector &, integer &, const ALGEB &);

std::map<int,void*> hashTable;
std::map<int, int> typeTable;


/*  End -- LBmaple back-end module "destructor".  Cleans up all external data
 *  when the LBmaple module goes out of scope and is garbage collected.
 *  Pre-Condition:  attached maple module has gone out of scope and is
 *         garbage collected.
 *  Post-Condition: all dynamically allocated data is freed
 */

extern "C"
{
  ALGEB End(MKernelVector kv, ALGEB* args)
  {
    /* Four types of objects stored in memory */

    TriplesBBi* BBip;
    TriplesBBI* BBIp;
    Vectorl* Vip;
    VectorI* VIp;

   std::map<int,int>::iterator f_i;
   std::map<int,void*>::iterator h_i;

   // For each object in the hashTable . . .
    for( h_i = hashTable.begin(); h_i != hashTable.end(); ++h_i) {

      // Find the objects type and switch appropriately. . .
      f_i = typeTable.find(h_i->first);
      switch( f_i->second ) {

	// In each case, delete the object pointed to
        case BlackBoxi:
	   BBip = (TriplesBBi*) h_i->second;
	   delete BBip;
	   break;

         case BlackBoxI:
	   BBIp = (TriplesBBI*) h_i->second;
	   delete BBIp;
	   break;

         case SmallV:
	   Vip = (std::vector<long> *) h_i->second;
	   delete Vip;
	   break;

         case LargeV:
	   VIp = ( std::vector<integer> *) h_i->second;
	   delete VIp;
	   break;

	   /* Of course there are more cases to follow */
      }

    }
    // When all objects are deleted, return
    return ToMapleNULL(kv);
  }
}

/* initBB - Backend "constructor" for Maple BlackBox objects. Takes as input 3 types of objects
 * from Maple:  A Maple NAG sparse Matrix, 3 vectors defining a word-size entry matrix, or
 * 3 vectors defining a multi-word size entry matrix
 * Declares the object in dynamic memory and hashes into hashTable and typeTable
 * using the key supplied by Maple.
 * Pre-Condition:  Called from Maple to create blackbox type
 * Post-Condition:  Blackbox object exists in dynamic memory, can be called upon with key
 * supplied by Maple
 */

extern "C"
{
  ALGEB initBB(MKernelVector kv, ALGEB* args)
  {

    // First get flag and key
    // IN this c
    int flag = MapleToInteger32(kv,args[1]), key = MapleToInteger32(kv,args[2]);
    size_t m, n, nonzeros, i;

    // perform an intial search to see if the key is stored in memory.  If so, return the key,
    // the object already exists

    // Otherwise, switch on the type
    switch(flag) {

      // The first type is the Maple NAG sparse matrix.  It is passed straight into the routine
      // the prime is passed first, then the object.  All pertinant data can be lifted using
      // Maple's extensive rtable API
       case 1: {
	 long p;
	 int *data;
	 NAG_INT *rowP, *colP;
	 TriplesBBi* In;

	 p = MapleToInteger32(kv, args[3]);
	 m = RTableUpperBound(kv, args[4], 1);
	 n = RTableUpperBound(kv, args[4], 2);
	 nonzeros = RTableNumElements(kv, args[4]);

	 // Get pointers to the data
	 rowP = RTableSparseIndexRow(kv,args[4],1);
	 colP = RTableSparseIndexRow(kv, args[4],2);
	 data = (int*) RTableDataBlock(kv,args[4]);

	 // Declare a new object on the heap
	 LinBox::Modular<long> modF(p);
	 In = new TriplesBBi(modF, m, n, nonzeros);

	 // Add each entry
	 for(i = 0; i < nonzeros; ++i) {
	   In->addEntry( data[i], rowP[i], colP[i]);
	 }

	 // hash the pointer into the hasTable, the type into the typeTable using the key
	 hashTable.insert(std::pair<int,void*>(key, In));
	 typeTable.insert(std::pair<int,int>( key, BlackBoxi ));
       }
       break;


       // In the second case, the user either A) Created a BlackBox from an existing sparse
       // matrix in Maple, or created a Matrix by passing in parameters and a procedure
       // In either case, the data is formatted in Maple, and a number of parameters are passed
       // to this routine.
       // First is the prime, then the row list, then the column list, the entry list, followed
       // by the rows, columns, and number of non-zero entries.
       case 2: {
	 long p;
	 TriplesBBi* In;

	 p = MapleToInteger32(kv, args[3]);
	 m = (size_t) MapleToInteger32(kv,args[7]);
	 n = (size_t) MapleToInteger32(kv,args[8]);
	 nonzeros = (size_t) MapleToInteger32(kv,args[9]);

	 // Declares field and blackbox
	 LinBox::Modular<long> modF(p);
	 In = new TriplesBBi(modF, m, n, nonzeros);
	 // Populates blackbox w/ entries
	 for(i = 1; i <= nonzeros; i++) {
	   In->addEntry( MapleToInteger32(kv, MapleListSelect(kv, args[4], i)), MapleToInteger32(kv, MapleListSelect(kv, args[5], i)), MapleToInteger32(kv, MapleListSelect(kv, args[6], i)));
	 }

	 hashTable.insert(std::pair<int,void*>(key, (void*) In));
	 typeTable.insert(std::pair<int,int>(key, BlackBoxi));
       }
       break;


       // In this case, a sparse Matrix was input, or a matrix was created w/ a procedure, but
       // the matrix is defined by a multi-word size prime.  Data is passed in the order
       // mentioned for case #2 above.  Notice that the function calls the MtoLI function
       // to convert into a gmp integer
      case 3:{
	integer blank, iPrime;
	VectorI Elements;
	TriplesBBI* In;

	iPrime = MtoLI(kv, iPrime, args[3]);
	m = (size_t) MapleToInteger32(kv,args[7]);
	n = (size_t) MapleToInteger32(kv,args[8]);
	nonzeros = (size_t) MapleToInteger32(kv,args[9]);

	// Declare Field and blackbox
	LinBox::Modular<integer> modF(iPrime);
	In = new TriplesBBI(modF, m, n, nonzeros);
	for(i = 1; i <= nonzeros; i++) {
	  In->addEntry( MtoLI(kv, blank, MapleListSelect(kv, args[4],i)), MapleToInteger32(kv, MapleListSelect(kv, args[5], i)), MapleToInteger32(kv, MapleListSelect(kv, args[6], i)));
	}

	hashTable.insert(std::pair<int,void*>(key, In));
	typeTable.insert(std::pair<int,int>(key,BlackBoxI));
      }
      break;


      // In case some weird operation is requested, just bail.
      default:
	MapleRaiseError(kv, "ERROR!  Confused by request.  No action performed.");
	break;
    }

    return ToMapleNULL(kv);
  }
}





/* initV - Backend "constructor" for LBmaple vector type.
 * Can create either a wordsize entry
 * vector, or a multi-precision vector (distinction made by entry size flag passed
 * from Maple).  Declares the object in dynamic memory and hashes into hashTable and typeTable
 * using the key supplied by Maple.
 * Pre-Condition:  Called from Maple with initialized data in a correct format
 * Post-Condition:  Blackbox object exists in dynamic memory, can be called upon with key
 * supplied by Maple
 */

extern "C" {
  ALGEB initV(MKernelVector kv, ALGEB* args)
  {

    // First get the flag and key to use.  The flag indicates what type of object was passed in
    // the key what key it should be hashed against
    int flag = MapleToInteger32(kv,args[1]), key = MapleToInteger32(kv,args[2]), length, index, i, j, L;
    std::vector<long>* vP;
    std::vector<integer>* VP;
    integer blank;


    // Performs a quick table check to see if the Vector has already been
    // created.  If so, simply return it
    std::map<int,void*>::iterator h_i = hashTable.find(key);
    if( h_i != hashTable.end() ) {
      return ToMapleInteger(kv, (long) key);
    }

    // A few good choices
    switch( flag ) {

      // For this case, a complete maple list of single-word size entries is passed
      // in that will be translated (1-1) into an STL vector and hased in.
       case 1:
	 // First get the vector length
	 length = MapleToInteger32(kv, args[4]);

	 // dynamically allocate the new object, and reserve space in memory for savings
	 vP = new std::vector<long>;
	 vP->reserve(length);

	 // initialize entries and hash into hashTable
	 for(i = 1; i <= length; i++) {
	   vP->push_back( MapleToInteger32(kv, MapleListSelect(kv, args[3], (M_INT) i)));
	 }
	 hashTable.insert(std::pair<int, void*>(key, vP));
	 break;

	 // For thsi case, a complete maple list of  multi-word sized entries is passed
	 // in that will be translate (1-1) into an STL vector and hashed in.
       case 2:
	 length = MapleToInteger32(kv, args[4]);
	 VP = new std::vector<integer>;
	 VP->reserve(length);
	 for(i = 1; i <= length; i++) {
	   VP->push_back( MtoLI(kv, blank, MapleListSelect(kv, args[3], (M_INT) i)));
	 }

	 hashTable.insert(std::pair<int,void*>(key, VP));
	 break;


      // In this case, the vector passed in was parsed into a two lists, the first containing
      // the data entries, the second containing the indice in which this entry goes.  Each of
      // these entries are single word integers
       case SmallV:

	 // Initialzie variables, get length and first indice
	 i = 1; j = 1;
	 length = MapleToInteger32(kv, args[5]);
	 L = MapleToInteger32(kv, args[6]);


	 // For speed, ensures only 1 memory management adjustment in creation
	 vP = new std::vector<long>;
	 vP->reserve(L);
	 for( ; i <= length; ++i) {

	   index = MapleToInteger32(kv, MapleListSelect(kv, args[3], i));

	   // All lists are passed back in sparse form, so just add zeros
	   // until you get to the correct index
	   while(j < index) {
	     vP->push_back(0L);
	     ++j;
	   }

	   // Add the nonzero entry
	   vP->push_back( (long) MapleToInteger32(kv,MapleListSelect(kv,args[4], i)));

	   // ++index
	   ++j;
	 }


	 // Now populate the rest of V with 0's
	 for(i = index + 1 ; i <= L; ++i)
	   vP->push_back( 0L );

	 hashTable.insert(std::pair<int,void*>(key, vP));
	 break;

    // Multi-word integer case.  As above, 2 sparse lists are passed in.  The first contains
    // all the non-zero entries.  The second contains the indice of each element.
    case LargeV:

      i = 1; j = 1;
      length = MapleToInteger32(kv,args[5]);
      L = MapleToInteger32(kv, args[6]);

      VP = new std::vector<integer>;
      VP->reserve(L );
      for( ; i <= length; ++i ) {

	index = MapleToInteger32(kv, MapleListSelect(kv, args[3], i));


	while(j < index) {
	  VP->push_back( integer(0) );
	  ++j;
	}

	// Add a nonzero entry
	VP->push_back( MtoLI(kv, blank, MapleListSelect(kv,args[4], i) ));

	++j;
	// ++index :-)
	index = MapleToInteger32(kv,MapleListSelect(kv,args[3],i));
      }

      for( ; i <= L; ++i)
	VP->push_back( integer(0) );


      hashTable.insert(std::pair<int,void*>(key, VP));
      break;

      // In case a weird action is performed, just bail out
      default:
	MapleRaiseError(kv,"ERROR! Confused by request.  Bailing out.");

    }

    // Corrects for list cases.  Ensures that the "flag" variable is either 3 or 4, so it will
    // properly correspond to the type settings in the type table
    // insert the type in the typeTable
    if(flag < 3) flag += 2;

    typeTable.insert(std::pair<int,int>(key,flag));
    return ToMapleNULL(kv);
  }
}






/* killMatrix - BackEnd "destructor" for a Maple blackbox element that goes out of scope.
 * Deletes the blackbox.
 * Pre-condition: The object hashed to by key is an active BlackBox matrix.
 * Post-condition:  The BlackBox matrix hashed to by key is freed
 */


extern "C"
{
  ALGEB killMatrix(MKernelVector kv, ALGEB* args)
  {
    char err[] = "ERROR!  Associated blackbox object does not exist!";
    int key = MapleToInteger32(kv,args[1]), flag;

    // Gets the type in question
    std::map<int,int>::iterator i = typeTable.find(key);

    // If the associated object code is not there, we have a problem.  Bail.
    if( i == typeTable.end() )
      MapleRaiseError(kv,err);

    // Otherwise we're good
    flag = i->second;

    // Erase the entry from the typeTable
    typeTable.erase(key);

    // hash out the pointer
    std::map<int, void*>::iterator h_i = hashTable.find(key);
    if( h_i != hashTable.end() ) {

      switch( flag ) { // Free the data, whatever it is
        case BlackBoxi:{
	  TriplesBBi* ptr = (TriplesBBi*) h_i->second;
	  delete ptr;
	}
	break;

        case BlackBoxI: {
	  TriplesBBI* ptr = (TriplesBBI*) h_i->second;
	  delete ptr;
	}
	break;

      }

      // remove the pointer and bail
      hashTable.erase(key);
    }

    return ToMapleNULL(kv);
  }
}






/* killVector - Backend "destructor" for a linbox vector object in Maple.  Takes as input a key to the
 * data.  Frees up dynamic memory and removes the object from the hash tables
 * Pre-Condition: key is a key to a valid vector object
 * Post-Condition:  Dynamic memory is freed and the object is removed from the hash tables.
 */

extern "C"
{
  ALGEB killVector(MKernelVector kv, ALGEB* args)
  {
    // Gets the key
    int key = MapleToInteger32(kv, args[1]), flag;
    char err[] = "ERROR! Associated Vector object does not exist!";

    // Checks to see if the object is there.  If not, bail
    std::map<int,int>::iterator f_i = typeTable.find(key);
    if( f_i == typeTable.end() )
      MapleRaiseError( kv, err);

    // We're good otherwise
    flag = f_i->second;

    // In case the flag is there but the data isn't
    typeTable.erase(key);

    // Get ahold of the data
    std::map<int,void*>::iterator h_i = hashTable.find(key);
    if( h_i != hashTable.end() ) {

      // Switch according to the vector type, small or large
      switch( flag) {

         case SmallV: {
	   Vectorl* ptr = (Vectorl*) h_i->second;
	   delete ptr;
	 }
	 break;

         case LargeV: {
	   VectorI* ptr = (VectorI*) h_i->second;
	   delete ptr;
	 }
	 break;

      }

      hashTable.erase(key);
    }

    // Do nothing if it isn't in there
    return ToMapleNULL(kv);
  }
}







/* getMatrix - Returns an rtable object initialized with the data from the blackbox.  This is a bit
 * of a hack, as the differences between Maple versions becomes an issue.  There are two methods:
 * 1)  (if using Maple 7) use string tools to build a string that embodies the proper command to
 * declare an rtable object, and invoke this command with the MapleEvalStatement() function, or
 * 2)  use Maple's RTable API to create a new RTable object (this works in Maple 8, and I'm
 * presuming Maple 6 and 5.  Special thanks to tech support at Maple Soft for informing me that in
 * maple 7, the RTableCreate() function has been disabled).
 * Pre-Condition:  The key provided by the call hashes to a valid BlackBox object
 * Post-Condition:  An RTable object initialized with the entries of the Black Box is returned to Maple
 */

extern "C"
{
  ALGEB getMatrix(MKernelVector kv, ALGEB* args)
  {
    // Get the key
    int key = MapleToInteger32(kv,args[1]), flag;
    char err[] = "ERROR!  The associated BlackBox object does not exist!";
    M_INT index[2], bound[4];
    RTableData d;
    ALGEB rtable, blank;
    RTableSettings s;
    std::vector<size_t> Row, Col;
    std::vector<size_t>::const_iterator r_i, c_i;
    char MapleStatement[100] = "rtable(1..";


    // Get the data type of the blackbox
    std::map<int,int>::iterator f_i = typeTable.find(key);
    if( f_i == typeTable.end() ) // In case the blackbox isn't there
      MapleRaiseError(kv,err);
    flag = f_i->second; // Otherwise, get the blackbox type

    // Check that the data is there
    std::map<int,void*>::iterator h_i = hashTable.find(key);
    if(h_i != hashTable.end() ) {

      // Switch according to mode - regular or "special fix" mode
      switch( MapleToInteger32(kv, args[3])) {

      case 1: // This is the Maple 7 case, "special fix" mode
	      // Use the EvalMapleStatement() to call the rtable constructor in the
	      // Maple environment

	   // Switch according to the type
	   switch(flag) {
	   case BlackBoxi:{ // For single word entry matrices

	     // Extract the necessary data
	     TriplesBBi* BB = (TriplesBBi*) h_i->second;
	     Vectorl Data = BB->getData();
	     Row = BB->getRows();
	     Col = BB->getCols();
	     Vectorl::const_iterator d_i;

	     // Builds the statement that will be used in the Maple 7 callback

	     sprintf(MapleStatement + strlen(MapleStatement), "%d", BB->rowdim() );
	     strcat(MapleStatement, ",1..");

	     sprintf(MapleStatement + strlen(MapleStatement), "%d", BB->coldim() );
	     strcat(MapleStatement, ", subtype=Matrix, storage=sparse);");

	     // Perform the callback
	     rtable = kv->evalMapleStatement(MapleStatement);

	     // Insert each non-zero entry
	     for(d_i = Data.begin(), r_i = Row.begin(), c_i = Col.begin(); r_i != Row.end(); ++d_i, ++c_i, ++r_i) {
	       index[0] = *r_i; index[1] = *c_i;
	       d.dag = ToMapleInteger(kv, *d_i); // d is a union, dag is the
	                                        // ALGEB union field
	       RTableAssign(kv, rtable, index, d);
	     }
	   }
	   break;

	   case BlackBoxI: { // For multi-word size matrix types
	     TriplesBBI* BB = (TriplesBBI*) h_i->second;
	     VectorI Data = BB->getData();
	     VectorI::const_iterator d_i;

	     // Build and execute the Maple callback
	     sprintf(MapleStatement + strlen(MapleStatement), "%d", BB->rowdim() );
	     strcat(MapleStatement, ", 1..");
	     sprintf(MapleStatement + strlen(MapleStatement), "%d", BB->coldim() );
	     strcat(MapleStatement, ", subtype=Matrix, storage=sparse);");
	     rtable = kv->evalMapleStatement(MapleStatement);

	     for(d_i = Data.begin(), r_i = Row.begin(), c_i = Col.begin(); r_i != Row.end(); ++d_i, ++r_i, ++c_i) {
	       index[0] = *r_i; index[1] = *c_i;

	       //    * Okay, here's how this line works.  Basically,
	       //    * in order to set the entries of this RTable to
	       //    * multi-precision integers, I have to first use my own conversion
	       //    * method, LiToM, to convert the integer entry to a ALGEB structure,
	       //    * then do a callback into Maple that calls the ExToM procedure,
	       //    * which converts the results of LiToM into a Maple multi-precision
	       //    * integer. At the moment, this is the best idea I've got as to
	       //    * how to convert a GMP integer into a Maple representation in one shot.
	       //    *

	       d.dag = EvalMapleProc(kv,args[2],1,LiToM(kv, *d_i, blank));
	       RTableAssign(kv, rtable, index, d);
             }
	   }
	   break;

	   // In this case the object is not a BlackBox type
	   default:
	     MapleRaiseError(kv,err);
	     break;
	   }
	break;

      case 2: // Okay, here is the normal case.
	      // Use RTableCreate to create a Maple rtable object

	    kv->rtableGetDefaults(&s);
	    // Get default settings - set datatype to Maple,
	    // DAGTAG to anything

	    s.subtype = RTABLE_MATRIX; // Subtype set to Matrix
	    s.storage = RTABLE_SPARSE; // Storage set to sparse
	    s.num_dimensions = 2; // What do you think this means :-)
	    bound[0] = bound[2] = 1; // Set the lower bounds of each dimension to 0, which for maple is 1

	    switch(flag) { // Switch on data type

	    case BlackBoxi:{ // word size entry Matrix
		TriplesBBi* BB = (TriplesBBi*) h_i->second;
		Vectorl Data = BB->getData();
		Row = BB->getRows();
		Col = BB->getCols();
		Vectorl::const_iterator d_i;

		bound[1] = BB->rowdim();
		bound[3] = BB->coldim();
		rtable = kv->rtableCreate(&s, NULL, bound); // This is the RTableCreate function, it's
		                                            // just the one that works

		// Assign all the non-zero rows
		for( d_i = Data.begin(), r_i = Row.begin(), c_i = Col.begin(); r_i != Row.end(); ++d_i, ++c_i, ++r_i) {
		  index[0] = *r_i; index[1] = *c_i;
		  d.dag = ToMapleInteger(kv, *d_i); // d is a union, dag is the
	                                        // ALGEB union field
		  RTableAssign(kv, rtable, index, d);
		}
	      }
	      break;

	    case BlackBoxI: { // For multi-word entry Matrices
	      TriplesBBI* BB = (TriplesBBI*) h_i->second;
	      VectorI Data = BB->getData();

	      // Setup the Create() call
	      VectorI::const_iterator d_i;
	      Row = BB->getRows();
	      Col = BB->getCols();
	      bound[1] = BB->rowdim();
	      bound[3] = BB->coldim();
	      rtable = kv->rtableCreate(&s, NULL, bound); // Create an empty RTable

	      // Populate the RTable using the callback method described below
	      for(d_i = Data.begin(), r_i = Row.begin(), c_i = Col.begin(); r_i != Row.end(); ++d_i, ++r_i, ++c_i) {
		index[0] = *r_i; index[1] = *c_i;

	    //    * Okay, here's how this line works.  Basically,
	   //    * in order to set the entries of this RTable to
	   //    * multi-precision integers, I have to first use my own conversion
	   //    * method, LiToM, to convert the integer entry to a ALGEB structure,
	   //    * then do a callback into Maple that calls the ExToM procedure,
	   //    * which converts the results of LiToM into a Maple multi-precision
	   //    * integer. At the moment, this is the best idea I've got as to
	   //    * how to convert a GMP integer into a Maple representation in one shot.

	      d.dag = EvalMapleProc(kv,args[2],1,LiToM(kv, *d_i, blank));
	      RTableAssign(kv, rtable, index, d);
	      }
	    }
	  break;
       }
      }
    }
    else
      MapleRaiseError(kv,err);

    return rtable;
  }
}





/* getVector - Returns a Maple vector initalized with the data from the LinBox vector.  As above,
 * is a hack in order to support two differing versions in Maple.  Normally, would use the RTableCreate()
 * function to properly create the Vector, but in Maple7 this isn't supported.  The code is essentially
 * the same as above, accept that we are building a vector rather than a Matrix
 * Pre-Condition:  the key provided hashes to a viable vector object in the hashTable
 * Post-Condition:  A properly initalized vector object is returned to maple
 */

extern "C"
{
  ALGEB getVector(MKernelVector kv, ALGEB* args)
  {
    // Get the key, declare variables
    int key = MapleToInteger32(kv,args[1]), flag;
    char err[] = "ERROR!  The associated Vector object does not exist!";
    M_INT index, bound[2];
    RTableData d;
    RTableSettings s;
    ALGEB rtable, blank;
    char MapleStatement[100] = "rtable(1..";


    // Check to see if the object pointed to by key is in the type table.  If not, panic
    std::map<int,int>::iterator f_i = typeTable.find(key);
    if(f_i == typeTable.end() ) {
      MapleRaiseError(kv, err);
    }

    // Otherwise, we have our object
    flag = f_i->second;

    // Get a pointer to the actual data
    std::map<int,void*>::iterator h_i = hashTable.find(key);
    if(h_i != hashTable.end() ) {

      // Diverge over whether we are using maple 7 or 8 ( and 5 & 6)
      // in Maple, arg 3 is a flag indicating which method to use
      switch( MapleToInteger32(kv, args[3])) {

	// In this case, Maple 7 is being used, we have to construct a call using "EvalMapleStatement()"
	// to call the RTable constructor
         case 1:

	   switch(flag) {
	   case SmallV:{
	     // Get the vector
	     Vectorl* V = (Vectorl*) h_i->second;
	     Vectorl::const_iterator V_i;

	     // Create the Maple object
	     sprintf(MapleStatement + strlen(MapleStatement), "%d", V->size() );
	     strcat(MapleStatement, ", subtype=Vector[column], storage=sparse)");
	     rtable = kv->evalMapleStatement(MapleStatement);

	     // populate the Maple vector w/ the entries from V above
	     for(index = 1, V_i = V->begin(); V_i != V->end(); ++V_i, ++index) {
	       d.dag = ToMapleInteger(kv, *V_i); // d is a union, dag is the
	                                         // ALGEB union field
	       RTableAssign(kv, rtable, &index, d);
	     }
	   }
	   break;

	   case LargeV: {
	     // This part works the same way as above
	     VectorI* V = (VectorI*) h_i->second;
	     VectorI::const_iterator V_i;
	     sprintf(MapleStatement + strlen(MapleStatement), "%d", V->size() );
	     strcat(MapleStatement, ",subtype=Vector[column], storage=sparse)");
	     rtable = kv->evalMapleStatement(MapleStatement);

	     // Use maple callback to call the procedure from Maple that translates a gmp integer
	     // into a large maple integer.  Then put this into the Maple vector
	     for(index = 1, V_i = V->begin(); V_i != V->end(); ++V_i, ++index) {

	       /* Okay, here's how this line works.  Basically,
		* in order to set the entries of this RTable to
		* multi-precision integers, I have to first use my own conversion
		* method, LiToM, to convert the integer entry to a ALGEB structure,
		* then do a callback into Maple that calls the ExToM procedure,
		* which converts the results of LiToM into a Maple multi-precision
		* integer. At the moment, this is the best idea I've got as to
		* how to convert a GMP integer into a Maple representation in one shot.
		*/

	       d.dag = EvalMapleProc(kv,args[2],1,LiToM(kv, *V_i, blank));
	       RTableAssign(kv, rtable, &index, d);
	     }
	   }
	   break;

	   default:
	     MapleRaiseError(kv, err);
	     break;
	   }
	   break;

	   // In this case, use the simpler RTableCreate function, rather than building a string
	   // that must be parsed by maple

	   case 2:

	     kv->rtableGetDefaults(&s); // Get default settings - set datatype to Maple,
                               // DAGTAG to anything
	     s.subtype = 2; // Subtype set to column vector
	     s.storage = 4; // Storage set to rectangular
	     s.num_dimensions = 1; // What do you think this means :-)
	     bound[0] = 1; // Set the lower bounds of each dimension to 0

	     switch(flag) {// Switch on data type of vector
	     case SmallV:{ // single word integer entry vector
	       Vectorl* V = (Vectorl*) h_i->second;
	       Vectorl::const_iterator V_i;
	       bound[1] = V->size();
	       rtable = kv->rtableCreate(&s, NULL, bound); // Create the Maple vector

	       for(index = 1, V_i = V->begin(); V_i != V->end(); ++V_i, ++index) {
		 d.dag = ToMapleInteger(kv, *V_i); // d is a union, dag is the
	                                        // ALGEB union field
		 RTableAssign(kv, rtable, &index, d);
	       }
	     }
	     break;

	     case LargeV: { // Same as above for multi-word integer entry vector
	       VectorI* V = (VectorI*) h_i->second;
	       VectorI::const_iterator V_i;
	       bound[1] = V->size();
	       rtable = kv->rtableCreate(&s, NULL, bound);

	       for(index = 1, V_i = V->begin(); V_i != V->end(); ++V_i, ++index) {


		 /* Okay, here's how this line works.  Basically,
		  * in order to set the entries of this RTable to
		  * multi-precision integers, I have to first use my own conversion
		  * method, LiToM, to convert the integer entry to a ALGEB structure,
		  * then do a callback into Maple that calls the ExToM procedure,
		  * which converts the results of LiToM into a Maple multi-precision
		  * integer. At the moment, this is the best idea I've got as to
		  * how to convert a GMP integer into a Maple representation in one shot.
		  */

		 d.dag = EvalMapleProc(kv,args[2],1,LiToM(kv, *V_i, blank));
		 RTableAssign(kv, rtable, &index, d);
	       }
	     }
	     break;

	     default:
	       MapleRaiseError(kv, err);
	       break;
	     }
	     break; // breaks case 2.
	     // This was causing a wicked error :-)


      default:
	MapleRaiseError(kv, err);
	break;

      }
    }
    else {
      MapleRaiseError(kv, err);
    }

    return rtable;
  }
}






/* apply - Application function,  Takes a Blackbox A and Vector X and performs the apply Ax.
 * There's alot of code here that essentially leads to simply calling the BlackBox's apply routine,
 * but there is alot that can go wrong with this call.
 * Pre-Condition:  BlackBox & Vector objects exist in the hashtable; Blackbox data-type can handle
 *      size of vector (ie - not trying to apply a large entry vector to a small entry Matrix
 * Post-Condition: New Vector is produced, put into hash table for future calls, access key returned
 *      to Maple
 */

extern "C"
{
  ALGEB apply(MKernelVector kv, ALGEB* args)
  {
    // There are probably a million better random number generators, but for the moment I use this one
    int BBKey = MapleToInteger32(kv,args[1]), VKey = MapleToInteger32(kv,args[2]), bflag, vflag, nKey;
    VectorI *tempIV, newV;
    Vectorl *tempiV, *vp;
    char MisMatchErr[] = "ERROR!  The Vector elements are not in the field of the Matrix!";
    char BBnoFind[] = "ERROR!  The associated Blackbox object does not exist!";
    char VectnoFind[] = "ERROR!  The associated Vector object does not exist!";

    std::map<int,int>::iterator f_i;
    std::map<int, void*>::iterator h_i;

    // Get the random key from Maple
    nKey = MapleToInteger32(kv, args[3]);

    // Hash out the blackbox type
    f_i = typeTable.find(BBKey);
    if( f_i == typeTable.end() )
      MapleRaiseError(kv, BBnoFind);
    else
      bflag = f_i->second;

    // Hash out the vector type
    f_i = typeTable.find(VKey);
    if( f_i == typeTable.end() )
      MapleRaiseError(kv, VectnoFind);
    else
      vflag = f_i->second;

    // check that there isn't a type mismatch
    if( bflag == BlackBoxi && vflag == LargeV)
      // If this comes up, it means we have a large int vector, and a small int blackbox.  Bad
      MapleRaiseError(kv, MisMatchErr);
    else if( bflag == BlackBoxI && vflag == SmallV) {

      // If this comes up, it means we have a small int vector, and a large int blackbox.  Not fatal.
      h_i = hashTable.find(BBKey);
      if( h_i == hashTable.end() )
	MapleRaiseError(kv, BBnoFind);
      TriplesBBI* BB = (TriplesBBI*) h_i->second;

      h_i = hashTable.find(VKey);
      if(h_i == hashTable.end() )
	MapleRaiseError(kv, VectnoFind);

      vp = (Vectorl*) h_i->second;

      // Converts the small vector into a large vector
      // Note, the reserve ensures there will only be at most one reallocation neccessary
      newV.reserve(vp->size() );
      for(Vectorl::iterator i = vp->begin(); i != vp->end(); ++i)
	newV.push_back(integer(*i));

      // Puts tempIV on the stack, so that it will persist.  Also, set the correct size
      tempIV = new VectorI;
      tempIV->resize( BB->coldim() );

      // Perform the apply
      BB->apply( *tempIV, newV);
      hashTable.insert(std::pair<int,void*>(nKey, tempIV));
      typeTable.insert(std::pair<int,int>(nKey, LargeV));

    }
    else { // In this case, the types match (single word integer matrix w/ single word integer vector,
           // or multi-word entry w/ multi-word entry.  This simplifies the code a good bit

      switch( bflag) {

         case BlackBoxi: {

	   // Get the BlackBox
	   h_i = hashTable.find(BBKey);
	   if( h_i == hashTable.end() )
	     MapleRaiseError(kv, BBnoFind);
	   TriplesBBi* BB = (TriplesBBi*) h_i->second;

	   // Get the Vector
	   h_i = hashTable.find(VKey);
	   if( h_i == hashTable.end() )
	     MapleRaiseError(kv, VectnoFind);
	   vp = (Vectorl*) h_i->second;

	   // Perform the apply.  Note, sets the temp vector's size to the columnar dimension
	   tempiV = new Vectorl;
	   tempiV->resize( BB->coldim() );
	   BB->apply(*tempiV, *vp);

	   // Hash the results
	   hashTable.insert(std::pair<int, void*>(nKey, tempiV));
	   typeTable.insert(std::pair<int,int>(nKey,SmallV));
	 }
	 break;

         case BlackBoxI: {

	   // Get the BlackBox
	   h_i = hashTable.find(BBKey);
	   if( h_i == hashTable.end() )
	     MapleRaiseError(kv, BBnoFind );
	   TriplesBBI* BB = (TriplesBBI*) h_i->second;

	   // Get the vector
	   h_i = hashTable.find(VKey);
	   if( h_i == hashTable.end() )
	     MapleRaiseError(kv, VectnoFind );
	   VectorI* Vp = (VectorI*) h_i->second;

	   // Perform the apply
	   tempIV = new VectorI;
	   tempIV->resize( BB->coldim() );
	   BB->apply( *tempIV, *Vp);

	   // Hash in the results
	   hashTable.insert(std::pair<int,void*>(nKey, tempIV));
	   typeTable.insert(std::pair<int,int>(nKey,LargeV));
	 }
	 break;
      }
    }

    // If you get to this point, something has been hashed in, so return the key
    return ToMapleInteger(kv, (long) nKey);
  }
}




/* applyT - Application function.  Takes a Blackbox A and Vector X and performs the transpose apply,
 * Atx.  There's alot of code here that essentially leads to simply calling the BlackBox's apply routine,
 * but there is alot that can go wrong with this call.
 * Pre-Condition:  BlackBox & Vector objects exist in the hashtable; Blackbox data-type can handle
 *      size of vector (ie - not trying to apply a large entry vector to a small entry Matrix
 * Post-Condition: New Vector is produced, put into hash table for future calls, access key returned
 *      to Maple
 */

extern "C"
{
  ALGEB applyT(MKernelVector kv, ALGEB* args)
  {

    int BBKey = MapleToInteger32(kv,args[1]), VKey = MapleToInteger32(kv,args[2]), bflag, vflag, nKey;
    char misMatchErr[] = "ERROR! Vector not in field of blackbox!";
    char BBnoFindErr[] = "ERROR! The associated blackbox object does not exist!";
    char VectNoFindErr[] = "ERROR! The associated vector object does not exist!";
    VectorI *tempIV, newV;
    Vectorl *tempiV, *vp;

    std::map<int,int>::iterator f_i;
    std::map<int, void*>::iterator h_i;


    // Creates the new key
    nKey = MapleToInteger32(kv, args[3]);

    // Hash out the blackbox type
    f_i = typeTable.find(BBKey);
    if( f_i == typeTable.end() )
      MapleRaiseError(kv,BBnoFindErr);
    else
      bflag = f_i->second;

    // Hash out the vector type
    f_i = typeTable.find(VKey);
    if( f_i == typeTable.end() )
      MapleRaiseError(kv,VectNoFindErr);
    else
      vflag = f_i->second;

    // check that there isn't a type mismatch
    if( bflag == BlackBoxi && vflag == LargeV)
      // If this comes up, it means we have a large int vector, and a small int blackbox.  Bad
      MapleRaiseError(kv, misMatchErr);
    else if( bflag == BlackBoxI && vflag == SmallV) {

      // If this comes up, it means we have a small int vector, and a large int blackbox.  Not fatal.


      // Bail if data can't be found
      h_i = hashTable.find(BBKey);
      if( h_i == hashTable.end() )
	MapleRaiseError(kv, BBnoFindErr);
      TriplesBBI* BB = (TriplesBBI*) h_i->second;

      h_i = hashTable.find(VKey);
      if(h_i == hashTable.end() )
	MapleRaiseError(kv, VectNoFindErr);

      // Build a large entry vector and populate it with the converted data from the small entry vector
      vp = (Vectorl*) h_i->second;
      newV.reserve( vp->size() );
      for(Vectorl::iterator i = vp->begin(); i != vp->end(); ++i)
	newV.push_back(integer(*i));

      // Puts tempIV on the stack, so that it will persist
      tempIV = new VectorI;
      tempIV->resize( BB->rowdim() );

      // Perform the apply
      BB->applyTranspose( *tempIV, newV);
      hashTable.insert(std::pair<int,void*>(nKey, tempIV));
      typeTable.insert(std::pair<int,int>(nKey, LargeV));

    }
    else { // For types that match (small Blackbox w/ small Vector or large blackbox w/ large vector)

      switch( bflag) {

         case BlackBoxi: {

	   // Get the BlackBox
	   h_i = hashTable.find(BBKey);
	   if( h_i == hashTable.end() )
	     MapleRaiseError(kv, BBnoFindErr);
	   TriplesBBi* BB = (TriplesBBi*) h_i->second;

	   // Get the Vector
	   h_i = hashTable.find(VKey);
	   if( h_i == hashTable.end() )
	     MapleRaiseError(kv, VectNoFindErr);
	   vp = (Vectorl*) h_i->second;

	   // Perform the apply
	   tempiV = new Vectorl;
	   tempiV->resize( BB->rowdim() );
	   BB->applyTranspose( *tempiV, *vp);

	   // Hash the results
	   hashTable.insert(std::pair<int, void*>(nKey, tempiV));
	   typeTable.insert(std::pair<int,int>(nKey,SmallV));
	 }
	 break;

         case BlackBoxI: {

	   // Get the BlackBox
	   h_i = hashTable.find(BBKey);
	   if( h_i == hashTable.end() )
	     MapleRaiseError(kv, BBnoFindErr);
	   TriplesBBI* BB = (TriplesBBI*) h_i->second;

	   // Get the vector
	   h_i = hashTable.find(VKey);
	   if( h_i == hashTable.end() )
	     MapleRaiseError(kv, VectNoFindErr);
	   VectorI* Vp = (VectorI*) h_i->second;

	   // Perform the apply
	   tempIV = new VectorI;
	   tempIV->resize( BB->rowdim() );
	   BB->applyTranspose( *tempIV, *Vp);

	   // Hash in the results
	   hashTable.insert(std::pair<int,void*>(nKey, tempIV));
	   typeTable.insert(std::pair<int,int>(nKey,LargeV));
	 }
	 break;
      }
    }

    // If you get to this point, something has been hashed in, so return the key
    return ToMapleInteger(kv, (long) nKey);
  }
}





/* rank - one of the three main "solution" functions, computes the rank of a BlackBox in the hash table
 * Rank is computed using the rank function from the linbox repository
 * Pre-Condition: key provided hashes to a valid blackbox object
 * Post-Condition: the rank is returned to maple as an int
 */

extern "C"
{
  ALGEB rank(MKernelVector kv, ALGEB* args)
  {
    int key = MapleToInteger32(kv,args[1]), flag;
    unsigned long result;
    char err[] = "ERROR!  Associated blackbox object does not exist!";
    TriplesBBi* BBi;
    TriplesBBI* BBI;

    std::map<int,int>::iterator f_i = typeTable.find(key);
    std::map<int,void*>::iterator h_i;

    if(f_i == typeTable.end() )
      MapleRaiseError(kv,err);
    else
      flag = f_i->second;


    h_i = hashTable.find(key);
    if( h_i != hashTable.end() ) {
      switch( flag ) {
	case BlackBoxi:
	  BBi = (TriplesBBi*) h_i->second;
	  return ToMapleInteger(kv, LinBox::rank(result, *BBi, BBi->field() ) ); // <- actual computation starts
	  break;                                                         //    here :-)

         case BlackBoxI:
	  BBI = (TriplesBBI*) h_i->second;
	  return ToMapleInteger(kv,LinBox::rank(result, *BBI, BBI->field() ));
	  break;
      }
    }
    else
      MapleRaiseError(kv,err);
  }
}


/* det - Another major "application" function, this computes the determinant of the blackbox mod the
 * prime using the LinBox determinant methods, of which I know very little.  Look at the
 * documentation in the solutions directory for more details.
 * Pre-condition: Key hashes to a valid BlackBox object
 * Post-condition: The modular determinant of the blackbox is returned as a maple integer
 */

extern "C"
{
  ALGEB det(MKernelVector kv, ALGEB* args)
  {
    char err[] = "ERROR!  Associated blackbox object does not exist!";
    int key = MapleToInteger32(kv,args[1]), flag;
    std::map<int,int>::iterator f_i;
    std::map<int,void*>::iterator h_i;
    long resulti;
    integer resultI;
    ALGEB blank;
    TriplesBBi *BBi;
    TriplesBBI *BBI;


    // Get the blackbox
    f_i = typeTable.find(key);
    if( f_i == typeTable.end() )
      MapleRaiseError(kv,err);
    else
      flag = f_i->second;

    // If it's there
    h_i = hashTable.find(key);
    if( h_i != hashTable.end() ) {

      switch( flag ) { // switch on the type
       case BlackBoxi:
	 BBi = (TriplesBBi*) h_i->second;
	 return ToMapleInteger(kv, LinBox::det(resulti, *BBi, BBi->field() ) );
	 break;

        case BlackBoxI:
	 BBI = (TriplesBBI*) h_i->second;
	 return LiToM(kv,  LinBox::det(resultI, *BBI, BBI->field() ), blank);

	break;
      }
    }
    else
      MapleRaiseError(kv, err);
  }
}





/* minpoly - Compute the minimal polynomial of the BlackBox matrix using linbox's minpoly solution.
 * I have no clue how it works, just that it does.  This one is slightly more complicated than the
 * other two above
 * Pre-condition:  Key maps to a valid BlackBox object in the hashTable
 * Post-condition: A Maple list of the polynomial coefficients, lowest degree first, is returned
 */

extern "C"
{
ALGEB minpoly(MKernelVector kv, ALGEB* args)
  {
    int i;
    ALGEB retlist, blank;
    char err[] = "ERROR!  Associated blackbox object does not exist!";
    int key = MapleToInteger32(kv,args[1]), flag;
    std::map<int,int>::iterator f_i;
    std::map<int,void*>::iterator h_i;

    // Get the data from the hash table
    f_i = typeTable.find(key);
    if( f_i == typeTable.end() )
      MapleRaiseError(kv,err);
    else
      flag = f_i->second;


    h_i = hashTable.find(key);
    if(h_i != hashTable.end() ) { // We've got data
      switch( flag ) {
	// Getting the minimal polynomial is rather complicated, so both instances of this code were
	// wrapped up inside a block, to cut down at the clutter at the top of this function.
	// First declares a vector of the proper type and casts the pointer.  Then computes the minimal
	// polynomial.  It then builds the proper Maple list structure for this application.


         case BlackBoxi: {
	   Vectorl mpreturn;
	    Vectorl::iterator mp_i;
	    TriplesBBi* BB = (TriplesBBi*) h_i->second;
	    LinBox::minpoly( mpreturn, *BB, BB->field() );
	    retlist = MapleListAlloc(kv, mpreturn.size() );
	    for(i = 1, mp_i = mpreturn.begin(); mp_i != mpreturn.end(); ++mp_i, ++i)
	      MapleListAssign(kv, retlist, i, ToMapleInteger(kv, *mp_i));
	 }
	 break;

         case BlackBoxI: {
	   VectorI mpreturn;
	   VectorI::iterator mp_i;
	   TriplesBBI* BB = (TriplesBBI*) h_i->second;
	   LinBox::minpoly( mpreturn, *BB, BB->field() );
	   retlist = MapleListAlloc(kv, mpreturn.size());
	   for(i = 1, mp_i = mpreturn.begin(); mp_i != mpreturn.end(); ++mp_i, ++i)
	     MapleListAssign(kv, retlist, i, LiToM(kv, *mp_i, blank));

	 }
	 break;
      }
    }
    else
      MapleRaiseError(kv,err);

    return retlist;

  }
}

/* LiToM - Converts from a GMP integer to a Maple list of word-size chunks
 * these chunks are passed back into Maple and converted by into maple represented integers
 * Pre-Condition:  In is an initialized GMP number, Out can be initialized (or re-initialized)
 * Post-Condition: Out is an initalized maple integer containing the value of In
 */

ALGEB & LiToM(MKernelVector & kv, const integer & In, ALGEB & Out)
{
  // If we get lucky, this thing fits in one word, which is good, so we just
  // straight convert it
  if(In.size() == 1)
    Out = ToMapleInteger(kv, Integer2long(In) );
  else {
    // Otherwise, create a list, and for each entry, add a word-sized chunk.
    Out = MapleListAlloc(kv, In.size() );
    for(size_t i = 1; i <= In.size(); ++i)
      MapleListAssign(kv,Out,i,ToMapleInteger(kv,In[i-1]));
  }

  // Return the result for redundancy
  return Out;
}

/* MtoLI - conversion between long integers passed from maple into gmp integers
 * Uses Horners method.
 * Pre-Condition:  Out is un-initialized (or can be re-initialized), In is either an int or a
 *   list
 * Post-Condition:  Out contains a GMP version of In.
 */


integer & MtoLI(MKernelVector & kv, integer & Out, const ALGEB &In)
{
  // We just need 1 base.  By default, Maple breaks all integers of a certain size up into
  // mod 10000 chunks.  The chunks are put into a maple list, and then passed into code
  static integer base(10000);
  int num, i;

  // This could be an int
  if( IsMapleInteger32(kv,In) )
    // if so, just convert it
    Out = integer(MapleToInteger32(kv, In));
  else {
    // else, get the number of chunks, and convert the highest chunk of In
    num = MapleToInteger32(kv, MapleListSelect(kv, In, 1) );
    Out = integer(MapleToInteger32(kv, MapleListSelect(kv,In, num) ) );

    for(i = num - 1; i > 1; --i) {
      // For each additional chunk, multiply Out by the base, and add in the new chunk
      Out *= base;
      Out += MapleToInteger32(kv, MapleListSelect( kv, In, i));
    }
  }

  // redundant
  return Out;
}
