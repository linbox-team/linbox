#include "maplec.h"
#include "linbox/field/modular.h"
#include "linbox/integer.h"
// #include "linbox/blackbox/integer.h" <- When linbox supports computations
//                                         over the whole integers
#include "MapleBB.h"
#include "linbox/solutions/rank.h"
#include "linbox/solutions/det.h"
#include "linbox/solutions/minpoly.h"
// #include "linbox/solutions/ssolve.h" <- When linbox supports System solver

#include <vector>
#include <map>
#include <cstdlib>
#include <ctime>

using LinBox::integer;

typedef std::vector<int> Vectori;
typedef std::vector<long> Vectorl;
typedef std::vector<integer> VectorI;
typedef LinBox::MapleBB<LinBox::Modular<long>, Vectorl> MapleBBi;
typedef LinBox::MapleBB<LinBox::Modular<integer>,VectorI> MapleBBI;

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

static std::map<int,void*> hashTable;
static std::map<int, int> typeTable;

extern "C"
{
  ALGEB End(MKernelVector kv, ALGEB* args)
  {
    MapleBBi* BBip;
    MapleBBI* BBIp;
    Vectorl* Vip;
    VectorI* VIp;

   std::map<int,int>::iterator f_i;
   std::map<int,void*>::iterator h_i;

    for( h_i = hashTable.begin(); h_i != hashTable.end(); ++h_i) {
      f_i = typeTable.find(h_i->first);
      switch( f_i->second ) {
         case BlackBoxi: 
	   BBip = (MapleBBi*) h_i->second;
	   delete BBip;
	   break;

         case BlackBoxI:
	   BBIp = (MapleBBI*) h_i->second;
	   delete BBIp;
	   break;

         case SmallV:
	   Vip = (Vectorl*) h_i->second;
	   delete Vip;
	   break;

         case LargeV:
	   VIp = (VectorI*) h_i->second;
	   delete VIp;
	   break;

	   /* Of course there are more cases to follow */
      }

    }    

    return ToMapleNULL(kv);
  }
}

extern "C"
{
  ALGEB initBB(MKernelVector kv, ALGEB* args)
  {
    int flag = MapleToInteger32(kv,args[1]), key = MapleToInteger32(kv,args[2]);
    long i;
    size_t m, n, nonzeros;
    Vectori Row, Col;
  
    // First perform an inital key search, if the key is already there, don't
    // do anything but return the pre-existing key
    std::map<int,void*>::iterator h_i = hashTable.find(key);
    if(h_i != hashTable.end() )
      return ToMapleInteger(kv,(long) key);

    switch(flag) {

       case 1: { 
	 long p;
	 int *data;
	 NAG_INT *rowP, *colP;
	 Vectorl Elements;
	 MapleBBi* In;

	 p = MapleToInteger32(kv, args[3]);
	 m = RTableUpperBound(kv, args[4], 1);
	 n = RTableUpperBound(kv, args[4], 2);
	 nonzeros = RTableNumElements(kv, args[4]);
	 rowP = RTableSparseIndexRow(kv,args[4],1);
	 colP = RTableSparseIndexRow(kv, args[4],2);
	 data = (int*) RTableDataBlock(kv,args[4]);

	 for(i = 0; i < nonzeros; ++i) {
	   Elements.push_back(data[i]);
	   Row.push_back(rowP[i]);
	   Col.push_back(colP[i]);
	 }
	 LinBox::Modular<long> modF(p);
	 In = new MapleBBi(modF, Elements, Row, Col, m, n, nonzeros);
	 hashTable.insert(std::pair<int,void*>(key, In));
	 typeTable.insert(std::pair<int,int>( key, BlackBoxi ));
       }
       break;
   
       case 2: {
	 long p;
	 Vectorl Element;
	 MapleBBi* In;

	 p = MapleToInteger32(kv, args[3]);
	 m = (size_t) MapleToInteger32(kv,args[7]);
	 n = (size_t) MapleToInteger32(kv,args[8]);
	 nonzeros = (size_t) MapleToInteger32(kv,args[9]);
	 for(i = 1; i <= nonzeros; i++) {
	   Element.push_back(MapleToInteger32(kv,MapleListSelect(kv, args[4],i)));
	   Row.push_back(MapleToInteger32(kv,MapleListSelect(kv,args[5],i)));
	   Col.push_back(MapleToInteger32(kv,MapleListSelect(kv,args[6],i)));
	 }
	 LinBox::Modular<long> modF(p);
	 In = new MapleBBi(modF, Element, Row, Col, m, n, nonzeros);
	 hashTable.insert(std::pair<int,void*>(key, In));
	 typeTable.insert(std::pair<int,int>(key, BlackBoxi));
       }
       break;

      case 3:{
	integer blank, iPrime;
	VectorI Elements;
	MapleBBI* In;

	iPrime = MtoLI(kv, iPrime, args[3]);
	m = (size_t) MapleToInteger32(kv,args[7]);
	n = (size_t) MapleToInteger32(kv,args[8]);
	nonzeros = (size_t) MapleToInteger32(kv,args[9]);

	for(i = 1; i <= nonzeros; i++) {
	  Elements.push_back( MtoLI( kv, blank, MapleListSelect(kv,args[4],i) ) );
	  Row.push_back(MapleToInteger32(kv,MapleListSelect(kv,args[5],i) ) );
	  Col.push_back(MapleToInteger32(kv,MapleListSelect(kv,args[6],i) ) );
	}
	LinBox::Modular<integer> modF(iPrime);
	In = new MapleBBI(modF, Elements, Row, Col, m, n, nonzeros);
	hashTable.insert(std::pair<int,void*>(key, In));
	typeTable.insert(std::pair<int,int>(key,BlackBoxI));
      }
      break;
      
      default:
	MapleRaiseError(kv, "ERROR!  Confused by request.  No action performed.");
	break;
    }

    return ToMapleNULL(kv);
  }
}

extern "C" {
  ALGEB initV(MKernelVector kv, ALGEB* args)
  {
    int flag = MapleToInteger32(kv,args[1]), key = MapleToInteger32(kv,args[2]), length, index, i, j;
    std::vector<long>* vP;
    std::vector<integer>* VP;
    integer blank;

    // Performs a quick table check to see if the Vector has already been
    // created.  If so, simply return it
    std::map<int,void*>::iterator h_i = hashTable.find(key);
    if( h_i == hashTable.end() )
      return ToMapleInteger(kv, (long) key);


    switch( flag ) {
      
      // Single word integer case
       case SmallV: 
	 i = 1; j = 1;
	 length = MapleToInteger32(kv, args[5]);
	 index = MapleToInteger32(kv,MapleListSelect(kv, args[3], i));
	 for( ; i <= length; ++i) {

	   // All lists are passed back in sparse form, so just add zeros
	   // until you get to the correct index
	   while(j < index) {
	     vP->push_back(0L);
	     ++j;
	   }
	   
	   // Add the nonzero entry
	   vP->push_back( (long) MapleToInteger32(kv,MapleListSelect(kv,args[4], i)));
	   
	   ++j;
	   index = MapleToInteger32(kv,MapleListSelect(kv, args[3],i));
	 }
	 
	 hashTable.insert(std::pair<int,void*>(key,vP));
	 break;

    // Multi-word integer case
    case LargeV:
      
      i = 1; j = 1;
      length = MapleToInteger32(kv,args[5]);
      index = MapleToInteger32(kv,MapleListSelect(kv,args[3],i));
      for( ; i <= length; ++i ) {

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

      hashTable.insert(std::pair<int,void*>(key,VP));
      break;

      default:
	MapleRaiseError(kv,"ERROR! Confused by request.  Bailing out.");

    }
    typeTable.insert(std::pair<int,int>(key,flag));
    return ToMapleNULL(kv);
  }
}

extern "C"
{
  ALGEB killMatrix(MKernelVector kv, ALGEB* args)
  {
    char err[] = "ERROR!  Associated blackbox object does not exist!";
    int key = MapleToInteger32(kv,args[1]), flag;
    std::map<int,int>::iterator i = typeTable.find(key);

    // If the associated object code is not there, we have a problem
    if( i == typeTable.end() )
      MapleRaiseError(kv,err);

    // Otherwise we're good
    flag = i->second;

    // In case the flag is here, but the data isn't
    typeTable.erase(key);

    std::map<int, void*>::iterator h_i = hashTable.find(key);
    if( h_i != hashTable.end() ) {

      switch( flag ) { // Free the data, whatever it is
        case BlackBoxi:{
	  MapleBBi* ptr = (MapleBBi*) h_i->second;
	  delete ptr;
	}
	break;
	
        case BlackBoxI: {
	  MapleBBI* ptr = (MapleBBI*) h_i->second;
	  delete ptr;
	}
	break;

      }

      hashTable.erase(key);
    }

    return ToMapleNULL(kv);
  }
}

extern "C" 
{
  ALGEB killVector(MKernelVector kv, ALGEB* args) 
  {
    int key = MapleToInteger32(kv, args[1]), flag;
    char err[] = "ERROR! Associated Vector object does not exist!";

    std::map<int,int>::iterator f_i = typeTable.find(key);
    if( f_i == typeTable.end() )
      MapleRaiseError( kv, err);
      
    // We're good otherwise
    flag = f_i->second;

    // In case the flag is there but the data isn't
    typeTable.erase(key);
    
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

    // Do nothing if it isn't in there, we do nothing
    return ToMapleNULL(kv);
  }
}

extern "C"
{
  ALGEB getMatrix(MKernelVector kv, ALGEB* args)
  {
    int key = MapleToInteger32(kv,args[1]), flag;
    char err[] = "ERROR!  The associated BlackBox object does not exist!";
    M_INT bound[4], index[2];
    RTableSettings s;
    RTableData d;
    ALGEB rtable, blank;
    Vectori Row, Col;
    Vectori::const_iterator r_i, c_i;
    
    kv->rtableGetDefaults(&s); // Get default settings - set datatype to Maple,
                               // DAGTAG to anything 
    s.subtype = RTABLE_MATRIX; // Subtype set to Matrix
    s.storage = RTABLE_SPARSE; // Storage set to sparse
    s.num_dimensions = 2; // What do you think this means :-)
    bound[0] = bound[2] = 1; // Set the lower bounds of each dimension to 0, which for maple is 1

    // Get the data type of the blackbox
    std::map<int,int>::iterator f_i = typeTable.find(key);
    if( f_i == typeTable.end() ) // In case the blackbox isn't there
      MapleRaiseError(kv,err);
    flag = f_i->second; // Otherwise, get the blackbox type

    std::map<int,void*>::iterator h_i = hashTable.find(key);
    if(h_i != hashTable.end() ) {
      // Switch according to the type
      switch(flag) {
          case BlackBoxi:{
	    MapleBBi* BB = (MapleBBi*) h_i->second;
	    Vectorl Data = BB->getData();
	    Vectorl::const_iterator d_i;
	    Row = BB->getRows();
	    Col = BB->getCols();
	    bound[1] = BB->rowdim();
	    bound[3] = BB->coldim();
	    rtable = kv->rtableCreate(&s, NULL, bound);

	    for(d_i = Data.begin(), r_i = Row.begin(), c_i = Col.begin(); r_i != Row.end(); ++d_i, ++c_i, ++r_i) {
	      index[0] = *r_i; index[1] = *c_i;
	      d.dag = ToMapleInteger(kv, *d_i); // d is a union, dag is the
	                                        // ALGEB union field
	      RTableAssign(kv, rtable, index, d);
	    }
	  }
	  break;

          case BlackBoxI: {
	    MapleBBI* BB = (MapleBBI*) h_i->second;
	    VectorI Data = BB->getData();
	    VectorI::const_iterator d_i;
	    Row = BB->getRows();
	    Col = BB->getCols();
	    bound[1] = BB->rowdim(); 
	    bound[3] = BB->coldim();
	    rtable = kv->rtableCreate(&s, NULL, bound);
	    
	    for(d_i = Data.begin(), r_i = Row.begin(), c_i = Col.begin(); r_i != Row.end(); ++d_i, ++r_i, ++c_i) {
	      index[0] = *r_i; index[1] = *c_i;
	  
	      /* Okay, here's how this line works.  Basically,
	       * in order to set the entries of this RTable to
	       * multi-precision integers, I have to first use my own conversion
	       * method, LiToM, to convert the integer entry to a ALGEB structure,
	       * then do a callback into Maple that calls the ExToM procedure, 
	       * which converts the results of LiToM into a Maple multi-precision
	       * integer. At the moment, this is the best idea I've got as to 
	       * how to convert a GMP integer into a Maple representation in one shot.
	       */

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
    }
    else
      MapleRaiseError(kv,err);
    
    return rtable;
  }
}

extern "C"
{
  ALGEB getVector(MKernelVector kv, ALGEB* args)
  {
    int key = MapleToInteger32(kv,args[1]), flag;
    char err[] = "ERROR!  The associated Vector object does not exist!";
    M_INT bound[2], index;
    RTableSettings s;
    RTableData d;
    ALGEB rtable, blank;
    
    kv->rtableGetDefaults(&s); // Get default settings - set datatype to Maple,
                               // DAGTAG to anything 
    s.subtype = 2; // Subtype set to column vector
    s.storage = 4; // Storage set to rectangular
    s.num_dimensions = 1; // What do you think this means :-)
    bound[0] = 1; // Set the lower bounds of each dimension to 0

    std::map<int,int>::iterator f_i = typeTable.find(key);
    if(f_i == typeTable.end() )
      MapleRaiseError(kv, err);
    
    flag = f_i->second;

    std::map<int,void*>::iterator h_i = hashTable.find(key);
    if(h_i != hashTable.end() ) {

      switch(flag) {
          case SmallV:{
	    Vectorl* V = (Vectorl*) h_i->second;
	    Vectorl::const_iterator V_i;
	    bound[1] = V->size();
	    rtable = kv->rtableCreate(&s, NULL, bound);

	    for(index = 1, V_i = V->begin(); V_i != V->end(); ++V_i, ++index) {
	      d.dag = ToMapleInteger(kv, *V_i); // d is a union, dag is the
	                                        // ALGEB union field
	      RTableAssign(kv, rtable, &index, d);
	    }
	  }
	  break;

          case BlackBoxI: {
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
    }
    else
      MapleRaiseError(kv, err );
    
    return rtable;
  }
}

extern "C"
{
  ALGEB apply(MKernelVector kv, ALGEB* args)
  {
    // There are probably a million better random number generators, but for the moment I use this one
    srand( time(0) );
    int BBKey = MapleToInteger32(kv,args[1]), VKey = MapleToInteger32(kv,args[2]), bflag, vflag, nKey;
    VectorI *tempIV, newV, *Vp;
    Vectorl *tempiV, *vp;
    char MisMatchErr[] = "ERROR!  The Vector elements are not in the field of the Matrix!";
    char BBnoFind[] = "ERROR!  The associated Blackbox object does not exist!";
    char VectnoFind[] = "ERROR!  The associated Vector object does not exist!";

    std::map<int,int>::iterator f_i;
    std::map<int, void*>::iterator h_i;

    // Creates the new key, once again I need a new random number generator
    nKey = rand();

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
      MapleBBI* BB = (MapleBBI*) h_i->second;
      
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
    else {

      switch( bflag) {
	
         case BlackBoxi: {

	   // Get the BlackBox
	   h_i = hashTable.find(BBKey);
	   if( h_i == hashTable.end() )
	     MapleRaiseError(kv, BBnoFind);
	   MapleBBi* BB = (MapleBBi*) h_i->second;

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
	   MapleBBI* BB = (MapleBBI*) h_i->second;

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


extern "C"
{
  ALGEB applyT(MKernelVector kv, ALGEB* args)
  {

    srand( time(0) );
    int BBKey = MapleToInteger32(kv,args[1]), VKey = MapleToInteger32(kv,args[2]), bflag, vflag, nKey;
    char misMatchErr[] = "ERROR! Vector not in field of blackbox!";
    char BBnoFindErr[] = "ERROR! The associated blackbox object does not exist!";
    char VectNoFindErr[] = "ERROR! The associated vector object does not exist!";
    VectorI *tempIV, newV, *Vp;
    Vectorl *tempiV, *vp;

    std::map<int,int>::iterator f_i;
    std::map<int, void*>::iterator h_i;


    // Creates the new key
    nKey = rand();

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
      h_i = hashTable.find(BBKey);
      if( h_i == hashTable.end() )
	MapleRaiseError(kv, BBnoFindErr);
      MapleBBI* BB = (MapleBBI*) h_i->second;
      
      h_i = hashTable.find(VKey);
      if(h_i == hashTable.end() )
	MapleRaiseError(kv, VectNoFindErr);

      vp = (Vectorl*) h_i->second;
      newV.reserve( vp->size() );
      // Converts the small vector into a large vector
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
    else {

      switch( bflag) {
	
         case BlackBoxi: {

	   // Get the BlackBox
	   h_i = hashTable.find(BBKey);
	   if( h_i == hashTable.end() )
	     MapleRaiseError(kv, BBnoFindErr);
	   MapleBBi* BB = (MapleBBi*) h_i->second;

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
	   MapleBBI* BB = (MapleBBI*) h_i->second;

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

extern "C"
{
  ALGEB rank(MKernelVector kv, ALGEB* args)
  {
    int key = MapleToInteger32(kv,args[1]), flag;
    unsigned long result;
    char err[] = "ERROR!  Associated blackbox object does not exist!";
    std::map<int,int>::iterator f_i = typeTable.find(key);
    std::map<int,void*>::iterator h_i;

    if(f_i == typeTable.end() )
      MapleRaiseError(kv,err);
    else
      flag = f_i->second;


    h_i = hashTable.find(key);
    if( h_i != hashTable.end() ) {
      switch( flag ) {
	case BlackBoxi: {
	  MapleBBi* BB = (MapleBBi*) h_i->second;
	  LinBox::Modular<long> field;
	  BB->getField(field);
	  return ToMapleInteger(kv, LinBox::rank(result, *BB, field));
	}
	break;

      case BlackBoxI: {
	  MapleBBI* BB = (MapleBBI*) h_i->second;
	  LinBox::Modular<integer> field;
	  BB->getField(field);
	  return ToMapleInteger(kv,LinBox::rank(result, *BB, field));
      }
      break;
      
      }
    }
    else 
      MapleRaiseError(kv,err);
  }
}

extern "C" 
{
  ALGEB det(MKernelVector kv, ALGEB* args)
  {
    char err[] = "ERROR!  Associated blackbox object does not exist!";
    int key = MapleToInteger32(kv,args[1]), flag;
    std::map<int,int>::iterator f_i;
    std::map<int,void*>::iterator h_i;

    f_i = typeTable.find(key);
    if( f_i == typeTable.end() )
      MapleRaiseError(kv,err);
    else
      flag = f_i->second;

    h_i = hashTable.find(key);
    if( h_i != hashTable.end() ) {

      switch( flag ) {
       case BlackBoxi:{
	 long result;
	 MapleBBi* BB = (MapleBBi*) h_i->second;
	 LinBox::Modular<long> field;
	 BB->getField(field);
	 return ToMapleInteger(kv, LinBox::det(result, *BB, field));
       }
       break;
	 
        case BlackBoxI: {
         integer result2;
	 ALGEB blank;
	 MapleBBI* BB = (MapleBBI*) h_i->second;
	 LinBox::Modular<integer> field;
	 BB->getField(field);
	 return LiToM(kv,  LinBox::det(result2, *BB, field), blank);
	}
	break;
      }
    }
    else
      MapleRaiseError(kv, err);
  }
}

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

    f_i = typeTable.find(key);
    if( f_i == typeTable.end() )
      MapleRaiseError(kv,err);
    else
      flag = f_i->second;


    h_i = hashTable.find(key);
    if(h_i != hashTable.end() ) {
      switch( flag ) {
         case BlackBoxi: {
	   Vectorl mpreturn;
	    Vectorl::iterator mp_i;
	    MapleBBi* BB = (MapleBBi*) h_i->second;
	    LinBox::Modular<long> field;
	    BB->getField(field);
	    LinBox::minpoly( mpreturn, *BB, field );
	    retlist = MapleListAlloc(kv, mpreturn.size() );
	    for(i = 1, mp_i = mpreturn.begin(); mp_i != mpreturn.end(); ++mp_i, ++i)
	      MapleListAssign(kv, retlist, i, ToMapleInteger(kv, *mp_i));
	 }
	 break;

         case BlackBoxI: {
	   VectorI mpreturn;
	   VectorI::iterator mp_i;
	   MapleBBI* BB = (MapleBBI*) h_i->second;
	   LinBox::Modular<integer> field;
	   BB->getField(field);
	   LinBox::minpoly( mpreturn, *BB, field );
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

ALGEB & LiToM(MKernelVector & kv, const integer & In, ALGEB & Out) 
{
  // If we get lucky, this thing fits in one word, which is good, so we just
  // straight convert it
  if(In.size() == 1)
    Out = ToMapleInteger(kv, Integer2long(In) );
  else {
    Out = MapleListAlloc(kv, In.size() );
    for(int i = 1; i <= In.size(); ++i)
      MapleListAssign(kv,Out,i,ToMapleInteger(kv,In[i-1]));
  }

  return Out;
}

integer & MtoLI(MKernelVector & kv, integer & Out, const ALGEB &In)
{
  static integer base(10000);
  integer b = integer(1); Out = integer(0);
  if( IsMapleInteger32(kv,In) )
    Out = integer(MapleToInteger32(kv, In));
  else
    for(int i = 2; i <= MapleToInteger32(kv,MapleListSelect(kv,In,1)); ++i) {
      Out += b * MapleToInteger32(kv,MapleListSelect(kv,In,i));
      b *= base;
    }

  return Out;
}
