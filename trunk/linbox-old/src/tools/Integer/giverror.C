// ==========================================================================
// $Source$
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Authors: T. Gautier
// $Id$
// ==========================================================================
// Description:
// - error exception 

#include "giverror.h"
#include <iostream.h>


ostream& GivError::print( ostream& o ) const
{ return o << strg ; }


void GivError::throw_error( const GivError& err )
{
  throw err;
}

ostream& operator<< (ostream& o, const GivError& E) 
{
   E.print(o) ; 
   return o ;
}

