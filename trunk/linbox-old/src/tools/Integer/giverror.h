#ifndef _GIVERROR_H_
#define _GIVERROR_H_
// ==========================================================================
// $Source$
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Authors: T. Gautier
// $Id$
// ==========================================================================
// Description:
// - error exception 

#include <iostream.h>

// ------------------------------- GivError
// - Base class for execption handling in Givaro
class GivError {
public:
  GivError(const char* msg =0 ) 
  : strg(msg) {};

  // -- virtual print of the error message
  virtual ostream& print( ostream& o )  const;
  
  // -- non virtual output operator
  friend ostream& operator<< (ostream& o, const GivError& E);

  // - useful to setup a break point on it
  static void throw_error( const GivError& err );

protected:
  const char* strg;  
};

class GivMathError : public GivError {
public:
  GivMathError(const char* msg ) : GivError(msg) {}
};

/*
class GivMathDivZero : public GivMathError {
public:
  GivMathDivZero(const char* msg ) : GivMathError(msg) {}
};
*/


// -- Exception thrown in input of data structure 
class GivBadFormat : public GivError {
public:
  GivBadFormat(const char* msg ) : GivError(msg) {}
};

#include <LinBox/giverror.C>

#endif
