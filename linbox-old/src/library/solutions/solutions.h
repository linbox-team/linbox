#ifndef __SOLUTIONS_H__
#define __SOLUTIONS_H__
#include <LinBox/integer.h>
/** Solutions

  The {/em solutions} are C++ functions implementing our ``standard" method 
  of solution for the various goal functions of LinBox.  The goal 
  functions are rank, determinant, linsolve(+variants), minpoly, charpoly, 
  frobenius form, etc.

  Often there are multiple methods which may be used for a specific goal function
  and/or settable parameters (such as a block size) which guide the computation.  
  The optional last argumment to each function is a MethodTrait.  It is used to 
  allow selection of method and/or selection of specific values for parameters to 
  the method.  
  */

/*
  Shorthand type names:
  
  typedef BlackBoxMatrix BB;
  typedef vector-or-dense-matrix-class Vector;
  typedef vector<element> Polynomial;  // degree is size.
  typedef MethodTrait MT;
*/
namespace LinBox
{
  typedef size_t index;
  typedef double Pbound;  // an upper bound for the probability of erroneous result.

  //  The basic function signatures are:

  /** r <- rank( A, p) 
   * Returns $r \le \rank(A)$, 
   * Pbound $p$ is set so that the $\P( r \lt \rank(A) ) \le p$.
  template <class BB, class MT>
  index rank(const BB& A, Pbound p, MT mt = 0)
  ;

  /// Return value = d = det(A).
  template <class Field, class BB, class MT>
  Field::element& determinant(Field::element& d, const BB& A, MT mt = 0)
  ;

  /** Returns true if $A$ is singular, false if $A$ is nonsingular.
   *  if true is returned, $A$ is singular.
   *  if false is returned. $\P$( $A$ is singular )~$\le p$.
   *  This isSingular() function is more efficient than isNonSingular().
   */
  template <class BB, class MT>
  bool isSingular(const BB& A, Pbound p, MT mt = 0)
  ;

  /// Returns false if A is singular, true if .A is nonsingular.
  template <class BB, class MT>
  bool isNonSingular(const BB& A, MT mt = 0)
  ;

  /** v <-- the solution to Ax = b.  A must be nonsingular.
   */
  template <class Vector, class BB, class MT>
  Vector& linsolve_nonsingular(Vector& v, const BB& A, const Vector& b)
  ;

  /** linsolve0.
   *  v <-- a nonzero solution to Ax = 0 if A is singular.
   *  v <-- 0 if A is nonsingular.
   */
  template <class Vector, class BB, class MT>
  Vector& linsolve(Vector& v, const BB& A, MT mt = 0)
  ;

  /** v <-- a random solution to Ax = 0.
   *  If the nullspace of $A$ has dimension $n$, no vector has probability 
   *  greater than $1/{r^n}$ of being returned.
   */
  template <class Vector, class BB, class MT>
  Vector& linsolve_random(Vector& v, const BB& A, integer r, MT mt = 0)
  ;

  /** v <-- a solution to Ax = b 
   *  If the Pbound $p$ returned is positive then $v$ is not a solution vector
   *  and the system is consistent with probability less than $p$.
   *  Separate calls are independent trials.
   */
  template <class Vector, class BB, class MT>
  Vector& linsolve(Vector& v, Pbound& p, const BB& A, const Vector& b, MT mt = 0)
  ;

  /** v <-- The Moore-Penrose least norm least squares solution to Ax = b 
   *  That is, $v$ satisfies 
   *  (1) $|Av - b|_2$ is minimal, and 
   *  (2) Among all such minimal solutions, $|v|_2$ is minimum.
   *  Equivalently, $v = A^\dagger b$, where $A^\dagger$ is the pseudoinverse of $A$.
   *  In particular, if $Ax = b$ is consistent, $v$ is the least norm solution.
   */
  template <class Vector, class BB, class MT>
  Vector& MPlnls(Vector& v, const BB& A, const Vector& b, MT mt = 0)

  /** v <-- a random solution to Ax = b.
   *  If the nullspace of $A$ has dimension n, no vector has probability 
   *  greater than $1/{r^n}$ of being returned.
   *  If the Pbound $p$ returned is positive then $v$ is not a solution vector
   *  and the system is consistent with probability less than $p$.
   *  Separate calls are independent trials.
   */
  template <class Vector, class BB, class MT>
  Vector& linsolve_random(Vector& v, Pbound& p, const BB& A, const Vector& b, integer r, MT mt = 0)
  ;

  /** p <-- minpoly(A), a factor of the minpoly of A.
   *  $\P ( p \neq \minpoly(A) ) \le p$.
   *  Separate calls are independent trials.
   */
  template <class Polynomial, class BB, class MT>
  Polynomial& minpoly(Polynomial& p, Pbound& p, const BB& A, MT mt = 0)
  ;

  /** return value = p = charpoly(A), the characteristic polynomial of A.
   */
  template <class Polynomial, class BB, class MT>
  Polynomial& charpoly(Polynomial& p, const BB& A, MT mt = 0)
  ;
} //namespace LinBox
#endif
