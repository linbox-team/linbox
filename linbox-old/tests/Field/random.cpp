/* File: src/examples/random.cpp
 * Author: William J Turner for the LinBox group
 */

#include "test_base.h"

#include "test_sparsemat.h"
#include "vector_utility.h"

class test_random : test_base
{
public:
  
  /** Constructor from int and array of C-style strings.
   * Creates input and output streams and calls \Ref{test_linbox}
   * to do actual test of LinBox.
   * @param argc number of arguments
   * @param argv array of input arguments:
   *         argv[0]  program name, 
   *         argv[1]  file from which to read input (default = cin), 
   *         argv[2]  file to which to print output (default = cout), 
   *         argv[3]  file to which to log messages (default = clog)
   */
  test_random(int argc, char* argv[]) : test_base(argc, argv) {}

  bool test(void) const;

}; // class test_random : test_random

bool test_random::test(void) const
{
  random_vector rand_vec(0, *in_ptr, *out_ptr, *log_ptr);
  rand_vec.create();

  random_sparsemat rand_sparse(0, 0, *in_ptr, *out_ptr, *log_ptr);
  rand_sparse.create();

  return true;
} // bool test_random::test(void) const
 
int main(int argc, char* argv[])
{
  test_random T(argc, argv);
  T.test();
} // int main(int argc, char* argv[])
