/* File: src/examples/field/unparameteric/fuzzy_envelope.cpp
 * Author: William Turner
 */

#include "test_fuzzy_envelope.h"

/** Test of LinBox system.
 * Creates input and output streams and calls \Ref{test_linbox}
 * to do actual test of LinBox.
 * @param argc number of arguments
 * @param argv array of input arguments:
 *             argv[0]  program name, 
 *             argv[1]  file from which to read input (default = cin), 
 *             argv[2]  file to which to print output (default = cout), 
 *             argv[3]  file to which to log messages (default = clog)
 */
int main(int argc, char* argv[])
{
  test_linbox T(argc, argv);
  T.test<test_linbox::field_categories::fuzzy_envelope_tag>();

} // int main(int argc, char* argv[])
