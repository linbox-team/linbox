/* File: include/LinBox/integer.h
 * Author: William J Turner for the LinBox group
 *
 * This file should not (probably) be in the final release, but
 * is only intended to be used until the integer type is *finally* nailled down.
 */

// Namespace in which all LinBox library code resides
namespace LinBox
{
  /** Integer type.
   * This type should be a widening of long in that it implements all 
   * operators the C++ basic type long does.  It should also be able to be 
   * converted to and from a long.
   */
  typedef long integer;

} // namespace LinBox

