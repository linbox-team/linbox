/* -*- mode: C++; style: linux -*- */

/* linbox/field/field-interface.h
 * Copyright (C) 2002 David Saunders
 *
 * For licensing information see COPYING
 */

#ifndef __FIELD_INTERFACE_H
#define __FIELD_INTERFACE_H

namespace LinBox
{
/** LinBox Field Interface
 * The LinBox {@link Fields field} common object {@link Interfaces interface}.
 * The field interface includes the following public members:
 *
 * Types: <tt>Element</tt> and <tt>RandIter</tt>.
 *
 * Object management member functions:
 *   null constructor, copy constructor, destructor, assignment operator, 
 *   <tt>convert(), init(), assign(), characteristic(), cardinality()</tt>.
 *
 * Predicates on field elements:
 *   <tt>areEqual(), isZero(), isOne()</tt>.
 *
 * Basic arithmetic functions:
 *   <tt>axpy(), add(), neg(), sub(), mul(), inv(), div()</tt>.
 *
 * Inplace arithmetic functions:
 *   <tt>axpyin(), addin(), negin(), subin(), mulin(), invin(), divin()</tt>.
 *
 * I/O functions:
 *   <tt>read()</tt> and <tt>write()</tt> for I/O of the field itself and for I/O of its elements.
 *
 * The field archetype class is is the reference instantiation of this 
 * interface and contains the generic specifications of the member functions.
 * Documentation of other field classes may only explain special properties 
 * of the member functions specific to the class and explain the constructors 
 * and any other functionality unique to the class.
 *
 *  @see Interfaces
*/
class FieldInterface { };// empty class so doc++ makes a nice hierarchy.

} // namespace LinBox

#endif // __FIELD_INTERFACE_H
