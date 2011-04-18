/* linbox/field/field-interface.h
 * Copyright (C) 2002 David Saunders
 *
 * For licensing information see COPYING
 */

#ifndef __LINBOX_field_interface_H
#define __LINBOX_field_interface_H

namespace LinBox
{
// LinBox Field Interface
///
/*
 * The LinBox {@link Fields field} common object {@link Interfaces interface}.
 * The field interface includes the following public members:
 *
 * Types: {\tt Element} and {\tt RandIter}.
 *
 * Object management member functions:
 *   null constructor, copy constructor, destructor, assignment operator, 
 *   {\tt convert(), init(), assign(), characteristic(), cardinality()}.
 *
 * Predicates on field elements:
 *   {\tt areEqual(), isZero(), isOne()}.
 *
 * Basic arithmetic functions:
 *   {\tt axpy(), add(), neg(), sub(), mul(), inv(), div()}.
 *
 * Inplace arithmetic functions:
 *   {\tt axpyin(), addin(), negin(), subin(), mulin(), invin(), divin()}.
 *
 * I/O functions:
 *   {\tt read()} and {\tt write()} for I/O of the field itself and for I/O of its elements.
 *
 * The field archetype class is is the reference instantiation of this 
 * interface and contains the generic specifications of the member functions.
 * Documentation in other field classes is more limited. It serves primarily to explain special properties 
 * specific to the class of the interface member functions and to explain any constructors 
 * or other functionality unique to the class.
 *
 *  @see Interfaces
*/
/** 
 * \brief This field base class exists solely to aid documentation organization.


 *  For the general field member function documentation consult the {@link FieldArchetype
 FieldArchetype}. For specific properties of individual representations consult the specific field classes.
 \ingourp field
 */
	class FieldInterface 
{
/*
    public:
	// this just demo's that some declarations could be here.
	typedef ElementArchetype Element; 
	virtual Element& mul(Element& c, const Element& a, const Element& b) const = 0;
*/
};// empty class so doc++ makes a nice hierarchy.

} // namespace LinBox

#endif // __LINBOX_field_interface_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
