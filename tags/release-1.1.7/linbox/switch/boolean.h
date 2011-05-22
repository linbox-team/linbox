/* linbox/switch/boolean.h
 * Copyright (C) 1999-2001 William J Turner
 *
 * Written by William J Turner <wjturner@math.ncsu.edu>,
 *
 * -----------------------------------------------------------
 * 2002-09-26  Bradford Hovinen  <bghovinen@math.uwaterloo.ca>
 *
 * Refactoring: The switch object now contains only the information for one 2x2
 * block. A vector of switches is maintained by the butterfly preconditioner
 * instead of the switch object. Since there will be many switch objects, they
 * should be kept very lightweight, so the field is not maintained in the object
 * itself, but instead passed to apply and applyTranspose. Since those methods
 * are inline, this does not create overhead. apply and applyTranspose now take
 * four field elements: the source elements and destination elements. This
 * eliminates the need to keep an additional temporary in the class, and
 * eliminates the need for copying in the butterfly.
 *
 * -----------------------------------------------------------
 * 2002-08-20  Bradford Hovinen  <hovinen@cis.udel.edu>
 *
 * Brought this file into the current Linbox framework:
 *   - Renamed file as boolean.h
 *   - Renamed class boolean_switch as BooleanSwitch
 *   - Reindent
 * -----------------------------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __LINBOX_boolean_H
#define __LINBOX_boolean_H

#include <vector>

namespace LinBox
{

class BooleanSwitchFactory;

/** Boolean switch object.
 * This is a switch predicate object that is applied
 * to two references to elements to switch them as needed 
 * by the \Ref{Butterfly Switching Network BlackBox Matrix Object}.
 */
class BooleanSwitch
{
    public:

	typedef BooleanSwitch Self_t;
	typedef BooleanSwitchFactory Factory;

	/** Constructor from an STL vector of booleans.
	 * The switch is applied using the vector of booleans.
	 * A true value means to swap the two elements, and a false
	 * value means not to.
	 * The apply function starts at the beginning of the vector moving 
	 * forward through it, and applyTranspose function starts at the end
	 * moving backwards.  Both repeat the vector after they pass through it.
	 * @param switches vector of switches
	 */
	BooleanSwitch (const bool s)
		: _s (s)
	{}

	/** Destructor.
	 */
	~BooleanSwitch () {}

	/** Apply switch function.
	 * Switches the elements in references according to current boolean
	 * value.  Swaps the elements if boolean is true, otherwise does nothing.
	 * It is templatized by the element type to be swapped.
	 * @return bool true if swapped, false otherwise
	 * @param x reference to first element to be switched
	 * @param y reference to second element to be switched
	 */
	template <class Field>
	bool apply (const Field             &F,
		    typename Field::Element &x,
		    typename Field::Element &y) const;

	/** Apply switch transpose function.
	 * Switches the elements in references according to current boolean
	 * value.  Swaps the elements if boolean is true, otherwise does nothing.
	 * It is templatized by the element type to be swapped.
	 * @return bool true if swapped, false otherwise
	 * @param x reference to first element to be switched
	 * @param y reference to second element to be switched
	 */
	template <class Field>
	bool applyTranspose (const Field             &F,
			     typename Field::Element &x,
			     typename Field::Element &y) const;

        template<typename _Tp1>
        struct rebind
        { 
            typedef BooleanSwitch other;

        };
    
    protected:

	bool _s;

}; // class boolean_switch

/** Boolean switch factory
 *
 * This class facilitates construction of boolean switch objects by the
 * butterfly matrix.
 */

class BooleanSwitchFactory 
{
    public:
	/** Constructor from an STL vector of bools
	 */
	BooleanSwitchFactory (const std::vector<bool> &switches)
		: _switches (switches), _iter (switches.begin ())
	{}

	/** Construct and return a boolean switch object
	 *
	 * This function walks through the switches object given in the
	 * constructor, advancing on each invocation. It wraps around to the
	 * beginning of the vector when it reaches the end.
	 */
	BooleanSwitch makeSwitch ()
	{
		if (_iter == _switches.end ())
			_iter = _switches.begin ();

		return BooleanSwitch (*_iter++);
	}

    private:

	const std::vector<bool> &_switches;
	std::vector<bool>::const_iterator _iter;
};

template <class Field> 
inline bool BooleanSwitch::apply (const Field             &F,
				  typename Field::Element &x,
				  typename Field::Element &y) const
{
	if (_s)
		std::swap (x, y);

	return _s;
}
  
template <class Field> 
inline bool BooleanSwitch::applyTranspose (const Field             &F,
					   typename Field::Element &x,
					   typename Field::Element &y) const
{
	if (_s)
		std::swap (x, y);

	return _s;
}

}

#endif // __LINBOX_boolean_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
