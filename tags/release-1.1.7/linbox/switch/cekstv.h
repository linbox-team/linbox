/* linbox/switch/cekstv.h
 * Copyright (C) 1999-2001 William J Turner
 *
 * Written by William J Turner <wjturner@math.ncsu.edu>
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
 *   - Renamed file as cekstv.h
 *   - Renamed class cekstv_switch as CekstvSwitch
 *   - Reindent
 * -----------------------------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __LINBOX_cekstv_H
#define __LINBOX_cekstv_H

#include <vector>

namespace LinBox
{

template <class Field>
class CekstvSwitchFactory;

/** Butterfly switch object from preconditioner paper.
 * This is a switch predicate object that is applied
 * to two references to elements to switch them as needed 
 * by the \Ref{Butterfly Switching Network BlackBox Matrix Object}
 * following the exchange matrix introduced in "Efficient Matrix
 * Preconditioners for Black Box Linear Algebra" by Chen, Eberly, 
 * Kaltofen, Saunders, Turner, and Villard.
 * This class is templatized by the field in which the arithmetic
 * is done.
 */
template <class Field>
class CekstvSwitch
{
    public:

	/// Typedef
	typedef typename Field::Element Element;
	typedef CekstvSwitch<Field> Self_t;
	typedef CekstvSwitchFactory<Field> Factory;

	CekstvSwitch () {}

	/** Constructor from a field and a field element.
	 * @param F field in which arithmetic is done
	 * @param switches vector of switches
	 */
	CekstvSwitch (const typename Field::Element &a)
		: _a (a) 
	{}

	/** Destructor.
	 */
	~CekstvSwitch () {}

	/** Apply switch function.
	 * Switches the elements in references according to the
	 * exchange matrix introduced in "Efficient Matrix
	 * Preconditioners for Black Box Linear Algebra" by Chen, Eberly, 
	 * Kaltofen, Saunders, Turner, and Villard and the current field element
	 * specified in the switch object.
	 * @return bool true if swapped, false otherwise
	 * @param x reference to first element to be switched
	 * @param y reference to second element to be switched
	 */
	bool apply (const Field &F, Element &x, Element &y) const;

	/** Apply switch transpose function.
	 * Switches the elements in references according to the
	 * transpose of the exchange matrix introduced in "Efficient Matrix
	 * Preconditioners for Black Box Linear Algebra" by Chen, Eberly, 
	 * Kaltofen, Saunders, Turner, and Villard and the current field element
	 * specified in the switch object.
	 * @return bool true if swapped, false otherwise
	 * @param x reference to first element to be switched
	 * @param y reference to second element to be switched
	 */
	bool applyTranspose (const Field &F, Element &x, Element &y) const;


        template<typename _Tp1>
        struct rebind
        { 
            typedef CekstvSwitch<_Tp1> other;

                // special rebind operator() with two fields, 
                // indeed local field is not stored in the switch
            void operator() (other & Ap, const Self_t& A, const _Tp1& T, const Field& F) {
                Hom<Field, _Tp1>(F,T).image(Ap.getData(), A.getData());
            }
        };
    
    typename Field::Element& getData() { return _a; }
    const typename Field::Element& getData() const { return _a; }


   private:

	// Parameter of this 2x2 block
	typename Field::Element _a;
};

/** Cekstv switch factory
 *
 * This class facilitates construction of cekstv switch objects by the
 * butterfly matrix.
 */

template <class Field>
class CekstvSwitchFactory 
{
    public:
	/** Constructor from an STL vector of bools
	 */
	CekstvSwitchFactory (typename Field::RandIter r)
		: _r (r)
	{}

	/** Construct and return a boolean switch object
	 */
	CekstvSwitch<Field> makeSwitch ()
		{ typename Field::Element a; return CekstvSwitch<Field> (_r.random (a)); }

    private:

	typename Field::RandIter _r;
};

template <class Field> 
inline bool CekstvSwitch<Field>::apply (const Field             &F,
					typename Field::Element &x,
					typename Field::Element &y) const
{
	F.axpyin (x, _a, y);
	F.addin (y, x);

	return true;
}
  
template <class Field> 
inline bool CekstvSwitch<Field>::applyTranspose (const Field             &F,
						 typename Field::Element &x,
						 typename Field::Element &y) const
{
	F.addin (x, y);
	F.axpyin (y, _a, x);

	return true;
}

}

#endif // __LINBOX_cekstv_H
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
