/* -*- mode: C++; style: linux -*- */

/* linbox/algorithms/blackbox-container.h
 * Copyright (C) 1999, 2001 Jean-Guillaume Dumas, Bradford Hovinen
 *
 * Written by Jean-Guillaume Dumas <Jean-Guillaume.Dumas@imag.fr>,
 *            Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __BLACKBOX_CONTAINER_H
#define __BLACKBOX_CONTAINER_H

#include <linbox/randiter/archetype.h>
#include <linbox/algorithms/blackbox-container-base.h>
#include <linbox/util/timer.h>

// Define to enable timing facilities
#undef INCLUDE_TIMING

namespace LinBox 
{

template<class Field, class Vector, class RandIter = typename Field::RandIter>
class BlackboxContainer : public BlackboxContainerBase<Field, Vector> {
    public:
	BlackboxContainer () {} 

	BlackboxContainer(const Blackbox * D, const Field &F, const Vector &u0) 
		: BlackboxContainerBase<Field, Vector> (D, F)
		{ init (u0, u0); w = u; }
    
	BlackboxContainer(const Blackbox * D, const Field &F, const Vector &u0, const Vector& v0) 
		: BlackboxContainerBase<Field, Vector> (D, F)
		{ init (u0, v0); w = u;}
    
	BlackboxContainer(const Blackbox * D, const Field &F, RandIter &g) 
		: BlackboxContainerBase<Field, Vector> (D, F)
		{ init (g); w = u; }

#ifdef INCLUDE_TIMING
	double applyTime () const { return _applyTime; }
	double dotTime   () const { return _dotTime; }
#endif // INCLUDE_TIMING

    protected:
	Vector w;

#ifdef INCLUDE_TIMING
	Timer _timer;
	double _applyTime = 0.0, _dotTime = 0.0;
#endif // INCLUDE_TIMING

	void _launch () {
		if (even) {
#ifdef INCLUDE_TIMING
			_timer.start ();
#endif // INCLUDE_TIMING
			_BB->apply (v, w);  // GV

#ifdef INCLUDE_TIMING
			_timer.stop ();
			_applyTime += _timer.realtime ();
			_timer.start ();
#endif // INCLUDE_TIMING

			_VD.dot (_value, u, v);  // GV 

#ifdef INCLUDE_TIMING
			_timer.stop ();
			_dotTime += _timer.realtime ();
#endif // INCLUDE_TIMING

			even = 0;
		} else {
#ifdef INCLUDE_TIMING
			_timer.start ();
#endif // INCLUDE_TIMING
			_BB->apply (w, v);  // GV

#ifdef INCLUDE_TIMING
			_timer.stop ();
			_applyTime += _timer.realtime ();
			_timer.start ();
#endif // INCLUDE_TIMING

			_VD.dot (_value, u, w);  // GV

#ifdef INCLUDE_TIMING
			_timer.stop ();
			_dotTime += _timer.realtime ();
#endif // INCLUDE_TIMING

			even = 1;
		}  
	}

	void _wait () {}
};
 
}

#endif // __BLACKBOX_CONTAINER_H

