/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

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

namespace LinBox 
{

/// \brief Limited doc so far.
template<class Field, class _Blackbox, class RandIter = typename Field::RandIter>
class BlackboxContainer : public BlackboxContainerBase<Field, _Blackbox> {
    public:
	typedef _Blackbox Blackbox;

	BlackboxContainer () {} 

	template<class Vector>

	BlackboxContainer(const Blackbox * D, const Field &F, const Vector &u0) 
		: BlackboxContainerBase<Field, Blackbox> (D, F)
	{
		init (u0, u0); w = this->u;
#ifdef INCLUDE_TIMING
		_applyTime = _dotTime = 0.0;
#endif
	}

	// Pascal Giorgi 16.02.2004
	template<class Vector>
	BlackboxContainer(const Blackbox * D, const Field &F, const Vector &u0, unsigned long size) 
		: BlackboxContainerBase<Field, Blackbox> (D, F,size)
	{
		init (u0, u0); w = this->u;
#ifdef INCLUDE_TIMING
		_applyTime = _dotTime = 0.0;
#endif
	}

    
	template<class Vector1, class Vector2>
	BlackboxContainer(const Blackbox * D, const Field &F, const Vector1 &u0, const Vector2& v0) 
		: BlackboxContainerBase<Field, Blackbox> (D, F)
	{
		init (u0, v0); w = this->v;
#ifdef INCLUDE_TIMING
		_applyTime = _dotTime = 0.0;
#endif
	}
    
	BlackboxContainer(const Blackbox * D, const Field &F, RandIter &g) 
		: BlackboxContainerBase<Field, Blackbox> (D, F)
	{
		init (g); w = this->u;
#ifdef INCLUDE_TIMING
		_applyTime = _dotTime = 0.0;
#endif
	}

#ifdef INCLUDE_TIMING
	double applyTime () const { return _applyTime; }
	double dotTime   () const { return _dotTime; }
#endif // INCLUDE_TIMING

    protected:
	std::vector<typename Field::Element> w;

#ifdef INCLUDE_TIMING
	Timer _timer;
	double _applyTime, _dotTime;
#endif // INCLUDE_TIMING

	void _launch () {
		if (this->casenumber) {
#ifdef INCLUDE_TIMING
			_timer.start ();
#endif // INCLUDE_TIMING
			this->_BB->apply (this->v, w);  // GV

#ifdef INCLUDE_TIMING
			_timer.stop ();
			_applyTime += _timer.realtime ();
			_timer.start ();
#endif // INCLUDE_TIMING

			this->_VD.dot (this->_value, this->u, this->v);  // GV 

#ifdef INCLUDE_TIMING
			_timer.stop ();
			_dotTime += _timer.realtime ();
#endif // INCLUDE_TIMING

			this->casenumber = 0;
		} else {
#ifdef INCLUDE_TIMING
			_timer.start ();
#endif // INCLUDE_TIMING
			this->_BB->apply (w, this->v);  // GV

#ifdef INCLUDE_TIMING
			_timer.stop ();
			_applyTime += _timer.realtime ();
			_timer.start ();
#endif // INCLUDE_TIMING

			this->_VD.dot (this->_value, this->u, w);  // GV

#ifdef INCLUDE_TIMING
			_timer.stop ();
			_dotTime += _timer.realtime ();
#endif // INCLUDE_TIMING

			this->casenumber = 1;
		}  
	}

	void _wait () {}
};
 
}

#endif // __BLACKBOX_CONTAINER_H

