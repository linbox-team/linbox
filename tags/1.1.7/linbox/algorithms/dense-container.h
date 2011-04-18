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

#ifndef __LINBOX_blackbox_container_H
#define __LINBOX_blackbox_container_H

#include <linbox/randiter/archetype.h>
#include <linbox/algorithms/blackbox-container-base.h>
#include <linbox/util/timer.h>

namespace LinBox 
{

/// \brief Limited doc so far.
template<class Field, class Vector, class RandIter = typename Field::RandIter>
class DenseContainer : public BlackboxContainerBase<Field, Vector> {
    public:
	typedef typename BlackboxContainerBase<Field, Vector>::Blackbox Blackbox;

	DenseContainer () {} 

	DenseContainer(const Blackbox * D, typename Field::Element * U, size_t ldu, 
		       const Field &F, const Vector &u0) 
		: BlackboxContainerBase<Field, Vector> (D, F)
	{
		init (u0, u0); w = u;
		Up = U;
		_ldu = ldu;
#ifdef INCLUDE_TIMING
		_applyTime = _dotTime = 0.0;
#endif
	}
    
	DenseContainer(const Blackbox * D, typename Field::Element * U, size_t ldu,
		       const Field &F, const Vector &u0, const Vector& v0) 
		: BlackboxContainerBase<Field, Vector> (D, F)
	{
                // JGD 22.03.03
// 		init (u0, v0); w = u;
		init (u0, v0); w = v;
		Up = U;
		_ldu = ldu;
#ifdef INCLUDE_TIMING
		_applyTime = _dotTime = 0.0;
#endif
	}
    
	DenseContainer(const Blackbox * D, typename Field::Element * U, size_t ldu, 
		       const Field &F, RandIter &g) 
		: BlackboxContainerBase<Field, Vector> (D, F)
	{
		init (g); w = u;
		Up = U;
		_ldu = ldu;
#ifdef INCLUDE_TIMING
		_applyTime = _dotTime = 0.0;
#endif
	}

#ifdef INCLUDE_TIMING
	double applyTime () const { return _applyTime; }
	double dotTime   () const { return _dotTime; }
#endif // INCLUDE_TIMING

    protected:
	Vector w;
	typename Field::Element * Up;

	size_t _ldu;

#ifdef INCLUDE_TIMING
	Timer _timer;
	double _applyTime, _dotTime;
#endif // INCLUDE_TIMING

	void _launch () {
		typename Vector::iterator it;
		Integer tmp;
		size_t i;
		if (casenumber) {
#ifdef INCLUDE_TIMING
			_timer.start ();
#endif // INCLUDE_TIMING
			_BB->apply (v, w);  // GV

#ifdef INCLUDE_TIMING
			_timer.stop ();
			_applyTime += _timer.realtime ();
			_timer.start ();
#endif // INCLUDE_TIMING

			// Copy of v into a row of U
			it = v.begin();
			i=0;
			for (; it!=v.end(); it++, i++){
				*(Up+i) = *it;
				_F.convert(tmp,*it);
				cerr<<" copie of "<<tmp;
			}
			cerr<<endl;
			Up += _ldu;
			_VD.dot (_value, u, v);  // GV 

#ifdef INCLUDE_TIMING
			_timer.stop ();
			_dotTime += _timer.realtime ();
#endif // INCLUDE_TIMING

			casenumber = 0;
		} else {
#ifdef INCLUDE_TIMING
			_timer.start ();
#endif // INCLUDE_TIMING
			_BB->apply (w, v);  // GV

			// Copy of v into a row of U
			it = w.begin();
  			i=0;
  			for (; it!=w.end(); it++, i++){
 				_F.convert(tmp,*it);
  				cerr<<" copie of "<<tmp;
  				*(Up+i) = *it;
  			}
  			cerr<<endl;
  			Up += _ldu;
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

			casenumber = 1;
		}  
	}

	void _wait () {}
};
 
}

#endif // __LINBOX_blackbox_container_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
