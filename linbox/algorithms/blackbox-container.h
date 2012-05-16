/* linbox/algorithms/blackbox-container.h
 * Copyright (C) 1999, 2001 Jean-Guillaume Dumas, Bradford Hovinen
 *
 * Written by Jean-Guillaume Dumas <Jean-Guillaume.Dumas@imag.fr>,
 *            Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * ------------------------------------
 *
 * 
 * ========LICENCE========
 * This file is part of the library LinBox.
 * 
 * LinBox is free software: you can redistribute it and/or modify
 * it under the terms of the  GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * 
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 *.
 */

#ifndef __LINBOX_blackbox_container_H
#define __LINBOX_blackbox_container_H

#include "linbox/randiter/archetype.h"
#include "linbox/algorithms/blackbox-container-base.h"
#include "linbox/util/timer.h"

namespace LinBox
{

	/// \brief Limited doc so far.
	template<class Field, class _Blackbox, class RandIter = typename Field::RandIter>
	class BlackboxContainer : public BlackboxContainerBase<Field, _Blackbox> {
	public:
		typedef _Blackbox Blackbox;

		BlackboxContainer () {}

		template<class Vector>

		BlackboxContainer(const Blackbox * D, const Field &F, const Vector &u0) :
			BlackboxContainerBase<Field, Blackbox> (D, F)
		{
			init (u0, u0); w = this->u;
#ifdef INCLUDE_TIMING
			_applyTime = _dotTime = 0.0;
#endif
		}

		// Pascal Giorgi 16.02.2004
		template<class Vector>
		BlackboxContainer(const Blackbox * D, const Field &F, const Vector &u0, unsigned long size) :
			BlackboxContainerBase<Field, Blackbox> (D, F,size)
		{
			init (u0, u0); w = this->u;
#ifdef INCLUDE_TIMING
			_applyTime = _dotTime = 0.0;
#endif
		}


		template<class Vector1, class Vector2>
		BlackboxContainer(const Blackbox * D, const Field &F, const Vector1 &u0, const Vector2& v0) :
			BlackboxContainerBase<Field, Blackbox> (D, F)
		{
			init (u0, v0); w = this->v;
#ifdef INCLUDE_TIMING
			_applyTime = _dotTime = 0.0;
#endif
		}

		BlackboxContainer(const Blackbox * D, const Field &F, RandIter &g) :
			BlackboxContainerBase<Field, Blackbox> (D, F)
		{
			this->casenumber = 1;
			this->u.resize (this->_BB->coldim ());
			for (long i = this->u.size (); i--;)
				g.random (this->u[i]);
			this->w.resize (this->_BB->coldim ());
			for (long i = this->w.size (); i--;)
				g.random (this->w[i]);
			this->v.resize (this->_BB->rowdim ());
			this->_VD.dot (this->_value, this->u, this->w);
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
			}
			else {
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

#endif // __LINBOX_blackbox_container_H


// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

