/* -*- mode: c; style: linux -*- */

/* linbox/src/algorithms/massey-domain.h
 * Copyright (C) 1999, 2001 Jean-Guillaume Dumas, Bradford Hovinen
 *
 * Written by Jean-Guillaume Dumas <Jean-Guillaume.Dumas@imag.fr>,
 *            Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

// ======================================================================= //
// Linbox project 1999
// Domain Massey
// - Computation is stopped when the polynomials remain the same
//   for more than EARLY_TERM_THRESOLD
// - When minimal polynomial equals characteristic polynomial,
//   2 additional iterations are needed to compute it 
//   (parameter DEFAULT_ADDITIONAL_ITERATION), but those
//   iterations are not needed for the rank
// Time-stamp: <27 Aug 01 18:18:12 Jean-Guillaume.Dumas@imag.fr> 
// ======================================================================= //

#ifndef __MASSEY_DOMAIN_H
#define __MASSEY_DOMAIN_H

#include <linbox/util/commentator.h>

namespace LinBox 
{

#ifndef MIN
#define MIN(a,b) ((a)<(b)?(a):(b))
#endif

#define DEFAULT_EARLY_TERM_THRESHOLD 20
#define DEFAULT_ADDITIONAL_ITERATION 2

template<class Field, class Sequence>
class MasseyDomain {
private:
	Sequence      *_container;
	Field          _field;
	Commentator    _Comm;
	unsigned long  EARLY_TERM_THRESHOLD;

public:
	typedef typename Field::Element Element;

        //-- Constructors
	MasseyDomain (unsigned long ett_default = DEFAULT_EARLY_TERM_THRESHOLD) 
		: _container           (), 
		  _field               (), 
		  _Comm                (PRINT_NOTHING, PRINT_NOTHING),
		  EARLY_TERM_THRESHOLD (ett_default)
		{}

	MasseyDomain (const MasseyDomain<Field, Sequence> &M, unsigned long ett_default = DEFAULT_EARLY_TERM_THRESHOLD) 
		: _container           (M._container), 
		  _field               (M._field), 
		  _Comm                (M._Comm),
		  EARLY_TERM_THRESHOLD (ett_default)
		{}

	MasseyDomain (Sequence *D, unsigned long ett_default = DEFAULT_EARLY_TERM_THRESHOLD) 
		: _container           (D), 
		  _field               (D->getField ()),
		  _Comm                (PRINT_NOTHING, PRINT_NOTHING),
		  EARLY_TERM_THRESHOLD (ett_default)
		{}
  
	MasseyDomain (Sequence *MD, const Field &F, unsigned long ett_default = DEFAULT_EARLY_TERM_THRESHOLD) 
		: _container           (MD), 
		  _field               (F), 
		  _Comm                (PRINT_NOTHING, PRINT_NOTHING),
		  EARLY_TERM_THRESHOLD (ett_default) 
		{}

	//-- with Commentator
	MasseyDomain (const Commentator &C, unsigned long ett_default = DEFAULT_EARLY_TERM_THRESHOLD) 
		: _container           (), 
		  _field               (), 
		  _Comm                (C),
		  EARLY_TERM_THRESHOLD (ett_default)
		{}

	MasseyDomain (const Commentator &C, Sequence *D, unsigned long ett_default = DEFAULT_EARLY_TERM_THRESHOLD) 
		: _container           (D),
		  _field               (D->getField ()),
		  _Comm                (C),
		  EARLY_TERM_THRESHOLD (ett_default)
		{}
  
	MasseyDomain (const Commentator &C, Sequence *MD, const Field &F, unsigned long ett_default = DEFAULT_EARLY_TERM_THRESHOLD) 
		: _container           (MD), 
		  _field               (F),
		  _Comm                (C),
		  EARLY_TERM_THRESHOLD (ett_default)
		{}

        //-- Principal method
	template<class Polynomial>
	void operator () (Polynomial &C, bool full_poly = 0) {
		massey (C, full_poly);
	};

        //-- Domains access
	const Field &getField () const { return _field; }
	Sequence *getSequence () const { return _container; }

private:
	// -----------------------------------------------
	// Polynomial emulation
	// Only container aspects of polynomials
	// AND degree and valuation are needed !
	// -----------------------------------------------
#define _DEGINFTY_ -1

	// Degree of v
	template <class V>
	long v_degree (V& v) {
		long i = v.size()-1;

		if (i == _DEGINFTY_)
			return _DEGINFTY_;
		else if (!_field.isZero (v[i]))
			return i;

		// We must re-compute the degree :
		for (long j=i-1; j>=0 ; j--) {
			if (!_field.isZero((v)[j]) ) {
				v.resize(j+1);
				return j ;
			}
		}

		return _DEGINFTY_ ;
	}

	// Valuation of v
	template <class V>
	long v_val(V& v) {
		long i = v.size()-1;

		if (i == _DEGINFTY_)
			return _DEGINFTY_ ;
		else if (!_field.isZero (v[0]) )
			return 0;

		// We must compute the valuation :
		for (long j = 1; j <= i; j++)
			if (!_field.isZero ((v)[j])) return j ;

		return _DEGINFTY_ ;
	}

	// -------------------------------------------------------------------
	// Berlekamp/Massey algorithm with Massey's Sequence generation
	// -------------------------------------------------------------------

	template<class Polynomial>
	void massey (Polynomial &C, bool full_poly = 0) { 
//              const long ni = _container->n_row (), nj = _container->n_col ();
//              const long n = MIN(ni,nj);
		const long END = _container->size () + (full_poly ? DEFAULT_ADDITIONAL_ITERATION:0);
		const long n = END >> 1;

		Integer card;

		_Comm.start ("Massey", LVL_NORMAL,INTERNAL_DESCRIPTION) 
			<< END << endl;

		// ====================================================
		// Sequence and iterator initialization
		//
		typename Sequence::const_iterator _iter (_container->begin ());
		Polynomial S (END + 1);

		// -----------------------------------------------
		// Preallocation. No further allocation.
		//
		C.reserve (n + 1);    C.resize (1); _field.init (C[0], 1);
		Polynomial B (n + 1); B.resize (1); _field.init (B[0], 1);

		long L = 0;
		Element b, d, Ds;
		long x = 1, b_deg = 0, c_deg = 0, l_deg;

		_field.init (b, 1);

		for (long N = 0; N < END && x < EARLY_TERM_THRESHOLD; ++N, ++_iter) {
			if (!(N % 1000)) 
				_Comm.progress ("m-v prods", LVL_IMP, N, END);

			// ====================================================
			// Next coefficient in the sequence
			// Discrepancy computation
			// 
			d = S[N] = *_iter; 
//			d = _container->next (S[N]);

			for (long i = MIN (L, c_deg); i; --i)
				_field.axpyin (d, C[i], S[N - i]);

			if (_field.isZero (d)) {
				++x;
			} else {
				if (L > (N >> 1)) {
					// -----------------------------------------------
					// C = C + (Polynome(X,x,-d/b) * B);
					// 
					_field.divin (_field.neg (Ds, d), b);
					long i = l_deg = (x + b_deg);
					if (l_deg > c_deg) {
						C.resize (l_deg + 1);
						if (x > c_deg) {
							for (; i >= x; --i)
								_field.mul (C[i], Ds, B[i-x]);
							for (; i > c_deg; --i)
								_field.init (C[i], 0);
						} else {
							for (; i > c_deg; --i)
								_field.mul (C[i], Ds, B[i-x]);
							for (; i >= x; --i)
								_field.axpyin (C[i], Ds, B[i-x]);
						}
					} else {
						for (; i >= x; --i)
							_field.axpyin (C[i], Ds, B[i-x]);
					}
					// -----------------------------------------------
					c_deg = v_degree(C);
					++x;
				} else {
					// -----------------------------------------------
					// C = C + (Polynome(X,x,-d/b) * B); // B = C;
					// 
					_field.divin (_field.neg (Ds, d), b);
					long i = l_deg = x + b_deg;
					B.resize (C.size ());
					if (l_deg > c_deg) {
						C.resize (l_deg+1);
						if (x > c_deg) {
							for (; i >= x; --i)
								_field.mul (C[i], Ds, B[i-x]);
							for (; i > c_deg; --i)
								_field.init (C[i], 0);
						} else {
							for (; i > c_deg; --i)
								_field.mul (C[i], Ds, B[i-x]);
							for (; i >= x; --i)
								_field.axpy (C[i], Ds, B[i-x], B[i] = C[i]);
						}
					} else {
						for (i = c_deg; i > l_deg; --i)
							B[i] = C[i];
						for (; i >= x; --i)
							_field.axpy (C[i], Ds, B[i-x], B[i] = C[i] );
					}

					for (; i >= 0; --i) B[i] = C[i];

					// -----------------------------------------------
					L = N+1-L;
					b_deg = c_deg;
					c_deg = v_degree (C);
					b = d;
					x = 1;
				}
			}
			// ====================================================
		}
		_Comm.stop(LVL_NORMAL,PARTIAL_RESULT) 
			<< "Degree : " << v_degree (C) - v_val (C)
			<< " over GF(" << _field.cardinality (card) << "), 0:" << x << endl;
	}

public:
	// ---------------------------------------------
	// Massey
	//
	void pseudo_rank (unsigned long &rank) {
		vector<Element> phi;
		massey (phi, 0);
		rank = v_degree (phi) - v_val (phi);
	}
 
	void valence (Element &valence, unsigned long &rank) {
		_Comm.start ("Valence",LVL_NORMAL,INTERNAL_DESCRIPTION) << endl;
		vector<Element> phi;
		massey (phi, 1);
		rank = v_degree (phi) - v_val (phi);
		valence = phi[v_degree (phi)];

		_field.write (_Comm.stop(LVL_NORMAL,PARTIAL_RESULT), valence)
			<< ", rank : " << rank
			<< endl;
	}

	template<class Polynomial>
	void pseudo_minpoly (Polynomial &phi, unsigned long &rank) {
		massey (phi, 1);
		rank = v_degree (phi) - v_val (phi);

		if (phi.size () > 0) {
			long dp = v_degree (phi);
			for (long i = dp >> 1; i > 0; --i) {
				phi[0] = phi[i];
				phi[i] = phi[dp-i];
				phi[dp-i] = phi[0];
			}
			phi[0] = phi[dp];
			_field.init (phi[dp], 1);
		}
	}

};
 
}
    
#endif // __MASSEY_DOMAIN_H
