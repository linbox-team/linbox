/* linbox/algorithms/massey-domain.h
 * Copyright (C) 1999, 2001 Jean-Guillaume Dumas, Bradford Hovinen
 *
 * Written by Jean-Guillaume Dumas <Jean-Guillaume.Dumas@imag.fr>,
 *            Bradford Hovinen <hovinen@cis.udel.edu>
 *            JGD 12.06.2002 : 
 *            			-- Put back domain.zero
 *            			-- Put back domain.one
 *            			-- Put back full_poly
 *            			-- Put back pseudo_minpoly as it was before
 *            			-- not yet fully checked since previous changes
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */
#ifndef __LINBOX_massey_domain_H
#define __LINBOX_massey_domain_H

// ======================================================================= 
// Linbox project 1999
// Domain Massey
// - Computation is stopped when the polynomials remain the same
//   for more than EARLY_TERM_THRESOLD
// - When minimal polynomial equals characteristic polynomial,
//   2 additional iterations are needed to compute it 
//   (parameter DEFAULT_ADDITIONAL_ITERATION), but those
//   iterations are not needed for the rank
// Time-stamp: <27 Aug 01 18:18:12 Jean-Guillaume.Dumas@imag.fr> 
// ======================================================================= 

#include "linbox/util/commentator.h"
#include "linbox/vector/reverse.h"
#include "linbox/vector/subvector.h"
#include "linbox/vector/vector-domain.h"
#include "linbox/util/timer.h"

namespace LinBox 
{

#ifndef MIN
#  define MIN(a,b) ((a)<(b)?(a):(b))
#endif

#define DEFAULT_EARLY_TERM_THRESHOLD 20
#ifndef DEFAULT_ADDITIONAL_ITERATION
#define DEFAULT_ADDITIONAL_ITERATION 2
#endif

const long _DEGINFTY_ = -1;

/** \brief Berlekamp/Massey algorithm. 

   Domain Massey
   - Computation is stopped when the polynomials remain the same
     for more than EARLY_TERM_THRESOLD
   - When minimal polynomial equals characteristic polynomial,
     2 additional iterations are needed to compute it 
     (parameter DEFAULT_ADDITIONAL_ITERATION), but those
     iterations are not needed for the rank
*/
template<class Field, class Sequence>
class MasseyDomain {
    private:
	Sequence            *_container;
	Field                _F;
	VectorDomain<Field>  _VD;
	unsigned long         EARLY_TERM_THRESHOLD;

#ifdef INCLUDE_TIMING
	// Timings
	double                _discrepencyTime;
	double                _fixTime;
#endif // INCLUDE_TIMING

    public:
	typedef typename Field::Element Element;

        //-- Constructors
	MasseyDomain (unsigned long ett_default = DEFAULT_EARLY_TERM_THRESHOLD) 
		: _container           (), 
		  _F                   (), 
		  _VD                  (_F),
		  EARLY_TERM_THRESHOLD (ett_default)
		{}

	MasseyDomain (const MasseyDomain<Field, Sequence> &M, unsigned long ett_default = DEFAULT_EARLY_TERM_THRESHOLD) 
		: _container           (M._container), 
		  _F                   (M._F), 
		  _VD                  (M._F),
		  EARLY_TERM_THRESHOLD (ett_default)
		{}

	MasseyDomain (Sequence *D, unsigned long ett_default = DEFAULT_EARLY_TERM_THRESHOLD) 
		: _container           (D), 
		  _F                   (D->getField ()),
		  _VD                  (D->getField ()),
		  EARLY_TERM_THRESHOLD (ett_default)
		{}
  
	MasseyDomain (Sequence *MD, const Field &F, unsigned long ett_default = DEFAULT_EARLY_TERM_THRESHOLD) 
		: _container           (MD), 
		  _F                   (F),
		  _VD                  (F),
		  EARLY_TERM_THRESHOLD (ett_default) 
		{}

        //-- Principal method
	template<class Polynomial>
	long operator () (Polynomial &C, bool full_poly = false) {
		return massey (C, full_poly);
	}

        //-- Domains access
	const Field &getField    () const { return _F; }
	Sequence    *getSequence () const { return _container; }

#ifdef INCLUDE_TIMING
	double       discrepencyTime () const { return _discrepencyTime; }
	double       fixTime         () const { return _fixTime; }
#endif // INCLUDE_TIMING

    private:
	// -----------------------------------------------
	// Polynomial emulation
	// Only container aspects of polynomials
	// AND degree and valuation are needed !
	// -----------------------------------------------

	// Degree of v
	template <class V>
	long v_degree (V& v) {
		long i = v.size()-1;

		if (i == _DEGINFTY_)
			return _DEGINFTY_;

		else if (!_F.isZero (v[i]))
			return i;

		// We must re-compute the degree :
		for (long j = i - 1; j >= 0; j--) {
			if (!_F.isZero (v[j])) {
				v.resize (j + 1);
				return j;
			}
		}

		return _DEGINFTY_ ;
	}

	// Valuation of v
	template <class V>
	long v_val(V& v) {
		long i = v.size() - 1;

		if (i == _DEGINFTY_)
			return _DEGINFTY_;

		else if (!_F.isZero (v[0]))
			return 0;

		// We must compute the valuation :
		for (long j = 1; j <= i; j++)
			if (!_F.isZero ((v)[j])) return j ;

		return _DEGINFTY_ ;
	}

	// -------------------------------------------------------------------
	// Berlekamp/Massey algorithm with Massey's Sequence generation
	// -------------------------------------------------------------------

	template<class Polynomial>
	long massey (Polynomial &C, bool full_poly = false) { 
//              const long ni = _container->n_row (), nj = _container->n_col ();
//              const long n = MIN(ni,nj);
		const long END = _container->size () + (full_poly ? DEFAULT_ADDITIONAL_ITERATION:0);
		const long n = END >> 1;

#ifdef INCLUDE_TIMING
		Timer timer;

		_discrepencyTime = _fixTime = 0.0;
#endif // INCLUDE_TIMING

		integer card;

		commentator.start ("Massey", "masseyd", END);

		// ====================================================
		// Sequence and iterator initialization
		//
		typename Sequence::const_iterator _iter (_container->begin ());
		Polynomial S (END + 1);
		
		Element Zero, One;
		_F.init(Zero, 0);
		_F.init(One, 1);

		// -----------------------------------------------
		// Preallocation. No further allocation.
		//
		C.reserve    (n + 1); C.resize (1); _F.assign (C[0], One);
		Polynomial B (n + 1); B.resize (1); _F.assign (B[0], One);

		long L = 0;
		Element b, d, Ds;
		long x = 1, b_deg = 0, c_deg = 0, l_deg;
		long COMMOD = (END > 40) ? (END / 20) : 2;                

		_F.assign (b, One);


		for (long N = 0; N < END && x < (long) EARLY_TERM_THRESHOLD; ++N, ++_iter) {

			if (!(N % COMMOD)) 
				commentator.progress (N);

			// ====================================================
			// Next coefficient in the sequence
			// Discrepancy computation
			// 
			S[N] = *_iter; 
                        
			// 
#ifdef INCLUDE_TIMING
			timer.start ();
#endif // INCLUDE_TIMING

			long poly_len = MIN (L, c_deg);
			Subvector<typename Polynomial::iterator> Cp (C.begin () + 1, C.begin () + poly_len + 1);
			Subvector<typename Polynomial::iterator> Sp (S.begin () + (N - poly_len), S.begin () + N);
			ReverseVector<Subvector<typename Polynomial::iterator> > Spp (Sp);
			_VD.dot (d, Cp, Spp);

			_F.addin (d, S[N]);

#ifdef INCLUDE_TIMING
			timer.stop ();

			_discrepencyTime += timer.realtime ();

			timer.start ();
#endif // INCLUDE_TIMING

			if (_F.isZero (d)) {
				++x;
			} else {
				if (L > (N >> 1)) {
					// -----------------------------------------------
					// C = C + (Polynome(X,x,-d/b) * B);
					// 
					_F.divin (_F.neg (Ds, d), b);
					long i = l_deg = (x + b_deg);
					if (l_deg > c_deg) {
						C.resize (l_deg + 1);
						if (x > c_deg) {
							for (; i >= x; --i)
								_F.mul (C[i], Ds, B[i-x]);
							for (; i > c_deg; --i)
								_F.assign (C[i], Zero);
						} else {
							for (; i > c_deg; --i)
								_F.mul (C[i], Ds, B[i-x]);
							for (; i >= x; --i)
								_F.axpyin (C[i], Ds, B[i-x]);
						}
					} else {
						for (; i >= x; --i)
							_F.axpyin (C[i], Ds, B[i-x]);
					}
					// -----------------------------------------------
					c_deg = v_degree(C);
					++x;
				} else {
					// -----------------------------------------------
					// C = C + (Polynome(X,x,-d/b) * B); 					// 
					_F.divin (_F.neg (Ds, d), b);
					long i = l_deg = x + b_deg;
					B.resize (C.size ());
					if (l_deg > c_deg) {
						C.resize (l_deg+1);
						if (x > c_deg) {
							for (; i >= x; --i)
								_F.mul (C[i], Ds, B[i-x]);
							for (; i > c_deg; --i)
								_F.assign (C[i], Zero);
						} else {
							for (; i > c_deg; --i)
								_F.mul (C[i], Ds, B[i-x]);
							for (; i >= x; --i)
								_F.axpy (C[i], Ds, B[i-x], B[i] = C[i]);
						}
					} else {
						for (i = c_deg; i > l_deg; --i)
							B[i] = C[i];
						for (; i >= x; --i)
							_F.axpy (C[i], Ds, B[i-x], B[i] = C[i] );
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

#ifdef INCLUDE_TIMING
			timer.stop ();

			_fixTime += timer.realtime ();
#endif // INCLUDE_TIMING
		}

                commentator.stop ("done", NULL, "masseyd");
//		commentator.stop ("Done", "Done", "LinBox::MasseyDomain::massey");

		return L;
	}

public:
	// ---------------------------------------------
	// Massey
	//
	void pseudo_rank (unsigned long &rank) {
		std::vector<Element> phi;
		massey (phi, 0);
		rank = v_degree (phi) - v_val (phi);
	}
 
	void valence (Element &valence, unsigned long &rank) {
		commentator.start ("Valence", "LinBox::MasseyDomain::valence");

		std::vector<Element> phi;
		massey (phi, 1);
		rank = v_degree (phi) - v_val (phi);
		valence = phi[v_degree (phi)];

		commentator.stop ("Done", "Done", "LinBox::MasseyDomain::valence");
	}

	template<class Polynomial>
	unsigned long pseudo_minpoly (Polynomial &phi, unsigned long &rank, bool full_poly = true) {
		unsigned long L = massey (phi, full_poly);
		long dp = v_degree(phi);
		rank = dp - v_val (phi);
        	if (phi.size()) {
			for(long i = dp >> 1;i > 0; --i) {
				phi[0] = phi[i];
				phi[i] = phi[dp-i];
				phi[dp-i] = phi[0];
			}
			phi[0] = phi[dp];
			_F.init (phi[dp], 1UL);
		}
                return L;
	}

	template<class Polynomial>
	void minpoly (Polynomial &phi, unsigned long &rank, bool full_poly = true) {
            long dp = massey (phi, full_poly);
            rank = v_degree(phi) - v_val (phi);
            if (phi.size () > 0) {
                phi.resize (dp+1);
                for (long i = dp >> 1; i > 0; --i)
                    std::swap (phi[i], phi[dp-i]);
                phi[0] = phi[dp];
                _F.init (phi[dp], 1UL);
            }
	}
};
 
}
    
#endif // __LINBOX_massey_domain_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
