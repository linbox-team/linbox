/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
/* linbox/algorithms/lifting-container-base.h
 * Copyright (C) 2004 Zhendong Wan, Pascal Giorgi
 *
 * Written by Zhendong Wan <wan@mail.eecis.udel.edu>
 * Modified by Pascal Giorgi <pascal.giorgi@ens-lyon.fr>
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

#ifndef __LINBOXX__RECONSTRUCTION_H__
#define __LINBOXX__RECONSTRUCTION_H__

#include <linbox/util/debug.h>

//#define DEBUG_RR

namespace LinBox {

	

	/// @memo Limited doc so far.  Used, for instance, after LiftingContainer.
template< class _LiftingContainer >
class RationalReconstruction {
		
public:
	typedef _LiftingContainer                  LiftingContainer;
	typedef typename LiftingContainer::Ring                Ring;
	typedef typename Ring::Element                      Integer;
	typedef typename LiftingContainer::IVector           Vector;		
	typedef typename LiftingContainer::Field              Field;		
	typedef typename Field::Element                     Element;
	  
		
	// data
protected:
		
	// pointer to digit generator
	const LiftingContainer& _lcontainer;
		
	// Ring
	Ring _r;
		
	// store early termination threshold.
	int _threshold;
		
public:
			
	/** @memo Constructor 
	 *  maybe use different ring than the ring in lcontainer
	 */
	RationalReconstruction (const LiftingContainer& lcontainer, const Ring& r = Ring(), int THRESHOLD =10) : 
		_lcontainer(lcontainer), _r(r), _threshold(THRESHOLD) {
			
		if ( THRESHOLD < 10) _threshold = 10;			
	}

	/** @memo Get the LiftingContainer
	 */
	const LiftingContainer& getContainer() const {
		return _lcontainer;
	}
		
	/** @memo Reconstruct a vector of rational numbers
	 *  from p-adic digit vector sequence.
	 * An early termination technique is used.
	 *  Answer is a vector of pair (num, den)
	 */
	template<class Vector1>
	bool getRational(Vector1& answer) const { 
 			
		linbox_check(answer.size() == (size_t)_lcontainer.size());

		// prime 
		Integer prime = _lcontainer.prime();


		// how answer we get so far
		int count =0;

		// frequency of same rational reconstructed
		std::vector<int> repeat(_lcontainer.size(),0);
		
		// store next digit
		Vector digit(_lcontainer.size()); 
			
		// store modulus
		Integer modulus;

		// store the integer number		
		_r. init(modulus,0);

		std::vector<Integer> zz(_lcontainer.size(), modulus);

		// initially set modulus to be 1
		_r.init(modulus, 1);

		// store den. upper bound
		Integer denbound;
			
		// store num.  upper bound
		Integer numbound;
			
		_r. init (denbound, 1);

		_r. init (numbound, 1);
	
		// some iterator
		typename Vector1::iterator answer_p;

		typename std::vector<Integer>::iterator zz_p;
			
		std::vector<int>::iterator repeat_p;
			
		typename Vector::iterator digit_p;
			
		long tmp;

		// store previous modulus
		Integer pmodulus;
			
		// temporary integer
		Integer tmp_i;
			
		/** due to the slow rational reconstruction
		 * Change it reconstruct the rational number at every 
		 * _threshold.
		 *  By z. wan
		 */
		int i = 0;
			
		typename LiftingContainer::const_iterator iter = _lcontainer.begin(); 

		// do until getting all answer
		while (count < _lcontainer.size() && iter != _lcontainer.end()) {
			
			++ i;		
#ifdef DEBUG_RR		
			cout<<"i: "<<i<<endl;
#endif
			// get next p-adic digit
			bool nextResult = iter.next(digit);
			if (!nextResult) {
				cout << "ERROR in lifting container. Are you using <double> ring with large norm?" << endl;
				return false;
			}
				
			// preserve the old modulus
			_r.assign (pmodulus, modulus);

			// upate _modulus *= _prime
			_r.mulin (modulus,  prime);

			// update den and num bound
			if (( i % 2) == 0) {
					
				_r. mulin (denbound, prime);

				_r. assign (numbound, denbound);
			}
				
			for ( digit_p = digit.begin(), repeat_p = repeat.begin(), zz_p = zz.begin();
			      digit_p != digit.end();
			      ++ digit_p, ++ repeat_p, ++ zz_p) {
					
				// already get answer
				if ( *repeat_p >= 2) continue;
					
				// update *zz_p += pmodulus * (*digit_p)
				_r.axpyin(*zz_p, pmodulus, *digit_p);
			}
				
#ifdef DEBUG_RR
			cerr<<"approximation mod p^"<<i<<" : \n";
			cerr<<"[";
			for (size_t j=0;j< zz.size( )-1;j++)
				cerr<<zz[j]<<",";
			cerr<<zz.back()<<"]\n";
			cerr<<"digit:\n";				
			for (size_t j=0;j< digit.size();++j)
				cerr<<digit[j]<<",";
			cerr<<endl;
#endif
			
			if ( i % _threshold && i < (int)_lcontainer.length() - _threshold) continue;
				
			for ( zz_p = zz.begin(), answer_p = answer.begin(), repeat_p = repeat.begin();
			      zz_p != zz. end();
			      ++ zz_p, ++ answer_p, ++ repeat_p) {

				if (*repeat_p >= 2) continue;
					
				// a possible answer exits
				if ( *repeat_p) {
					
					_r. mul (tmp_i, answer_p -> second, *zz_p);
					_r. subin (tmp_i, answer_p -> first);
					_r. remin (tmp_i, modulus);
					//cerr<<tmp_i<<endl;
					// if the rational number works for _zz_p mod _modulus
					if (_r.isZero (tmp_i)) {
							
						++ *repeat_p;
							
						++ count;
							
					}
						
					// previus result is fake
					else {
							
						// try to reconstruct a rational number
						tmp = _r.reconstructRational(answer_p -> first,
									     answer_p -> second,
									     *zz_p, modulus,
									     numbound, denbound);
							
						// there exists a possible answer
						if (tmp) *repeat_p = 1;
							
						// no answer
						else *repeat_p = 0;
					}
						
				}
					
				// not previous result
				else {

					tmp = _r.reconstructRational(answer_p -> first, 
								     answer_p -> second,
								     *zz_p, modulus,
								     denbound, numbound);
						
					if (tmp) *repeat_p = 1; 
 
				}
			}

		}
		return true; //lifted ok
	}

	/** @memo Reconstruct a vector of rational numbers
	 *  from p-adic digit vector sequence.
	 *  An early termination technique is used.
	 *  Answer is a vector of pair (num, den)
	 */
	template<class Vector1>
	bool getRational2(Vector1& answer) const { 
 			
		linbox_check(answer.size() == (size_t)_lcontainer.size());

		Integer prime = _lcontainer.prime(); 		        // prime used for lifting
		std::vector<int> accuracy(_lcontainer.size(), 0); 	// accuracy (in powers of prime) of each component of answer so far
		Vector digit(_lcontainer.size());  		        // to store next digit
		Integer modulus;			                // store modulus (power of prime)
		Integer prev_modulus;       	                        // store previous modulus
		bool gotAll;                                            // are we done?
		int numConfirmed;                                       // number of reconstructions which passed at least twice
		_r.init(modulus, 0);
		std::vector<Integer> zz(_lcontainer.size(), modulus);   // stores the truncated p-adic approximation of each component
		_r.init(modulus, 1);

		Integer numbound, denbound;
		_r. init (denbound, 1);
		_r. init (numbound, 1);
	
		typename Vector1::iterator answer_p;
		typename std::vector<Integer>::iterator zz_p;
		std::vector<int>::iterator accuracy_p;
		typename Vector::iterator digit_p;
			
		long tmp;
		Integer tmp_i;

		int i = 0;
			
		typename LiftingContainer::const_iterator iter = _lcontainer.begin(); 

		// do until getting all answer
		do {
			++ i;		
#ifdef DEBUG_RR		
			cout<<"i: "<<i<<endl;
#endif
			// get next p-adic digit
			bool nextResult = iter.next(digit);
			if (!nextResult) {
				cout << "ERROR in lifting container. Are you using <double> ring with large norm?" << endl;
				return false;
			}
				
			// preserve the old modulus
			_r.assign (prev_modulus, modulus);

			// upate _modulus *= _prime
			_r.mulin (modulus,  prime);

			// update den and num bound
			if (( i % 2) == 0) {
				_r. mulin (denbound, prime);
				_r. assign (numbound, denbound);
			}
				
			// update truncated p-adic approximation
			for ( digit_p = digit.begin(), zz_p = zz.begin(); digit_p != digit.end(); ++ digit_p, ++ zz_p) 
				_r.axpyin(*zz_p, prev_modulus, *digit_p);
#ifdef DEBUG_RR
			cerr<<"approximation mod p^"<<i<<" : \n";
			cerr<<"[";
			for (size_t j=0;j< zz.size( )-1;j++)
				cerr<<zz[j]<<",";
			cerr<<zz.back()<<"]\n";
			cerr<<"digit:\n";				
			for (size_t j=0;j< digit.size();++j)
				cerr<<digit[j]<<",";
			cerr<<endl;
#endif			
			numConfirmed = 0;
			if ( i % _threshold && i < (int)_lcontainer.length() - _threshold) continue;			
			// change to geometric skips

			gotAll = true;
			bool justConfirming = true;
			// try to construct all unconstructed numbers
			for ( zz_p = zz.begin(), answer_p = answer.begin(), accuracy_p = accuracy.begin();
			      zz_p != zz.end();  ++ zz_p, ++ answer_p, ++ accuracy_p) {

				if ( *accuracy_p == 0) {
					justConfirming = false;
					// if no answer yet (or last answer became invalid)
					// try to reconstruct a rational number
					tmp = _r.reconstructRational(answer_p -> first,
								     answer_p -> second,
								     *zz_p, modulus,
								     numbound, denbound);
					// update 'accuracy' according to whether it worked or not
					if (tmp) {
						*accuracy_p = i;
						linbox_check (!_r.isZero(answer_p->second));
					}
					else {
						*accuracy_p = 0;
						gotAll = false;
					}
				}
			}

			// if all were already done, check old approximations and try to reconstruct broken ones
			if (justConfirming)
			for ( zz_p = zz.begin(), answer_p = answer.begin(), accuracy_p = accuracy.begin();
			      gotAll && zz_p != zz.end(); ++ zz_p, ++ answer_p, ++ accuracy_p) {
				
				if ( *accuracy_p < i ) {
					// check if the rational number works for _zz_p mod _modulus
					_r. mul (tmp_i, answer_p -> second, *zz_p);
					_r. subin (tmp_i, answer_p -> first);
					_r. remin (tmp_i, modulus);
					if (_r.isZero (tmp_i)) {
						*accuracy_p = i;
						numConfirmed++;
					}
					else {
						// previous result is fake, reconstruct new answer
						tmp = _r.reconstructRational(answer_p -> first,
									     answer_p -> second,
									     *zz_p, modulus,
									     numbound, denbound);
						if (tmp) {
							*accuracy_p = i;
							linbox_check (!_r.isZero(answer_p->second));
						}
						else {
							*accuracy_p = 0;
							gotAll = false;
						}
					}
				}
			}
		}
		while (numConfirmed < _lcontainer.size() && iter != _lcontainer.end());
		//still probabilstic, but much less so

		//while (iter != _lcontainer.end());
		return true; //lifted ok, assuming norm was correct
	}
};
		
}

#endif
