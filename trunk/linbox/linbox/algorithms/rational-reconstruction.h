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
//#define DEBUG_RR_BOUNDACCURACY
#define DEF_THRESH 50

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
	  
#ifdef RSTIMING
	mutable Timer tRecon, ttRecon;	
	mutable int _num_rec;
#endif	
		
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
	RationalReconstruction (const LiftingContainer& lcontainer, const Ring& r = Ring(), int THRESHOLD =DEF_THRESH) : 
		_lcontainer(lcontainer), _r(r), _threshold(THRESHOLD) {
			
		//if ( THRESHOLD < DEF_THRESH) _threshold = DEF_THRESH;
	}

	/** @memo Get the LiftingContainer
	 */
	const LiftingContainer& getContainer() const {
		return _lcontainer;
	}
		

	/** @memo  Handler to switch between different rational reconstruction strategy.
	 *  Allow  early termination and direct fast method
	 *  Switch is made by using a threshold as the third argument 
	 *  (default is set to that of constructor THRESHOLD
	 * 0    -> direct method
	 * > 0  -> early termination with 
	 */ 
	template <class Vector>
	bool getRational(Vector& num, Integer& den, int switcher) const { 
		if ( switcher > 0)
			return getRational2(num,den);
		else
			return getRational3(num,den);
	}

	template <class Vector>
	bool getRational(Vector& num, Integer& den) const { 
		if ( _threshold > 0)
			return getRational2(num,den);
		else
			return getRational3(num,den);
	}
	
	/** @memo Reconstruct a vector of rational numbers
	 *  from p-adic digit vector sequence.
	 *  An early termination technique is used.
	 *  Answer is a pair (numerator, common denominator)
	 *  The trick to reconstruct the raitonal solution (V. Pan) is implemented.
	 */
	template<class Vector>
	bool getRational1(Vector& num, Integer& den) const { 
 			
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
		// store den. upper bound
		Integer denbound;
		// store num.  upper bound
		Integer numbound;

		_r. init (modulus, 1);
		_r. init (denbound, 1);
		_r. init (numbound, 1);
	
		// some iterator
		typename Vector::iterator answer_p;
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
			
		_r. init (den, 1);
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
			cout<<"approximation mod p^"<<i<<" : \n";
			cout<<"[";
			for (size_t j=0;j< zz.size( )-1;j++)
				cout<<zz[j]<<",";
			cout<<zz.back()<<"]\n";
			cout<<"digit:\n";				
			for (size_t j=0;j< digit.size();++j)
				cout<<digit[j]<<",";
			cout<<endl;
#endif
			
			if ( i % _threshold && i < (int)_lcontainer.length() - _threshold) continue;
				
			for ( zz_p = zz.begin(), num_p = num.begin(), repeat_p = repeat.begin();
			      zz_p != zz. end();
			      ++ zz_p, ++ num_p, ++ repeat_p) {

				if (*repeat_p >= 2) continue;
					
				// a possible answer exits
				if ( *repeat_p) {
					
					_r. mul (tmp_i, den, *zz_p);
					_r. subin (tmp_i, *num_p);
					_r. remin (tmp_i, modulus);
					//cout<<tmp_i<<endl;
					// if the rational number works for _zz_p mod _modulus
					if (_r.isZero (tmp_i)) {
							
						++ *repeat_p;
						++ count;
					}
						
					// previus result is fake
					else {

						Integer tmp_nu, tmp_den, lcm, gcd;
							
						// try to reconstruct a rational number
						tmp = _r.reconstructRational(answer_p -> first,
									     answer_p -> second,
									     *zz_p, modulus,
									     numbound, denbound);
							
						// there exists a possible answer
						if (tmp) {
							*repeat_p = 1;
							_r. assin (*num_p, tmp_num);

							if (! _r. equal (tmp_den, den)) {
								Integer lcm, t1, t2;
								_r. lcm (lcm, tmp_den, den);
								_r. div (t1, lcm, tmp_den);
								_r. mulin (*num_p, t1);
								_r. div (t2, lcm, den);
								_r. assign (den, lcm);
								for (typename Vector::iterator tmp_p = num. begin (); tmp_p != num_p; ++ tmp_p)
									_r. mulin (*tmp_p, t2);
								}
							// no answer
						else *repeat_p = 0;
						}
						
					}
				}
					
				// not previous result
				else {

					// try if den if a multiple of the denominator of rational.
					Integer tmp_num, tmp_den;
					_r. assign (tmp_den, den);
					_r. mul (tmp_num, co_den, *zz_p);
					_r. remin (tmp_num, modulus);

					// assign tmp_num = one of tmp_num and tmp_num - modulus with smallest absolute value.
					Integer n_num;
					_r. sub (n_num, tmp_num, modulus);
					Integer abs_n, abs_nn;
					_r. abs (abs_n, tmp_num);
					_r. abs (abs_nn, n_num);
					if (_r. compare (abs_n, abs_nn) > 0)
						_r. assign (tmp_num, n_num);
					Integer g;
					_r. gcd (g, tmp_num, tmp_den);
					if (!_r. isUnit (g)) {
						_r. divin (tmp_num, g);
						_r. divin (tmp_den, g);
					}

					// check if (tmp_num, tmp_den) is an answer
					_r. abs (abs_n, tmp_num);
					_r. abs (abs_nn, tmp_den);
					//yes
					if (_r. compare (abs_n, numbound) < 0 && _r. compare (abs_nn, denbound) < 0) {
						*repeat_p = 1;
						continue;
					}

					// no
					tmp = _r.reconstructRational(tmp_num, tmp_den, *zz_p, modulus, denbound, numbound);

					if (tmp) {
						*repeat_p = 1;
						_r. assin (*num_p, tmp_num);
						if (! _r. equal (tmp_den, den)) {
							Integer lcm, t1, t2;
							_r. lcm (lcm, tmp_den, den);
							_r. div (t1, lcm, tmp_den);
							_r. mulin (*num_p, t1);
							_r. div (t2, lcm, den);
							_r. assign (den, lcm);
							for (typename Vector::iterator tmp_p = num. begin (); tmp_p != num_p; ++ tmp_p)
								_r. mulin (*tmp_p, t2);
						}
					}
				
					else *repeat_p = 0;
				}
			}

		}
		Integer g;
		_r. init (g, 1);
		_r. gcdin (g, den);
		for (num_p = num. begin(); num_p != num. end(); ++ num_p)
			_r. gcdin (g, *num_p);

		if (!_r. isOne (g) && !_r. isZero(g)) {
			for (num_p = num. begin(); num_p != num. end(); ++ num_p)
				_r. divin (*num_p, g);
			_r. divin (den, g);
		}
		return true; //lifted ok
	}
	
	/** @memo Reconstruct a vector of rational numbers
	 *  from p-adic digit vector sequence.
	 *  An early termination technique is used.
	 *  Answer is a vector of pair (num, den)
	 *
	 * Note, this may fail.
	 * Generically, the probability of failure should be 1/p^n where n is the number of elements being constructed
	 * since p is usually quite large this should be ok
	 */
	template<class Vector>
	bool getRational2(Vector& num, Integer& den) const { 
#ifdef RSTIMING
		ttRecon.clear();
		tRecon.start();
#endif
		linbox_check(num.size() == (size_t)_lcontainer.size());

		_r. init (den, 1);
		Integer prime = _lcontainer.prime(); 		        // prime used for lifting
		std::vector<size_t> accuracy(_lcontainer.size(), 0); 	// accuracy (in powers of p) of each answer so far
		Vector digit(_lcontainer.size());  		        // to store next digit
		Integer modulus;        		                // store modulus (power of prime)
		Integer prev_modulus;       	                        // store previous modulus
		Integer numbound, denbound;                             // current num/den bound for early termination
		size_t numConfirmed;                                       // number of reconstructions which passed twice
		_r.init(modulus, 0);
		std::vector<Integer> zz(_lcontainer.size(), modulus);   // stores each truncated p-adic approximation
		_r.init(modulus, 1);

		size_t len = _lcontainer.length(); // should be ceil(log(2*numbound*denbound)/log(prime))

		// should grow in rough proportion to overall num/denbound, 
		// but MUST have product less than p^i / 2
		// The heuristic used here is  
		// -bound when i=1  : N[0]=D[0]=sqrt(1/2)
		// -bound when i=len: N[len] = N*sqrt(p^len /2/N/D) , D[len] = D*sqrt(p^len /2/N/D)
		// -with a geometric series in between

		// 2 different ways to compute growing num bound (den is always picked as p^len/2/num)
		// when log( prime ^ len) is not too big (< 150 digits) we use logs stored as double; this way for very
		// small primes we don't accumulate multiplicative losses in precision
		
		// when log( prime ^ len) is very big we keep multiplying a factor into the num/den bound
		// at each level of reconstruction

		// note: currently it is usually the case that numbound > denbound. If it ever
		// happens that denbound >>> numbound, then denbound should be computed first at each step
		// and then numbound set to p^i / 2 / denbound, for less accumulated precision loss

		double dtmp;
		double half_log_p = 0.5 * log(_r.convert(dtmp, prime)); 
		const double half_log_2 = 0.34657359027997265472;
		double multy = 0;
		Integer numFactor;

		bool verybig = half_log_p * len > 75 * log(10.0);

		if (len > 1) {
			if (!verybig)
				multy = 0.5 * (log(_r.convert(dtmp, _lcontainer.numbound())) 
					       -log(_r.convert(dtmp, _lcontainer.denbound()))) / (len - 1);
			else {
				LinBox::integer iD, iN, pPower, tmp;
				_r.convert(iD, _lcontainer.numbound());
				_r.convert(iN, _lcontainer.denbound());
				_r.convert(pPower, prime);
				pPower = pow(pPower, len-1);

				tmp = pPower * iN;
				tmp /= iD;
				tmp = root(tmp, 2*len);
				_r.init(numFactor, tmp);

				// inital numbound is numFactor/sqrt(2)
				tmp = (tmp * 29) / 41;
				_r.init(numbound, tmp);
			}
		}
			
#ifdef DEBUG_RR
		cout << "nbound, dbound:" << _lcontainer.numbound() << ",  " << _lcontainer.denbound() << endl;
#endif

		typename Vector::iterator num_p;
		typename std::vector<Integer>::iterator zz_p;
		std::vector<size_t>::iterator accuracy_p;
		typename Vector::iterator digit_p;
			
		long tmp;
		Integer tmp_i;

		size_t i = 0;
			
		typename LiftingContainer::const_iterator iter = _lcontainer.begin(); 

		bool gotAll = false; //set to true if all values are reconstructed on a particular step

		// do until getting all answer
		do {
			++ i;		
#ifdef DEBUG_RR		
			cout<<"i: "<<i<<endl;
#endif
#ifdef RSTIMING
			tRecon.stop();
			ttRecon += tRecon;
#endif
			// get next p-adic digit
			bool nextResult = iter.next(digit);
			if (!nextResult) {
				cout << "ERROR in lifting container. Are you using <double> ring with large norm?" << endl;
				return false;
			}
#ifdef RSTIMING
			tRecon.start();
#endif			
			// preserve the old modulus
			_r.assign (prev_modulus, modulus);

			// upate _modulus *= _prime
			_r.mulin (modulus,  prime);

			// update truncated p-adic approximation
			for ( digit_p = digit.begin(), zz_p = zz.begin(); digit_p != digit.end(); ++ digit_p, ++ zz_p) 
				_r.axpyin(*zz_p, prev_modulus, *digit_p);
#ifdef DEBUG_RR
			cout<<"approximation mod p^"<<i<<" : \n";
			cout<<"[";
			for (size_t j=0;j< zz.size( )-1;j++)
				cout<<zz[j]<<",";
			cout<<zz.back()<<"]\n";
			cout<<"digit:\n";				
			for (size_t j=0;j< digit.size();++j)
				cout<<digit[j]<<",";
			cout<<endl;
#endif			
			if (verybig && i>1 && i<len)
				_r.mulin(numbound, numFactor);

			numConfirmed = 0;
			if ( !gotAll && (i % _threshold != 0) && (i + _threshold < len)) continue;
			// change to geometric skips?

			// update den and num bound, see above for details
			if (i == len) { //in this case we set bound to exact values
				_r. assign(numbound, _lcontainer.numbound());
				_r. assign(denbound, _lcontainer.denbound());
			}
			else {
				if (!verybig) 
					_r.init(numbound, exp(i*half_log_p - half_log_2 + multy * (i - 1)));

				Integer tmp;
				_r.init(tmp, 2);
				_r.mulin(tmp, numbound);
				_r.quo(denbound, modulus, tmp);
			}
#ifdef DEBUG_RR
			cout << "i, N, D bounds: " << i << ", " << numbound << ", " << denbound << endl;
#endif
			gotAll = true;
			bool justConfirming = true;
			int index = 0;
			// try to construct all unconstructed numbers
			for ( zz_p = zz.begin(), num_p = num.begin(), accuracy_p = accuracy.begin();
			      zz_p != zz.end();  ++ zz_p, ++ num_p, ++ accuracy_p, ++index) {

				if ( *accuracy_p == 0) {
					justConfirming = false;
					// if no answer yet (or last answer became invalid)
					// try to reconstruct a rational number
					// try if co_den if a multiple of the denominator of rational.
					Integer tmp_num, tmp_den;
					_r. assign (tmp_den, den);
					_r. mul (tmp_num, den, *zz_p);
					_r. remin (tmp_num, modulus);

					// assign tmp_num = one of tmp_num and tmp_num - modulus with smallest absolute value.
					Integer n_num;
					_r. sub (n_num, tmp_num, modulus);
					Integer abs_n, abs_nn;
					_r. abs (abs_n, tmp_num);
					_r. abs (abs_nn, n_num);
					if (_r. compare (abs_n, abs_nn) > 0)
						_r. assign (tmp_num, n_num);
				
					Integer g;
					_r. gcd (g, tmp_num, tmp_den);
					if (!_r. isUnit (g)) {
						_r. divin (tmp_num, g);
						_r. divin (tmp_den, g);
					}
					// check if (tmp_num, tmp_den) is an answer
					_r. abs (abs_n, tmp_num);
					_r. abs (abs_nn, tmp_den);
					// yes
					if (_r. compare (abs_n, numbound) < 0 && _r. compare (abs_nn, denbound) < 0) {
						*accuracy_p = i;
						continue;
					}
					//no
				  	justConfirming = false;
					// if no answer yet (or last answer became invalid)
					// try to reconstruct a rational number
					tmp = _r.reconstructRational(*num_p, tmp_den, *zz_p, modulus, numbound, denbound);
					// update 'accuracy' according to whether it worked or not
					if (tmp) {
						linbox_check (!_r.isZero(tmp_den));
						if (! _r. areEqual (tmp_den, den)) {
							Integer lcm, t1, t2;
							_r. lcm (lcm, tmp_den, den);
							_r. div (t1, lcm, tmp_den);
							_r. mulin (*num_p, t1);
							_r. div (t2, lcm, den);
							_r. assign (den, lcm);
							for (typename Vector::iterator tmp_p = num. begin (); tmp_p != num_p; ++ tmp_p)
								_r. mulin (*tmp_p, t2);
						}
					*accuracy_p = i;
					}
					else {
						*accuracy_p = 0;
						gotAll = false;
					}
				}
			}

			// if all were already done, check old approximations and try to reconstruct broken ones
			// also need to do this when we're on last iteration
			index = 0;
			if (justConfirming || i == len)
			for ( zz_p = zz.begin(), num_p = num.begin(), accuracy_p = accuracy.begin();
			      gotAll && zz_p != zz.end(); ++ zz_p, ++ num_p, ++ accuracy_p, index++) {
				
				if ( *accuracy_p < i ) {
					// check if the rational number works for _zz_p mod _modulus
					_r. mul (tmp_i, den, *zz_p);
					_r. subin (tmp_i, *num_p);
					_r. remin (tmp_i, modulus);
					if (_r.isZero (tmp_i)) {
						*accuracy_p = i;
						numConfirmed++;
					}
					else {
						// previous result is fake, reconstruct new answer
						Integer tmp_den;
						tmp = _r.reconstructRational(*num_p, tmp_den, *zz_p, modulus, numbound, denbound);
						if (tmp) {
							linbox_check (!_r.isZero(den));
							if (! _r. areEqual (tmp_den, den)) {
							Integer lcm, t1, t2;
							_r. lcm (lcm, tmp_den, den);
							_r. div (t1, lcm, tmp_den);
							_r. mulin (*num_p, t1);
							_r. div (t2, lcm, den);
							_r. assign (den, lcm);
							for (typename Vector::iterator tmp_p = num. begin (); tmp_p != num_p; ++ tmp_p)
								_r. mulin (*tmp_p, t2);
							}
							*accuracy_p = i;
						}
						else {
							*accuracy_p = 0;
							gotAll = false;
						}
					}
				}
			}
		}
		while (numConfirmed < _lcontainer.size() && i < len);
		//still probabilstic, but much less so
#ifdef RSTIMING
		tRecon.stop();
		ttRecon += tRecon;
#endif	
#ifdef DEBUG_RR_BOUNDACCURACY
		cout << "Computed " << i << " digits out of estimated " << len << endl;
#endif


		Integer g;
		_r. init (g, 0);
		_r. gcdin (g, den);
		for (num_p = num. begin(); num_p != num. end(); ++ num_p)
			_r. gcdin (g, *num_p);

		if (!_r. isOne (g) && !_r. isZero(g)) {
			for (num_p = num. begin(); num_p != num. end(); ++ num_p)
				_r. divin (*num_p, g);
			_r. divin (den, g);
		}
		return true; //lifted ok, assuming norm was correct
	}



        void PolEval(Vector& y, std::vector<Vector>::const_iterator &Pol, size_t deg, Integer &x) const {
		
		
		if (deg == 1){
			for (size_t i=0;i<y.size();++i)
				_r.assign(y[i],(*Pol)[i]);			
		}
		else{
			size_t deg_low, deg_high;
			deg_high = deg/2;
			deg_low  = deg - deg_high;
			Integer zero;
			_r.init(zero,0);
			Vector y1(y.size(),zero), y2(y.size(),zero);
			Integer x1=x, x2=x;

			PolEval(y1, Pol, deg_low, x1);

			std::vector<Vector>::const_iterator Pol_high= Pol+deg_low;
			PolEval(y2, Pol_high, deg_high, x2);
						
			for (size_t i=0;i< y.size();++i)
				_r.axpy(y[i],x1,y2[i],y1[i]);			
			_r.mul(x,x1,x2);
		}
	}


	/** @memo Reconstruct a vector of rational numbers
	 *  from p-adic digit vector sequence.
	 *  compute all digits and reconstruct rationals only once
	 *  Result is a vector of numerators and one common denominator
	 */
	template<class Vector1>
	bool getRational3(Vector1& num, Integer& den) const { 
		
#ifdef RSTIMING
		ttRecon.clear();
		tRecon.start();
#endif
		linbox_check(num.size() == (size_t)_lcontainer.size());
		
		// prime 
		Integer prime = _lcontainer.prime();
		
		// length of whole approximation
		size_t length=_lcontainer.length();
		
		// size of solution
		size_t size= _lcontainer.size();
		 
		
		Integer zero;
		_r.init(zero,0);
		Vector zero_digit(_lcontainer.size(),zero);	
		
		// store approximation as a polynomial and evaluate by baby step giant step
		std::vector<Vector>  digit_approximation(length,zero_digit); 

		// store real approximation
		Vector real_approximation(size,zero);


		// store modulus (intially set to 1)
		Integer modulus;
		_r.init(modulus, 1);
		
		// denominator upper bound
		Integer denbound;
		_r.assign(denbound,_lcontainer.denbound());
		
		// numerator  upper bound
		Integer numbound;
		_r.assign(numbound,_lcontainer.numbound());
	
#ifdef RSTIMING
		tRecon.stop();
		ttRecon += tRecon;
#endif		
#ifdef LIFTING_PROGRESS
		Commentator lifting_commentator;
		lifting_commentator.start("Padic Lifting","LinBox::LiftingContainer",_lcontainer.length());
#endif
	
		//Timer eval_horner,eval_horn;
		//eval_horner.clear();
		// Compute all the approximation using liftingcontainer
		typename LiftingContainer::const_iterator iter = _lcontainer.begin();
		for (size_t i=0 ; iter != _lcontainer.end() && iter.next(digit_approximation[i]);++i) {

#ifdef LIFTING_PROGRESS			
			lifting_commentator.progress(i);
#endif
			//eval_horn.start();
			//for (size_t j=0;j<size;++j)
			//	_r.axpyin(real_approximation[j],modulus, digit_approximation[i][j]);							
			//eval_horn.stop();
			//eval_horner+=eval_horn;

			_r.mulin(modulus,prime); 
		}

#ifdef LIFTING_PROGRESS			
		lifting_commentator.stop ("Done", "Done", "LinBox::LinBox::LiftingContainer");	 
#endif
	
		// problem occured during lifting
		if (iter!= _lcontainer.end()){
			cout << "ERROR in lifting container. Are you using <double> ring with large norm?" << endl;
			return false;
		}
		
#ifdef RSTIMING
		tRecon.start();
#endif
	
		//Timer eval_dac, eval_bsgs;
		
		//eval_bsgs.start();
		// sqrt of approximation's length
		//int sqrt_length= (int) sqrt((double) length);				
		/*
		 * Baby-Step/ Giant-Step Polynomial evaluation of digit approximation
		 */
		/*
		{
			// store intermediate baby-step/ giant-step polynomial evaluation of the approximation in prime
			std::vector<Vector> baby_approx (sqrt_length+1,zero_digit);
					
			// perform baby-step
			int skip=-sqrt_length;
			for (int k=0;k<sqrt_length;++k){
				skip+=sqrt_length;
				for (int i= sqrt_length-1; i>=0; --i) 
					for (size_t j=0;j<size;++j) {					
						_r.mulin(baby_approx[k][j] , prime);
						_r.addin(baby_approx[k][j], digit_approximation[skip+i][j]);	
					}
			}
			
			for (int i= length -1; i>= skip+sqrt_length; --i)
				for (size_t j=0;j<size;++j) {					
					_r.mulin(baby_approx[sqrt_length][j] , prime);
					_r.addin(baby_approx[sqrt_length][j], digit_approximation[i][j]);				
				}
			
			LinBox::integer p_to_sqrt, p;
			_r.convert(p,prime);
			p_to_sqrt= pow(p,sqrt_length);
			Integer prime_to_sqrt;
			_r.init(prime_to_sqrt, p_to_sqrt);
			
					
			// perform giant step
			for (int i= sqrt_length; i>= 0; --i)
				for (size_t j=0;j<size;j++) {
					_r.mulin(real_approximation[j] , prime_to_sqrt );
					_r.addin(real_approximation[j], baby_approx[i][j]);
				}
		}
		eval_bsgs.stop();
		*/

		//eval_dac.start();
		Integer xeval=prime;
		typename std::vector<Vector>::const_iterator poly_digit= digit_approximation.begin();
		PolEval(real_approximation, poly_digit, length, xeval);

		//eval_dac.stop();
		
// 		integer modulus_size;
// 		_r.convert(modulus_size,modulus);
// 		cout<<"number of bit : "<< modulus_size.bitsize()<<endl;
// 		cout<<"length        : "<< length<<endl;
// 		cout<<"prime         : "<< prime<<endl;
// 		cout<<"evaluation divide&conquer  : "<<eval_dac<<endl;
// 		cout<<"evaluation baby/giant step : "<<eval_bsgs<<endl;
// 		cout<<"evaluation horner method   : "<<eval_horner<<endl;


		/* 
		 * dumb rational reconstruction (this is just for timing comparison)
		 */
// 		{
// 			Timer dumb_ratrecon;
// 			dumb_ratrecon.start();
// 			Vector den_r(num.size()), num_r(num.size());			
// 			typename Vector::iterator   iter_a  = real_approximation.begin();
// 			typename Vector::iterator   iter_n  = num_r.begin();
// 			typename Vector::iterator   iter_d  = den_r.begin();
						
// 			for (size_t i=0; iter_a != real_approximation.end(); ++iter_a, ++ iter_n, ++iter_d, ++i){				
// 				if (!_r.reconstructRational(*iter_n, *iter_d,
// 						    *iter_a, modulus, numbound, denbound))
// 				{					
// 					cout << "ERROR in reconstruction ?\n" << endl;
// 				}
// 			}
// 			dumb_ratrecon.stop();
// 			cout<<"full rational reconstruction : "<<dumb_ratrecon.usertime()<<endl;
// 		}
		

		/*
		 * Rational Reconstruction of each coefficient according to a common denominator
		 */
		
		//Timer ratrecon;
		//ratrecon.start();
		Integer common_den, common_den_mod_prod, bound,two,tmp;
		_r.init(common_den,1);
		_r.init(common_den_mod_prod,1);
		_r.init(two,2);		
		
		Vector denominator(num.size());

		int counter=0;
		typename Vector::iterator   iter_approx = real_approximation.begin();
		typename Vector1::iterator  iter_num    = num.begin();
		typename Vector::iterator   iter_denom  = denominator.begin();
		
		//numbound=denbound;
		
		for (size_t i=0; iter_approx != real_approximation.end(); ++iter_approx, ++ iter_num, ++iter_denom, ++i){
			_r.mulin( *iter_approx , common_den_mod_prod);
			_r.remin( *iter_approx , modulus);
			if (!_r.reconstructRational(*iter_num, *iter_denom,
						    *iter_approx, modulus, numbound, denbound))
				{					
					cout << "ERROR in reconstruction ?\n" << endl;
#ifdef DEBUG_RR
					cout<<" try to reconstruct :\n";
					cout<<"approximation: "<<*iter_approx<<endl;
					cout<<"modulus: "<<modulus<<endl;
					cout<<"numbound: "<<numbound<<endl;
					cout<<"denbound: "<<denbound<<endl;
#endif
					return false;
				}			
			
			_r.mulin(common_den, *iter_denom);
			if (i != size-1){
				if (! _r.isUnit(*iter_denom)) {counter++;

				_r.quoin(denbound , *iter_denom);
				_r.mul(bound, denbound,numbound);
				_r.mulin(bound,two);
				_r.div(tmp,modulus,prime);
				while(tmp > bound) {
					_r.assign(modulus,tmp);
					_r.div(tmp,modulus,prime);
				}
				_r.rem(tmp , *iter_denom , modulus);
				_r.remin(common_den_mod_prod , modulus);
				_r.mulin(common_den_mod_prod , tmp);
				_r.remin(common_den_mod_prod , modulus);
				}	
			}
			
		}
		       
		typename Vector1::reverse_iterator rev_iter_num   = num.rbegin();
		typename Vector::reverse_iterator  rev_iter_denom = denominator.rbegin();
		_r.init(tmp,1);
		for (; rev_iter_num != num.rend(); ++rev_iter_num, ++rev_iter_denom){
			_r.mulin(*rev_iter_num,tmp);
			_r.mulin(tmp, *rev_iter_denom);
		}

		den = common_den;
		
		//ratrecon.stop();
		//cout<<"partial rational reconstruction : "<<ratrecon.usertime()<<endl;
#ifdef RSTIMING
		tRecon.stop();
		ttRecon += tRecon;
		_num_rec=counter;
#endif
		
		return true;

	} // end of getRational3
	
};
	
}

#endif
