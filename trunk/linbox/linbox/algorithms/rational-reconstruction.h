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

#ifndef __LINBOX_reconstruction_H
#define __LINBOX_reconstruction_H

#include <linbox/algorithms/rational-reconstruction-base.h>
#include <linbox/algorithms/classic-rational-reconstruction.h>
//#include <linbox/algorithms/fast-rational-reconstruction.h>

#include <linbox/linbox-config.h>
#include <linbox/util/debug.h>

//#define DEBUG_RR
//#define DEBUG_RR_BOUNDACCURACY
#define DEF_THRESH 50

#ifdef __LINBOX_HAVE_NTL
#include <NTL/LLL.h>
#endif


//#define __LINBOX_HAVE_FPLLL
//!@bug BB: demander FPLLL comme dépendance optionnelle...
#ifdef __LINBOX_HAVE_FPLLL
extern "C" {
#include </home/pgiorgi/Library/fplll-1.1/myheuristic.h>
#include </home/pgiorgi/Library/fplll-1.1/myproved.h>
}
#include <linbox/algorithms/short-vector.h>
#endif

namespace LinBox {	

	/// \brief Limited doc so far.  Used, for instance, after LiftingContainer.
template< class _LiftingContainer,
	  class RatRecon = RReconstruction<typename _LiftingContainer::Ring, ClassicMaxQRationalReconstruction<typename _LiftingContainer::Ring> > >
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
	RatRecon RR;

	/** \brief Constructor 
	 *  maybe use different ring than the ring in lcontainer
	 */
	RationalReconstruction (const LiftingContainer& lcontainer, const Ring& r = Ring(), int THRESHOLD =DEF_THRESH) : 
		_lcontainer(lcontainer), _r(r), _threshold(THRESHOLD), RR(_r) {
			
		//if ( THRESHOLD < DEF_THRESH) _threshold = DEF_THRESH;
	}

	/** \brief Get the LiftingContainer
	 */
	const LiftingContainer& getContainer() const {
		return _lcontainer;
	}
		

	/** \brief  Handler to switch between different rational reconstruction strategy.
	 *  Allow  early termination and direct fast method
	 *  Switch is made by using a threshold as the third argument 
	 *  (default is set to that of constructor THRESHOLD
	 * 0    -> direct method
	 * > 0  -> early termination with 
	 */ 
	template <class Vector>
	bool getRational(Vector& num, Integer& den, int switcher) const { 
		if ( switcher == 0)
			return getRational3 (num, den);
			//{getRational1(num,den); print (num); std::cout << "Denominator: " << den << "\n";
			//getRational3(num, den);print (num); std::cout << "Denominator: " << den << "\n";}
			
		else
			//return getRational6(num,den, switcher);
			return getRational1 (num, den);
			//{getRational1(num,den); print (num); std::cout << "Denominator: " << den << "\n";
			//getRational3(num, den);print (num); std::cout << "Denominator: " << den << "\n";}
		return 1;
	}

	template <class Vector>
	bool getRational(Vector& num, Integer& den) const { 
		if ( _threshold == 0)
			return getRational3 (num, den);
			//{getRational1(num,den); print (num); std::cout << "Denominator: " << den << "\n";
			//getRational3(num, den);print (num); std::cout << "Denominator: "  << den << "\n";}
		else
			return getRational1 (num, den);
			//{getRational1(num,den); print (num); std::cout << "Denominator: " << den << "\n";
			//getRational3(num, den);print (num); std::cout << "Denominator: "  << den << "\n"; }

		return 1;
	}
	
	template <class InVect1, class InVect2>
	Integer& dot (Integer& d, const InVect1& v1, const InVect2& v2) const {
		typename InVect1::const_iterator v1_p;
		typename InVect2::const_iterator v2_p;
		_r. init (d, 0);
		for (v1_p = v1. begin(), v2_p = v2. begin(); v1_p != v1. end(); ++ v1_p, ++ v2_p)
			_r. axpyin (d, *v1_p, *v2_p);

		return d;
	}
	/** \brief Reconstruct a vector of rational numbers
	 *  from p-adic digit vector sequence.
	 *  An early termination technique is used.
	 *  Answer is a pair (numerator, common denominator)
	 *  The trick to reconstruct the raitonal solution (V. Pan) is implemented.
	 *  Implement the certificate idea, preprint submitted to ISSAC'05
	 */
	
	/*
	template <class Vector>
	void print (const Vector& v) const {
		typename Vector::const_iterator v_p;
		std::cout << "[";
		for (v_p = v. begin(); v_p != v. end(); ++ v_p)
			std::cout << *v_p << ", ";
		std::cout << "]\n";
	}
	*/
	template<class Vector>
	bool getRational1(Vector& num, Integer& den) const { 
 		
#ifdef RSTIMING
		ttRecon.clear();
		tRecon.start();
#endif			
		linbox_check(num. size() == (size_t)_lcontainer.size());
		typedef Vector IVector;
		typedef std::vector<IVector> LVector;
		int n = num. size();
		int len = _lcontainer. length();
		Integer prime = _lcontainer.prime();//prime
		LVector digits; //Store all p-adic digits
		digits. resize (len); //reserve space for all digits
		Integer modulus; //store current modulus
		Integer denbound; // store current bound for den
		Integer numbound; //store current bound for num

		_r. init (modulus, 1);
		_r. init (denbound, 1);
		_r. init (numbound, 1);
		Integer c1, c2, c1_den, c1_num, c2_den, c2_num;
		IVector r1(num.size()), r2(num.size());
		_r. init (c1, 0); _r. init(c1_den, 1); _r. init (c1_num, 0);
		_r. init (c2, 0); _r. init(c2_den, 1); _r. init (c2_num, 0);
	
		typename IVector::iterator r_p;
		for (r_p = r1. begin(); r_p != r1. end(); ++ r_p)
			_r. init (*r_p, rand());
		for (r_p = r2. begin(); r_p != r2. end(); ++ r_p)
			_r. init (*r_p, rand());

		//std::cout << "Random vecotor1: " ;
		//print (r1);
		//std::cout << "Random vecotor2: " ;
		//print (r2); 

		Integer pmodulus; //store previous modulus
		Integer tmp; //temprary integer
			
		_r. init (den, 1);
		typename LiftingContainer::const_iterator iter = _lcontainer.begin(); 

		Integer tmp_den, tmp_num, rem1, rem2;
		// Do until it meets early termination conditions.
		int step = 0;
		//std::cout << "length:= " << len << '\n';
		//std::cout << "threshold is: "<< _threshold<<std::endl;
		typename LVector::iterator digits_p = digits. begin();

		
#ifdef RSTIMING
		tRecon.stop();
		ttRecon+=tRecon;
#endif
		while (step < len) {
			
			//std::cout << "In " << step << "th step:\n";
			
			IVector& dig = *digits_p;
			++step; ++ digits_p;

			dig. resize (n);

			// get next p-adic digit
			bool nextResult = iter.next(dig);
			if (!nextResult) {
				std::cout << "ERROR in lifting container. Are you using <double> ring with large norm?" << std::endl;
				return false;
			}
			//std::cout << "New digits:\n";
			//print (dig); 

			// preserve the old modulus
			_r.assign (pmodulus, modulus);
			// upate _modulus *= _prime
			_r.mulin (modulus,  prime);
			// update den and num bound
			if (( step % 2) == 0) {
				_r. mulin (denbound, prime);
				_r. assign (numbound, denbound);
			}
			
			if ((step % _threshold) == 0) {
				
				//std::cout << "Previous (Current) modulus: " << pmodulus << "( " << modulus << ")\n";
				dot (tmp, r1, dig); _r. remin (tmp, prime); _r. axpyin (c1, tmp, pmodulus);
				//std::cout << "r1 * digit: " << tmp << '\n';
				dot (tmp, r2, dig); _r. remin (tmp, prime); _r. axpyin (c2, tmp, pmodulus);
				//std::cout << "r2 * digit: " << tmp << '\n';
				//std::cout << "c1, c2: " << c1 << ", " << c2 << "\n";
				
				_r. mul (rem1, c1, c1_den); _r. subin (rem1, c1_num); _r. remin (rem1, modulus);
				_r. mul (rem2, c2, c2_den); _r. subin (rem2, c2_num); _r. remin (rem2, modulus);
				
				//Early termination condition is met.
				
				if(_r. isZero (rem1) && _r. isZero (rem2)) {
					//std::cout << "Early termination happens:\n";					
					break;
				}
				
				if (!_r. isZero (rem1)) {
					int status;
					status = _r.reconstructRational(tmp_num, tmp_den, c1, modulus, numbound, denbound);
					if(status) {
						_r. assign (c1_den, tmp_den); _r. assign (c1_num, tmp_num);
					}					
				}
				
				if (!_r. isZero (rem2)) {
					int status;
					status = _r.reconstructRational(tmp_num, tmp_den, c2, modulus, numbound, denbound);
					if(status) {
						_r. assign (c2_den, tmp_den); _r. assign (c2_num, tmp_num);
					}
				}
			}
		}
		IVector res (n);
		typename LVector::const_iterator digit_begin = digits. begin();
		PolEval (res, digit_begin, (size_t)step, prime);
		if(step < len) _r. lcm (den, c1_den, c2_den);
		else {
			_r. sqrt(denbound, modulus);
			_r. assign (numbound, denbound);
		}
		
		//std::cout << "Numbound (Denbound): " << numbound << ", " << denbound << '\n';
		//std::cout << "Answer mod(" << modulus << "): ";// print (res);

#ifdef RSTIMING	
		tRecon.start();
#endif	
		std::cout << "Start rational reconstruction:\n";
		typename Vector::iterator num_p; typename IVector::iterator res_p;
		Integer tmp_res, neg_res, abs_neg, l, g;
		_r. init (den, 1);
		int counter=0;
		for (num_p = num. begin(), res_p = res. begin(); num_p != num. end(); ++ num_p, ++ res_p) {
			_r. mul (tmp_res, *res_p, den);
			_r. remin (tmp_res, modulus);
			_r. sub (neg_res, tmp_res, modulus);
			_r. abs (abs_neg, neg_res);

			if (_r. compare(tmp_res, numbound) < 0) 
				_r. assign (*num_p, tmp_res);
			else if (_r. compare(abs_neg, numbound) < 0)
				_r. assign (*num_p, neg_res);
			else {
				int status;
				status = _r. reconstructRational(tmp_num, tmp_den, *res_p, modulus, numbound, denbound);
				if (!status) {
					std::cout << "ERROR in reconstruction ?\n" << std::endl;
#ifdef DEBUG_RR
					std::cout<<" try to reconstruct :\n";
					//	std::cout<<"approximation: "<<*iter_approx<<std::endl;
					std::cout<<"modulus: "<<modulus<<std::endl;
					std::cout<<"numbound: "<<numbound<<std::endl;
					std::cout<<"denbound: "<<denbound<<std::endl;
#endif
					return false;
				}
				counter++;								
				_r. lcm (l, den, tmp_den);
				_r. div (g, l, den);

				if (!_r. isOne (g)) {
					typename Vector::iterator num_p1;
					for (num_p1 = num. begin(); num_p1 != num_p; ++ num_p1)
						_r. mulin (*num_p1, g);
				}

				_r. div (g, l, tmp_den);
				_r. mul(*num_p, g, tmp_num);
				_r. assign (den, l);
			}
		}
		
#ifdef RSTIMING
		tRecon.stop();
		ttRecon+=tRecon;
		_num_rec=counter;
#endif
		return true; //lifted ok
	} // end of getRational1
	
	/** \brief Reconstruct a vector of rational numbers
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
				root(tmp, tmp, 2*len);
				_r.init(numFactor, tmp);

				// inital numbound is numFactor/sqrt(2)
				tmp = (tmp * 29) / 41;
				_r.init(numbound, tmp);
			}
		}
			
#ifdef DEBUG_RR
		std::cout << "nbound, dbound:" << _lcontainer.numbound() << ",  " << _lcontainer.denbound() << std::endl;
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
			std::cout<<"i: "<<i<<std::endl;
#endif
#ifdef RSTIMING
			tRecon.stop();
			ttRecon += tRecon;
#endif
			// get next p-adic digit
			bool nextResult = iter.next(digit);
			if (!nextResult) {
				std::cout << "ERROR in lifting container. Are you using <double> ring with large norm?" << std::endl;
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
			std::cout<<"approximation mod p^"<<i<<" : \n";
			std::cout<<"[";
			for (size_t j=0;j< zz.size( )-1;j++)
				std::cout<<zz[j]<<",";
			std::cout<<zz.back()<<"]\n";
			std::cout<<"digit:\n";				
			for (size_t j=0;j< digit.size();++j)
				std::cout<<digit[j]<<",";
			std::cout<<std::endl;
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
			std::cout << "i, N, D bounds: " << i << ", " << numbound << ", " << denbound << std::endl;
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
		std::cout << "Computed " << i << " digits out of estimated " << len << std::endl;
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



	template <class ConstIterator>
        void PolEval(Vector& y, ConstIterator& Pol, size_t deg, Integer &x) const {
		
		
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

			ConstIterator Pol_high= Pol+deg_low;
			PolEval(y2, Pol_high, deg_high, x2);
			

			for (size_t i=0;i< y.size();++i){
				_r.assign(y[i],y1[i]);
				_r.axpyin(y[i],x1,y2[i]);
				//_r.axpy(y[i],x1,y2[i],y1[i]);			
			}
					
			_r.mul(x,x1,x2);
		}
	} // end of getRational2


	/** \brief Reconstruct a vector of rational numbers
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
			std::cout << "ERROR in lifting container. Are you using <double> ring with large norm?" << std::endl;
			return false;
		}
		
#ifdef RSTIMING
		tRecon.start();
#endif
	
		Timer eval_dac;//, eval_bsgs;
		
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

		eval_dac.start();
		Integer xeval=prime;
		typename std::vector<Vector>::const_iterator poly_digit= digit_approximation.begin();
		PolEval(real_approximation, poly_digit, length, xeval);

		//std::std::cout << "Another way get answer mod(" << modulus << "): "; print(real_approximation);

		eval_dac.stop();
		
// 		integer modulus_size;
// 		_r.convert(modulus_size,modulus);
// 		std::cout<<"number of bit : "<< modulus_size.bitsize()<<std::endl;
// 		std::cout<<"length        : "<< length<<std::endl;
// 		std::cout<<"prime         : "<< prime<<std::endl;
//		std::cout<<"evaluation divide&conquer  : "<<eval_dac<<std::endl;
// 		std::cout<<"evaluation baby/giant step : "<<eval_bsgs<<std::endl;
// 		std::cout<<"evaluation horner method   : "<<eval_horner<<std::endl;


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
// 					std::cout << "ERROR in reconstruction ?\n" << std::endl;
// 				}
// 			}
// 			dumb_ratrecon.stop();
// 			std::cout<<"full rational reconstruction : "<<dumb_ratrecon.usertime()<<std::endl;
// 		}
		
		
		/*
		 * Rational Reconstruction of each coefficient according to a common denominator
		 */
		
		Timer ratrecon;
		ratrecon.start();
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
		Integer neg_approx, abs_approx;
		int idx_last_den=0;

		for (size_t i=0; iter_approx != real_approximation.end(); ++iter_approx, ++ iter_num, ++iter_denom, ++i){
			//_r.mulin( *iter_approx , common_den_mod_prod);
			_r.mulin( *iter_approx , common_den);
			_r.remin( *iter_approx , modulus);
			_r. sub (neg_approx, *iter_approx, modulus);
			_r. abs (abs_approx, neg_approx);

			if ( _r.compare(*iter_approx, numbound) < 0){
				_r.assign(*iter_num, *iter_approx);
				_r.init(*iter_denom, 1);
			}
			else if (_r.compare(abs_approx, numbound) <0){
				_r.assign(*iter_num, neg_approx);
				_r.init(*iter_denom, 1);
			}
			else {
				if  (!_r.reconstructRational(*iter_num, *iter_denom, *iter_approx, modulus, numbound, denbound))
					{					
						std::cout << "ERROR in reconstruction ?\n" << std::endl;
#ifdef DEBUG_RR
						std::cout<<" try to reconstruct :\n";
						std::cout<<"approximation: "<<*iter_approx<<std::endl;
						std::cout<<"modulus: "<<modulus<<std::endl;
						std::cout<<"numbound: "<<numbound<<std::endl;
						std::cout<<"denbound: "<<denbound<<std::endl;
#endif
						return false;
					}			
			
				_r.mulin(common_den, *iter_denom);
				idx_last_den=(int)i;
				counter++;
				/*
				if (i != size-1){
					//if (! _r.isUnit(*iter_denom)) {counter++;

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
				*/
				
			}
			
		}
		
		_r.init(tmp,1);
		for (int i= idx_last_den ; i>=0;--i){
			_r.mulin(num[i],tmp);
			_r.mulin(tmp,denominator[i]);
		}
			

		/*
		typename Vector1::reverse_iterator rev_iter_num   = num.rbegin();
		typename Vector::reverse_iterator  rev_iter_denom = denominator.rbegin();
		_r.init(tmp,1);
		for (; rev_iter_num != num.rend(); ++rev_iter_num, ++rev_iter_denom){
		
			_r.mulin(*rev_iter_num,tmp);
			_r.mulin(tmp, *rev_iter_denom);
		}
		*/
		den = common_den;
		
		ratrecon.stop();
		//std::cout<<"partial rational reconstruction : "<<ratrecon.usertime()<<std::endl;
#ifdef RSTIMING
		tRecon.stop();
		ttRecon += tRecon;
		_num_rec=counter;
#endif
		
		return true;

	} // end of getRational3

/*
 * early terminated analog of getRational3
 */ 

	template<class Vector1>
        bool getRationalET(Vector1& num, Integer& den, const Integer& den_app =1) const {
                //cout << "ET p ading lifting using ClassicMaxQRationalReconstruction by default or given RReconstruction\n";
#ifdef RSTIMING
		ttRecon.clear();
		tRecon.start();
		_num_rec = 0;
#endif
		
		linbox_check(num.size() == (size_t)_lcontainer.size());
	
		Integer init_den = den_app;
		if (den > 0) lcm(init_den,den,den_app);
			_r. init (den, 1);
#ifdef DEBUG_RR
		cout << "debug: den " << den;
#endif
		Integer prime = _lcontainer.prime();                    // prime used for lifting
		Vector digit(_lcontainer.size());                       // to store next digit
		Integer modulus;                                        // store modulus (power of prime)
		Integer prev_modulus;                                   // store previous modulus
		_r.init(modulus, 0);
		std::vector<Integer> zz(_lcontainer.size(), modulus);   // stores each truncated p-adic approximation
		_r.init(modulus, 1);
		
		size_t len = _lcontainer.length(); // should be ceil(log(2*numbound*denbound)/log(prime))
#ifdef DEBUG_RR
		std::cout << "nbound, dbound:" << _lcontainer.numbound() << ",  " << _lcontainer.denbound() << std::endl;
#endif
		
		typename Vector::iterator num_p;
		typename std::vector<Integer>::iterator zz_p;
		typename Vector::iterator digit_p;
		
		size_t i = 0;
		typename LiftingContainer::const_iterator iter = _lcontainer.begin();	
		
		bool gotAll = false; //set to true if all values are reconstructed on a particular step
		bool terminated = false; // set to true if same values are reconstructed and confirmed (reconstructed twice)
#ifdef RSTIMING		
		long counter=0;//counts number of RR
#endif	
		// do until getting all answera
		while ((i < len) && (!terminated)) {
			++ i;		
#ifdef DEBUG_RR		
			std::cout<<"i: "<<i<<std::endl;
#endif
#ifdef RSTIMING
			tRecon.stop();
			ttRecon += tRecon;
			_num_rec +=counter;
			counter = 0;
#endif
			// get next p-adic digit
			bool nextResult = iter.next(digit);
			if (!nextResult) {
				std::cout << "ERROR in lifting container. Are you using <double> ring with large norm?" << std::endl;
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
			std::cout<<"approximation mod p^"<<i<<" : \n";
			std::cout<<"[";
			for (size_t j=0;j< zz.size( )-1;j++)
				std::cout<<zz[j]<<",";
			std::cout<<zz.back()<<"]\n";
			std::cout<<"digit:\n";				
			for (size_t j=0;j< digit.size();++j)
				std::cout<<digit[j]<<",";
			std::cout<<std::endl;
#endif			
			if ((!gotAll) && (i % _threshold != 0) && (i + _threshold < len)) continue;
			if ((!gotAll) && (!RR.scheduled(i-1)) && (i + _threshold < len)) continue;

			if (gotAll) {
				terminated = true;
				for ( zz_p = zz.begin(), num_p = num.begin(); zz_p != zz.end();  ++ zz_p, ++ num_p) {
					Integer a = *num_p;   
					Integer bx= *zz_p; _r. mulin (bx, den); _r. remin(bx, modulus);
					Integer _bx = bx-modulus;
					if (!_r.areEqual(a,bx) && !_r.areEqual(a,_bx)) {
						terminated = false; 
						break;
					}
				}
				if (terminated) break;
			}

			gotAll = true;
			_r. assign (den, init_den);
			for ( zz_p = zz.begin(), num_p = num.begin();
			      zz_p != zz.end();  ++ zz_p, ++ num_p) {

				Integer tmp_den=0;
				Integer zz_p_den (*zz_p);
				_r. mulin (zz_p_den,den);
				_r. remin (zz_p_den,modulus);
				bool tmp = RR.reconstructRational(*num_p, tmp_den, zz_p_den, modulus); 
#ifdef RSTIMING
++counter;
#endif
				if (tmp) {
					linbox_check (!_r.isZero(tmp_den));
					if (! _r. isOne (tmp_den)) {
						_r. mulin (den, tmp_den);
						for (typename Vector::iterator tmp_p = num. begin (); tmp_p != num_p; ++ tmp_p)
							_r. mulin (*tmp_p, tmp_den);	
					}
				}
				else {
					gotAll = false;
					break;
				}	
			}
		}	

		// if last iteration - reconstruct the result using num and den bounds
		if (i == len) {
			den = 1;
			for ( zz_p = zz.begin(), num_p = num.begin();
			      zz_p != zz.end(); ++ zz_p, ++ num_p) {
			
				Integer tmp_den;
                                Integer zz_p_den (*zz_p);
                                _r. mulin (zz_p_den,den);
                                _r. remin (zz_p_den,modulus);

				bool tmp = RR.reconstructRational(*num_p, tmp_den, zz_p_den, modulus, _lcontainer.numbound(), _lcontainer.denbound()); 
#ifdef RSTIMING
++counter;
#endif
				if (tmp) {
					linbox_check (!_r.isZero(tmp_den));
					if (! _r. isOne (tmp_den)) {
						_r. mulin (den, tmp_den);
						for (typename Vector::iterator tmp_p = num. begin (); tmp_p != num_p; ++ tmp_p)
							_r. mulin (*tmp_p, tmp_den);
					}
				}
				else {
					cout << "ERROR in reconstruction ?\n" << flush;
				}
		
			}
		} 		
#ifdef RSTIMING
		tRecon.stop();
		ttRecon += tRecon;
		_num_rec +=counter;
		counter = 0;
#endif	
#ifdef DEBUG_RR_BOUNDACCURACY
		//std::cout << "Computed " << i << " digits out of estimated " << len << std::endl;
#endif
		Integer g;
		Integer abs_num_p(0);
		_r. init (g, 0);
		_r. gcdin (g, den);
		for (num_p = num. begin(); num_p != num. end(); ++ num_p) { 
			_r. gcdin (g, *num_p);
		}	

		if (!_r. isOne (g) && !_r. isZero(g)) {
			for (num_p = num. begin(); num_p != num. end(); ++ num_p)
				_r. divin (*num_p, g);
			_r. divin (den, g);
		}
		//std::cerr << "Computed num, den of size " << sizeN << ", " << sizeD << "\n By " << i << " digits out of estimated " << len << std::endl;
		return true; //lifted ok, assuming size was correct
 		
	} // end of getRationalET


#ifdef __LINBOX_HAVE_NTL
	/* 
	 * Rational reconstruction using Lattice base reduction	 
	 */
	template<class Vector1>
	bool getRational4(Vector1& num, Integer& den, size_t thresh) const { 

#ifdef RSTIMING
		ttRecon.clear();
		tRecon.start();
#endif
	
		linbox_check(num.size() == (size_t)_lcontainer.size());
		
		// prime 
		Integer prime = _lcontainer.prime();
		
		// length of whole approximation
		size_t length=_lcontainer.length();
		
		// size of the solution
		size_t size= _lcontainer.size();

		// numerator  upper bound
		Integer numbound;
		_r.assign(numbound,_lcontainer.numbound());
		 
		// parameter used for the lattice dimension
		size_t k = (5> size)? size:5 ;

		// number of padic steps to perform first (use of LinBox::integer)
		LinBox::integer N,D, bound, mod;
		_r.convert(N, _lcontainer.numbound());
		_r.convert(D, _lcontainer.denbound());
		_r.convert(mod, prime);
		::root (D, D, k); D+=1;
		bound=2*N*D;
		std::cout<<"size in bit of the bound : "<<bound.bitsize()<<std::endl;
		size_t minsteps = logp(bound, mod)+1;

		Timer magn;
		magn.start();
		// magnitude of A and b
		integer maxValue=0,value, MagnA, Magnb;
		typename LiftingContainer::IMatrix::ConstRawIterator it = _lcontainer.getMatrix().rawBegin();
		for (; it != _lcontainer.getMatrix().rawEnd(); ++it) {
			_r.convert(value,*it);
			if (value<0) value=-value;
			if (value> maxValue)
				maxValue= value;						    
		}
		MagnA=maxValue;

		maxValue=0;
		typename LiftingContainer::IVector::const_iterator it_b = _lcontainer.getVector().begin();
		for (;it_b!= _lcontainer.getVector().end();++it_b){
			_r.convert(value,*it_b);
			if (value<0) value=-value;
			if (value> maxValue)
				maxValue= value;	
		}
		Magnb=maxValue;
		magn.stop();
		std::cout<<"magnitude time:                 "<<magn<<"\n";
		
		// some constants
		Integer zero;
		_r.init(zero,0);
		Vector zero_digit(_lcontainer.size(),zero);	

		// store approximation as a polynomial and evaluate by baby step giant step
		std::vector<Vector>  digit_approximation(length,zero_digit); 

		// store real approximation
		Vector real_approximation(size,zero); 
		Vector last_real_approximation;

		// store modulus (intially set to 1)
		Integer modulus, last_modulus;
		_r.init(modulus, 1);
		
		typename LiftingContainer::const_iterator iter = _lcontainer.begin();

		
		size_t moresteps    = thresh;
		size_t startingsteps =0;
		size_t endingsteps  =minsteps;

		// switchers
		bool latticeOK=false, numeratorOK=false, domoresteps=true, domorelattice=true;
			
		// common denominator
		Integer common_denom;
		_r.init(common_denom,1);

		// bad numerator index
		size_t bad_num_index=0;		


		bool neg_denom=false;

#ifdef RSTIMING
		tRecon.stop();
		ttRecon += tRecon;	
		std::cout<<"\ninitialization time :           "<<tRecon<<"\n";
#endif	
		



		do {// main loop
			
			// keep track on the last power of the approximation
			_r.assign(last_modulus,modulus);
			
			if (domoresteps){
			
				// compute the padic digits
				for (size_t i = startingsteps ;  (i< endingsteps) && (iter.next(digit_approximation[i]));++i) {
					_r.mulin(modulus,prime); 
				}
					
			
#ifdef RSTIMING
				tRecon.clear();
				tRecon.start();
#endif
				// evaluate the padic digit into an integer approximation	
				Integer xeval=prime;
				typename std::vector<Vector>::const_iterator poly_digit= digit_approximation.begin()+startingsteps;
				PolEval(real_approximation, poly_digit, endingsteps - startingsteps, xeval);				
								
				if (startingsteps != 0){
					for (size_t i=0;i<size;++i){
						_r.axpyin(last_real_approximation[i],real_approximation[i], last_modulus); 
					}
					real_approximation=last_real_approximation;
				}
				else 
					last_real_approximation = real_approximation;	
#ifdef RSTIMING
		tRecon.stop();
		ttRecon += tRecon;		
		std::cout<<"evaluation time :               "<<tRecon<<"\n";
#endif	
			}
#ifdef RSTIMING
			tRecon.clear();
			tRecon.start();
#endif	
		
			
			// construct the lattice
			NTL::mat_ZZ Lattice;
			NTL::ZZ m, tmp, det;
			integer tmp_int;
			_r.convert(mod, modulus);
			m=NTL::to_ZZ((std::string(mod)).c_str());
			

			Lattice.SetDims(k+1, k+1);
			NTL::clear(Lattice);
			Lattice[0][0]=1;
			for (size_t i= bad_num_index+1;i< bad_num_index+k+1;++i){
				Lattice[i][i]=m;//not working when bad index <> 0
				_r.convert(tmp_int, real_approximation[i-1]);
				tmp=NTL::to_ZZ((std::string(tmp_int)).c_str());
				Lattice[0][i]=tmp;//not working when bad index <> 0
			}
		
			// ratio to check the validity of the denominator compare to the entries in the reduced lattice
			NTL::ZZ ratio;
			ratio=NTL::to_ZZ(100L);
			
			// reduce the lattice using LLL algorithm
			Timer chrono;
			chrono.start();
			//NTL::LLL(det, Lattice);
			NTL::LLL_XD(Lattice);
			chrono.stop();
			std::cout<<"lattice reduction time :        "<<chrono<<std::endl;

		
			// check if the 1st row is the short vector
			latticeOK=true;
			tmp=abs(Lattice[0][0])*ratio;	
			for (size_t i=1;i<k+1;++i){
				for (size_t j=0;j<k+1;++j)
					if (tmp > abs(Lattice[i][j])){
						latticeOK=false;
						break;
					}
			}
		
			
			

			if (latticeOK) {// lattice ok
				Timer  checknum;
				checknum.start();
				bool neg=false;
				// get the denominator from the lattice
				tmp =Lattice[0][0];
				if (sign(tmp) <0) neg=true;
				long b = NumBytes(tmp);
				unsigned char* byteArray;
				byteArray = new unsigned char[(size_t)b ];
				BytesFromZZ(byteArray, tmp, b);
				integer base(256);
				integer dd= integer(0);	    
				for(long i = b - 1; i >= 0; --i) {
					dd *= base;
					dd += integer(byteArray[i]);
				}
				delete [] byteArray;
				Integer denom;
				_r.init(denom,dd);
				if (neg) _r.negin(denom);							

				neg_denom= neg_denom^neg;

				// compute the lcm of the denomintator and the last denominator
				_r.lcmin(common_denom, denom);				

				Integer neg_approx, abs_approx;
				
				numeratorOK=true;
				// compute the numerators and check their validity according to the numerator  bound
				for (size_t i=0;i<size;++i){
					_r.mulin(real_approximation[i], denom);
					_r.remin(real_approximation[i], modulus);
					_r. sub (neg_approx, real_approximation[i], modulus);
					_r. abs (abs_approx, neg_approx);

					if ( _r.compare(real_approximation[i], numbound) < 0)
						_r.assign(num[i], real_approximation[i]);									
					else if (_r.compare(abs_approx, numbound) <0)
						_r.assign(num[i], neg_approx);		
					else {
						bad_num_index= std::min(i, size-k);
						numeratorOK=false;
						break;
					}
				}
				checknum.stop();
				std::cout<<"checking numerator time :       "<<checknum<<"\n";

				if (numeratorOK) {//numerator ok					
					Timer checksol;
					checksol.start();
					// compute the magnitude of the numerator
					integer maxnum=0;
					typename Vector::const_iterator it_num=num.begin();
					for (; it_num != num.end(); ++it_num) {
						_r.convert(value,*it_num);
						if (value<0) value=-value;
						if (value> maxnum)
							maxnum= value;						    
					}
				
					// check the validity of the solution according to n.||A||.||num||+ d.||b|| < modulus
					integer check= size*MagnA*maxnum+dd*Magnb;
					
					checksol.stop();
					std::cout<<"checking solution time :        "<<checksol<<"\n\n";
					
					domorelattice=false;

					if (check < mod)
						domoresteps=false;
					else
						domoresteps=true;
			
				}
				else{
					domorelattice=true;
					domoresteps=false;
				}
				
			} 
			else{
				std::cout<<"lattice failed\n";
				domoresteps=true;
				domorelattice=false;
			}
			
			if (domoresteps) std::cout<<"do more steps\n";
			if (domorelattice) std::cout<<"do more lattice\n";

			startingsteps = endingsteps;
			if (domoresteps){
				bad_num_index=0;
				endingsteps+= moresteps;
				if (endingsteps>length)
					endingsteps=length;
			}
#ifdef RSTIMING
			tRecon.stop();
			ttRecon += tRecon;		
#endif
		}
		while (domoresteps||domorelattice);
	
#ifdef RSTIMING
		tRecon.clear();
		tRecon.start();
#endif	
		_r.assign(den, common_denom);
		
		if (neg_denom){
			for (size_t i=0;i<size;++i)
				_r.negin(num[i]);			
		}
#ifdef RSTIMING
		tRecon.stop();
		ttRecon += tRecon;		
#endif	
		return true;

	} // end of getRational4

#endif // end of __LINBOX_HAVE_NTL	





#ifdef __LINBOX_HAVE_FPLLL

	/* 
	 * Rational reconstruction using Lattice base reduction	 
	 */
	template<class Vector1>
	bool getRational5(Vector1& num, Integer& den, size_t thresh) const { 

#ifdef RSTIMING
		ttRecon.clear();
		tRecon.start();
#endif
	
		linbox_check(num.size() == (size_t)_lcontainer.size());
		
		// prime 
		Integer prime = _lcontainer.prime();
		
		// length of whole approximation
		size_t length=_lcontainer.length();
		
		// size of the solution
		size_t size= _lcontainer.size();

		// numerator  upper bound
		Integer numbound;
		_r.assign(numbound,_lcontainer.numbound());
		 
		// parameter used for the lattice dimension
		size_t k = (2> size)? size:2 ;

		// number of padic steps to perform first (use of LinBox::integer)
		LinBox::integer N,D, bound, mod;
		_r.convert(N, _lcontainer.numbound());
		_r.convert(D, _lcontainer.denbound());
		_r.convert(mod, prime);
		::root (D, D, k); D+=1;
		bound=2*N*D;
		std::cout<<"size in bit of the bound : "<<bound.bitsize()<<std::endl;
		size_t minsteps = logp(bound, mod)+1;

		Timer magn;
		magn.start();
		// magnitude of A and b
		integer maxValue=0,value, MagnA, Magnb;
		typename LiftingContainer::IMatrix::ConstRawIterator it = _lcontainer.getMatrix().rawBegin();
		for (; it != _lcontainer.getMatrix().rawEnd(); ++it) {
			_r.convert(value,*it);
			if (value<0) value=-value;
			if (value> maxValue)
				maxValue= value;						    
		}
		MagnA=maxValue;

		maxValue=0;
		typename LiftingContainer::IVector::const_iterator it_b = _lcontainer.getVector().begin();
		for (;it_b!= _lcontainer.getVector().end();++it_b){
			_r.convert(value,*it_b);
			if (value<0) value=-value;
			if (value> maxValue)
				maxValue= value;	
		}
		Magnb=maxValue;
		magn.stop();
		std::cout<<"magnitude time:                 "<<magn<<"\n";
		
		// some constants
		Integer zero;
		_r.init(zero,0);
		Vector zero_digit(_lcontainer.size(),zero);	

		// store approximation as a polynomial and evaluate by baby step giant step
		std::vector<Vector>  digit_approximation(length,zero_digit); 

		// store real approximation
		Vector real_approximation(size,zero); 
		Vector last_real_approximation;

		// store modulus (intially set to 1)
		Integer modulus, last_modulus;
		_r.init(modulus, 1);
		
		typename LiftingContainer::const_iterator iter = _lcontainer.begin();

		
		size_t moresteps    = thresh;
		size_t startingsteps =0;
		size_t endingsteps  =minsteps;

		// switchers
		bool latticeOK=false, numeratorOK=false, domoresteps=true, domorelattice=true;
			
		// common denominator
		Integer common_denom;
		_r.init(common_denom,1);

		// bad numerator index
		size_t bad_num_index=0;		


		bool neg_denom=false;

#ifdef RSTIMING
		tRecon.stop();
		ttRecon += tRecon;	
		std::cout<<"\ninitialization time :           "<<tRecon<<"\n";
#endif	
		



		do {// main loop
			
			// keep track on the last power of the approximation
			_r.assign(last_modulus,modulus);
			
			if (domoresteps){
			
				linbox_check(startingsteps != endingsteps);
				
				// compute the padic digits
				for (size_t i = startingsteps ;  (i< endingsteps) && (iter.next(digit_approximation[i]));++i) {
					_r.mulin(modulus,prime); 
				}
					
			
#ifdef RSTIMING
				tRecon.clear();
				tRecon.start();
#endif
				// evaluate the padic digit into an integer approximation	
				Integer xeval=prime;
				typename std::vector<Vector>::const_iterator poly_digit= digit_approximation.begin()+startingsteps;
				PolEval(real_approximation, poly_digit, endingsteps - startingsteps, xeval);				
								
				if (startingsteps != 0){
					for (size_t i=0;i<size;++i){
						_r.axpyin(last_real_approximation[i],real_approximation[i], last_modulus); 
					}
					real_approximation=last_real_approximation;
				}
				else 
					last_real_approximation = real_approximation;	
#ifdef RSTIMING
		tRecon.stop();
		ttRecon += tRecon;		
		std::cout<<"evaluation time :               "<<tRecon<<"\n";
#endif	
			}
#ifdef RSTIMING
			tRecon.clear();
			tRecon.start();
#endif	
		
			
			// construct the lattice
			mpz_t **Lattice;
			Lattice= new mpz_t*[k+2];
			for (size_t i=0;i<k+2;++i){
				Lattice[i]= new mpz_t[k+2];
				for (size_t j=0;j<k+2;++j)
					mpz_init(Lattice[i][j]);
			}
			
			integer tmp=1;
			_r.convert(mod, modulus);
					       
			mpz_set(Lattice[1][1], tmp.get_mpz());
			for (size_t i=2;i< k+2;++i){
				mpz_set(Lattice[i][i],mod.get_mpz());
				_r.convert(tmp, real_approximation[bad_num_index+i-2]);
				mpz_set(Lattice[1][i],tmp.get_mpz());
			}
		
			// ratio to check the validity of the denominator compare to the entries in the reduced lattice
			integer ratio;
			ratio=100L;
			
			// reduce the lattice using LLL algorithm
			Timer chrono;
			chrono.start();
			myLLLproved(Lattice, k+1,k+1);
			chrono.stop();
			std::cout<<"lattice reduction time :        "<<chrono<<std::endl;

		
			// check if the 1st row is the short vector
			latticeOK=true;
			mpz_mul(tmp.get_mpz(), Lattice[1][1],ratio.get_mpz());	
			for (size_t i=2;i<k+2;++i){
				for (size_t j=1;j<k+2;++j)
					if (mpz_cmpabs(tmp.get_mpz() , Lattice[i][j])> 0){
						latticeOK=false;
						break;
					}
			}
		
			integer dd;
			// get the denominator from the lattice
			mpz_set(dd.get_mpz(),Lattice[1][1]);
			
			//delete the lattice
			for (size_t i=0;i<k+2;++i){
				for (size_t j=0;j<k+2;++j)
					mpz_clear(Lattice[i][j]);
				delete[] Lattice[i];
			}
			delete[] Lattice;


			if (latticeOK) {// lattice ok
				Timer  checknum;
				checknum.start();
			

				Integer denom;
				_r.init(denom,tmp);
				
				bool neg=true;
				if (dd < 0){
					neg=true;std::cout<<"negative det\n";}
				
				neg_denom= neg_denom^neg;

				// compute the lcm of the denomintator and the last denominator
				_r.lcmin(common_denom, denom);				

				Integer neg_approx, abs_approx;
				
				numeratorOK=true;
				// compute the numerators and check their validity according to the numerator  bound
				for (size_t i=0;i<size;++i){
					_r.mulin(real_approximation[i], common_denom);
					_r.remin(real_approximation[i], modulus);
					_r. sub (neg_approx, real_approximation[i], modulus);
					_r. abs (abs_approx, neg_approx);

					if ( _r.compare(real_approximation[i], numbound) < 0)
						_r.assign(num[i], real_approximation[i]);									
					else if (_r.compare(abs_approx, numbound) <0)
						_r.assign(num[i], neg_approx);		
					else {
						bad_num_index= std::min(i, size-k);
						numeratorOK=false;
						break;
					}
				}
				checknum.stop();
				std::cout<<"checking numerator time :       "<<checknum<<"\n";

				if (numeratorOK) {//numerator ok					
					Timer checksol;
					checksol.start();
					// compute the magnitude of the numerator
					integer maxnum=0;
					typename Vector::const_iterator it_num=num.begin();
					for (; it_num != num.end(); ++it_num) {
						_r.convert(value,*it_num);
						if (value<0) value=-value;
						if (value> maxnum)
							maxnum= value;						    
					}
				
					// check the validity of the solution according to n.||A||.||num||+ d.||b|| < modulus
					integer check= size*MagnA*maxnum+dd*Magnb;
					
					checksol.stop();
					std::cout<<"checking solution time :        "<<checksol<<"\n\n";
					
					domorelattice=false;

					if (check < mod)
						domoresteps=false;
					else
						domoresteps=true;
			
				}
				else{
					domorelattice=true;
					domoresteps=false;
				}
				
			} 
			else{
				std::cout<<"lattice failed\n";
				domoresteps=true;
				domorelattice=false;
			}
			
			if (domoresteps) std::cout<<"do more steps\n";
			if (domorelattice) std::cout<<"do more lattice\n";

			startingsteps = endingsteps;
			if (domoresteps){
				bad_num_index=0;
				endingsteps+= moresteps;
				if (endingsteps>length)
					endingsteps=length;
			}
#ifdef RSTIMING
			tRecon.stop();
			ttRecon += tRecon;		
#endif
		}
		while (domoresteps||domorelattice);
	
#ifdef RSTIMING
		tRecon.clear();
		tRecon.start();
#endif	
		_r.assign(den, common_denom);
		
		if (neg_denom){
			for (size_t i=0;i<size;++i)
				_r.negin(num[i]);			
		}
#ifdef RSTIMING
		tRecon.stop();
		ttRecon += tRecon;		
#endif	
		return true;

	} // end of getRational5


	/* 
	 * Rational reconstruction using Lattice base reduction	 
	 */
	template<class Vector1>
	bool getRational6(Vector1& num, Integer& den, size_t thresh) const { 

#ifdef RSTIMING
		ttRecon.clear();
		tRecon.start();
#endif
	
		linbox_check(num.size() == (size_t)_lcontainer.size());
		
		// prime 
		Integer prime = _lcontainer.prime();
		
		// length of whole approximation
		size_t length=_lcontainer.length();
		
		// size of the solution
		size_t size= _lcontainer.size();

		// numerator  upper bound
		Integer numbound;
		_r.assign(numbound,_lcontainer.numbound());
		 
		// parameter used for the lattice dimension
		size_t k = (2> size)? size:2 ;

		// number of padic steps to perform first (use of LinBox::integer)
		LinBox::integer N,D, bound, mod;
		_r.convert(N, _lcontainer.numbound());
		_r.convert(D, _lcontainer.denbound());
		_r.convert(mod, prime);
		::root (D, D, k); D+=1;
		bound=2*N*D;
		std::cout<<"size in bit of the bound : "<<bound.bitsize()<<std::endl;
		size_t minsteps = logp(bound, mod)+1;

		Timer magn;
		magn.start();
		// magnitude of A and b
		integer maxValue=0,value, MagnA, Magnb;
		typename LiftingContainer::IMatrix::ConstRawIterator it = _lcontainer.getMatrix().rawBegin();
		for (; it != _lcontainer.getMatrix().rawEnd(); ++it) {
			_r.convert(value,*it);
			if (value<0) value=-value;
			if (value> maxValue)
				maxValue= value;						    
		}
		MagnA=maxValue;

		maxValue=0;
		typename LiftingContainer::IVector::const_iterator it_b = _lcontainer.getVector().begin();
		for (;it_b!= _lcontainer.getVector().end();++it_b){
			_r.convert(value,*it_b);
			if (value<0) value=-value;
			if (value> maxValue)
				maxValue= value;	
		}
		Magnb=maxValue;
		magn.stop();
		std::cout<<"magnitude time:                 "<<magn<<"\n";
		
		// some constants
		Integer zero;
		_r.init(zero,0);
		Vector zero_digit(_lcontainer.size(),zero);	

		// store approximation as a polynomial and evaluate by baby step giant step
		std::vector<Vector>  digit_approximation(length,zero_digit); 

		// store real approximation
		Vector real_approximation(size,zero); 
		Vector last_real_approximation;

		// store modulus (intially set to 1)
		Integer modulus, last_modulus;
		_r.init(modulus, 1);
		
		typename LiftingContainer::const_iterator iter = _lcontainer.begin();

		
		size_t moresteps    = thresh;
		size_t startingsteps =0;
		size_t endingsteps  =minsteps;

		// switchers
		bool latticeOK=false, numeratorOK=false, domoresteps=true, domorelattice=true;
 			
		// common denominator
		Integer common_denom;
		_r.init(common_denom,1);

		// bad numerator index
		size_t bad_num_index=0;		


		bool neg_denom=false;

#ifdef RSTIMING
		tRecon.stop();
		ttRecon += tRecon;	
		std::cout<<"\ninitialization time :           "<<tRecon<<"\n";
#endif	
		



		do {// main loop
			
			// keep track on the last power of the approximation
			_r.assign(last_modulus,modulus);
			
			if (domoresteps){
			
				linbox_check(startingsteps != endingsteps);
				
				// compute the padic digits
				for (size_t i = startingsteps ;  (i< endingsteps) && (iter.next(digit_approximation[i]));++i) {
					_r.mulin(modulus,prime); 
				}
					
			
#ifdef RSTIMING
				tRecon.clear();
				tRecon.start();
#endif
				// evaluate the padic digit into an integer approximation	
				Integer xeval=prime;
				typename std::vector<Vector>::const_iterator poly_digit= digit_approximation.begin()+startingsteps;
				PolEval(real_approximation, poly_digit, endingsteps - startingsteps, xeval);				
								
				if (startingsteps != 0){
					for (size_t i=0;i<size;++i){
						_r.axpyin(last_real_approximation[i],real_approximation[i], last_modulus); 
					}
					real_approximation=last_real_approximation;
				}
				else 
					last_real_approximation = real_approximation;	
#ifdef RSTIMING
		tRecon.stop();
		ttRecon += tRecon;		
		std::cout<<"evaluation time :               "<<tRecon<<"\n";
#endif	
			}
#ifdef RSTIMING
			tRecon.clear();
			tRecon.start();
#endif	
		
			
			// construct the lattice			
			std::vector<integer> Lattice(9);
			
			integer tmp=1;
			_r.convert(mod, modulus);
					       
			Lattice[0]=tmp;
			for (size_t i=1;i< k+1;++i){
				Lattice[i*(k+1)+i]=mod;				
				_r.convert(tmp, real_approximation[bad_num_index+i-1]);
				Lattice[i]=tmp;
			}
			
			// ratio to check the validity of the denominator compare to the entries in the reduced lattice
			integer ratio;
			ratio=100L;
			
			// reduce the lattice using LLL algorithm
			Timer chrono;
			chrono.start();
			TernaryLattice L3(Lattice);
			L3.reduce();
			chrono.stop();
			std::cout<<"lattice reduction time :        "<<chrono<<std::endl;

		
			// check if the 1st row is the short vector
			latticeOK=true;
			tmp = ::abs(L3[0][0]*ratio);
			for (size_t i=1;i<k+1;++i){
				for (size_t j=0;j<k+1;++j)
					if (tmp > ::abs(L3[i][j])){
						latticeOK=false;
						break;
					}
			}
		
			integer dd;
			// get the denominator from the lattice
			dd= L3[0][0];
			//L3.print();

			if (latticeOK) {// lattice ok
				Timer  checknum;
				checknum.start();
			

				Integer denom;
				_r.init(denom,tmp);
				
				bool neg=true;
				if (dd < 0){
					neg=true;std::cout<<"negative det\n";}
				
				neg_denom= neg_denom^neg;

				// compute the lcm of the denominator and the last denominator
				_r.lcmin(common_denom, denom);				

				Integer neg_approx, abs_approx;
				
				numeratorOK=true;
				// compute the numerators and check their validity according to the numerator  bound
				for (size_t i=0;i<size;++i){
					_r.mulin(real_approximation[i], common_denom);
					_r.remin(real_approximation[i], modulus);
					_r. sub (neg_approx, real_approximation[i], modulus);
					_r. abs (abs_approx, neg_approx);

					if ( _r.compare(real_approximation[i], numbound) < 0)
						_r.assign(num[i], real_approximation[i]);									
					else if (_r.compare(abs_approx, numbound) <0)
						_r.assign(num[i], neg_approx);		
					else {
						bad_num_index= std::min(i, size-k);
						numeratorOK=false;
						break;
					}
				}
				checknum.stop();
				std::cout<<"checking numerator time :       "<<checknum<<"\n";

				if (numeratorOK) {//numerator ok					
					Timer checksol;
					checksol.start();
					// compute the magnitude of the numerator
					integer maxnum=0;
					typename Vector::const_iterator it_num=num.begin();
					for (; it_num != num.end(); ++it_num) {
						_r.convert(value,*it_num);
						if (value<0) value=-value;
						if (value> maxnum)
							maxnum= value;						    
					}
				
					// check the validity of the solution according to n.||A||.||num||+ d.||b|| < modulus
					integer check= size*MagnA*maxnum+dd*Magnb;
					
					checksol.stop();
					std::cout<<"checking solution time :        "<<checksol<<"\n\n";
					
					domorelattice=false;

					if (check < mod)
						domoresteps=false;
					else
						domoresteps=true;
			
				}
				else{
					domorelattice=true;
					domoresteps=false;
				}
				
			} 
			else{
				std::cout<<"lattice failed\n";
				domoresteps=true;
				domorelattice=false;
			}
			
			if (domoresteps) std::cout<<"do more steps\n";
			if (domorelattice) std::cout<<"do more lattice\n";

			startingsteps = endingsteps;
			if (domoresteps){
				bad_num_index=0;
				endingsteps+= moresteps;
				if (endingsteps>length)
					endingsteps=length;
			}
#ifdef RSTIMING
			tRecon.stop();
			ttRecon += tRecon;		
#endif
		}
		while (domoresteps||domorelattice);
	
#ifdef RSTIMING
		tRecon.clear();
		tRecon.start();
#endif	
		_r.assign(den, common_denom);
		
		if (neg_denom){
			for (size_t i=0;i<size;++i)
				_r.negin(num[i]);			
		}
#ifdef RSTIMING
		tRecon.stop();
		ttRecon += tRecon;		
#endif	
		return true;

	} // end of getRational6



#endif // end of __LINBOX_HAVE_FPLLL	






};
	
}

#undef DEF_THRESH
#endif //__LINBOX_reconstruction_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
