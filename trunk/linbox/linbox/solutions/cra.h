/* -*- mode: C++; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*- */
/* author: B. David Saunders and Zhendong Wan*/
#ifndef __CRA_H
#define __CRA_H

#include <vector>
//#include <linbox/util/timer.h>

namespace LinBox {

/* Warning, we won't detect repeat primes */
/* Warning, we won't detect bad primes */

template<class _Integer>
class CRA{
	public:
	typedef _Integer Integer;
	typedef std::vector<Integer> Vector;
	typedef std::list<Vector> LVector;
	typedef std::list<Integer> List;
	typedef std::vector<int> IntVector;

	protected:

	Integer m;	// current modulus
	Integer cert1, cert2;
	std::vector<int> randv1, randv2;

	LVector holdres;
	List holdprime;

	size_t isTerminating;
	bool first_call;

	int k; // step counter

	//UserTimer timer;
	//double cra_time;
	//double mod_time;
	public:

	CRA(): isTerminating(0), first_call(true), k(0) {
		m = 1;
		//cra_time = 0;
		//mod_time = 0;
	}

	// should also allow a bound to be given.
	/* possible set sub related to d. size() */
	template <class V>
	void step(Integer p, const V& d){

		++ k;
		//std::clog << '\r' << k << " " << p;
		/*
		std::cout << k << "th input: " << p << ' ';
		printVect (d);
		std::cout << "endl\n";
		*/

		if (first_call) {

			//timer. start();
			// set up certificate system
			//   fixme: use custom (simpler) certs for very short vectors.  
			first_call = false;
			std::vector<int>::iterator int_p;
			randv1. resize (d. size());
			randv2. resize (d. size());
			for (int_p = randv1. begin(); 
				 int_p != randv1. end(); ++ int_p) 
				*int_p = rand() % 10000;
			for (int_p = randv2. begin(); 
				 int_p != randv2. end(); ++ int_p) 
				*int_p = rand() % 10000;
            cert1 = cert2 = 0;

			// take partial answer
			/*
			holdprime. push_back (p);
			holdres. push_back (Vector());
			Vector& cur = holdres. back();
			cur. resize (d. size());
			std::copy (d. begin(), d. end(), cur. begin());

			dot (cert1, cur, randv1);
			dot (cert2, cur, randv2);
			normalize (cert1, p);
			normalize (cert2, p);
			m = p;
			*/
			//std::cout <<"Random Vector1: ";
			//printVect (randv1);
			//std::cout << std::endl;
			//std::cout <<"Random Vector2: ";
			//printVect (randv2);
			//std::cout << std::endl;
			return ;
		}
		/*
		else if (res.size() != d.size()) {
				std::clog << "cra.step: vector length mismatch" << std::endl; 
				return;
		}
		*/
		// take partial answer
		holdres. push_back (Vector());
		Vector& cur = holdres. back();
		cur. resize (d. size());
		std::copy (d. begin(), d. end(), cur. begin());
		holdprime. push_back (p);

		// make relatively prime 
		Integer g, u, v;
		Integer& cur_p = holdprime. back();
		gcd (g, m, cur_p, u, v);  
		while (g != 1) { // take gcd out of p and d
			cur_p /= g; 
			gcd (g, m, cur_p, u, v);  
		}
		if (cur_p == 1) { // nothing new contributed
			holdres. pop_back (); 
			holdprime. pop_back (); 
			return; 
		}
		else if (cur_p != p) // not a full contribution
			for (typename Vector::iterator i = cur. begin(); 
			                               i != cur. end(); ++i)
				*i %= cur_p;
		// endif

		// check certificates for termination
		u *= m; // u = 0 mod m, u = 1 mod cur_p
		v *= cur_p; // v = 1 mod m, v = 0 mod cur_p
		m *= cur_p;
		Integer tmp1, tmp2;
		dot (tmp1, cur, randv1);
		dot (tmp2, cur, randv2);

		Integer cert1_p, cert2_p;
		cert1_p = (cert1  - tmp1) % cur_p; cert2_p = (cert2 - tmp2) % cur_p;

		/*
		std::cout << "Previous certificates: " << cert1 << ' ' 
				  << cert2 << ' ' << m << '\n';
		std::cout << "Current  certificates: " << tmp1 << ' '
				  << tmp2 << ' ' << p << '\n'; 
		std::cout << "Difference " << cert1_p << ' '
				  << cert2_p << " mod " << p << '\n';
		*/
		
		if (sign(cert1_p)  || sign(cert2_p)) {
		/*
			holdres. push_back (Vector());
			Vector& cur = holdres. back();
			cur. resize (d. size());
			std::copy (d. begin(), d. end(), cur. begin());
			holdprime. push_back (p);
			Integer g, u, v;
			gcd (g, m, p, u, v);
			u *= m; // u = 0 mod m, u = 1 mod p
			v *= p; // v = 1 mod m, v = 0 mod p
			m *= p;
		*/
			cert1 = cert1 * v + tmp1 * u;
			cert2 = cert2 * v + tmp2 * u;
			normalize (cert1, m);
			normalize (cert2, m);
			isTerminating = 0;
		}
		else {
			++isTerminating; 
		}
	}

	int steps() {return k;}

	bool terminated() { return isTerminating > 0;}
	size_t changelessness() { return isTerminating;}

	template<class V>
	void result(V& w) { 
		//timer. stop();
		//mod_time += timer. time();

		//timer. start();
		buildTree (holdprime, holdres);
		normalize (holdres. front(), holdprime. front());
		//timer. stop();
		//cra_time += timer. time();
/*
		if (check (holdres. front(), d, p)) {
			std::cout << "Number of primes needed: " << k << '\n';
			//std::cout << "Modulus prime time: " << mod_time << std::endl;
			//std::cout << "CRA time: " << cra_time << std::endl;
			isTerminated = true;
		}
		else {
			std::clog << "Fake position:\n";
		}
		*/
		std::copy (holdres.front(). begin(), holdres.front().end(), w. begin());
	}
	 
    Integer& modulus() { return m; }

	protected:
	void buildTree (List& holdprime, LVector& holdres) {

		//std::cout << "In building tree:\n";
		//debug();
		typename LVector::iterator half_res;
		typename List::iterator half_p;
		while (holdres. size () > 1) {
			combine (holdprime, holdres);
			half_res = holdres. begin();
			half_p = holdprime. begin();
			for (size_t i = 0; i < ((holdres. size() + 1) / 2); ++ i) {
				++ half_res;
				++ half_p;
			}
			holdres. erase (half_res, holdres. end());
			holdprime. erase (half_p, holdprime. end());
			
		}
		//std::cout << "Answer: " << std::endl;
		//debug();
	}

	void combine (List& holdprime, LVector& holdres){
		//std::cout << "In building subtree:\n";
		//debug();
		typename LVector::iterator first_res, half_res;
		typename List::iterator first_p, half_p;
		first_res = half_res = holdres. begin();
		first_p = half_p = holdprime. begin();
		for ( size_t i = 0; i < ((holdres. size() + 1) / 2); ++ i) {
			++ half_res;
			++ half_p;
		}
		
		Integer g, u, v;
		for (;half_res != holdres. end(); 
			 ++ half_res, ++ first_res, ++ half_p, ++ first_p) {

			Integer& pk = *first_p;
			Integer& pn = *half_p;
			Vector& rk = *first_res;
			Vector& rn = *half_res;
			gcd (g, pk, pn, u, v);
			u *= pk; // u = 0 mod pk, u = 1 mod pn;
			v *= pn; // v = 1 mod pk, v = 0 mod pn;
			pk *= pn;
			typename Vector::iterator rk_p, rn_p;
			for (rk_p = rk. begin(), rn_p = rn. begin();
				 rk_p != rk. end(); ++ rk_p, ++ rn_p) {

				 	*rk_p = (*rk_p) * v + (*rn_p) * u;
					//normailze here?
					// *rk_p %= pk;
			}
		}
		//std::cout << "After calling:\n";
		//debug();
	}

	void normalize (Integer& res, const Integer& m) {
	 	res %= m ;
		Integer tmp;
		if (sign (res) > 0)
			tmp = res - m;
		else
			tmp = res + m;

		if (absCompare (res, tmp) > 0)
			res = tmp;
	}

	void normalize (Vector& res, const Integer& m) {
	 	typename Vector::iterator res_p;
		Integer tmp;
		for (res_p = res. begin(); res_p != res. end(); ++ res_p) 
			normalize (*res_p, m);
	}

	template <class Vect>
	void printVect (const Vect& v) {
		typename Vect::const_iterator v_p;
		std::cout << "[ ";
		for (v_p = v. begin(); v_p != v. end(); ++ v_p) {
			std::cout << *v_p << " ";
		}
		std::cout << "]" << std::endl;
	}

	void debug () {
		std::cout << "Data:\n";
		typename List::iterator prime_p;
		typename LVector::iterator res_p;
		for (prime_p = holdprime. begin(), res_p = holdres. begin();
			 prime_p != holdprime. end(); ++ prime_p, ++ res_p) {
			 	std::cout << "prime: " << *prime_p << std::endl;
			std::cout << "residue: ";
			printVect (*res_p);
		}
	}
	

	template <class Vect1, class Vect2>
	bool check (const Vect1& v1, const Vect2& v2, const Integer& p) {
		
		typename Vect1::const_iterator v1_p;
		typename Vect2::const_iterator v2_p;

		integer tmp;
		for (v1_p = v1. begin(), v2_p = v2. begin();
			 v1_p != v1. end(); ++ v1_p, ++ v2_p) {

			 tmp = (*v1_p - *v2_p) % p;
			 if (sign(tmp)) 
			 	return false;
		}
	
		return true;
	}

	template <class Vect1, class Vect2>
	void dot (Integer& d, const Vect1& v1, const Vect2& v2){
		d = 0;
		typename Vect1::const_iterator v1_p;
		typename Vect2::const_iterator v2_p;
		for (v1_p = v1. begin(), v2_p = v2. begin(); 
			 v1_p != v1. end(); ++ v1_p, ++ v2_p)
			 	
			d += (*v1_p) * (*v2_p);
	}

};

}

#endif
