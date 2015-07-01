/* -*- mode: C++; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*- */
/* author: B. David Saunders and Zhendong Wan*/
#ifndef __CRA_H
#define __CRA_H

#include <vector>
//#include <linbox/util/timer.h>

namespace LinBox {

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
	double lm;
	Integer cert;
	std::vector<int> randv;
	unsigned int EARLY_TERM_THRESHOLD;
	double UPPER_BOUND;

	LVector   holdres;
	List    holdvalue;
	List    holdprime;

	size_t occurency;

	int k; // step counter

            //UserTimer timer;
            //double cra_time;
            //double mod_time;
    public:

	CRA(size_t n=0, unsigned int EARLY=1, const double BOUND=0) 
                : m(1), lm(0), occurency(0), k(0) {
            initialize(n,EARLY, BOUND);
        }
	
	CRA(size_t n, unsigned int EARLY, const integer BOUND) 
                : m(1), lm(0), occurency(0), k(0) {
            initialize(n,EARLY, BOUND);
        }


	void initialize(size_t n, unsigned int EARLY, const integer BOUND) {
            initialize( n, EARLY, log( (double)BOUND ) );
	}

	void initialize(size_t n, unsigned int EARLY=1, const double BOUND=0) {
            k = 0;
            EARLY_TERM_THRESHOLD = EARLY;
            UPPER_BOUND = BOUND;
            if (EARLY_TERM_THRESHOLD > 0) {
                std::vector<int>::iterator int_p;
                randv. resize (n);
                for (int_p = randv. begin(); 
                     int_p != randv. end(); ++ int_p) 
                    *int_p = rand() % 10000;
            }
            cert = 0;
	}

	// Works also if Vect is an Integer
	/** \brief The CRA loop

	Given a function to generate residues mod a single prime, this loop produces the residues 
	resulting from the Chinese remainder process on sufficiently many primes to meet the 
	termination condition.

	\param F - Function object of two arguments, F(r, p), given prime p it outputs residue(s) r.
	This loop may be parallelized.  F must be reentrant, thread safe.
	For example, F may be returning the coefficients of the minimal polynomial of a matrix mod p.
	Warning - we won't detect bad primes.

	\param genprime - RandIter object for generating primes.
	\result res - an integer or object of a class meeting the FixedVector interface (with iterator
	requirement relaxed to forward iterators).  Vectors, SubVectors, STL lists may be used.
	*/
	template<class Vect, class Function, class RandPrime>
	Vect & operator() (Vect& res, const Function& F, RandPrime& genprime) {
		Integer p;
		Vect r;
		while( ! this->terminated() ) {
			genprime.randomPrime(p);
			this->progress( p, F(r, p) );
		}
		return this->result(res);
	}

        
    /** \brief Function for adding a new prime and it's residue to the CRA process.
	\param p - A modulus.  Process is most efficient if it is relatively prime to 
	all other moduli used.
	\param d - A residue, image  mod p of the desired value.
	*/
	void progress (const Integer& p, const Integer& d) {

            ++ k;	
		// make relatively prime 
            Integer g, u, v, dp;
            Integer cur_p = p;
            gcd (g, m, cur_p, u, v);  
            dp = d;
            while (g != 1) { // take gcd out of p and d
                cur_p /= g; 
                gcd (g, m, cur_p, u, v);  
            }		
            if (cur_p == 1) { // nothing new contributed
                return; 
            }
            else if (cur_p != p) // not a full contribution
                dp %= cur_p;
		
            if (EARLY_TERM_THRESHOLD == 0){
                holdprime. push_back (p);
                holdvalue. push_back (dp);			
            }		
		// nothing inside
            else if (m == 1) {
                m = p; cert = dp;
            }	   
            else{			
                    // compute the new result by CRA
                    // new_result = old_result * M * (
                Integer tmp, g, s, t;

                gcd(g, p, m, s, t);
                s *= p; // s == 1 mod m, == 0 mod p
                t *= m; // t == 0 mod m, == 1 mod p;
                tmp = s * cert + t * dp;
                    /*
                      s = (s * sert) % p;
                      t = (t * dp) % m;
                      tmp = t*m + s*p;
                    */

                m *=p;
                normalize (tmp, m);
                if (sign(tmp - cert)) {
                    cert = tmp;
                    occurency = 0;
                }
                else {
                    ++occurency;
                }
            }
	}


            // should also allow a bound to be given.
            /* possible set sub related to d. size() */
    /** \brief Taking a step when the CRA process is being applied ot a vector or list of values..
	\param p - A modulus.  Process is most efficient if it is relatively prime to 
	all other moduli used.
	\param d - A residue sequence: images  mod p of the desired value.  May be a list, vector,
	SubVector.
	*/
	template <class Vect>
	void progress(const Integer& p, const Vect& d){

		//		linbox_check ( d.size() == randv.size());
            ++ k;

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

            if (EARLY_TERM_THRESHOLD > 0) {
                    // check certificates for termination
                u *= m; // u = 0 mod m, u = 1 mod cur_p
                v *= cur_p; // v = 1 mod m, v = 0 mod cur_p
                m *= cur_p;
                Integer tmp;
                dot (tmp, cur, randv);
			
                Integer cert_p;
                cert_p = (cert  - tmp) % cur_p; 
			
			
                    /*
                      std::cout << "Previous certificates: " << cert1 << ' ' 
                      << cert2 << ' ' << m << '\n';
                      std::cout << "Current  certificates: " << tmp1 << ' '
                      << tmp2 << ' ' << p << '\n'; 
                      std::cout << "Difference " << cert1_p << ' '
                      << cert2_p << " mod " << p << '\n';
                    */
		
                if (sign(cert_p)) {
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
                    cert = cert * v + tmp * u;
                    normalize (cert, m);
                    occurency = 0;
                }
                else {
                    ++occurency;
                }
            }
            else{
                lm += log(double(cur_p))*1.442695040;
                m*= cur_p;
            }
					
	}

	int steps() {return k;}

	bool terminated() { 
            return ((EARLY_TERM_THRESHOLD && (occurency > EARLY_TERM_THRESHOLD)) || ((lm > UPPER_BOUND) && (UPPER_BOUND > 0)));
	}

    /** \brief Number of progress steps without change in the combined residue. 
	
	Allows flexibility in deciding early termination.
	(earlier early termination).
	*/
	size_t stableSteps() { return occurency;}


    /** \brief result mod the lcm of the moduli.
	
	A value mod the lcm of the progress step moduli
	which agrees with each residue mod the corresponding modulus. 
	*/
	integer& result (integer &d){
            if (EARLY_TERM_THRESHOLD>0)
                return d = cert; 
            else {
                buildTree (holdprime, holdvalue);
                normalize (holdvalue. front(), holdprime. front());
                return d = cert = holdvalue. front();
            }
	}

    /** \brief results mod the lcm of the moduli.
	
	A sequence of values mod the lcm of the progress step moduli 
	each entry of which agrees with the corresponding (sequence position) residue mod the corresponding 
	(progress step) modulus. 
	*/
	template<class Vect>
	Vect& result (Vect& w) { 
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
            return w;
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

//                 gcd (g, pk, pn, u, v);
//                 u *= pk; // u = 0 mod pk, u = 1 mod pn;
//                 v *= pn; // v = 1 mod pk, v = 0 mod pn;
//                 pk *= pn;
//                 typename Vector::iterator rk_p, rn_p;
//                 for (rk_p = rk. begin(), rn_p = rn. begin();
//                      rk_p != rk. end(); ++ rk_p, ++ rn_p) {

//                     *rk_p = (*rk_p) * v + (*rn_p) * u;
//                         //normailze here?
//                         // *rk_p %= pk;
//                 }

                inv(u, pk, pn);
                u *= pk;	// u = 0 mod pk, u = 1 mod pn;
                typename Vector::iterator rk_p, rn_p;
                for (rk_p = rk. begin(), rn_p = rn. begin();
                     rk_p != rk. end(); ++ rk_p, ++ rn_p) {
                    v = *rn_p;
                    v -= *rk_p;        
                    v *= u;         // v = 0 mod pk, v = (rn-rk) mod pn;
                    *rk_p += v;
                }
                pk *= pn;
            }
		//std::cout << "After calling:\n";
		//debug();
	}

	void buildTree (List& holdp, List& holdv) {

            typename List::iterator half_res;
            typename List::iterator half_p;
            while (holdres. size () > 1) {
                combine (holdp, holdv);
                half_res = holdv. begin();
                half_p = holdp. begin();
                for (size_t i = 0; i < ((holdv. size() + 1) / 2); ++ i) {
                    ++ half_res;
                    ++ half_p;
                }
                holdv. erase (half_res, holdv. end());
                holdp. erase (half_p, holdp. end());
			
            }
	}

	void combine (List& holdp, List& holdv){
            typename List::iterator first_res, half_res;
            typename List::iterator first_p, half_p;
            first_res = half_res = holdv. begin();
            first_p = half_p = holdp. begin();
            for ( size_t i = 0; i < ((holdv. size() + 1) / 2); ++ i) {
                ++ half_res;
                ++ half_p;
            }
		
            Integer g, u, v;
		//for (;half_res != holdres. end(); 
            for (;half_res != holdv. end(); 
                 ++ half_res, ++ first_res, ++ half_p, ++ first_p) {

                Integer& pk = *first_p;
                Integer& pn = *half_p;
                Integer& rk = *first_res;
                Integer& rn = *half_res;
//                 gcd (g, pk, pn, u, v);
//                 u *= pk; // u = 0 mod pk, u = 1 mod pn;
//                 v *= pn; // v = 1 mod pk, v = 0 mod pn;
//                 pk *= pn;

//                 rk = u*rn  +  v*rk ;

                inv(u, pk, pn);
                u *= pk;	// u = 0 mod pk, u = 1 mod pn;
                v = rn;
                v -= rk;        
                v *= u;         // v = 0 mod pk, v = (rn-rk) mod pn;
                rk += v;
                pk *= pn;
            }
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



/* KEEP TRACK OF BRADFORD CODE */

// template <class Vector>
// integer &cra (integer      &res,
// 	      const Vector &residues,
// 	      const Vector &moduli)
// {
// 	linbox_check (residues.size () == moduli.size ());

// 	commentator.start ("Chinese remainder algorithm", "cra", residues.size ());

// 	integer pi = 1L, pi_m_j, s;
// 	typename Vector::const_iterator i, j;

// 	for (i = moduli.begin (); i != moduli.end (); ++i)
// 		integer::mulin (pi, (unsigned long) *i);

// 	commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
// 		<< "Product of moduli: " << pi << endl;

// 	res = 0L;

// 	for (i = residues.begin (), j = moduli.begin (); j != moduli.end (); ++i, ++j) {
// 		integer::div (pi_m_j, pi, (unsigned long) *j);

// 		// Dan Roche 8-7-04 The function invmod was changed to inv at
// 		// some point.  We should use this instead of making a new
// 		// Modular<integer> field every time.
// 		inv (s, pi_m_j, integer (*j));

// 		integer::mulin (s, (unsigned long) *i);
// 		integer::modin (s, (unsigned long) *j);
// 		integer::mulin (s, pi_m_j);
// 		integer::addin (res, s);
// 		commentator.progress ();
// 	}

// 	integer::modin (res, pi);
// 	integer pio2 = pi>>1;
// 	if ( res > pio2)
// 		integer::subin(res, pi);
// 	if ( res < -pio2)
// 		integer::addin(res, pi);
// 	commentator.stop ("done", NULL, "cra");

// 	return res;
// }		
