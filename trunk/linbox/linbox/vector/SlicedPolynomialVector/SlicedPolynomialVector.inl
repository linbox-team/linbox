#ifndef __LINBOX_matrix_SlicedPolynomialVector_SlicedPolynomialVector_INL
#define __LINBOX_matrix_SlicedPolynomialVector_SlicedPolynomialVector_INL

namespace Linbox
{
						//////////////////////////
		        			//irreducible polynomial//
						//////////////////////////

	template < class _Field, class _Rep, class _VectorElement >
	polynomial& SlicedPolynomialVector< _Field, _Rep, _VectorElement >::modulo(polynomial& g, polynomial&h, polynomial& f)
	{
		polynomial w1;
		Poly1Dom<Domain,Dense>::div(w1, h, f);
		polynomial w2;
		Poly1Dom<Domain,Dense>::mul(w2, w1, f);
		Poly1Dom<Domain,Dense>::sub(g, h, w2);
		return g;
	}
	
	/*
	algorithm description
	http://en.wikipedia.org/wiki/Factorization_of_polynomials_over_finite_fields#Rabin.27s_test_of_irreducibility
	Algorithm Rabin Irreducibility Test
	 Input: A monic polynomial f in Fq[x] of degree n, 
	        p1, ..., pk all distinct prime divisors of n.
	 Output: Either "f is irreducible" or "f is reducible".
	 Begin
	     for j = 1 to k do 
	        n_j=n/p_j;
	     for i = 1 to k do 
	        h:=x^{q^{n_i}}-x \bmod f;
	        g := gcd(f, h);
	        if g ? 1, then return 'f is reducible' and STOP;
	     end for;
	     h:= x^{q^{n}}-x \bmod f;
	     if h = 0, then return "f is irreducible", 
	         else return "f is reducible"
	 end.
	 */
	/*
	function description
	 Fq[x], f = x^n + bx + a
	 returns true if f is irreducible and false if f is reducible
	 */
	template < class _Field, class _Rep, class _VectorElement >
	bool SlicedPolynomialVector< _Field, _Rep, _VectorElement >::rabin(int n, int q, int a, int b)
	{
		polynomial f(n + 1);
		f[0] = a;
		f[1] = b;
		for (int i = 2; i < n; i++)
		{
			f[i] = 0;
		}
		f[n] = 1;
		std::vector<int> pd;
		factorize(n, pd); //factorizes n into primes 
		//n = pd[0]^alpha[0] * ... * pd[k - 1]^alpha[k - 1]
		int k = pd.size();
		std::vector<int> nd(k);
		for (int j = 0; j < k; j++)
		{
			nd[j] = n / pd[j];
		}
		polynomial vector_zero();
		assign(vector_zero, 0);
		polynomial vector_one();
		assign(vector_one, 1);
		for (int j = 0; j < k; j++)
		{
			int h_size = pow(q, nd[j]) + 1;
			polynomial h(h_size);
			h[0] = 0;
			h[1] = -1;
			for (int i = 2; i < h_size - 1; i++)
			{
				h[i] = 0;
			}
			h[h_size - 1] = 1;
			polynomial g;
			Poly1Dom<IntField, Dense>::gcd (g, f, h);
			bool g_equals_1 = g.areEqual(vector_one);
			if (! g_equals_1)
			{
				return true; //f is irreducible
			}
		}
		int h_size = pow(q, n) + 1;
		polynomial h(h_size);
		h[0] = 0;
		h[1] = -1;
		for (int i = 2; i < h_size - 1; i++)
		{
			h[i] = 0;
		}
		h[h_size - 1] = 1;
		polynomial g;
		modulo(g, h, f);//!!!
		bool g_equals_0 = g.areEqual(vector_zero);
		return g_equals_0; //if g = 0, then return "f is irreducible", else return "f is reducible" 
	}
	
	template < class _Field, class _Rep, class _VectorElement >
	void SlicedPolynomialVector< _Field, _Rep, _VectorElement >::setIrreduciblePlynomial(int max_steps)
	{
		bool found = false;
		int w = F.characteristic();
		for (int step = 0; (step < max_steps) && (!found); step++)
		{
			int a = rand() % w;
			int b = rand() % w;
			if (a + b > 0)
			{
				found = rabin(n, w, a, b);
			}
		}
		irreducible.push_back(a);
		irreducible.push_back(b);
		for (int i = 1; i < n; i++)
		{
			irreducible.push_back(0);
		}
		irreducible.push_back(1);
		return;
	}
	
						////////////////
		        			//Constructors//
						////////////////

	template < class _Field, class _Rep, class _VectorElement >
	SlicedPolynomialVector< _Field, _Rep, _VectorElement >::SlicedPolynomialVector (const _Field &BF)
	{
		GF = BF;
		F_temp = IntField(GF.characteristic());
		F = F_temp;
		n = GF.exponent(); //GF = GF(p^n)
		for (size_t r = 0; r < n; r++)
		{
			V.emplace_back(BlasVector<IntField>(F));
		}
		setIrreduciblePolynomial();
	}

	template < class _Field, class _Rep, class _VectorElement >
	SlicedPolynomialVector< _Field, _Rep, _VectorElement >::SlicedPolynomialVector (const _Field &BF, const size_t &m)
	{
		GF = BF;
		F_temp = IntField(GF.characteristic()); //public function to set characteristic?
		F = F_temp;
		n = GF.exponent(); //GF = GF(p^n)
		for (size_t r = 0; r < n; r++)
		{
			V.emplace_back(BlasVector<IntField>(F, m));
		}
		setIrreduciblePolynomial();
	}
	
	template < class _Field, class _Rep, class _VectorElement >
	SlicedPolynomialVector< _Field, _Rep, _VectorElement >::SlicedPolynomialVector (const _Field &BF, polynomial& pp)
	{
		GF = BF;
		F_temp = IntField(GF.characteristic()); //public function to set characteristic?
		F = F_temp;
		n = GF.exponent(); //GF = GF(p^n)
		for (size_t r = 0; r < n; r++)
		{
			V.emplace_back(BlasVector<IntField>(F));
		}
		irreducible = pp;
	}

	template < class _Field, class _Rep, class _VectorElement >
	SlicedPolynomialVector< _Field, _Rep, _VectorElement >::SlicedPolynomialVector (const _Field &BF, const size_t &m, polynomial& pp)
	{
		GF = BF;
		F_temp = IntField(GF.characteristic()); //public function to set characteristic?
		F = F_temp;
		n = GF.exponent(); //GF = GF(p^n)
		for (size_t r = 0; r < n; r++)
		{
			V.emplace_back(BlasVector<IntField>(F, m));
		}
		irreducible = pp;
	}

						///////////////
						// Destructor//
						///////////////

	template < class _Field, class _Rep, class _VectorElement >
	SlicedPolynomialVector< _Field, _Rep, _VectorElement >::~SlicedPolynomialVector()
	{
		//LidiaGfq has a destructor, GivaroGfq doesn't, so currently field type members GF and F are not destroyed
		V.~vector();
		//if some members are added, delete them			
	}
						////////////////////////
		        			//dimensions of vector//
						////////////////////////

        template < class _Field, class _Rep, class _VectorElement >
	size_t SlicedPolynomialVector< _Field, _Rep, _VectorElement >::length() const
	{
		return V.size();				
	}

	template < class _Field, class _Rep, class _VectorElement >
	size_t SlicedPolynomialVector< _Field, _Rep, _VectorElement >::rowdim() const
	{
		return V[0].size();				
	}
	
	                    			/////////////////
	                    			//return fields//
	                    			/////////////////

	template < class _Field, class _Rep, class _VectorElement >
	const Field& SlicedPolynomialVector< _Field, _Rep, _VectorElement >::fieldGF() const
	{
		return GF;
	}

	template < class _Field, class _Rep, class _VectorElement >
	const IntField& SlicedPolynomialVector< _Field, _Rep, _VectorElement >::fieldF() const
	{
		return F;
	}

						/////////////////////////
		        			//functions for entries//
						/////////////////////////
		
        template < class _Field, class _Rep, class _VectorElement >
	void SlicedPolynomialVector< _Field, _Rep, _VectorElement >::setEntry (size_t m, size_t k, const VectorElement &a_mk)
	{
		V[m].setEntry(k, a_mk);
	}

	template < class _Field, class _Rep, class _VectorElement >
	_VectorElement & SlicedPolynomialVector< _Field, _Rep, _VectorElement >::refEntry (size_t m, size_t k)
	{
		return V[m].refEntry(k);

	}

	template < class _Field, class _Rep, class _VectorElement >
	_VectorElement & SlicedPolynomialVector< _Field, _Rep, _VectorElement >::getEntry (size_t m, size_t k)
	{
		return V[m].getEntry(k);

	}
	
						/////////////////////////////////////
		                		//functions for matrix-coefficients//
						/////////////////////////////////////

	template < class _Field, class _Rep, class _VectorElement >
	void SlicedPolynomialVector< _Field, _Rep, _VectorElement >::setMatrixCoefficient (size_t m, const BlasVector<IntField> &V_m)
	{
		V[m] = V_m;
	}

	template < class _Field, class _Rep, class _VectorElement >
	BlasVector<IntField> &SlicedPolynomialVector< _Field, _Rep, _VectorElement >::refMatrixCoefficient (size_t m)
	{
		return V[m];
	}

	template < class _Field, class _Rep, class _VectorElement >
	const BlasVector<IntField> &SlicedPolynomialVector< _Field, _Rep, _VectorElement >::getMatrixCoefficient (size_t m) const
	{
		return V[m];
	}

						/////////
		                		//swaps//
						/////////

	template < class _Field, class _Rep, class _VectorElement >
	void SlicedPolynomialVector< _Field, _Rep, _VectorElement >::swapRows(size_t k1, size_t k2)
	{
		for (size_t m = 0; m < this->length(); m++)
		{
			VectorElement c = this->getEntry(m, k1);
			this->setEntry(m, k1, this->getEntry(m, k2));
			this->setEntry(m, k2, c);
		}
	}
	
						//////////////////
		                		//input / output//
						//////////////////	
	template < class _Field, class _Rep, class _VectorElement >
	std::istream& SlicedPolynomialVector< _Field, _Rep, _VectorElement >::read (std::istream &file)
	{
		int M = this->length();
		int K = this->rowdim();
		VectorElement c;
		for (int m = 0; m < M; m++)
		{
			for (int k = 0; k < K; k++)
			{
				file >> c;
				this->setEntry(m, k, c);
			}
		}
		return file;
	}
	
	template < class _Field, class _Rep, class _VectorElement >
	std::ostream& SlicedPolynomialVector< _Field, _Rep, _VectorElement >::write (std::ostream &file)
	{
		int M = this->length();
		int K = this->rowdim();
		for (int m = 0; m < M; m++)
		{
			for (int k = 0; k < K; k++)
			{
					file << this->getEntry(m, k) << std::endl;
			}
			file << std::endl;
		}
		return file;
	}
}

#endif