#ifndef __LINBOX_matrix_SlicedPolynomialMatrix_SlicedPolynomialMatrix_INL
#define __LINBOX_matrix_SlicedPolynomialMatrix_SlicedPolynomialMatrix_INL

namespace Linbox
{
						//////////////////////////
		        			//irreducible polynomial//
						//////////////////////////

	template < class _Field, class _Rep, class _MatrixElement >
	polynomial& SlicedPolynomialMatrix< _Field, _Rep, _MatrixElement >::modulo(polynomial& g, polynomial&h, polynomial& f)
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
	template < class _Field, class _Rep, class _MatrixElement >
	bool SlicedPolynomialMatrix< _Field, _Rep, _MatrixElement >::rabin(int n, int q, int a, int b)
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
	
	template < class _Field, class _Rep, class _MatrixElement >
	void SlicedPolynomialMatrix< _Field, _Rep, _MatrixElement >::setIrreduciblePlynomial(int max_steps)
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

	template < class _Field, class _Rep, class _MatrixElement >
	SlicedPolynomialMatrix< _Field, _Rep, _MatrixElement >::SlicedPolynomialMatrix (const _Field &BF)
	{
		GF = BF;
		F_temp = IntField(GF.characteristic()); //public function to set characteristic?
		F = F_temp;
		int n = GF.exponent(); //GF = GF(p^n)
		for (size_t m = 0; m < n; m++)
		{
			V.emplace_back(BlasMatrix<IntField>(F));
		}
		setIrreduciblePolynomial();
	}

	template < class _Field, class _Rep, class _MatrixElement >
	SlicedPolynomialMatrix< _Field, _Rep, _MatrixElement >::SlicedPolynomialMatrix (const _Field &BF, const size_t & m1, const size_t &m2)
	{
		GF = BF;
		F_temp = IntField(GF.characteristic()); //public function to set characteristic?
		F = F_temp;
		int n = GF.exponent(); //GF = GF(p^n)
		for (size_t m = 0; m < n; m++)
		{
			V.emplace_back(BlasMatrix<IntField>(F, m1, m2));
		}
		setIrreduciblePolynomial();
	}
	
	template < class _Field, class _Rep, class _MatrixElement >
	SlicedPolynomialMatrix< _Field, _Rep, _MatrixElement >::SlicedPolynomialMatrix (const _Field &BF, polynomial& pp)
	{
		GF = BF;
		F_temp = IntField(GF.characteristic()); //public function to set characteristic?
		F = F_temp;
		int n = GF.exponent(); //GF = GF(p^n)
		for (size_t m = 0; m < n; m++)
		{
			V.emplace_back(BlasMatrix<IntField>(F));
		}
		irreducible = pp;
	}

	template < class _Field, class _Rep, class _MatrixElement >
	SlicedPolynomialMatrix< _Field, _Rep, _MatrixElement >::SlicedPolynomialMatrix (const _Field &BF, const size_t & m1, const size_t &m2, polynomial& pp)
	{
		GF = BF;
		F_temp = IntField(GF.characteristic()); //public function to set characteristic?
		F = F_temp;
		int n = GF.exponent(); //GF = GF(p^n)
		for (size_t m = 0; m < n; m++)
		{
			V.emplace_back(BlasMatrix<IntField>(F, m1, m2));
		}
		irreducible = pp;
	}

						///////////////
						// Destructor//
						///////////////

	template < class _Field, class _Rep, class _MatrixElement >
	SlicedPolynomialMatrix< _Field, _Rep, _MatrixElement >::~SlicedPolynomialMatrix()
	{
		//LidiaGfq has a destructor, GivaroGfq doesn't, so currently field type members GF and F are not destroyed
		V.~vector();
		//if some members are added, delete them			
	}
						////////////////////////
		        			//dimensions of vector//
						////////////////////////

        template < class _Field, class _Rep, class _MatrixElement >
	size_t SlicedPolynomialMatrix< _Field, _Rep, _MatrixElement >::length() const
	{
		return V.size();				
	}

	template < class _Field, class _Rep, class _MatrixElement >
	size_t SlicedPolynomialMatrix< _Field, _Rep, _MatrixElement >::rowdim() const
	{
		return V[0].rowdim();				
	}

	template < class _Field, class _Rep, class _MatrixElement >
	size_t SlicedPolynomialMatrix< _Field, _Rep, _MatrixElement >::coldim() const
	{
		return V[0].coldim();				
	}
	
	                    			/////////////////
	                    			//return fields//
	                    			/////////////////

	template < class _Field, class _Rep, class _MatrixElement >
	const Field& SlicedPolynomialMatrix< _Field, _Rep, _MatrixElement >::fieldGF() const
	{
		return GF;
	}

	template < class _Field, class _Rep, class _MatrixElement >
	const IntField& SlicedPolynomialMatrix< _Field, _Rep, _MatrixElement >::fieldF() const
	{
		return F;
	}

						/////////////////////////
		        			//functions for entries//
						/////////////////////////
		
        template < class _Field, class _Rep, class _MatrixElement >
	const MatrixElement& SlicedPolynomialMatrix< _Field, _Rep, _MatrixElement >::setEntry (size_t m, size_t i, size_t j, const MatrixElement &a_mij)
	{
		V[m].setEntry(i, j, a_mij);
        return a_mij;
	}

	template < class _Field, class _Rep, class _MatrixElement >
	_MatrixElement & SlicedPolynomialMatrix< _Field, _Rep, _MatrixElement >::refEntry (size_t m, size_t i, size_t j)
	{
		return V[m].refEntry(i, j);

	}

	template < class _Field, class _Rep, class _MatrixElement >
	_MatrixElement & SlicedPolynomialMatrix< _Field, _Rep, _MatrixElement >::getEntry (size_t m, size_t i, size_t j)
	{
		return V[m].getEntry(i, j);

	}
	
						/////////////////////////////////////
		                		//functions for matrix-coefficients//
						/////////////////////////////////////

	template < class _Field, class _Rep, class _MatrixElement >
	void SlicedPolynomialMatrix< _Field, _Rep, _MatrixElement >::setMatrixCoefficient (size_t m, const BlasMatrix<IntField> &V_m)
	{
		V[m] = V_m;
	}

	template < class _Field, class _Rep, class _MatrixElement >
	BlasMatrix<IntField> &SlicedPolynomialMatrix< _Field, _Rep, _MatrixElement >::refMatrixCoefficient (size_t m)
	{
		return V[m];
	}

	template < class _Field, class _Rep, class _MatrixElement >
	const BlasMatrix<IntField> &SlicedPolynomialMatrix< _Field, _Rep, _MatrixElement >::getMatrixCoefficient (size_t m) const
	{
		return V[m];
	}

						/////////
		                		//swaps//
						/////////

	template < class _Field, class _Rep, class _MatrixElement >
	void SlicedPolynomialMatrix< _Field, _Rep, _MatrixElement >::swapRows(size_t i1, size_t i2)
	{
		for (size_t m = 0; m < this->length(); m++)
		{
			for (size_t j = 0; j < this->coldim(); j++)
			{
				MatrixElement c = this->getEntry(m, i1, j);
				this->setEntry(m, i1, j, this->getEntry(m, i2, j));
				this->setEntry(m, i2, j, c);
			}
		}
	}

	template < class _Field, class _Rep, class _MatrixElement >
	void SlicedPolynomialMatrix< _Field, _Rep, _MatrixElement >::swapCols(size_t j1, size_t j2)
	{
		for (size_t m = 0; m < this->length(); m++)
		{
			for (size_t i = 0; i < this->colrow(); i++)
			{
				MatrixElement c = this->getEntry(m, i, j1);
				this->setEntry(m, i, j1, this->getEntry(m, i, j2));
				this->setEntry(m, i, j2, c);
			}
		}
	}
	
						/////////////
		                		//transpose//
						/////////////

	template < class _Field, class _Rep, class _MatrixElement >
	SlicedPolynomialMatrix< _Field, _Rep, _MatrixElement > SlicedPolynomialMatrix< _Field, _Rep, _MatrixElement >::transpose(SlicedPolynomialMatrix< _Field, _Rep, _MatrixElement > & tV) const
	{
		//check dimensions
		for (size_t m = 0; m < this->length(); m++)
		{
			this->getMatrixCoefficent(m).transpose(tV.refMatrixCoefficent(m));
		} 
		return tV;
	}

						//////////////////
		                		//input / output//
						//////////////////	
	template < class _Field, class _Rep, class _MatrixElement >
	std::istream& SlicedPolynomialMatrix< _Field, _Rep, _MatrixElement >::read (std::istream &file)
	{
		int K = this->length();
		int I = this->rowdim();
		int J = this->coldim();
		MatrixElement c;
		for (int k = 0; k < K; k++)
		{
			for (int i = 0; i < I; i++)
			{
				for (int j = 0; j < J; j++)
				{
					file >> c;
					this->setEntry(k, i, j, c);
				}
			}
		}
		return file;
	}
	
	template < class _Field, class _Rep, class _MatrixElement >
	std::ostream& SlicedPolynomialMatrix< _Field, _Rep, _MatrixElement >::write (std::ostream &file)
	{
		int K = this->length();
		int I = this->rowdim();
		int J = this->coldim();
		for (int k = 0; k < K; k++)
		{
			for (int i = 0; i < I; i++)
			{
				for (int j = 0; j < J; j++)
				{
					file << this->getEntry(k, i, j) << " ";
				}
				file << std::endl;
			}
			file << std::endl;
		}
		return file;
	}
}

#endif

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
