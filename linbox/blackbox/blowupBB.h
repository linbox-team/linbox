/* linbox/blackbox/blowup.h
 * Copyright (C) 2023 BDS
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

/*! @file blackbox/blowup.h
 * @ingroup blackbox
 * @brief vis a vis irreducible degree d poly f(x), blowup of A \in (F[x]/<f^e>)^nxn to a 
 * a matrix coefficient field, specifically in F^ndexnde.
 * This class is especially meant for use in Frobenius form algorithms.
 * Based on ElSheik et al, "Fast Computation for Smith Forms of Sparse Matrices Over Local Rings", ISSAC 2012. 
 */

#ifndef __LINBOX_blowup_H
#define __LINBOX_blowup_H

#include <vector>
#include "linbox/linbox-config.h"
#include "linbox/matrix/dense-matrix.h"
#include "linbox/vector/vector-domain.h"
#include "linbox/vector/blas-vector.h"
#include "bb.h"

namespace LinBox
{

	/**
    * @brief Given blackbox matrix A in F^nxn and prime power polynomial f(x) of degree d, 
    * this is a blowup of xI-A modulo f to a dimension nd blackbox matrix over 
    * the coefficient field F.
    *
    * Future extension might represent arbitrary polynomial matrices modulo f in this way.
    *
	 * \ingroup blackbox
    * This class is especially meant for use in Frobenius form algorithms.
    * Based on ElSheik et al, "Fast Computation for Smith Forms of Sparse Matrices Over Local Rings", ISSAC 2012. 
	 *
	 */


template <class _Field>
class Blowup
: public BB<_Field> {

public:
	typedef Blowup<_Field> Self_t;

   using Father_t = BB<_Field>;
   using Field = typename Father_t::Field;
   using Element = typename Father_t::Element;
   using Matrix = typename Father_t::Matrix;
   using Vector_t = BlasVector<Field>;

	/// \brief cstor setting this to the blowup of xI-A modulo f
	Blowup(BB<Field> &baseBB, Vector_t & modulusPoly) 
   : _A(&baseBB), _f(modulusPoly) {}
         
	~Blowup(){}

   BBType bbTag() const { return other; } // for now

   BB<Field> & baseBB() const { return &_A; }

   Vector_t modulusPoly() const { return _f; }

   void mywrite(const Matrix & V) const {
   typename Field::Element tmp;
      for (size_t i = 0; i < V.rowdim(); ++i) {
         for (size_t j = 0; j < V.coldim(); ++j)
            std::cout << V.getEntry(tmp, i, j) << ", ";
         std::cout << std::endl;
      }
      std::cout << std::endl;
   }

   // strawman implementation
	template <class OutVector, class InVector>
	OutVector &apply (OutVector &v, const InVector &u) const {
      return applyBase(v, u, false);
   }

	template <class OutVector, class InVector>
	OutVector &applyTranspose (OutVector &v, const InVector &u) const {
      return applyBase(v, u, true);
   }

	template <class OutVector, class InVector>
	OutVector &applyBase (OutVector &v, const InVector &u, bool Tr) const {
      size_t n = _A->rowdim();
      size_t d = _f.size() - 1; // poly degree
      typename Field::Element tmp; field().init(tmp);
      VectorDomain<Field> VD(field());
      Vector_t w(field(), n*d);

   // mul by x modulo _f
#if 0
      //shiftMod(v, u); 
      typename Field::Element a;   field().init(a);
      for (size_t i = 0; i < n*d; i += d) {
         // In row i, let g(x) be poly whose j-th coeff is u[i+j].
         a = u[i + d-1]; // leading coeff of g(x).
         field().negin(a); 
         // h(x) = xg(x) + a f(x) is xg(x) modulo f and has degree < d.
         field().mul(v[i + 0], a, _f[0]); 

         // v[i + j] = a*_f[j] + u[i + j - 1].
         for (size_t j = 1; j < d; ++j ) {
            field().axpy(v[i + j], a, _f[j], u[i + j - 1]); // h_j
         }
      } //now v = phi(xI mod f)*u or u*phi(xI mod f).
#endif
#if 1
   // mul by _A 
      Matrix U(field(), (Tr ? d : n), (Tr ? n : d));
      Matrix V(field(), (Tr ? d : n), (Tr ? n : d));

      //blockify(U, u, Tr);
      // For each d-block as poly g(x), make h(x) = xg(x) - g(d-1)f(x).
      // Then each d-block w of v becomes w = h(x) - w.
      for(size_t i = 0; i < n; ++i) 
         for(size_t j = 0; j < d; ++j)
            U.setEntry((Tr ? j : i), (Tr ? i : j), u[i*d + j]);

      if (Tr) _A->applyLeft (V, U); 
      else    _A->applyRight(V, U); 

      //unblockify(w, V, Tr); 
      for(size_t i = 0; i < n; ++i) 
         for(size_t j = 0; j < d; ++j){
            V.getEntry(tmp, (Tr ? j : i), (Tr ? i : j));
            field().assign(w[i*d + j], tmp);
      // now w = phi(A mod f)*u or u*phi(A^T mod f).
      }
//  std::cout << "U in apply:" << Tr << std::endl; mywrite(U);
//  std::cout << "V in apply:" << Tr << std::endl; mywrite(V);
      VD.subin(v, w);
      // now v = phi((xI-A) mod f)*u or u*phi((xI-A)^T mod f).
#endif
      return v;
   }

		Matrix& applyRight(Matrix& Y, const Matrix& X) const // Y = AX
		{   MatrixDomain<Field> MD(field());
		    return MD.mul(Y, *this, X);
		}

		Matrix& applyLeft(Matrix& Y, const Matrix& X) const // Y = XA
		{   MatrixDomain<Field> MD(field());
		    return MD.mul(Y, X, *this);
		}
      
		size_t rowdim(void) const { 
         return _A->rowdim()*(_f.size()-1);
      }
		size_t coldim(void) const { 
         return _A->coldim()*(_f.size()-1);
      }
		const Field& field() const { return _A->field(); }

		std::ostream& write(std::ostream& os) const
		{
			return os;
		}

		std::istream& read(std::istream& is)
		{
			return is;
		}

		template<typename _Tp1>
		struct rebind {
			typedef Blowup<_Tp1> other;

			void operator() (other & Ap, const Self_t& A)
			{}
      };
#if 0
		template<typename _Tp1>
		struct rebind {
			typedef Blowup<_Tp1> other;

			void operator() (other & Ap, const Self_t& A)
			{

            rebind( &Ap._A, &A._A);

				Hom<typename Self_t::Field, _Tp1> hom(A.field(), Ap.field());

            Ap._f.resize(A._f.size());
				typename BlasVector<_Tp1>::iterator nit = Ap._f.begin();
				typename Vector_t::const_iterator oit = A._f.begin();
				for( ; oit != A._f.end() ; ++nit, ++oit)
					hom.image (*nit, *oit);
			}

		};

#endif
		template<typename _Tp1>
		Blowup(const Blowup<_Tp1>& B, const Field& F) : 
         _f(F)
		{
			typename Blowup<_Tp1>::template rebind<Field>() (*this, B);
		}
	protected:

      BB<Field> * _A;
      // later: BB<PolyRing> & _B;
      //PolyRing *_R;
      Vector_t _f;

	}; // template <Field, Vector> class Blowup

} // namespace LinBox

#endif // __LINBOX_blowup_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
