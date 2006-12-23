/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/blackbox/direct-sum.h
 * David Saunders
 *
 * See COPYING about license
 */

#ifndef __DIRECT_SUM_H
#define __DIRECT_SUM_H

#include <linbox/blackbox/null-matrix.h>
#include <linbox/vector/vector-traits.h>
#include <linbox/blackbox/blackbox-interface.h>
#include <linbox/vector/subvector.h>
#include <linbox/matrix/matrix-domain.h>

namespace LinBox
{

	template <class Blackbox1, class Blackbox2 = Blackbox1>
	class DirectSum;

	/** \brief If C = DirectSum(A, B) and y = xA and z = wB, then (y,z) = (x,w)C.

	 * And similarly for apply. 
	\ingroup blackbox
	 */
	template <class Blackbox1, class Blackbox2>
	class DirectSum : public BlackboxInterface
	{
	    public:
            typedef DirectSum<Blackbox1, Blackbox2> Self_t;
            typedef Blackbox1 Blackbox1_t;
            typedef Blackbox2 Blackbox2_t;

		typedef typename Blackbox1::Field Field;
		typedef typename Blackbox1::Element Element;

		/** Constructor from two black box matrices.
		 * This becomes direct sum of A and B.
		 * They may be rectangular.  
		 * @param A, B:  black box matrices over a common field.
		 */
		DirectSum(const	Blackbox1& A, const Blackbox2& B)
			: _Ap(&A), _Bp(&B)
		{}

		/** Constructor from two black box matrix pointers.
		 * This becomes direct sum of A and B.
		 * They may be rectangular.  They must be over the same field (or ring). 
		 * @param A_ptr pointer to black box matrix A.
		 * @param B_ptr pointer to black box matrix B.
		 */
	// 	DirectSum(): _Ap(& DirectSum<Blackbox1, Blackbox2>::_NullMatrix),
// 			     _Bp(& DirectSum<Blackbox1, Blacbox2>::_NullMatrix) 
// 		{}



		DirectSum(const	Blackbox1* Ap, const Blackbox2* Bp)
			: _Ap(Ap), _Bp(Bp)
		{}

		/// Copy constructor.
		DirectSum (const DirectSum<Blackbox1, Blackbox2>& M) 
			: _Ap (M._Ap), _Bp (M._Bp)
		{}


		/// Destructor
		~DirectSum (void)
		{}

		template<class OutVector, class InVector>
		OutVector& apply (OutVector& y, const InVector& x) const
		{
			linbox_check(x.size() == coldim());
			linbox_check(y.size() == rowdim());
		// FIXED: I want to use subvectors to avoid copying and memory allocation, but problems with it...

			const Subvector<typename InVector::const_iterator> x1(x.begin(), x.begin() + _Ap->coldim());
			const Subvector<typename InVector::const_iterator> x2(x.begin() + _Ap->coldim(), x.end());
			Subvector<typename OutVector::iterator> y1(y.begin(), y.begin() + _Ap->rowdim());
			Subvector<typename OutVector::iterator> y2(y.begin() + _Ap->rowdim(), y.end());
			_Ap -> apply(y1,x1);
			_Bp -> apply(y2,x2);
		/*
			if (x.size() == 0) return y;  // Null matrix

			InVector xA(_Ap->coldim());
			OutVector yA(_Ap->rowdim());
			std::copy (x.begin(), x.begin() + _Ap->coldim(), xA.begin());
			_Ap->apply (yA, xA);
			std::copy (yA.begin(), yA.end(), y.begin());

			InVector xB(_Bp->coldim());
			OutVector yB(_Bp->rowdim());
			std::copy (x.begin() + _Ap->coldim(), x.end(), xB.begin());
			_Bp->apply (yB, xB);
			std::copy (yB.begin(), yB.end(), y.begin() + _Ap->rowdim());

			y.resize(_Ap->rowdim());
			x.resize(_Ap->coldim());
			_Ap->apply(y,x);
		*/
			return y;
		}


		template<class OutVector, class InVector>
		OutVector& applyTranspose (OutVector& y, const InVector& x) const
		{
			linbox_check(x.size() == rowdim());
			linbox_check(y.size() == coldim());
		// FIXED: I want to use subvectors to avoid copying and memory allocation, but problems with it...
			const Subvector<typename InVector::const_iterator> x1(x.begin(), x.begin() + _Ap->rowdim());
			const Subvector<typename InVector::const_iterator> x2(x.begin() + _Ap->rowdim(), x.end());
			Subvector<typename OutVector::iterator> y1(y.begin(), y.begin() + _Ap->coldim());
			Subvector<typename OutVector::iterator> y2(y.begin() + _Ap->coldim(), y.end());
			_Ap -> applyTranspose(y1,x1);
			_Bp -> applyTranspose(y2,x2);

		/*
			Vector local_x(x1.size());
			Vector local_y(y1.size());
			copy (x1.begin(), x1.end(), local_x.begin());
			_Ap->applyTranspose (local_y, local_x);
			copy (local_y.begin(), local_y.end(), y1.begin());

			local_x.resize(x1.size());
			local_y.resize(y2.size());
			copy (x2.begin(), x2.end(), local_x.begin());
			_Bp->applyTranspose (local_y, local_x);
			copy (local_y.begin(), local_y.end(), y2.begin());

			//_Ap->applyTranspose (y1, x1);
			//_Bp->applyTranspose (y2, x2);

			if (x.size() == 0 ) return y;
			InVector xAT(_Ap->rowdim());
			OutVector yAT(_Ap->coldim());
			std::copy (x.begin(), x.begin() + _Ap->rowdim(), xAT.begin());
			_Ap->apply (yAT, xAT);
			std::copy (yAT.begin(), yAT.end(), y.begin());

			InVector xBT(_Bp->rowdim());
			OutVector yBT(_Bp->coldim());
			std::copy (x.begin() + _Ap->rowdim(), x.end(), xBT.begin());
			_Bp->apply (yBT, xBT);
			std::copy (yBT.begin(), yBT.end(), y.begin() + _Ap->rowdim());
		*/


			return y;
		}

            template<typename _Tp1, typename _Tp2 = _Tp1>
            struct rebind
            { 
                typedef DirectSum<typename Blackbox1::template rebind<_Tp1>::other, typename Blackbox2::template rebind<_Tp2>::other> other; 

		void operator() (other *& Ap, const Self_t& A, const _Tp1& F) {
                    typename other::Blackbox1_t * A1;
                    typename Blackbox1_t::template rebind<_Tp1> () ( A1, *(A._Ap), F);
                    typename other::Blackbox2_t * A2;
                    typename Blackbox2_t::template rebind<_Tp1> () ( A2, *(A._Bp), F);
                    Ap = new other(*A1, *A2);
                }

            };

		inline size_t rowdim (void) const
		{

			return _Ap->rowdim () + _Bp->rowdim ();

		}
    
		inline size_t coldim(void) const 
		{
			return _Ap->coldim () + _Bp->coldim ();
		}


		const Field& field() const { return _Ap -> field();}
	    protected:
		// the direct summands
		const Blackbox1* _Ap;
		const Blackbox2* _Bp; 
		

	}; // template <Vector> class DirectSum

	

        template <class Blackbox>
        class DirectSum<Blackbox, Blackbox>
        {
            public:
            typedef DirectSum<Blackbox, Blackbox> Self_t;
            typedef Blackbox Blackbox_t;

                typedef typename Blackbox::Field Field;
                typedef typename Blackbox::Element Element;
		DirectSum() : m(0),n(0){}
                DirectSum(const Blackbox& A, const Blackbox& B) {
                	_VB.push_back(&A);
			_VB.push_back(&B);
			m = A.rowdim() + B.rowdim();
			n = A.coldim() + B.coldim();
		}

                DirectSum(const Blackbox* Ap, const Blackbox* Bp) {
			_VB.push_back(Ap);
			_VB.push_back(Bp);
			m = Ap -> rowdim() + Bp -> rowdim();
			n = Ap -> coldim() + Bp -> coldim();
		}
		template<class Vector>
		DirectSum(const Vector& v) :_VB(v.begin(),v.end()) {
			m = 0; n = 0;
			typename std::vector<const Blackbox* >::iterator bp;
			for( bp = _VB.begin(); bp != _VB.end(); ++bp) {
				m += (*bp) -> rowdim();
				n += (*bp) -> coldim();
			}
		}
		
                /// Copy constructor.
                DirectSum (const DirectSum<Blackbox, Blackbox>& M)
                        :_VB(M._VB),m(M.m),n(M.n)
                {}


                /// Destructor
                ~DirectSum (void)
                {}

            template<typename _Tp1>
            struct rebind
            { 
                typedef DirectSum<typename Blackbox::template rebind<_Tp1>::other, typename Blackbox::template rebind<_Tp1>::other> other; 

                void operator() (other *& Ap, const Self_t& A, const _Tp1& F) {
                    std::vector<typename other::Blackbox_t *> newPtrV;                  
                    typename std::vector<typename other::Blackbox_t *>::iterator np;
                    typename std::vector<const Blackbox_t* >::const_iterator bp;
                    for( bp = A._VB.begin(), np = newPtrV.begin(); 
                         bp != A._VB.end(); ++bp, ++np) {
                        typename Blackbox_t::template rebind<_Tp1> () (*np, *(*bp), F);
                    }
                    Ap = new other(newPtrV);
                }
            };

               template<class OutVector, class InVector>
                OutVector& apply (OutVector& y, const InVector& x) const
                {
                        linbox_check(y.size() == rowdim());
                        linbox_check(x.size() == coldim());
			int offset_x = 0;
                        int offset_y = 0;
                        typename std::vector<const Blackbox* >::const_iterator bp;
                        for(bp = _VB.begin(); bp != _VB.end(); ++bp){
                                const Subvector<typename InVector::const_iterator> x1(x.begin() + offset_x, x.begin() + (offset_x + (*bp)->coldim()));
                                Subvector<typename OutVector::iterator> y1(y.begin() + offset_y, y.begin() + (offset_y + (*bp)->rowdim()));
                                (*bp) -> apply(y1,x1);
                                offset_x += (*bp) -> coldim();
                                offset_y += (*bp) -> rowdim();
                        }
                        return y;

                }


                template<class OutVector, class InVector>
                OutVector& applyTranspose (OutVector& y, const InVector& x) const
                {
                        linbox_check(y.size() == coldim());
                        linbox_check(x.size() == rowdim());
			int offset_x = 0;
			int offset_y = 0;
			typename std::vector<const Blackbox* >::const_iterator bp;
			for(bp = _VB.begin(); bp != _VB.end(); ++bp){
                        	const Subvector<typename InVector::const_iterator> x1(x.begin() + offset_x, x.begin() + (offset_x + (*bp)->rowdim()));
                        	Subvector<typename OutVector::iterator> y1(y.begin() + offset_y, y.begin() + (offset_y + (*bp)->coldim()));
                        	(*bp) -> applyTranspose(y1,x1);
				offset_x += (*bp) -> rowdim();
				offset_y += (*bp) -> coldim();
			}
                        return y;
                }

                inline size_t rowdim (void) const
                {

                        return m;

                }

                inline size_t coldim(void) const
                {
                        return n;
                }

                const Field& field() const { return _VB.front() -> field();}
            protected:
		std::vector<const Blackbox* > _VB;
		size_t m;
		size_t n;
        }; 

		template <class Matrix> struct MatrixTraits;
		template<>
		template <class BB1, class BB2>
		struct MatrixTraits< DirectSum<BB1, BB2> >
		{
    		typedef DirectSum<BB1, BB2> MatrixType;
	    	typedef MatrixCategories::BlackboxTag MatrixCategory;
		};
		                                                                                                
																										template<>
		template <class BB1, class BB2>
		struct MatrixTraits< const DirectSum<BB1, BB2> >
		{
			typedef const DirectSum<BB1, BB2> MatrixType;
			typedef MatrixCategories::BlackboxTag MatrixCategory;
		};
																												                                                                                                

}; // namespace LinBox

#endif // __DIRECT_SUM_H
