/* linbox/blackbox/sum.h
 * Copyright (C) 2002 The LinBox group
 *
 * Time-stamp: <06 Jul 10 18:38:02 Jean-Guillaume.Dumas@imag.fr> 
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __LINBOX_sum_H
#define __LINBOX_sum_H

#include "linbox/vector/vector-domain.h"
#include "linbox/util/debug.h"
#include <linbox/blackbox/blackbox-interface.h>

namespace LinBox
{
    template <class _Blackbox1, class _Blackbox2 = _Blackbox1>
    class Sum;
    
    template <class _Blackbox1, class _Blackbox2 = _Blackbox1>
    class SumOwner;
}


// Namespace in which all LinBox library code resides
namespace LinBox
{

	/** \brief blackbox of a matrix sum without copying.  

\ingroup blackbox
         * Adds only at apply time.
	 * Given two black boxes A and B of the same dimensions, form a black
	 * box representing A+B, i.e., Sum(A,B)x=(A+B)x=Ax+Bx
	 * @param Vector \ref{LinBox} dense or sparse vector of field elements
	 */
	template <class _Blackbox1, class _Blackbox2>
	class Sum : public BlackboxInterface 
	{
                typedef Sum<_Blackbox1, _Blackbox2> Self_t;
	    public:
		typedef _Blackbox1 Blackbox1;
		typedef _Blackbox2 Blackbox2;

		typedef typename Blackbox1::Field Field;
		typedef typename Blackbox1::Element Element;


		/** Constructor from black box matrices.
		 * This constructor creates a matrix that is the sum,
		 * A + B, of black box matrices.
		 * @param A, B:  black box matrices.
		 */
		Sum (const Blackbox1 &A, const Blackbox2 &B)
			: _A_ptr(&A), _B_ptr(&B), VD( field() )
		{
			linbox_check (A.coldim () == B.coldim ());
			linbox_check (A.rowdim () == B.rowdim ());

			VectorWrapper::ensureDim (_z1, A.rowdim ());
			VectorWrapper::ensureDim (_z2, A.coldim ());
		}

		/** Constructor from black box pointers.
		 * This constructor creates a matrix that is the sum,
		 * A + B, of black box matrices.
		 * @param A_ptr, B_ptr:  pointers to black box matrices.
		 */
		Sum (const Blackbox1 *A_ptr, const Blackbox2 *B_ptr)
			: _A_ptr(A_ptr), _B_ptr(B_ptr), VD( field() )
		{
			// create new copies of matrices in dynamic memory
			linbox_check (A_ptr != 0);
			linbox_check (B_ptr != 0);
			linbox_check (A_ptr->coldim () == B_ptr->coldim ());
			linbox_check (A_ptr->rowdim () == B_ptr->rowdim ());

			VectorWrapper::ensureDim (_z1, A_ptr->rowdim ());
			VectorWrapper::ensureDim (_z2, A_ptr->coldim ());
		}

		/** Copy constructor.
		 * Creates new black box objects in dynamic memory.
		 * @param M constant reference to compose black box matrix
		 */
		Sum (const Sum<Blackbox1, Blackbox2> &M)
			: _A_ptr (M._A_ptr), _B_ptr (M._B_ptr), VD(M.VD)
		{
			VectorWrapper::ensureDim (_z1, _A_ptr->rowdim ());
			VectorWrapper::ensureDim (_z2, _A_ptr->coldim ());
		}

		/// Destructor
		~Sum (void)
		{
		}

		/** Application of BlackBox matrix.
		 * y= (A+B)*x.
		 * Requires one vector conforming to the \ref{LinBox}
		 * vector {@link Archetypes archetype}.
		 * Required by abstract base class.
		 * @return reference to vector y containing output.
		 * @param  x constant reference to vector to contain input
		 */
		template<class OutVector, class InVector>
		inline OutVector &apply (OutVector &y, const InVector &x) const
		{
                    _A_ptr->apply (y, x);
                    _B_ptr->apply (_z1, x);
                    VD.addin (y, _z1);

                    return y;
		}

		/** Application of BlackBox matrix transpose.
		 * y= transpose(A+B)*x.
		 * Requires one vector conforming to the \ref{LinBox}
		 * vector {@link Archetypes archetype}.
		 * Required by abstract base class.
		 * @return reference to vector y containing output.
		 * @param  x constant reference to vector to contain input
		 */
		template<class OutVector, class InVector>
		inline OutVector &applyTranspose (OutVector &y, const InVector &x) const
		{
                    _A_ptr->applyTranspose (y, x);
                    _B_ptr->applyTranspose (_z2, x);
                    VD.addin (y, _z2);

			return y;
		}

            template<typename _Tp1, typename _Tp2 = _Tp1> 
            struct rebind                           
            { typedef SumOwner<
                  typename Blackbox1::template rebind<_Tp1>::other, 
                  typename Blackbox2::template rebind<_Tp2>::other
              > other; 

    		void operator() (other & Ap, const Self_t& A, const _Tp1& F) {
                    typename Blackbox1::template rebind<_Tp1> () ( Ap.getLeftData(), *(A.getLeftPtr()), F);
                    typename Blackbox2::template rebind<_Tp2> () ( Ap.getRightData(), *(A.getRightPtr()), F);
                }


            };
                /** Retreive row dimensions of BlackBox matrix.
		 * This may be needed for applying preconditioners.
		 * Required by abstract base class.
		 * @return integer number of rows of black box matrix.
		 */
		size_t rowdim (void) const
			{ return _A_ptr->rowdim (); }
    
		/** Retreive column dimensions of BlackBox matrix.
		 * Required by abstract base class.
		 * @return integer number of columns of black box matrix.
		 */
		size_t coldim (void) const 
			{ return _A_ptr->coldim (); }


		const Field& field() const { return _A_ptr -> field(); }

		const Blackbox1* getLeftPtr() const {return  _A_ptr;}
		
	        const Blackbox2* getRightPtr() const {return  _B_ptr;}

	    protected:

		// use a copy of the input field for faster performance (no pointer dereference).

		const Blackbox1       *_A_ptr;
		const Blackbox2       *_B_ptr;

		mutable std::vector<Element>  _z1;
		mutable std::vector<Element>  _z2;

    		VectorDomain<Field> VD;
	}; // template <Field, Vector> class Sum

} // namespace LinBox

// Namespace in which all LinBox library code resides
namespace LinBox
{

	/** \brief blackbox of a matrix sum without copying.  

\ingroup blackbox
         * Adds only at apply time.
	 * Given two black boxes A and B of the same dimensions, form a black
	 * box representing A+B, i.e., SumOwner(A,B)x=(A+B)x=Ax+Bx
	 * @param Vector \ref{LinBox} dense or sparse vector of field elements
	 */
	template <class _Blackbox1, class _Blackbox2>
	class SumOwner : public BlackboxInterface 
	{
                typedef SumOwner<_Blackbox1, _Blackbox2> Self_t;
	    public:
		typedef _Blackbox1 Blackbox1;
		typedef _Blackbox2 Blackbox2;

		typedef typename Blackbox1::Field Field;
		typedef typename Blackbox1::Element Element;


		/** Constructor from black box matrices.
		 * This constructor creates a matrix that is the sum,
		 * A + B, of black box matrices.
		 * @param A, B:  black box matrices.
		 */
		SumOwner (const Blackbox1 &A, const Blackbox2 &B)
			: _A_data(&A), _B_data(&B), VD( field() )
		{
			linbox_check (A.coldim () == B.coldim ());
			linbox_check (A.rowdim () == B.rowdim ());

			VectorWrapper::ensureDim (_z1, A.rowdim ());
			VectorWrapper::ensureDim (_z2, A.coldim ());
		}

		/** Constructor from black box pointers.
		 * This constructor creates a matrix that is the sum,
		 * A + B, of black box matrices.
		 * @param A_data, B_data:  pointers to black box matrices.
		 */
		SumOwner (const Blackbox1 *A_data, const Blackbox2 *B_data)
			: _A_data(A_data), _B_data(B_data), VD( field() )
		{
			// create new copies of matrices in dynamic memory
			linbox_check (A_data != 0);
			linbox_check (B_data != 0);
			linbox_check (A_data.coldim () == B_data.coldim ());
			linbox_check (A_data.rowdim () == B_data.rowdim ());

			VectorWrapper::ensureDim (_z1, A_data.rowdim ());
			VectorWrapper::ensureDim (_z2, A_data.coldim ());
		}

		/** Copy constructor.
		 * Creates new black box objects in dynamic memory.
		 * @param M constant reference to compose black box matrix
		 */
		SumOwner (const SumOwner<Blackbox1, Blackbox2> &M)
			: _A_data (M._A_data), _B_data (M._B_data), VD(M.VD)
		{
			VectorWrapper::ensureDim (_z1, _A_data.rowdim ());
			VectorWrapper::ensureDim (_z2, _A_data.coldim ());
		}

		/// Destructor
		~SumOwner (void)
		{
		}

		/** Application of BlackBox matrix.
		 * y= (A+B)*x.
		 * Requires one vector conforming to the \ref{LinBox}
		 * vector {@link Archetypes archetype}.
		 * Required by abstract base class.
		 * @return reference to vector y containing output.
		 * @param  x constant reference to vector to contain input
		 */
		template<class OutVector, class InVector>
		inline OutVector &apply (OutVector &y, const InVector &x) const
		{
                    _A_data.apply (y, x);
                    _B_data.apply (_z1, x);
                    VD.addin (y, _z1);
                    return y;
		}

		/** Application of BlackBox matrix transpose.
		 * y= transpose(A+B)*x.
		 * Requires one vector conforming to the \ref{LinBox}
		 * vector {@link Archetypes archetype}.
		 * Required by abstract base class.
		 * @return reference to vector y containing output.
		 * @param  x constant reference to vector to contain input
		 */
		template<class OutVector, class InVector>
		inline OutVector &applyTranspose (OutVector &y, const InVector &x) const
		{
                    _A_data.applyTranspose (y, x);
                    _B_data.applyTranspose (_z2, x);
                    VD.addin (y, _z2);

			return y;
		}

            template<typename _Tp1, typename _Tp2 = _Tp1> 
            struct rebind                           
            { typedef SumOwner<typename Blackbox1::template rebind<_Tp1>::other, typename Blackbox2::template rebind<_Tp2>::other> other; 

    		void operator() (other & Ap, const Self_t& A, const _Tp1& F) {
                    typename Blackbox1::template rebind<_Tp1> () ( Ap.getLeftData(), A.getLeftData(), F);
                    typename Blackbox2::template rebind<_Tp1> () ( Ap.getRightData(), A.getRightData(), F);
                 }


            };
            template<typename _BBt1, typename _BBt2, typename Field>
            SumOwner (const Sum<_BBt1, _BBt2> &M, const Field& F)
                    : _A_data(*(M.getLeftPtr()), F),
                      _B_data(*(M.getRightPtr()), F),
                      _z1(_A_data.rowdim()),
                      _z2(_A_data.coldim()),
                      VD(F)
                {
                    typename Sum<_BBt1, _BBt2>::template rebind<Field>()(*this,M,F);
                }

            template<typename _BBt1, typename _BBt2, typename Field>
            SumOwner (const SumOwner<_BBt1, _BBt2> &M, const Field& F)
                    : _A_data(M.getLeftData(), F),
                      _B_data(M.getRightData(), F) ,
                      _z1(_A_data.rowdim()),
                      _z2(_A_data.coldim()) ,
                      VD(F)     
            	{
                    typename SumOwner<_BBt1, _BBt2>::template rebind<Field>()(*this,M,F);
                }



	
		/** Retreive row dimensions of BlackBox matrix.
		 * This may be needed for applying preconditioners.
		 * Required by abstract base class.
		 * @return integer number of rows of black box matrix.
		 */
		size_t rowdim (void) const
			{ return _A_data.rowdim (); }
    
		/** Retreive column dimensions of BlackBox matrix.
		 * Required by abstract base class.
		 * @return integer number of columns of black box matrix.
		 */
		size_t coldim (void) const 
			{ return _A_data.coldim (); }


		const Field& field() const { return _A_data . field(); }

		// accessors to the blackboxes without ownership
		const Blackbox1& getLeftData() const {return  _A_data;}
		Blackbox1& getLeftData() {return  _A_data;}
		
	        const Blackbox2& getRightData() const {return  _B_data;}
	        Blackbox2& getRightData() {return  _B_data;}

	    protected:

		// use a copy of the input field for faster performance (no pointer dereference).

		Blackbox1       _A_data;
		Blackbox2       _B_data;

		mutable std::vector<Element>  _z1;
		mutable std::vector<Element>  _z2;

    		VectorDomain<Field> VD;
	}; // template <Field, Vector> class SumOwner

} // namespace LinBox

#endif // __LINBOX_sum_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
