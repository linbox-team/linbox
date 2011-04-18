/* Copyright (C)  LinBox
 * Written by 
 *            David Saunders
 *
 *
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

#ifndef __LINBOX_direct_sum_H
#define __LINBOX_direct_sum_H

#include <linbox/blackbox/null-matrix.h>
#include <linbox/vector/vector-traits.h>
#include <linbox/blackbox/blackbox-interface.h>
#include <linbox/vector/subvector.h>
#include <linbox/matrix/matrix-domain.h>
#include <linbox/vector/light_container.h>

namespace LinBox
{

    template <class Blackbox1, class Blackbox2 = Blackbox1>
    class DirectSum;

    template <class Blackbox1, class Blackbox2 = Blackbox1>
    class DirectSumOwner;
    
}


namespace LinBox
{

	/** \brief If C = DirectSum(A, B) and y = xA and z = wB, then (y,z) = (x,w)C.

        * And similarly for apply. 
	\ingroup blackbox
        */
    template <class _Blackbox1, class _Blackbox2>
    class DirectSum : public BlackboxInterface
    {
        typedef DirectSum<_Blackbox1, _Blackbox2> Self_t;
    public:
        typedef _Blackbox1 Blackbox1;
        typedef _Blackbox2 Blackbox2;

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
                _Ap->apply(y1,x1);
                _Bp->apply(y2,x2);
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
                _Ap->applyTranspose(y1,x1);
                _Bp->applyTranspose(y2,x2);

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
            typedef DirectSumOwner<
                typename Blackbox1::template rebind<_Tp1>::other,
                typename Blackbox2::template rebind<_Tp2>::other
            > other; 

            void operator() (other & Ap, const Self_t& A, const _Tp1& F) {
                typename Blackbox1::template rebind<_Tp1> () ( Ap.getLeftData(), *(A.getLeftPtr()), F);
                typename Blackbox2::template rebind<_Tp2> () ( Ap.getRightData(), *(A.getRightPtr()), F);
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


        const Field& field() const { return _Ap->field();}

            // accessors to the blackboxes
        const Blackbox1* getLeftPtr() const {return  _Ap;}
        const Blackbox2* getRightPtr() const {return  _Bp;}

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
        typedef std::vector<const Blackbox* > ListBB_t;


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
            m = Ap->rowdim() + Bp->rowdim();
            n = Ap->coldim() + Bp->coldim();
        }
        template<class Vector>
        DirectSum(const Vector& v) :_VB(v.begin(),v.end()) {
            m = 0; n = 0;
            typename std::vector<const Blackbox* >::iterator bp;
            for( bp = _VB.begin(); bp != _VB.end(); ++bp) {
                m += (*bp)->rowdim();
                n += (*bp)->coldim();
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
            typedef DirectSumOwner<
                typename Blackbox::template rebind<_Tp1>::other, 
                typename Blackbox::template rebind<_Tp1>::other
            > other; 


            void operator() (other & Ap, const Self_t& A, const _Tp1& F) {
                typename other::ListBB_t::iterator itp = Ap.getDataSum().begin();
                typename Self_t::ListBB_t::const_iterator it = A.getSum().begin();
                  
                for( ; it != A.getSum().end(); ++itp,++it) 
                    typename Blackbox::template rebind<_Tp1>()( *itp, *(*it), F);
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
                    (*bp)->apply(y1,x1);
                    offset_x += (*bp)->coldim();
                    offset_y += (*bp)->rowdim();
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
                    (*bp)->applyTranspose(y1,x1);
                    offset_x += (*bp)->rowdim();
                    offset_y += (*bp)->coldim();
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

        const Field& field() const { return _VB.front()->field();}

        std::vector<const Blackbox* >& getSum() { return _VB; }
        const std::vector<const Blackbox* >& getSum() const { return _VB; }

        size_t size() const {
            return this->_VB.size();
        }
    protected:
        std::vector<const Blackbox* > _VB;
        size_t m;
        size_t n;
    }; 



    template <class Matrix> struct MatrixTraits;

#ifndef __INTEL_COMPILER
    template<>
#endif
    template <class BB1, class BB2>
    struct MatrixTraits< DirectSum<BB1, BB2> >
    {
        typedef DirectSum<BB1, BB2> MatrixType;
        typedef MatrixCategories::BlackboxTag MatrixCategory;
    };
		                                                                                                
#ifndef __INTEL_COMPILER
    template<>
#endif
    template <class BB1, class BB2>
    struct MatrixTraits< const DirectSum<BB1, BB2> >
    {
        typedef const DirectSum<BB1, BB2> MatrixType;
        typedef MatrixCategories::BlackboxTag MatrixCategory;
    };
																												                                                                                                

}; // namespace LinBox


namespace LinBox
{

    template <class _Blackbox1, class _Blackbox2>
    class DirectSumOwner : public BlackboxInterface
    {
        typedef DirectSumOwner<_Blackbox1, _Blackbox2> Self_t;
    public:
        typedef _Blackbox1 Blackbox1;
        typedef _Blackbox2 Blackbox2;

        typedef typename Blackbox1::Field Field;
        typedef typename Blackbox1::Element Element;

            /** Constructor from two black box matrices.
             * This becomes direct sum of A and B.
             * They may be rectangular.  
             * @param A, B:  black box matrices over a common field.
             */
        DirectSumOwner(const	Blackbox1& A, const Blackbox2& B)
                : _A_data(A), _B_data(B)
            {}

            /** Constructor from two black box matrix pointers.
             * This becomes direct sum of A and B.
             * They may be rectangular.  They must be over the same field (or ring). 
             * @param A_ptr pointer to black box matrix A.
             * @param B_ptr pointer to black box matrix B.
             */
            // 	DirectSumOwner(): _Ap(& DirectSumOwner<Blackbox1, Blackbox2>::_NullMatrix),
// 			     _Bp(& DirectSumOwner<Blackbox1, Blacbox2>::_NullMatrix) 
// 		{}



        DirectSumOwner(const	Blackbox1* Ap, const Blackbox2* Bp)
                : _A_data(*Ap), _B_data(*Bp)
            {}

            /// Copy constructor.
        DirectSumOwner (const DirectSumOwner<Blackbox1, Blackbox2>& M) 
                : _A_data (M._A_data), _B_data (M._B_data)
            {}


            /// Destructor
        ~DirectSumOwner (void)
            {}

        template<class OutVector, class InVector>
        OutVector& apply (OutVector& y, const InVector& x) const
            {
                linbox_check(x.size() == coldim());
                linbox_check(y.size() == rowdim());
                    // FIXED: I want to use subvectors to avoid copying and memory allocation, but problems with it...

                const Subvector<typename InVector::const_iterator> x1(x.begin(), x.begin() + _A_data.coldim());
                const Subvector<typename InVector::const_iterator> x2(x.begin() + _A_data.coldim(), x.end());
                Subvector<typename OutVector::iterator> y1(y.begin(), y.begin() + _A_data.rowdim());
                Subvector<typename OutVector::iterator> y2(y.begin() + _A_data.rowdim(), y.end());
                _A_data->apply(y1,x1);
                _B_data->apply(y2,x2);
                return y;
            }


        template<class OutVector, class InVector>
        OutVector& applyTranspose (OutVector& y, const InVector& x) const
            {
                linbox_check(x.size() == rowdim());
                linbox_check(y.size() == coldim());
                const Subvector<typename InVector::const_iterator> x1(x.begin(), x.begin() + _A_data.rowdim());
                const Subvector<typename InVector::const_iterator> x2(x.begin() + _A_data.rowdim(), x.end());
                Subvector<typename OutVector::iterator> y1(y.begin(), y.begin() + _A_data.coldim());
                Subvector<typename OutVector::iterator> y2(y.begin() + _A_data.coldim(), y.end());
                _A_data->applyTranspose(y1,x1);
                _B_data->applyTranspose(y2,x2);


                return y;
            }

        template<typename _Tp1, typename _Tp2 = _Tp1>
        struct rebind
        { 
            typedef DirectSumOwner<
                typename Blackbox1::template rebind<_Tp1>::other,
                typename Blackbox2::template rebind<_Tp2>::other
            > other; 

            void operator() (other & Ap, const Self_t& A, const _Tp1& F) {
                typename Blackbox1::template rebind<_Tp1> () ( Ap.getLeftData(), A.getLeftData(), F);
                typename Blackbox2::template rebind<_Tp2> () ( Ap.getRightData(), A.getRightData(), F);
            }

        };


        template<typename _BBt1, typename _BBt2, typename Field>
        DirectSumOwner (const DirectSum<_BBt1, _BBt2> &M, const Field& F)
                : _A_data(*(M.getLeftPtr()), F),
                  _B_data(*(M.getRightPtr()), F) 
            {
                typename DirectSum<_BBt1, _BBt2>::template rebind<Field,Field>()(*this, M, F);
            }

        template<typename _BBt1, typename _BBt2, typename Field>
        DirectSumOwner (const DirectSumOwner<_BBt1, _BBt2> &M, const Field& F)
                : _A_data(M.getLeftData(), F),
                  _B_data(M.getRightData(), F) 
            {
                typename DirectSumOwner<_BBt1, _BBt2>::template rebind<Field,Field>()(*this, M, F);
            }



        inline size_t rowdim (void) const
            {

                return _A_data.rowdim () + _B_data.rowdim ();

            }
    
        inline size_t coldim(void) const 
            {
                return _A_data.coldim () + _B_data.coldim ();
            }


        const Field& field() const { return _A_data->field();}

            // accessors to the blackboxes without ownership
        const Blackbox1& getLeftData() const {return  _A_data;}
		
        const Blackbox2& getRightData() const {return  _B_data;}

    protected:

            // the direct summands
        Blackbox1 _A_data;
        Blackbox2 _B_data;		

    }; // template <Vector> class DirectSumOwner



    template <class Blackbox>
    class DirectSumOwner<Blackbox, Blackbox>
    {
    public:
        typedef DirectSumOwner<Blackbox, Blackbox> Self_t;
        typedef Blackbox Blackbox_t;
        typedef LightContainer<Blackbox> ListBB_t;
        
        typedef typename Blackbox::Field Field;
        typedef typename Blackbox::Element Element;

        DirectSumOwner() : m(0),n(0){}
        DirectSumOwner(const Blackbox& A, const Blackbox& B) {
            _VB_data.push_back(A);
            _VB_data.push_back(B);
            m = A.rowdim() + B.rowdim();
            n = A.coldim() + B.coldim();
        }

        DirectSumOwner(const Blackbox* Ap, const Blackbox* Bp) {
            _VB_data.push_back(*Ap);
            _VB_data.push_back(*Bp);
            m = Ap->rowdim() + Bp->rowdim();
            n = Ap->coldim() + Bp->coldim();
        }
        template<class Vector>
        DirectSumOwner(const Vector& v) :_VB_data(v.begin(),v.end()) {
            m = 0; n = 0;
            typename ListBB_t::iterator bp;
            for( bp = _VB_data.begin(); bp != _VB_data.end(); ++bp) {
                m += bp->rowdim();
                n += bp->coldim();
            }
        }
		
            /// Copy constructor.
        DirectSumOwner (const DirectSumOwner<Blackbox, Blackbox>& M)
                :_VB_data(M._VB_data),m(M.m),n(M.n)
            {}


            /// Destructor
        ~DirectSumOwner (void)
            {}

        template<typename _Tp1>
        struct rebind
        { 
            typedef DirectSumOwner<
                typename Blackbox::template rebind<_Tp1>::other, 
                typename Blackbox::template rebind<_Tp1>::other
            > other; 


            void operator() (other & Ap, const Self_t& A, const _Tp1& F) {
                typename other::ListBB_t::iterator itp = Ap.getDataSum().begin();
                typename Self_t::ListBB_t::const_iterator it = A.getDataSum().begin();

                  
                for( ; it != A.getDataSum().end(); ++itp,++it)
                    typename Blackbox::template rebind<_Tp1>()( *(itp->get()), *(*it), F);
            }

        };

        template<typename _BBt, typename Field>
        DirectSumOwner (const DirectSum<_BBt> &M, const Field& F)
                : _VB_data( M.size() ), m( M.rowdim() ), n( M.coldim())
            {
                typename DirectSum<_BBt>::template rebind<Field>()(*this, M, F);
            }

        template<typename _BBt, typename Field>
        DirectSumOwner (const DirectSumOwner<_BBt> &M, const Field& F)
                : _VB_data( M.size() ), m( M.rowdim() ), n( M.coldim())
            {
                typename DirectSumOwner<_BBt>::template rebind<Field>()(*this, M, F);
            }


        template<class OutVector, class InVector>
        OutVector& apply (OutVector& y, const InVector& x) const
            {
                linbox_check(y.size() == rowdim());
                linbox_check(x.size() == coldim());
                int offset_x = 0;
                int offset_y = 0;
                typename ListBB_t::const_iterator bp;
                for(bp = _VB_data.begin(); bp != _VB_data.end(); ++bp){
                    const Subvector<typename InVector::const_iterator> x1(x.begin() + offset_x, x.begin() + (offset_x + bp->coldim()));
                    Subvector<typename OutVector::iterator> y1(y.begin() + offset_y, y.begin() + (offset_y + bp->rowdim()));
                    bp->apply(y1,x1);
                    offset_x += bp->coldim();
                    offset_y += bp->rowdim();
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
                typename ListBB_t::const_iterator bp;
                for(bp = _VB_data.begin(); bp != _VB_data.end(); ++bp){
                    const Subvector<typename InVector::const_iterator> x1(x.begin() + offset_x, x.begin() + (offset_x + bp->rowdim()));
                    Subvector<typename OutVector::iterator> y1(y.begin() + offset_y, y.begin() + (offset_y + bp->coldim()));
                    bp->applyTranspose(y1,x1);
                    offset_x += bp->rowdim();
                    offset_y += bp->coldim();
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

        const Field& field() const { return _VB_data.front().field();}
        ListBB_t& getDataSum() { return _VB_data; }
        const ListBB_t& getDataSum() const { return _VB_data; }

        size_t size() const {
            return this->_VB_data.size();
        }
            
                    

    protected:
        ListBB_t _VB_data;
        size_t m;
        size_t n;
    }; 


	
}; // namespace LinBox


#endif // __LINBOX_direct_sum_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
