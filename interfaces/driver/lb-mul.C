/* lb-mul.C
 * Copyright (C) 2017 Pascal Giorgi
 *
 * Written by Pascal Giorgi <pascal.giorgi@lirmm.fr>
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 */

#include "linbox/linbox-config.h"
#include "linbox/matrix/matrix-domain.h"

#include "lb-blackbox-function.h"
#include "lb-blackbox.h"

extern BlackboxTable blackbox_hashtable;
extern DomainTable   domain_hashtable;


/***********************
 * Matrix Mul  Functor *
 ***********************/

template< class Field1, class Field2, class Field3, class Mtrait>
class LaunchMulFunctor {
public:
        template <class Blackbox, class MTrait>
        void operator()(Blackbox *C,Blackbox *A, Blackbox *B){
                throw lb_runtime_error("LinBox ERROR: incompatible domain in matrix mul");// throw an exception for incompatible data type
        }
        template <class Blackbox1,class Blackbox2,class Blackbox3>
        void operator()(Blackbox1 *C,Blackbox2 *A, Blackbox3 *B){
                throw lb_runtime_error("LinBox ERROR: incompatible blackbox type in matrix mul");// throw an exception for incompatible data type
        }
};


template< class Field, class MTrait>
class LaunchMulFunctor<Field,Field,Field, MTrait> {
public:
        template <class Blackbox>
        void operator()(Blackbox *C,Blackbox *A, Blackbox *B){
                LinBox::MatrixDomain<Field> BMD(A->field());
                BMD.mul(*C,*A,*B);                
        }
        template <class Blackbox1,class Blackbox2,class Blackbox3>
        void operator()(Blackbox1 *C,Blackbox2 *A, Blackbox3 *B){
                throw lb_runtime_error("LinBox ERROR: incompatible blackbox type in matrix mul");// throw an exception for incompatible data type
        }

};
template< class Field>
class LaunchMulFunctor<Field,Field,Field, LinBox::MatrixContainerCategory::BlasContainer> {
public:
        template <class Blackbox>
        void operator()(Blackbox *C,Blackbox *A, Blackbox *B){
                LinBox::BlasMatrixDomain<Field> BMD(A->field());
                BMD.mul(*C,*A,*B);                
        }
        template <class Blackbox1,class Blackbox2,class Blackbox3>
        void operator()(Blackbox1 *C,Blackbox2 *A, Blackbox3 *B){
                throw lb_runtime_error("LinBox ERROR: incompatible blackbox type in matrix mul");// throw an exception for incompatible data type
        }

};

template<>
class LaunchMulFunctor<Givaro::QField<Givaro::Rational>,Givaro::QField<Givaro::Rational>,Givaro::QField<Givaro::Rational>, LinBox::MatrixContainerCategory::BlasContainer> {
public:
        template <class Blackbox>
        void operator()(Blackbox *C,Blackbox *A, Blackbox *B){
                LinBox::MatrixDomain<Givaro::QField<Givaro::Rational> > BMD(A->field());
                BMD.mul(*C,*A,*B);                        
        }
        template <class Blackbox1,class Blackbox2,class Blackbox3>
        void operator()(Blackbox1 *C,Blackbox2 *A, Blackbox3 *B){
                throw lb_runtime_error("LinBox ERROR: incompatible blackbox type in matrix mul");// throw an exception for incompatible data type
        }

};


template<class Blackbox>
class FinalMatrixMulFunctor{
private:
        Blackbox  *_B;
public:

        FinalMatrixMulFunctor(Blackbox* B) :  _B(B) {}

        template <class Blackbox1,class Blackbox2>
        void operator() (Blackbox1 *C,  Blackbox2 *A) const {
                LaunchMulFunctor<typename Blackbox1::Field,typename Blackbox2::Field,typename Blackbox::Field, typename LinBox::MatrixContainerTrait<Blackbox>::Type>()(C,A,_B);
        }
};

class MatrixMulFunctor{
private:
        const BlackboxKey &_Akey;
public:
        MatrixMulFunctor(const BlackboxKey& k) : _Akey(k) {}

        template<class Blackbox>
        inline void operator()(const BlackboxKey &Bkey, Blackbox *C) const {
                MatrixMulFunctor fct(_Akey);
                BlackboxFunction::call(C, Bkey, fct);
        }
        template<class Blackbox1, class Blackbox2>
        inline void operator()(Blackbox1 *C, Blackbox2 *B) const {
                FinalMatrixMulFunctor<Blackbox2> fct(B);
                BlackboxFunction::call(C, _Akey, fct);
        }
        
};

/*******************************************************
 * API for matrix multiplication                       *
 * matrix result is given through a blackbox key       *
 *******************************************************/
void lb_mul(const BlackboxKey &Ckey, const BlackboxKey &Akey, const BlackboxKey &Bkey){
        MatrixMulFunctor fct(Akey);
        BlackboxFunction::call(Bkey,Ckey,fct);
}

/*******************************************************
 * API for matrix multiplication                       *
 * matrix result is given through a blackbox key       *
 *******************************************************/
const BlackboxKey & lb_mul(const BlackboxKey &Akey, const BlackboxKey &Bkey){
        BlackboxTable::iterator ita = blackbox_hashtable.find(Akey);
        BlackboxTable::iterator itb = blackbox_hashtable.find(Bkey);
        
        if (ita == blackbox_hashtable.end() || itb == blackbox_hashtable.end())
                throw lb_runtime_error("LinBox ERROR: blackbox does not exist (multiplying impossible)\n");

        const DomainKey *Dkeya = &ita->second->getDomainKey();
        const DomainKey *Dkeyb = &itb->second->getDomainKey();
        if (*Dkeya != *Dkeyb)
                throw lb_runtime_error("LinBox ERROR: trying to multiply matrices from different domains\n");
	
	
        std::pair<size_t,size_t> dima, dimb;
        dima = getBlackboxDimension(Akey);
        dimb = getBlackboxDimension(Bkey);
        if (dima.second != dimb.second)
                throw lb_runtime_error("LinBox ERROR: trying to multiply matrices with incompatible sizes\n");

        size_t rowdim = dima.first;
        size_t coldim = dimb.second;

        const BlackboxKey *res = &createBlackbox(*Dkeya, rowdim, coldim);

        lb_mul(*res, Akey, Bkey);
        return *res;
}






// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

