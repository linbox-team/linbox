/* linbox/algorithms/sigma-basis.h
 * Copyright (C) 2005,2013 Pascal Giorgi, Romain Lebreton
 *
 * Written by Pascal Giorgi pascal.giorgi@lirmm.fr
 *            Romain Lebreton lebreton@lirmm.fr
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

#include "linbox/matrix/polynomial-matrix.h"
#include "linbox/algorithms/polynomial-matrix/polynomial-matrix-domain.h"
#include <vector>
#include <algorithm>
using namespace std;

#define MBASIS_THRESHOLD_LOG 5
#define MBASIS_THRESHOLD (1<<MBASIS_THRESHOLD_LOG)

namespace LinBox {
        
        template< size_t K>
        struct EarlyTerm {
                size_t _count;
                size_t _val;
                
                EarlyTerm():  _count(0),_val(0){}
                
                void update(size_t r, const vector<size_t>& u){
                        vector<size_t> v(u);
                        sort(v.begin(),v.end());
                        size_t x=0;
                        for (size_t i=0;i<r;i++)
                                x+=v[i];
                        if (x==_val)
                                _count++;
                        else{
                                _val=x;
                                _count=0;
                        }                    
                }
                
                bool terminated() const {return _count>=K;}
                
                void reset() {_count=0;_val=0;}
        };

        template<class Field, class ET=EarlyTerm<4> >
        class OrderBasis {
        public:
                typedef PolynomialMatrix<PMType::polfirst,PMStorage::plain,Field> MatrixP;
                typedef PolynomialMatrix<PMType::matfirst,PMStorage::plain,Field> PMatrix;
        private:
                const Field*                     _field;
                PolynomialMatrixMulDomain<Field>   _PMD;
                BlasMatrixDomain<Field>            _BMD;
                ET                           _EarlyStop;                  
        public:
                double mul;
                Timer tmul;

                OrderBasis(const Field& f) : _field(&f), _PMD(f), _BMD(f) {
                        mul=0.;
                }

                inline const Field& field() const {return *_field;}

                
                // serie must have exactly order elements (i.e. its degree = order-1)
                // sigma can have at most order+1 elements (i.e. its degree = order)
                template<typename PMatrix1, typename PMatrix2>
                size_t PM_Basis(PMatrix1             &sigma,
                                const PMatrix2       &serie,
                                size_t                order,
                                vector<size_t>       &shift)
                {
                        if (order <= MBASIS_THRESHOLD) {                        
                                return M_Basis(sigma, serie, order, shift);
                        }
                        else {                                
                                size_t ord1,ord2,d1,d2;
                                ord1 = order>>1; 
                                ord2 = order-ord1; // ord1+ord2=order                        
                                size_t m,n,k;
                                m=sigma.rowdim();
                                n=sigma.coldim();
                                k=serie.coldim();

                                PMatrix sigma1(field(),m,n,ord1+1);
                                //typename PMatrix1::view       sigma1=sigma.at(0,ord1);
                                //PMatrix serie1(field(),n,k,ord1); serie1.copy(serie,0,ord1-1);
                                typename PMatrix2::const_view serie1=serie.at(0,ord1-1);
                                PMatrix sigma2(field(),m,n,ord2+1);

                                // first recursive call
                                d1 = PM_Basis(sigma1, serie1, ord1, shift);

                                // compute the serie update
                                //PMatrix serie2(field(),n,k,ord2+(ord2==ord1));//serie2 size=ord1+1 -> midproduct)
                                PMatrix serie2(field(),n,k,ord2);//serie2 size=ord1+1 -> midproduct)
                                // TODO: for Block Wiedemann, this step can use only the first column of sigma
                                _PMD.midproductgen(serie2, sigma1, serie, true, ord1+1,ord1+ord2);
                                // cout<<"updated serie:"<<endl<<serie<<endl;

                                // second recursive call
                                d2 = PM_Basis(sigma2, serie2, ord2, shift);

                                // compute the result
                                _PMD.mul(sigma, sigma2, sigma1);

                                sigma.resize(d1+d2+1);

                                return d1+d2;
                        }
                }

                // serie must have exactly order elements (i.e. its degree = order-1)
                template<typename PMatrix1, typename PMatrix2>
                size_t M_Basis(PMatrix1          &sigma,
                               const PMatrix2    &serie,
                               size_t            order,
                               vector<size_t>   &shift)
                {        
                        //cout<<"------------- mba : "<<order<<endl;
                        //cout<<serie<<endl;
                        size_t m=serie.rowdim();
                        size_t n=serie.coldim();
                        size_t rank=0;
                        BlasMatrix<Field> delta(field(),m,n);
                        size_t max_degree=0;        // a bound on the row degree of sigma
                        vector<size_t> degree(m,0); // a bound on each row degree of sigma

                        // set sigma to identity
                        for(size_t i=0;i<m*m;i++)
                                sigma.ref(i,0)=0;
                        for(size_t i=0;i<m;i++)
                                sigma.ref(i,i,0)=1;

                        BlasPermutation<size_t> Qt(m), P(n);
                        typedef BlasSubmatrix<BlasMatrix<Field> > View;

                        for (size_t k=0; k<order && !_EarlyStop.terminated(); k++){                           
                                // sort the shift in ascending order (-> minimize the shifted row-degree of sigma)
                                // -> store the permutation in Bperm
                                // -> permute the row degree at the same time
                                vector<size_t> perm(m);
                                for (size_t i=0;i<m;i++) perm[i]=i;
                                for (size_t i=0;i<m;++i) {
                                        size_t idx_min=i;
                                        for (size_t j=i+1;j<m;++j)
                                                if (shift[j]< shift[idx_min])
                                                        idx_min=j;
                                        std::swap( shift[i],  shift[idx_min]);
                                        std::swap(degree[i], degree[idx_min]);
                                        perm[i]=idx_min;
                                }                        
                                BlasPermutation<size_t> Bperm(perm);

                                // check if Qt is identity
                                size_t lll=0;
                                while(lll<m && Qt[lll]==lll) lll++;

                                // permute the row of the current sigma basis
                                // and compute the new discrepancy
                                //if (true){
                                if (k==0 || lll!=m){  
                                        _BMD.mulin_right(Bperm, sigma[0]);
                                        _BMD.mul(delta,sigma[0],serie[k]);                        
                                        for(size_t i=1;i<=min(k,max_degree);i++){
                                                _BMD.mulin_right(Bperm, sigma[i]);
                                                _BMD.axpyin(delta,sigma[i],serie[k-i]);
                                        }
                                }
                                else{         
                                        // if Qt is identity then the first rank rows of delta remain unchanged
                                        View delta1(delta,   rank,0,m-rank,n);
                                        View sigma1(sigma[0],rank,0,m-rank,m);
                                        _BMD.mul(delta1,sigma1,serie[k]);                        
                                        _BMD.mulin_right(Bperm, sigma[0]);
                                        for(size_t i=1;i<=min(k,max_degree);i++){
                                                View sigmak(sigma[i],rank,0,m-rank,m);
                                                _BMD.axpyin(delta1,sigmak,serie[k-i]);  
                                                _BMD.mulin_right(Bperm, sigma[i]);
                                        }
                                        _BMD.mulin_right(Bperm, delta);                                
                                }

                                //delta.write(cout,Tag::FileFormat::Plain);
                                //cout<<endl;

                                // Compute LQUP of Discrepancy
                                BlasMatrix<Field> delta_copy(delta);
                                LQUPMatrix<Field> LQUP(delta_copy,P,Qt);

                                // Get L from LQUP
                                TriangularBlasMatrix<Field> L(field(), m, m, Tag::Shape::Lower, Tag::Diag::Unit);
                                LQUP.getL(L);
                                rank=LQUP.getRank();  // the first rank entries of Qt give the pivot row

                                // inverse L in-place (COULD BE IMPROVED -> only compute the left kernel by trsm)
                                //FFPACK::ftrtri(field(),FFLAS::FflasLower,FFLAS::FflasUnit,m,L.getPointer(),L.getStride());

                                View L1(L,0   ,0,  rank,rank);
                                View L2(L,rank,0,m-rank,rank);
                                FFLAS::ftrsm(field(),FFLAS::FflasRight,FFLAS::FflasLower,
                                             FFLAS::FflasNoTrans,FFLAS::FflasUnit, 
                                             m-rank,rank, field().mOne, L1.getPointer(),m, L2.getWritePointer(),m);

                                // update sigma by L^(-1) (rank sensitive -> use only the left kernel basis)
                                for(size_t i=0;i<=max_degree;i++){
                                        // NEED TO APPLY Qt to sigma[i]
                                        View S1(sigma[i],0,0,rank,m);
                                        View S2(sigma[i],rank,0,m-rank,m);
                                        _BMD.axpyin(S2,L2,S1); 
                                        //_BMD.mulin_right(L,sigma[i]);

                                }

                                size_t dmax=0, smax=0;                                                
                                // update: the row-degree, the shifted row-degree,
                                //         the max pivot degree and the maximum row degree
                                for (size_t i=0;i<rank;++i) {
                                        dmax=max(dmax, degree[Qt[i]]);
                                        smax=max(smax, shift [Qt[i]]);
                                        degree[Qt[i]]++;
                                        shift [Qt[i]]++;                                                        
                                }                        
                                max_degree=max(max_degree,dmax+1);
                                for (size_t i=rank;i<m;i++){
                                        degree[Qt[i]]=max(dmax, degree[Qt[i]]);
                                        //THIS LINE IS NOT NEEDED -> shift [Qt[i]]=max(smax, shift [Qt[i]]);
                                }

                                // shift the pivot row of sigma by x 
                                for (int l=max_degree-1;l>=0;l--)
                                        for (size_t i=0;i<rank;i++)
                                                for (size_t j=0;j<m;j++)
                                                        sigma.ref(Qt[i],j,l+1)=sigma.ref(Qt[i],j,l);
                                for (size_t i=0;i<rank;i++)
                                        for (size_t j=0;j<m;j++)
                                                sigma.ref(Qt[i],j,0)=field().zero;
                                //cout<<"max degree="<<max_degree<<endl<<endl;
                                //cout<<sigma<<endl<<"******************"<<endl;

                                // update Early Termination
                                _EarlyStop.update(rank,shift);
                        }
                        sigma.resize(max_degree+1);
                        return max_degree;
                }



                inline size_t twoValuation(size_t x){
                        size_t i=0;
                        while (x!=0 && !(x&0x1)){
                                i++;
                                x>>=1;
                        }
                        return i;
                }

                template<class Polynomial1>
                void update_sigma(size_t m, size_t n, list<Polynomial1*>& L, size_t k){
                        Polynomial1 *P1,*P2,*P3;
                        for(size_t i=0;i<k;i++){
                                P2=L.back();L.pop_back();
                                P3=L.back();L.pop_back();
                                P1 = new Polynomial1(field(),m,n,P2->size()+P3->size()-1);
                                _PMD.mul(*P1,*P2,*P3);
                                L.push_back(P1);
                                delete P2; delete P3;
                        }
                }

                // serie must have exactly order elements (i.e. its degree = order-1)
                // sigma can have at most order+1 elements (i.e. its degree = order)
                // Algorithm from [Giorgi, Lebreton ISSAC'2014]
                template<typename PMatrix1, typename PMatrix2>
                void oPM_Basis(PMatrix1             &sigma,
                               const PMatrix2       &serie,
                               size_t                order,
                               vector<size_t>       &shift)
                {
                        size_t m,n,k,l,lp;
                        m=sigma.rowdim();
                        n=sigma.coldim();
                        k=serie.coldim();

                        // log of the order
                        size_t log_order=integer(order).bitsize();

                        //  leaf size of the recursive PM_Basis algorithm (must be a power of 2) 
                        size_t log_ord = MBASIS_THRESHOLD_LOG;
                        size_t ord     = min(1UL<<log_ord ,order);

                        // prepare the storage for each serie update
                        vector<PMatrix2*> L_serie(log_order+1);
                        for (size_t i=log_ord;i<log_order;i++)
                                L_serie[i]= new PMatrix2(field(),n,k,(1<<i));
                        L_serie[log_order]=const_cast<PMatrix2*>(&serie);

                        typedef typename PMatrix2::const_view cview;
                        list<PMatrix1*>  L_sigma;                
                        PMatrix1*         sigmak; 
                        cview             seriek;
                        typedef HalflineMPDomain<Field,PMatrix2,PMatrix1,PMatrix2> HFMPD;
                        list<HFMPD*> L_mp;
                        typename list<HFMPD*>::iterator iter, t_iter;
                        
                        // // Reset Early Termination
                        _EarlyStop.reset();
                        
                        sigmak = new PMatrix1(field(),m,n,ord+1);
                        seriek = serie.at(0,ord-1);
                        M_Basis(*sigmak, seriek, ord, shift);                  
                        L_sigma.push_back(sigmak);
                        size_t sss=0;
                        for(size_t k=ord;k<order &&  !_EarlyStop.terminated();k+=ord,sss++){  
                                // cout<<"------------ order="<<k<<endl;
                                l  = twoValuation(k);
                                lp = twoValuation(k-(1<<l)); lp=(lp==0?log_order:lp); 
                                // compute next element in the original serie and 
                                // update all subsequent computed series
                                for(iter=L_mp.begin(); iter!=L_mp.end(); ){
                                        (*iter)->update(ord);
                                        t_iter=iter;
                                        ++iter;    
                                
                                        if ((*t_iter)->terminated()){
                                                delete *t_iter; 
                                                L_mp.erase(t_iter);  
                                        }                                  
                                }

                                // compute the serie update
                                //seriek = const_cast<const PMatrix2*>(L_serie[lp])->at(0,min(order,(1UL<<(l+1)))-1);

                                //size_t update_max=min(order,1UL<<(l+1));
                                //_PMD.midproductgen(*L_serie[l],*L_sigma.back(), seriek, true, (1UL<<l)+1,update_max);

                        
                                /*
                                  cout<<"---------------"<<endl;
                                  cout<<"MP "<<(1UL<<l)<<"x"<<(1UL<<(l+1))<<endl;
                                  cout<<*(L_sigma.back())<<endl;
                                  cout<<lp<<" ----- "<<endl<<*(L_serie[lp])<<endl;
                                  cout<<"---------------"<<endl;
                                */

                                //L_mp.push_back(new HFMPD(field(),*(L_serie[l]),*(L_sigma.back()), seriek, 1UL<<l));
                                L_mp.push_back(new HFMPD(field(),*(L_serie[l]),*(L_sigma.back()),*(L_serie[lp]), 1UL<<l));
                   
                                // compute the new sigma
                                sigmak = new PMatrix(field(),m,n,ord+1);
                                size_t step= min(ord,order-k); // needed if order%ord <> 0
                                seriek = const_cast<const PMatrix2*>(L_serie[l])->at(0,step-1);
                        
                                // compute the next "step" elements of the serie L_serie[l]
                                L_mp.back()->update(ord);
                        
                                M_Basis(*sigmak, seriek, step, shift);                        
                                L_sigma.push_back(sigmak);

                                // update the sigma list
                                update_sigma(m, n, L_sigma, twoValuation(k/ord+1));
                        }
                        // get the product of all sigma                
                        if (L_sigma.size()>1){
                                update_sigma(m, n, L_sigma, L_sigma.size()-2);                
                                PMatrix1 *s1,*s2;
                                s1= L_sigma.back();L_sigma.pop_back();
                                s2= L_sigma.back();L_sigma.pop_back();
                                _PMD.mul(sigma,*s1,*s2);
                                delete s1; delete s2;
                        }
                        else {
                                sigma.resize(L_sigma.back()->size());
                                sigma.copy(*L_sigma.back(),0,L_sigma.back()->size()-1);
                                delete L_sigma.back();
                        }
                        for (size_t i=0;i<log_order;i++)
                                delete L_serie[i];

                        // Info about early termination
                        //if (_EarlyStop.terminated())
                        //        cout<<"Early termination at order "<<sss<<" ("<<order<<")"<<endl;
                }
        
        };

} // end of namespace LinBox

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
