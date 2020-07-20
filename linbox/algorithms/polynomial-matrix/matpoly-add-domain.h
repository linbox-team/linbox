/*
 * Copyright (C) 2013  Pascal Giorgi
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

#ifndef __LINBOX_matpoly_add_domain_H
#define __LINBOX_matpoly_add_domain_H

#include <algorithm>
#include "linbox/matrix/matrix-domain.h"

namespace LinBox {

	template <class Field>
	class PolynomialMatrixAddDomain {
	protected:
		MatrixDomain<Field>      _BMD;

	public:
		PolynomialMatrixAddDomain(const Field& F)
			: _BMD(F) {}

		// add function (c must be allocated with the right size)
		template<typename PMatrix1,typename PMatrix2,typename PMatrix3>
		void add(PMatrix1 &c, const PMatrix2 &a, const PMatrix3 &b) const {
            // std::clog<<"Add matpoly:"<<std::endl;
            // std::clog<<"A="<<a<<std::endl;
            // std::clog<<"B="<<b<<std::endl;
			size_t i=0;
			for(;i<std::min(a.size(),b.size());i++){
                auto c_tmp=c[i];
                _BMD.add(c_tmp,a[i],b[i]);
                c.setMatrix(c_tmp,i);
            }
			if (a.size()>b.size()){
				for(;i<a.size();i++)
                    c.setMatrix(a[i],i);                    
			}
			else{
				for(;i<c.size();i++)
                    c.setMatrix(b[i],i);                    
			}
            // std::clog<<"C="<<c<<std::endl;
            // std::clog<<"-----------"<<std::endl;
            
		}
  
		// addin function (a must be allocated with the right size)
		template<typename PMatrix1,typename PMatrix2>
		void addin(PMatrix1 &a, const PMatrix2 &b) const {
            // std::clog<<"Addin matpoly:"<<std::endl;
            // std::clog<<"A="<<a<<std::endl;
            // std::clog<<"B="<<b<<std::endl;

			for(size_t i=0;i<b.size();i++){
                auto a_tmp=a[i];
				_BMD.addin(a_tmp,b[i]);
                a.setMatrix(a_tmp,i);                    
            }
            // std::clog<<"A="<<a<<std::endl;
            // std::clog<<"-----------"<<std::endl;       
		}
  
		// sub function (c must be allocated with the right size)
		template<typename PMatrix1,typename PMatrix2,typename PMatrix3>
		void sub(PMatrix1 &c, const PMatrix2 &a, const PMatrix3 &b) const {
            size_t i=0;
			for(;i<std::min(a.size(),b.size());i++){
                auto c_tmp=c[i];
                _BMD.sub(c_tmp,a[i],b[i]);
                c.setMatrix(c_tmp,i);
            }
			if (a.size()>b.size()){
				for(;i<a.size();i++)
                    c.setMatrix(a[i],i);                    
			}
			else{
				for(;i<b.size();i++){
                    auto c_tmp=c[i];
                    _BMD.neg(c_tmp,b[i]);                                
                    c.setMatrix(c_tmp,i);
                }
			}
            std::clog<<"C="<<c<<std::endl;
            std::clog<<"-----------"<<std::endl;            
		}

		// subin function (a must be allocated with the right size)
		template<typename PMatrix1,typename PMatrix2>
		void subin(PMatrix1 &a, const PMatrix2 &b) const {
            // std::clog<<"Subin matpoly:"<<std::endl;
            // std::clog<<"A="<<a<<std::endl;
            // std::clog<<"B="<<b<<std::endl;


			for(size_t i=0;i<b.size();i++){
                auto a_tmp=a[i];
                _BMD.subin(a_tmp,b[i]);
                a.setMatrix(a_tmp,i);           
            }
            // std::clog<<"A="<<a<<std::endl;
            // std::clog<<"-----------"<<std::endl;           
		}
	};
}
#endif

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
