/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
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

		// add function (a must be allocated with the right size)
		template<typename PMatrix1,typename PMatrix2,typename PMatrix3>
		void add(PMatrix1 &a, const PMatrix2 &b, const PMatrix3 &c){
			size_t i=0;
			for(;i<std::min(b.size(),c.size());i++)
				_BMD.add(a[i],b[i],c[i]);
			if (b.size()>c.size()){
				for(;i<b.size();i++)
					a[i]=b[i];                                
			}
			else{
				for(;i<c.size();i++)
					a[i]=c[i];                                
			}
		}
  
		// addin function (a must be allocated with the right size)
		template<typename PMatrix1,typename PMatrix2>
		void addin(PMatrix1 &a, const PMatrix2 &b){	
			for(size_t i=0;i<b.size();i++)
				_BMD.addin(a[i],b[i]);
		}
  
		// sub function (a must be allocated with the right size)
		template<typename PMatrix1,typename PMatrix2,typename PMatrix3>
		void sub(PMatrix1 &a, const PMatrix2 &b, const PMatrix3 &c){
			size_t i=0;
			for(;i<std::min(b.size(),c.size());i++)
				_BMD.sub(a[i],b[i],c[i]);
			if (b.size()>c.size()){
				for(;i<b.size();i++)
					a[i]=b[i];                                
			}
			else{ 
				for(;i<c.size();i++)
					_BMD.neg(a[i],c[i]);                                
			}
		}

		// subin function (a must be allocated with the right size)
		template<typename PMatrix1,typename PMatrix2>
		void subin(PMatrix1 &a, const PMatrix2 &b){
			for(size_t i=0;i<b.size();i++)
				_BMD.subin(a[i],b[i]);
		}
	};
}
#endif

