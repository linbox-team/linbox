/* -*- mode: c; style: linux -*- */

/* linbox/field/vector-domain.C
 * Copyright (C) 2001-2002 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
 * Modified by Zhendong Wan <wan@mail.eecis.udel.edu>
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

#ifndef DENSE_VECTOR_DOMAIN_C
#define DENSE_VECTOR_DOMAIN_C
#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include "linbox/field/dense-vector-domain.h"
#include "linbox/util/debug.h"

namespace LinBox
{
                template<class Field>
  		template<class Vect>
		ostream& DenseVectorDomain<Field>::write (ostream &os, const Vect &x) const
		{
		  typename Vect::const_iterator p;
		  for(p=x.begin();p!=x.end();++p)
		    {
		      _F.write(os,*p);
		      os<<" ";
		    }
		  return os;
		}
                
                template<class Field>
		template<class Vect>
		istream& DenseVectorDomain<Field>::read (istream &is, Vect &x) const
		{
		  typename Vect::iterator p;
		  for(p=x.begin();p!=x.end();++p)
		    _F.read(is,*p);
		  return is;
		}

                template<class Field>
		template<class Vector1, class Vector2>
		DenseVectorDomain<Field>::Element& DenseVectorDomain<Field>::dotproduct (DenseVectorDomain<Field>::Element& res, const Vector1& v1, const Vector2& v2) const
		{
		  linbox_check(v1.size() == v2.size());
		  _F.init(res,0);
		  typename Vector1::const_iterator p1;
		  typename Vector2::const_iterator p2;
		  for(p1=v1.begin(),p2=v2.begin();p1!=v1.end();++p1,++p2)
		    _F.axpyin(res,*p1,*p2);
		  return res;
		}

                template<class Field>
		template<class Vector1, class Vector2>
		Vector1& DenseVectorDomain<Field>::mul (Vector1 &res, const Vector2 &x, const DenseVectorDomain<Field>::Element &a) const
		{
		  linbox_check(res.size() == x.size());
		  typename Vector1::iterator p1;
		  typename Vector2::const_iterator p2;
		  for(p1=res.begin(),p2=x.begin();p1!=res.end();++p1,++p2)
		    _F.mul(*p1,*p2,a);
		  return res;
		}

                template<class Field>
		template<class Vector1>
		Vector1& DenseVectorDomain<Field>::mulin (Vector1 &x, const Element &a) const
		{
		  typename Vector1::iterator p;
		  for(p=x.begin();p!=x.end();++p)
		    _F.mulin(*p,a);
		  return x;
		}
  
                template<class Field>
		template<class Vector1, class Vector2, class Vector3>
		Vector1& DenseVectorDomain<Field>::axpy (Vector1 &res, const DenseVectorDomain<Field>::Element &a, const Vector2 &x, const Vector3 &y) const
		{
		  linbox_check((res.size()==y.size())&&(res.size()==x.size()));
		  typename Vector1::iterator p1;
		  typename Vector2::const_iterator p2;
		  typename Vector3::const_iterator p3;
		  for(p1=res.begin(),p2=x.begin(),p3=y.begin();p1!=res.end();++p1,++p2,++p3)
		    _F.axpy(*p1,a,*p2,*p3);
		  return res;
		}
		  

                template<class Field>
		template<class Vector1, class Vector2>
		Vector1& DenseVectorDomain<Field>::axpyin (Vector1 &y, const Element &a, const Vector2 &x) const
		{
		  linbox_check(y.size()==x.size());
		  typename Vector1::iterator p1;
		  typename Vector2::const_iterator p2;
		  for(p1=y.begin(),p2=x.begin();p1!=y.end();++p1,++p2)
		    _F.axpyin(*p1,a,*p2);
		  return y;
		}
    
}
#endif
