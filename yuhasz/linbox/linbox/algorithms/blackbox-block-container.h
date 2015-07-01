/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/algorithms/blackbox-block-container.h
 * Copyright (C) 2002 Pascal Giorgi
 *
 * Written by Pascal Giorgi pascal.giorgi@ens-lyon.fr
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

#ifndef __BLACKBOX_BLOCK_CONTAINER_H
#define __BLACKBOX_BLOCK_CONTAINER_H

#include <linbox/randiter/archetype.h>

#include <linbox/algorithms/blackbox-block-container-base.h>

namespace LinBox 
{

  template<class Field, class Vector,class Matrix>
    class BlackboxBlockContainer : public BlackboxBlockContainerBase<Field,Vector,Matrix > {
    public:
    typedef typename Field::RandIter RandIter;
	
    BlackboxBlockContainer () {} 

    BlackboxBlockContainer(const Blackbox *D, const Field &F, const Block  &u0) 
      : BlackboxBlockContainerBase<Field, Vector,Matrix > (D, F) , w(u0.coldim(),Vector (u0.rowdim(),Zero))
      { init (u0, u0); }
    
    BlackboxBlockContainer(const Blackbox *D, const Field &F, const Block &u0, const Block& v0) 
      : BlackboxBlockContainerBase<Field, Vector,Matrix > (D, F) , w(v0.coldim(),Vector (u0.coldim(),Zero))
      { init (u0, v0);}
    
    BlackboxBlockContainer(const Blackbox *D, const Field &F,int n,int m) 
      : BlackboxBlockContainerBase<Field, Vector,Matrix > (D, F) , w(m,Vector (D->rowdim(),Zero))
      { init (n,m); }
    
    protected: 
    Block w;

    void _launch () {
      if (casenumber) {	
	Mul(w,*_BB,v);
	Mul(_value,u,w);
	casenumber = 0;
      } else { 
	Mul(v,*_BB,w);
	Mul(_value,u,v);
	casenumber = 1;
      }  
    }

    void _wait () {}
  };
 
}

#endif // __BLACKBOX_BLOCK_CONTAINER_H

