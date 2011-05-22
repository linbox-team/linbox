/* Copyright (C) 2010 LinBox
 * Written by 
 * zhendong wan 
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


#ifndef __LINBOX_bb_nullmatrix_H
#define __LINBOX_bb_nullmatrix_H

#include <linbox/util/debug.h>
#include <linbox/blackbox/blackbox-interface.h>
 
namespace LinBox
{
  
  // couldn't a null instance of one of the other classes serve as well?

  /** \brief  This is a representation of the 0 by 0 empty matrix which does not occupy memory.  
   * It has it's uses!
   * \ingroup blackbox
   */
  
	class NullMatrix : public  BlackboxInterface{
	public:
		NullMatrix() {}//cout << "NullMatrix default cstor" << endl;}
		NullMatrix(const NullMatrix& n) {}
	        virtual ~NullMatrix() {}
		
	public:
		
		template<class OutVector, class InVector>
		inline OutVector& apply(OutVector& y, const InVector& x) const {
			linbox_check(y.size()==0);
			linbox_check(x.size()==0);
			return y;
		}

		/* applyIn is depreciated.  If you have a desire to use it, please tell me about that.  -bds
		   virtual inline Vector& applyIn(Vector& x) const {
		   linbox_check(x.size()==0);
		   return x;
		   }
    */
		

		template<class OutVector, class InVector>
		inline OutVector& applyTranspose(OutVector& y, const InVector& x) const {
			linbox_check(y.size()==0);
			linbox_check(x.size()==0);
			return y;
		}

		/* applyIn is depreciated.  If you have a desire to use it, please tell me about that.  -bds
		   virtual inline Vector& applyTransposeIn(Vector& x) const {
		   linbox_check(x.size()==0);
		   return x;
		   }
		*/
		
		virtual inline size_t rowdim() const {return 0;}
		
		virtual inline size_t coldim() const {return 0;}
	
		template<typename _Tp1>
		struct rebind
		{ typedef NullMatrix other; };


	};
	
} 

#endif //__LINBOX_bb_nullmatrix_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
