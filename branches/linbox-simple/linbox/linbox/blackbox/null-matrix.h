/*-* mode:C++;tab-with=8;c-basic-offset:8 -*- */
/* zhendong wan */
#ifndef __NULLMATRIX_H
#define __NULLMATRIX_H

#include <linbox/util/debug.h>
#include <linbox/blackbox/blackbox-interface.h>
 
namespace LinBox{
  
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

#endif
