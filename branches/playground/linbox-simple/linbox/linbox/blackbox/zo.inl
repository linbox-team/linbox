/* -*- mode: C++; style: linux -*- */

/* linbox/blackbox/nag-sparse.h
 * Copyright (C) 2002 Rich Seagraves
 *
 * Written by Rich Seagraves <seagrave@cis.udel.edu>
 * Modified by Zhendong Wan, -bds
 * ------------------------------------
 *
 * See COPYING for license information.
 */

namespace LinBox
{
  
  template<class Field>
  template<class OutVector, class InVector>
  OutVector & ZeroOne<Field>::apply(OutVector & y, const InVector & x) const
  //OutVector & ZeroOne<Field>::apply(OutVector & y, const InVector & x)
  {
    //std::cout << "\n\napply\n\n";

    linbox_check((y.size()==rowdim())&&(x.size()==coldim()));
    
    FieldAXPY<Field> accum (_F);

    typename OutVector::iterator yp;
    typename InVector::const_iterator xp;
    PointerVector::const_iterator ip;
    IndexVector::const_iterator jp;
    
    if( sorted == sortedByCol ) switch_sort();

    xp=x.begin();
    yp=y.begin();
    accum.reset();

    for(ip = _indexP.begin(); ip < _indexP.end()-1; ++ip, ++yp)
    //for(ip = _indexP.begin(); ip < _indexP.end()-2; ++ip, ++yp)  //zigzag way
      {	
	for(jp = *ip; jp <*(ip + 1); ++jp)
	  accum.accumulate_special( *(xp + *jp) );
        accum.get(*yp);
        accum.reset();

	/*
	++ip;++yp;
	
	for(jp = *(ip + 1); jp >*ip; --jp)
	  accum.accumulate_special( *(xp + *(jp-1)) );
	accum.get(*yp);
	accum.reset();
	*/

      }

    return y;
  }

  /* if you want to keep two copies for the matrix, one of which is sorted by row, 
   * the other by column, then you want to use this applyTranspose function. In 
   * this case, un-comment this one and comment out the applyTranspose further down
   */
  /*
  template<class Field>
  template<class OutVector, class InVector>
  OutVector & ZeroOne<Field>::applyTranspose(OutVector & y, const InVector & x) const
  {
    linbox_check((y.size()==coldim())&&(x.size()==rowdim()));
    
    FieldAXPY<Field> accum (_F);
    
    typename OutVector::iterator yp;
    typename InVector::const_iterator xp;
    PointerVector::const_iterator ip; 
    IndexVector::const_iterator jp;      

        
    xp=x.begin();
    yp=y.begin();
    accum.reset();

    for(ip = _colP.begin(); ip < _colP.end()-1; ++ip, ++yp)
      {       
	for(jp = *ip; jp <*(ip + 1); ++jp)
	  accum.accumulate_special( *(xp + *jp) );
        accum.get(*yp);
        accum.reset();
      }

    return y;
    
  }
  */

  template<class Field>
  template<class OutVector, class InVector>
  OutVector & ZeroOne<Field>::applyTranspose(OutVector & y, const InVector & x) const
  {
    //std::cout << "\n\napplyTranspose\n\n";
    
    linbox_check((y.size()==coldim())&&(x.size()==rowdim()));
    
    FieldAXPY<Field> accum (_F);
    
    typename OutVector::iterator yp;
    typename InVector::const_iterator xp;
    PointerVector::const_iterator ip; 
    IndexVector::const_iterator jp;      
    
        
    xp=x.begin();
    yp=y.begin();
    accum.reset();

    if( sorted == sortedByRow ) switch_sort();

    for(ip = _indexP.begin(); ip < _indexP.end()-1; ++ip, ++yp)
      {       
	for(jp = *ip; jp <*(ip + 1); ++jp)
	  accum.accumulate_special( *(xp + *jp) );
        accum.get(*yp);
        accum.reset();
      }

    return y;
    
  }
       
}//End of LinBox
