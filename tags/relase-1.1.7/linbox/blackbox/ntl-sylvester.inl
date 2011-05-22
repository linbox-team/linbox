/*-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
 *    ntl-sylvester.h
 *    Copyright (C) 2003 Austin Lobo, B. David Saunders
 *
 *    Member functions for the sylvester matrix in one variable
 *    for polynomials in one variable.
 *    Linbox version 2003
 *
 *    See COPYING for licence information
 *-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+*/

#ifndef __LINBOX_bb_ntl_sylvester_INL
#define __LINBOX_bb_ntl_sylvester_INL

#include <iostream>

namespace LinBox
{

    /*----------------------------------------------------------------------
     *    Destructor
     *---------------------------------------------------------------------*/
    template <class Field>
    inline Sylvester<Field>::~Sylvester() 
    {
	  
#ifdef DBGMSGS
	  std::cout << "Sylvester::~Sylvester():\tDestroyed a " 
			<< rowDim << "x"<< colDim << " Sylvester matrix "<< std::endl;
#endif
	  
    }// End, ~Sylvester()
    
    
    
    /*----------------------------------------------------------------------
     *    Default constructor.
     *---------------------------------------------------------------------*/
    template <class Field>
    Sylvester<Field>::Sylvester()  
    {
	  sysDim =               // Default dimension is 0
	  rowDim =               // Default row dim is 0
	  colDim = 0;            // Default col dim is 0
	  
#ifdef DBGMSGS
	  std::cout << "Sylvester::Sylvester():\tCreated a " << rowDim << "x"<< colDim<<
		" Sylvester matrix "<< std::endl;
#endif

    }



    /*----------------------------------------------------------------------
     *    Constructor. Builds two rectangular Toeplitz blocks, the top from 
     *                 px with deg(qx) rows, and the bottom out qx, with deg(px)
     *                 rows.
     *---------------------------------------------------------------------*/
    template <class Field>
    Sylvester<Field>::Sylvester( 
		    const Field F,
		    const std::vector<typename Field::Element> &vp,
		    const std::vector<typename Field::Element> &vq 
		    ) :	  pdata(vp), qdata(vq)
    {
	  
	  // Set the row, col, and system dimensions
	  rowDim = 
	  colDim =
	  sysDim = vp.size() + vq.size()-2;  // square system with dimension deg(px)+deg(qx) 
	  
	  // Copy the vector to a polynomial representation, for px
	  pxdata.SetMaxLength( vp.size() );
	  for (size_t ip=0; ip< vp.size(); ip++ )
	  SetCoeff( pxdata, ip, vp[ip]);
	  
	  // Copy the vector to a polynomial representation, for qx
	  qxdata.SetMaxLength( vq.size() );
	  for ( size_t iq=0; iq< vq.size(); iq++ ) 
	  SetCoeff( qxdata, iq, vq[iq] );
	  
	  
#ifdef DBGMSGS
	  std::cout << "Sylvester::Sylvester(F,v,w):\tCreated a " << rowDim << 
		"x"<< colDim<< " Sylvester matrix "<< std::endl;
#endif
	  
	  /* Inserted for a timing calculation with the gcd
	   *  double start_time = GetTime();
	   *	  ZZ_pX gx;
	   *	  GCD( gx, pxdata, qxdata);
	   *
	   *	  double end_time = GetTime();
	   *	  cout << gx << endl;
	   *	  std::cout << "Time for gcd = " << end_time-start_time << " sec\n";
	   */
    }// Sylvester()



    /*----------------------------------------------------------------------
     *    Prints to an ouput stream or to stdout by default -- useful for
     *    debugging
     *---------------------------------------------------------------------*/
    template <class Field>
    void Sylvester<Field>::print(std::ostream& os) const 
    {
	  if ( sysDim < 20 ) 
		{
		    os << rowDim << " " << colDim << std::endl;
		    for (size_t irow=0; irow < qxdeg(); irow++ ) 
			  {
				os << "[  ";	
				for (size_t icoeff=0; icoeff < irow; icoeff++ ) 
				    os << "0  ";
				
				for (size_t icoeff = pxdeg(); icoeff > 0; icoeff-- ) 
				    os << pdata[ icoeff ] << " ";
				
				os << pdata[ 0 ] << " ";
				
				for (size_t icoeff = sysDim-1-pxdeg()-irow; icoeff > 0; icoeff-- ) 
				    os << "0  ";
				os << " ]\n";
			  }
		    for (size_t irow=0; irow < pxdeg(); irow++ ) 
			  {
				os << "[ "; 
				
				for (size_t icoeff=0; icoeff < irow; icoeff++ )
				    os << "0  "; 		
				
				for (size_t icoeff = qxdeg(); icoeff > 0; icoeff-- )
				    os << qdata[ icoeff ] << " "; 
				os << qdata[ 0 ] << " ";
				
				for (size_t icoeff = sysDim-1-qxdeg()-irow; icoeff > 0; icoeff-- ) 
				    os << "0  ";
				
				os << " ]\n";
			  }
		} 
	  else 
		{
		    os << rowDim << " " << colDim << std::endl;
		    os << pxdata;
		    os << "\n";
		    os << qxdata;
		    os << "\n";
		}
	  
    }// Print()
    



    /*----------------------------------------------------------------------
     *    Print to a file. By default print to stdout
     *---------------------------------------------------------------------*/

    template <class Field>
    void Sylvester<Field>::print(char *outFileName) const 
    {
	  if ( outFileName == NULL )
		print();
	  else {
		std::ofstream o_fp(outFileName, std::ios::out);
		o_fp << rowDim << " " << colDim << std::endl;
		o_fp << "[";
		for (size_t icoeff = 0; icoeff <= pxdeg() ; icoeff++ )
		    o_fp << pdata[ icoeff ] << " ";
		o_fp << "]\n";
		o_fp << "[";
		for (size_t icoeff = 0; icoeff <= qxdeg(); icoeff++ )
		    o_fp << qdata[ icoeff ] << " ";
		o_fp << "]\n";
		o_fp.close();
	  }


    }// Print



    /*----------------------------------------------------------------------
     *    Print to a file in the format specified by the French Linbox group
     *    for LUDivine. Format is <nrow> <ncol> M
     *                            <row> <col> <entry>
     *                              .....
     *                             0 0 0
     *--------------------------------------------------------------------- */

    template <class Field>
    void Sylvester<Field>::printcp(char *outFileName) const 
    {
	  if ( outFileName == NULL )
		print();
	  else {
		std::ofstream os(outFileName, std::ios::out);
		os << rowDim << " " << colDim << " M" << std::endl;
		size_t irow;
		size_t icoeff;
		for (irow=0; irow < qxdeg(); irow++ ) 
		    {
			  
			  for (icoeff = pxdeg(); icoeff > 0; icoeff-- ) 
				os << (irow+1) << " " << (irow + pxdeg()- icoeff +1) << " " << 
				    pdata[ icoeff ] << "\n";

			  os << (irow+1) << " " << (irow + pxdeg()+1) << " " << 
				    pdata[ icoeff ] << "\n";

		    }
		
		for (irow=0; irow < pxdeg(); irow++ ) 
		    {

			  for (icoeff = qxdeg(); icoeff > 0; icoeff-- )
				os << (irow+qxdeg()+1) << " " << (irow + pxdeg()- icoeff +1) << " " << 
				    qdata[ icoeff ] << "\n";
			  
			  os << (irow+qxdeg()+1) << " " << (irow + pxdeg()+1) << " " << 
				qdata[ icoeff ] << "\n";
		    }

		os << "0 0 0\n";

		os.close();
	  }


    }



    /* *----------------------------------------------------------------------
     *    Applytranspose. Does 2 FFT's each of degree (deg(qx)+deg(px))
     *                    Speed can be improved to 1 FFT in reverse direction
     *                    by saving pre-computed FFT values
     *---------------------------------------------------------------------*/

    template <class Field>
    template <class OutVector, class InVector>
    OutVector& Sylvester<Field>::apply( OutVector &v_out, 
							   const InVector& v_in) const
    {  
	  /* uncomment the following lines here and in apply() and swap method names
	   * to swap the apply() for the faster applyTranspose()
	   *
	   *	  static int flipper = 1;
	   *	  if (flipper == 1) 
	   *		{ cout << "\t sylvester applyTranspose is flipped apply\n"; 
	   *		  flipper++;
	   *		}
	   */
	  
	  /*--------------- Check the size of the output vector ----------*/
	  if ( v_out.size() != sysdim() )
		std::cout << "\tSylvester::apply()\t output vector not correct size, at "
			    << v_out.size() << ". System rowdim is" <<  sysdim() << std::endl;
	  NTL::ZZ_pX txOut, txIn;
	  
	  /*--------------- Convert input vector to a polynomial ---------*/
	  txIn.SetMaxLength( v_in.size() -1 );
	  for ( size_t i=0; i < v_in.size(); i++ )
		SetCoeff( txIn, i, v_in[i] );

	  /*--------------  Poly multiply the upper Sylvester poly by input ------*/
	  mul( txOut, txIn, pxdata);
	  int Nq = qxdeg();
	  int m  = pxdeg();
	  /*--------------  vout[0..deg(q)-1] <--- txout[deg(qx)...2deg(qx)-1] --- */
	  for ( int i=0; i < Nq; ++i ) 
		GetCoeff(v_out[i], txOut, m+i);

	  /*--------------  Poly multiply the lower Sylvester poly by input -----*/
	  mul( txOut, txIn, qxdata);

	  int Np = pxdeg();
	  int n  = qxdeg();

	  /*--------------  vout[deg(qx)..deg(qx)+deg(px)-1] <--- 
	   *                              txout[deg(qx)...deg(qx)+deg(px)-1] --- */
	  for ( int i=0; i < Np; ++i ) 
		GetCoeff(v_out[Nq+i], txOut, n+i );


	  return v_out;
    }// Apply -- Tested and works 3/17/03


    /*----------------------------------------------------------------------
     *    Apply Transpose: Does 2 polymults of degree deg(qx), and deg(px) 
     *                     respectively. It takes approx 55% of the time of
     *                     apply, in apparent violation of Tellegen's theorem
     *                     Actually, both apply and applyT can be done with
     *                     just one FFT per call, by saving previous FFTs
     *                     we use deg(qx) = m; and deg(px) = n; N=m+n below
     *---------------------------------------------------------------------*/
    template <class Field>
    template <class OutVector, class InVector>
    OutVector& Sylvester<Field>::applyTranspose( OutVector &v_out, 
									const InVector& v_in) const
    {  
	  /* uncomment the following lines here and in apply() and swap method names
	   * to swap the apply() for the faster applyTranspose()
	   *
	   *	  static int flipper=1;
	   *	  if (flipper == 1) 
	   *		{ cout << "\tSylvester apply is flipped applyTranspose\n"; flipper++;}
	   */	  

	  /*--------- Check for size-compatibility of output vectors ------------*/
	  if ( v_out.size() != sysdim() )
		std::cout << "\tSylvester::apply()\t output vector not correct size, at "
			    << v_out.size() << ". System rowdim is" <<  sysdim() << std::endl;
	  
	  NTL::ZZ_pX txOut, txIn;
	  NTL::ZZ_p tval;

	  /*--------- We need to reverse the input vector           -----------*/
	  /*--------- txIn[0...m-1] <--- vin[m-1...0]               ----------*/
	  txIn.SetMaxLength( qxdeg()-1 );
	  for ( size_t i=0; i < qxdeg(); i++ )
		SetCoeff( txIn, qxdeg()-1-i, v_in[i] );
	  
	  mul( txOut, txIn, pxdata);     // Do the poly multiply
	  
	  /*--------- We need to reverse the output vector           -----------*/
	  /*-------  v_out[0..N-1] <--- txOut[N-1...0]                  --------*/
	  for (size_t i=0; i < v_out.size(); i++) 
		GetCoeff(v_out[i], txOut, sysdim()-1-i); // Extract the coeffs

	  /*--------- We need to reverse the input vector           -----------*/
	  /*--------- txIn[m ... m+n-1] <--- vin[m+n-1 ... m]       -----------*/
	  txIn.SetMaxLength( pxdeg()-1 );
	  for ( size_t i=0; i < pxdeg(); i++ )
		SetCoeff( txIn,pxdeg()-1-i, v_in[qxdeg()+i] );

	  mul( txOut, txIn, qxdata);    // do the poly multiply

	  /*--------- We need to reverse the input vector           -----------*/
	  /*-------  v_out[0..N-1] <--- txOut[N-1...0]                  --------*/

	  for (size_t i=0; i < v_out.size(); i++) {
		GetCoeff(tval, txOut, sysdim()-1-i);
		add( v_out[i], v_out[i], tval ); // add to the accumulated output
	  }



	  return v_out;
    }



}//End, namespace LinBox


#endif //__LINBOX_bb_ntl_sylvester_INL
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
