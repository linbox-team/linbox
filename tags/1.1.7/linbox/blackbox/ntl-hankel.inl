/*-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+ 
 *    ntl-hankel.inl     NTL_Hankel.cpp file 
 *    Copyright (C) 2003 Austin Lobo, B. David Saunders
 *
 *    Author: Austin Lobo 
 *    Linbox version 2001 and 2002 
 *
 *    This file is included in the template description of ntl-Hankel.h
 *    it contains the implementations of templatized member functions in the 
 *    partial template  specialization for hankel matrices that
 *    are manipulated in fields and rings according to the arithmetic
 *    in the ntl package from V. Shoup
 *
 *    see COPYING for license information
 *
 *    Everything is in the Linbox namespace by virtue of the #include
 *    in ntl-Hankel.h
 *-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+*/

#ifndef __LINBOX_bb_ntl_hankel_INL
#define __LINBOX_bb_ntl_hankel_INL

#include <iostream>
#include <fstream>
#include <NTL/ZZ_pX.h>

namespace LinBox 
{
	/*-----------------------------------------------------------------
	 *----    Destructor
	 *----------------------------------------------------------------*/
	template <class Field>
	inline Hankel<Field>::~Hankel()
	{
#ifdef DBGMSGS
		std::cout << "Hankel::~Hankel():\tDestroyed a " << this->rowDim << "x"<< this->colDim<<
			" Hankel matrix "<< std::endl;
#endif
	}//---- Destructor ---- 
	
	
	
	/*-----------------------------------------------------------------
	 *----    Default Constructor    
	 *----------------------------------------------------------------*/
	template <class Field>
	Hankel<Field>::Hankel() 
	{
		this->shape.shape(BlackboxSpecifier::HANKEL);
#ifdef DBGMSGS
		std::cout << "Hankel::Hankel():\tCreated a " << this->rowDim << "x"<< this->colDim<<
			" Hankel matrix "<< std::endl;
#endif
		
	}//----- Zero Param Constructor ---- [Tested 6/14/02 -- Works]
	
	
	
	
	/*-----------------------------------------------------------------
	 *----- Constructor With User-Supplied First Row And Column
	 *----------------------------------------------------------------*/
	template <class Field>
	Hankel<Field>::Hankel( const Field F,
					   const std::vector<typename Field::Element>&v) 
	{
		// Assumes that the input is a vector of ZZ_p else things will FAIL
		if ( (1 & v.size()) == 0) 
		{
			std::cout << "There must be an ODD number of entries in the input vector " <<
				"The length given is " << v.size();
		}
		assert( (1 & v.size()) == 1);
		
		this->rowDim = (1+v.size())/2; // The vector is 0..2n-2;
		this->colDim = (1+v.size())/2;
		this->sysDim = (1+v.size())/2;
		
		this->pdata.SetMaxLength( v.size());
		//		rpdata.SetMaxLength( v.size());
		for (unsigned int i=0; i< v.size(); i++) 
		{
			this->P.setCoeff( this->pdata, i, v[i]);
			//SetCoeff( rpdata, i, v[v.size()-1-i]);
		}
		
#ifdef DBGMSGS
		std::cout << "Hankel::Hankel(F,V):\tCreated a " << this->rowDim << "x"<< this->colDim<<
			" Hankel matrix "<< std::endl;
#endif
		
	}//----- Constructor given a vector---- 
	
	
	
	/*-----------------------------------------------------------------
	 *-----    Print The Matrix To Screen
	 *----------------------------------------------------------------*/
	template <class Field>
	void Hankel<Field>::print(std::ostream& os) const 
	{
		register size_t i, N, j;
		
		os<< this->rowDim << " " << this->colDim << " " << this->shape.shape() << std::endl;
		N = (this->rowDim-1)<<1;
		
		if ( N < 20 ) {            // Print small matrices in dense format
                    Element tmp;
			for (i = N ; i >= this->colDim-1; i--) {
				for ( j = 0; j < this->colDim ; j++)
					os << " " << this->P.getCoeff(tmp, this->pdata,i-j) ;
				os << std::endl;
			}
		} 
		else {
			// Print large matrices' first row and col
			os << "<Hankel<";
			this->P.write(os, this->pdata) << ">>\n";
		} //[v(2n-2),....,v(0)]; where v(0) is the top right entry of the matrix
		
		return;
	} //---- print()----- 
	
	
	
	
	/*-----------------------------------------------------------------
	 *----    The infamous clone has been created here 
	 *----------------------------------------------------------------*/
	//template <class Field, class Vector>
	//BlackboxArchetype<Vector>* Hankel<Field, Vector>::clone() const 
	//{ 
		//return new Hankel(*this); 
	//}// ------ This is not tested. 
	
	
	/*-----------------------------------------------------------------
	 *----    Save To File, Given Destination Filename
	 *----------------------------------------------------------------*/
	template <class Field>
	void Hankel<Field>::print( char *outFileName) const
	{
		int i, j, N;
		
		std::cout << "Printing hankel matrix to " << outFileName << std::endl;
		
		if ( outFileName == NULL ) 
			print();    // Print to stdout if no file is specified
		else { 
			std::ofstream o_fp(outFileName, std::ios::out);
			o_fp << this->rowDim << " " << this->colDim << " " << this->shape.shape() << std::endl ;
			o_fp << "<Hankel<";
			this->P.write(o_fp, this->pdata) << ">>\n";
			
			o_fp.close();
		}
		return;
	} // print(char *) 
	
	
	
	/*-----------------------------------------------------------------
	 *    Make the matrix LOWER triangular with determinant 1.
	 *    i.e. clear the last this->coldim-1 elements in the this->data vector
	 *----------------------------------------------------------------*/
	template <class Field>
	void Hankel<Field>::setToUniModLT()
	{
                int L = (this->rowDim-1)<<1;
		this->shape.shape(BlackboxSpecifier::UNIMOD_LT);

                Element one,zero;
                this->K.init(one,1);
                this->K.init(zero,0);
		for (int i=this->rowDim-1; i <= L; i++ ) {
			// zero out the below-diagonal entries 
                    this->P.setCoeff(this->pdata,i,zero);
		}
                    // set the antidiagonal to 1
		this->P.setCoeff( this->pdata, this->rowDim-1, one);       // update the corresponding coeff of this->pdata
		//reverse(rpdata,this->pdata);        // no need to construct the transpose
		return;
	}// 
	
	
	
	/*-----------------------------------------------------------------
	 *    Make matrix a unimodular UPPER Triangular with det 1
	 *    i.e. clear the first N-1 elements in the this->data vector
	 *    and make the elements below the anti-diagonal all zero
	 *----------------------------------------------------------------*/
	template <class Field>
	void Hankel<Field>::setToUniModUT()
	{
		this->shape.shape(BlackboxSpecifier::UNIMOD_UT);
		
                Element one,zero;
                this->K.init(one,1);
                this->K.init(zero,0);

		for (size_t i=0; i < this->rowDim-1; i++ ) {
			// zero out the below-antidiagonal entries 
                    this->P.setCoeff(this->pdata, i , zero);
		}

                    // set antidiagonal to 1
		this->P.setCoeff(this->pdata,this->rowDim-1, one);      // update the corresponding coeff of this->pdata
		//reverse(rpdata,this->pdata);    // no need to construct the transpose
		
		return;
	}// 
	
	
	
	/*-----------------------------------------------------------------
	 *    Apply the matrix to a vector Here the input and output 
	 *    vectors are both over the SAME prime ZZ_p field as the 
	 *    Hankel matrix itself.
	 *----------------------------------------------------------------*/
	template <class Field>
	template<class OutVector, class InVector>
	OutVector& Hankel<Field>::apply( OutVector &v_out, 
										  const InVector& v_in) const
	{  
		if (v_out.size() != this->rowdim())
			std::cout << "\tHankel::apply()\t output vector not correct size, at "
					  << v_out.size() << ". System rowdim is" <<  this->rowdim() << std::endl;
		if ( v_out.size() != v_in.size())
			std::cout << "\tHankel::apply()\t input vector not correct size at " 
					  << v_in.size() << ". System coldim is" <<  this->coldim() << std::endl;
		assert((v_out.size() == this->rowdim()) && 
			   (v_in.size() == this->coldim()))  ;
		
		NTL::ZZ_pX pxOut, pxIn;
		pxIn.SetMaxLength( v_in.size()-1);
		for (unsigned int i=0; i< v_in.size(); i++)
			this->P.setCoeff( pxIn, i, v_in[i]);
		
#ifdef DBGMSGS
		std::cout << "\npX in is " << pxIn << std::endl;
		std::cout << "multiplied by " << this->pdata << std::endl;
#endif
		mul(pxOut,pxIn,this->pdata);
		
#ifdef DBGMSGS
		std::cout <<"pxOut is " << pxOut << std::endl;
#endif
		int N = this->rowdim();
		for ( int i= 0; i < N; i++) 
			this->P.getCoeff(v_out[N-1-i], pxOut, N-1+i);
		
		return v_out;
		
	}
	
	
	/*-----------------------------------------------------------------
	 *    Apply the transposed matrix to a vector Here the input and output 
	 *    vectors are both over the SAME prime ZZ_p field as the 
	 *    Hankel matrix itself. Calls the multiply from the Toeplitz matrix
	 *    Since Hankel is symmetric, this is the same as apply
	 *----------------------------------------------------------------*/
	template <class Field>
	template <class OutVector, class InVector>
	OutVector& Hankel<Field>::applyTranspose( OutVector &v_out, 
													 const InVector& v_in) const
	{  
		return(v_out = apply(v_out,v_in));

	}
	
	
	
} // namespace LinBox

#endif //__LINBOX_bb_ntl_hankel_INL
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
