/* -*- mode: C++; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*- */

/*-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+ 
 *    ntl-toeplitz.inl     NTL_Toeplitz.cpp file 
 *
 *    Copyright (C) 2002 Austin Lobo, B. David Saunders
 *    Author: Austin Lobo 
 *    Linbox version 2001 and 2002 
 *
 *    This file is included in the template description of ntl-Toeplitz.h
 *    it contains the implementations of templatized member functions in the 
 *    partial template  specialization for toeplitz matrices that
 *    are manipulated in fields and rings according to the arithmetic
 *    in the ntl package from V. Shoup
 *
 *    Everything is in the Linbox namespace by virtue of the #include
 *    in ntl-Toeplitz.h
 *-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+*/

#include <iostream>
#include <fstream>
#include <NTL/ZZ_pX.h>
using namespace NTL;
namespace LinBox 
{
	/*-----------------------------------------------------------------
	 *----    Destructor
	 *----------------------------------------------------------------*/
	template <class Field, class Vector>
	inline Toeplitz<Field, Vector>::~Toeplitz()
	{
#ifdef DBGMSGS
		std::cout << "Toeplitz::~Toeplitz():\tDestroyed a " << rowDim << "x"<< colDim<<
			" Toeplitz matrix "<< std::endl;
#endif
	}//---- Destructor ---- [Tested 6/14/02 -- Works]
	
	
	
	/*-----------------------------------------------------------------
	 *----    Zero Parameter Constructor    
	 *----------------------------------------------------------------*/
	template <class Field, class Vector>
	Toeplitz<Field, Vector>::Toeplitz()
	{
		shape  =
		sysDim =               // Default dimension is 0
		rowDim =               // Default row dim is 0
		colDim = 0;            // Default col dim is 0
#ifdef DBGMSGS
		std::cout << "Toeplitz::Toeplitz():\tCreated a " << rowDim << "x"<< colDim<<
			" Toeplitz matrix "<< std::endl;
#endif
		
	}//----- Zero Param Constructor ---- [Tested 6/14/02 -- Works]
	
	
	
	
	/*-----------------------------------------------------------------
	 *----- Constructor With User-Supplied First Row And Column
	 *----------------------------------------------------------------*/
	template <class Field, class Vector>
	Toeplitz<Field, Vector>::Toeplitz( const Field F, 
									   const std::vector<typename Field::Element>&v)
	{
		// Assumes that the input is a vector of ZZ_p else things will FAIL
		if ( (1 & v.size()) == 0) 
			{
				std::cout << "There must be an ODD number of entries in the input vector " <<
					"The length given is " << v.size();
			}
		assert( (1 & v.size()) == 1);
		
		rowDim = (1+v.size())/2; // The vector is 0..2n-2;
		colDim = (1+v.size())/2;
		sysDim = (1+v.size())/2;
		
		data = v;
		// bds //pdata.SetMaxLength((long) v.size());
		pdata.SetMaxLength( v.size());
		// bds //rpdata.SetMaxLength((long) v.size());
		rpdata.SetMaxLength( v.size());
		for ( size_t i=0; i< v.size(); i++) 
			{
				SetCoeff( pdata, i, v[i]);
				SetCoeff( rpdata, i, v[v.size()-1-i]);
			}
		
#ifdef DBGMSGS
		std::cout << "Toeplitz::Toeplitz(F,V):\tCreated a " << rowDim << "x"<< colDim<<
			" Toeplitz matrix "<< std::endl;
#endif
		
	}//----- Constructor given a vector---- [Tested 6/14/02 -- Works]
	
	
	
	/*-----------------------------------------------------------------
	 *-----    Print The Matrix To Screen
	 *----------------------------------------------------------------*/
	template <class Field, class Vector>
	void Toeplitz<Field, Vector>::print(std::ostream& os) const 
	{
		
		register int i, N;
		register unsigned int j;
		
		os<< rowDim << " " << colDim << " " << shape << std::endl;
		N = rowDim + colDim -1;
		
		if ( N < 20 )             // Print small matrices in dense format
			{
				for (i = colDim-1; i < N; i++) 
					{
						for ( j = 0; j < colDim ; j++)
							os << " " << data[i-j] ;
						os << std::endl;
					}
			} 
		else 
			{                    // Print large matrices' first row and col
				os << rowDim << " " << colDim << " " << shape << std::endl ;
				os << "[";
				for (int i=data.size()-1; i>= 0;i--)
					os << data[i] << " ";
				os << "]\n";
				os << pdata << std::endl;
			} //[v(2n-2),....,v(0)]; where v(0) is the top right entry of the matrix
		
		return;
	} //---- print()----- [Tested 6/14/02 -- Works]
	
	
	
	
	/*-----------------------------------------------------------------
	 *----    The infamous clone has been created here 
	 *----------------------------------------------------------------*/
	template <class Field, class Vector>
	BlackboxArchetype<Vector>* Toeplitz<Field, Vector>::clone() const 
	{ 
		return new Toeplitz(*this); 
	}// ------ This is not tested. 
	
	
	/*-----------------------------------------------------------------
	 *----    Save To File, Given Destination Filename
	 *----------------------------------------------------------------*/
	template <class Field, class Vector>
	void Toeplitz<Field, Vector>::print( char *outFileName) const
	{
		int i, j, N;
		
		std::cout << "Printing toeplitz matrix to " << outFileName << std::endl;
		
		if ( outFileName == NULL ) 
			print();    // Print to stdout if no file is specified
		else 
			{
				std::ofstream o_fp(outFileName, std::ios::out);
				o_fp << rowDim << " " << colDim << " " << shape << std::endl ;
				o_fp << "[";
				for (i=data.size()-1; i>= 0;i--) o_fp << data[i] << " ";
				o_fp << "]\n";
				
				o_fp.close();
			}
		return;
	} // print(char *) [Tested 6/14/02 -- Works]
	
	
	
	/*-----------------------------------------------------------------
	 *    Make the matrix upper triangular with determinant 1.
	 *    i.e. clear the last N-1 elements in the data vector
	 *----------------------------------------------------------------*/
	template <class Field, class Vector>
	void Toeplitz<Field, Vector>::setToUniModUT()
	{
		int L = data.size();
		int N = sysDim;
		shape = UnimodUT;
		
		for (int i=N; i<L; i++ )
			K.init(data[i],0);     // zero out the below-diagonal entries 
		K.init(data[N-1],1);
		// AAL : change here to zero out the higher degree terms in poly
		//       and the lower degree terms in rpoly
		return;
	}// [UNCOMMENTED PART Tested 6/14/02 -- Works]
	
	
	
	/*-----------------------------------------------------------------
	 *    Make matrix a unimodular Lower Triangular with det 1
	 *    i.e. clear the first N-1 elements in the data vector
	 *----------------------------------------------------------------*/
	template <class Field, class Vector>
	void Toeplitz<Field, Vector>::setToUniModLT()
	{
		int L = data.size();
		int N = sysDim;
		shape = UnimodUT;
		
		for (int i=0; i<N; i++ )
			K.init(data[i],0);     // zero out the ABOVE-diagonal entries 
		K.init(data[N-1],1);
		// AAL : change here to zero out the lower degree terms in poly
		//       and the lower higher terms in rpoly
		
		return;
	}// [UNCOMMENTED PART Tested 6/14/02 -- Works]
	
	
	
	/*-----------------------------------------------------------------
	 *    Apply the matrix to a vector Here the input and output 
	 *    vectors are both over the SAME prime ZZ_p field as the 
	 *    Toeplitz matrix itself.
	 *----------------------------------------------------------------*/
	template <class Field, class Vector>
	Vector& Toeplitz<Field, Vector>::apply( Vector &v_out, 
											const Vector& v_in) const
	{  
		
		if (v_out.size() != rowdim())
			std::cout << "\tToeplitz::apply()\t output vector not correct size, at "
					  << v_out.size() << ". System rowdim is" <<  rowdim() << std::endl;
		if ( v_in.size() != coldim() )
			std::cout << "\tToeplitz::apply()\t input vector not correct size at " 
					  << v_in.size() << ". System coldim is" <<  coldim() << std::endl;
		assert((v_out.size() == rowdim()) && 
			   (v_in.size() == coldim()))  ;
		
		NTL::ZZ_pX pxOut, pxIn;
		// bds // pxIn.SetMaxLength( (long) v_in.size()-1);
		pxIn.SetMaxLength( v_in.size()-1);
		for ( size_t  i=0; i< v_in.size(); i++ )
			SetCoeff( pxIn, i, v_in[i]);
		
#ifdef DBGMSGS
		std::cout << "\npX in is " << pxIn << std::endl;
		std::cout << "multiplied by " << pdata << std::endl;
#endif
		mul(pxOut,pxIn,pdata);
		
#ifdef DBGMSGS
		std::cout <<"pxOut is " << pxOut << std::endl;
#endif
		int N = rowdim();
		for ( int i= 0; i < N; i++) 
			GetCoeff(v_out[i], pxOut, N-1+i);
		
		return v_out;
		
	}
	
	
	
	
	/*-----------------------------------------------------------------
	 *    Apply the transposed matrix to a vector Here the input and output 
	 *    vectors are both over the SAME prime ZZ_p field as the 
	 *    Toeplitz matrix itself.
	 *----------------------------------------------------------------*/
	template <class Field, class Vector>
	Vector& Toeplitz<Field, Vector>::applyTranspose( Vector &v_out, 
													 const Vector& v_in) const
	{  
		
		if (v_out.size() != coldim())
			std::cout << "\tToeplitz::applyT()\t output vector not correct size, at "
					  << v_out.size() << ". System coldim is" <<  coldim() << std::endl;
		if ( v_in.size() != rowdim())
			std::cout << "\tToeplitz::applyT()\t input vector not correct size at " 
					  << v_in.size() << ". System rowdim is" <<  rowdim() << std::endl;
		assert((v_out.size() == coldim()) || 
			   (v_in.size() == rowdim()))  ;
		
		NTL::ZZ_pX pxOut, pxIn;
		// bds //pxIn.SetMaxLength( (long) v_in.size()-1);
		pxIn.SetMaxLength( v_in.size()-1);
		
		for (unsigned int i=0; i< v_in.size(); i++)
			SetCoeff( pxIn, i, v_in[i]);
		
#ifdef DBGMSGS
		std::cout << "\npX in is " << pxIn << std::endl;
		std::cout << "multiplied by " << rpdata << std::endl;
#endif
		mul(pxOut,pxIn,rpdata);
		
#ifdef DBGMSGS
		std::cout <<"pxOut is " << pxOut << std::endl;
#endif
		int N = rowdim();
		for ( int i= 0; i < N; i++) 
			GetCoeff(v_out[i], pxOut, N-1+i);
		
		return v_out;
		
		
	}
	
	
	/*-----------------------------------------------------------------
	 *----    Return the number of rows
	 *----------------------------------------------------------------*/
	template <class Field, class Vector>
	inline size_t Toeplitz<Field, Vector>::rowdim() const
	{
		return rowDim;
	}
	
	
	
	
	/*-----------------------------------------------------------------
	 *----    Return the number of columns
	 *----------------------------------------------------------------*/
	template <class Field, class Vector>
	inline size_t Toeplitz<Field, Vector>::coldim() const
	{
		return colDim;
	}
	
	
	
	
	/*-----------------------------------------------------------------
	 *----    Return the max of {rows, cols} as the system dimension of
	 *        a square matrix
	 *----------------------------------------------------------------*/
	
	template <class Field, class Vector>
	inline size_t Toeplitz<Field, Vector>::sysdim() const
	{
		return sysDim;
	}
	
} // namespace LinBox
