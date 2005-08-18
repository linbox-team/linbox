/* -*- mode: C++; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*- */

/*-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+ 
 *    toeplitz.inl     NTL_Toeplitz.cpp file 
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
#include <assert.h> // JGD 26.09.2003
#include <linbox/algorithms/toeplitz-det.h>

namespace LinBox 
{
	/*-----------------------------------------------------------------
	 *----    Destructor
	 *----------------------------------------------------------------*/
	template <class _CField, class _PField>
	inline ToeplitzBase<_CField,_PField>::~ToeplitzBase()
	{
#ifdef DBGMSGS
		std::cout << "Toeplitz::~Toeplitz():\tDestroyed a " << rowDim << "x"<< colDim<<
			" Toeplitz matrix "<< std::endl;
#endif
	}//---- Destructor ---- [Tested 6/14/02 -- Works]
	
	
	
	/*-----------------------------------------------------------------
	 *----    Field-only Constructor    
	 *----------------------------------------------------------------*/
	template <class _CField, class _PField>
	ToeplitzBase<_CField, _PField>::ToeplitzBase(const Field& F)
	:P(F), K(F)
	{
		shape  =
		sysDim =               // Default dimension is 0
		rowDim =               // Default row dim is 0
		colDim = 0;            // Default col dim is 0
#ifdef DBGMSGS
		std::cout << "Toeplitz::Toeplitz():\tCreated a " << rowDim << "x"<< colDim<<
			" Toeplitz matrix "<< std::endl;
#endif
		
	}//----- Field-only Constructor


	/*-----------------------------------------------------------------
	 *----    Zero Parameter Constructor    
	 *----------------------------------------------------------------*/
	template <class _CField, class _PField>
	ToeplitzBase<_CField, _PField>::ToeplitzBase()
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
	 *------ Polynomial Field constructor
	 *-----------------------------------------------------------------*/
	template< class _CField, class _PField >
	ToeplitzBase<_CField,_PField>::ToeplitzBase( const PField& PF )
		:P(PF), K(PF.getCoeffField())
	{
		shape = sysDim = rowDim = colDim = 0;

	}//------ Polynomial Field constructor

	/*-----------------------------------------------------------------
	 *------ Polynomial constructor
	 *-----------------------------------------------------------------*/
	template< class _CField, class _PField >
	ToeplitzBase<_CField,_PField>::ToeplitzBase
		( const PField& PF, const Poly& p, size_t m, size_t n )
		:P(PF), K(PF.getCoeffField()), rowDim(m), colDim(n), pdata(p)
	{
		shape = 0;
		if( n == 0 ) colDim = rowDim;
		if( rowDim >= colDim ) sysDim = rowDim;
		else sysDim = colDim;

		assert( P.deg(p) <= rowDim + colDim - 2 );
		
		P.rev(rpdata, pdata);

		// Account for possible trailing zeroes
		if( P.deg(pdata) < rowDim + colDim - 2 ) {
			Poly x;
			P.init(x,0);
			Element one;
			K.init(one,1);
			P.setCoeff(x, (rowDim + colDim - 2 - P.deg(pdata)), one);
			P.mulin( rpdata, x );
		}


	}//------ Polynomial constructor

	/*-----------------------------------------------------------------
	 *----- Constructor With User-Supplied First Row And Column
	 *----------------------------------------------------------------*/
	template <class _PField>
        void Toeplitz<typename _PField::CoeffField, _PField>::init_vector( const std::vector<Element>&v)	
        {
		if ( (1 & v.size()) == 0) 
			{
				std::cout << "There must be an ODD number of entries in the input vector " <<
					"The length given is " << v.size();
			}
		assert( (1 & v.size()) == 1);

		P.init(pdata, v);
		P.rev(rpdata, pdata);

		// Account for possible trailing zeroes
		if( P.deg(pdata) < v.size() - 1 ) {
			Poly x;
			P.init(x,0);
			Element one;
			K.init(one,1);
			P.setCoeff(x, (v.size() - 1 - P.deg(pdata)), one);
			P.mulin( rpdata, x );
		}
		
		rowDim = colDim = sysDim = (v.size()+1)/2;
		
		//data = v;
		
#ifdef DBGMSGS
		std::cout << "Toeplitz::Toeplitz(F,V):\tCreated a " << rowDim << "x"<< colDim<<
			" Toeplitz matrix "<< std::endl;
#endif
		
	}//----- Constructor given a vector---- [Tested 6/14/02 -- Works]
	
	

	/*-----------------------------------------------------------------
	 *-----    Print The Matrix To Screen
	 *----------------------------------------------------------------*/
	template <class _PField>
	void Toeplitz<typename _PField::CoeffField,_PField>::print(std::ostream& os) const 
	{
		
		register int i, N;
		register unsigned int j;
		Element temp;
		
		os<< rowDim << " " << colDim << " " << shape << std::endl;
		N = rowDim + colDim -1;

		if ( N < 20 )             // Print small matrices in dense format
			{
				for (i = colDim-1; i < N; i++) 
					{
						for ( j = 0; j < colDim ; j++)
							os << " " ;
							K.write(os,P.getCoeff(temp,pdata,static_cast<size_t>(i-j))) ;
						os << std::endl;
					}
			} 
		else 
			{                    // Print large matrices' first row and col
				os << "[";
				for (size_t i = rowDim + colDim - 2; i> 0;i--)
					K.write(os, P.getCoeff(temp,pdata,i) ) << " ";
				K.write(os,P.getCoeff(temp,pdata,0)) << "]\n";
				P.write(os, pdata) << std::endl;
			} //[v(2n-2),....,v(0)]; where v(0) is the top right entry of the matrix
		
		return;
	} //---- print()----- [Tested 6/14/02 -- Works]
	
	
// 	/*-----------------------------------------------------------------
// 	 *----    The infamous clone has been created here 
// 	 *----------------------------------------------------------------*/
// 	template <class Field, class Vector>
// 	BlackboxArchetype<Vector>* Toeplitz<Field, Vector>::clone() const 
// 	{ 
// 		return new Toeplitz(*this); 
// 	}// ------ This is not tested. 
	
	/*-----------------------------------------------------------------
	 *----    Save To File, Given Destination Filename
	 *----------------------------------------------------------------*/
	template <class _PField>
	void Toeplitz<typename _PField::CoeffField, _PField>::print( char *outFileName) const
	{
		Element temp;

		std::cout << "Printing toeplitz matrix to " << outFileName << std::endl;
		
		if ( outFileName == NULL ) 
			print();    // Print to stdout if no file is specified
		else 
			{
				std::ofstream o_fp(outFileName, std::ios::out);
				o_fp << rowDim << " " << colDim << " " << shape << std::endl ;
				o_fp << "[";
				for (size_t i = rowDim + colDim - 2; i>= 0;i--) 
					K.write(o_fp,P.getCoeff(temp,pdata,i))
					    << " ";
				o_fp << "]\n";
				
				o_fp.close();
			}
		return;
	} // print(char *) [Tested 6/14/02 -- Works]
	
	
	/*-----------------------------------------------------------------
	 *    Make the matrix upper triangular with determinant 1.
	 *    i.e. clear the last N-1 elements in the data vector
	 *----------------------------------------------------------------*/
	template <class _CField, class _PField>
	void ToeplitzBase<_CField, _PField>::setToUniModUT()
	{
		typename PField::Coeff zero, one;
		P.getCoeffField().init(zero,0);
		P.getCoeffField().init(one,1);

		for( size_t i = sysDim; i <= P.deg(pdata); ++i )
		    	P.setCoeff(pdata,i,zero);

		for( size_t i = 0; i < sysDim - 1; ++i )
			P.setCoeff(rpdata,i,zero);

		P.setCoeff(pdata,sysDim-1,one);
		P.setCoeff(rpdata,sysDim-1,one);

		shape = UnimodUT;
		return;
	}// [UNCOMMENTED PART Tested 6/14/02 -- Works]
	
	
	
	/*-----------------------------------------------------------------
	 *    Make matrix a unimodular Lower Triangular with det 1
	 *    i.e. clear the first N-1 elements in the data vector
	 *----------------------------------------------------------------*/
	template <class _CField, class _PField>
	void ToeplitzBase<_CField, _PField>::setToUniModLT()
	{
		typename PField::Coeff zero, one;
		P.getCoeffField().init(zero,0);
		P.getCoeffField().init(one,1);

		for( size_t i = sysDim; i <= P.deg(rpdata); ++i )
		    	P.setCoeff(rpdata,i,zero);

		for( size_t i = 0; i < sysDim - 1; ++i )
			P.setCoeff(pdata,i,zero);

		P.setCoeff(pdata,sysDim-1,one);
		P.setCoeff(rpdata,sysDim-1,one);

		shape = UnimodUT;
		return;
	}// [UNCOMMENTED PART Tested 6/14/02 -- Works]


	/*-----------------------------------------------------------------
	 *     Compute the determinant of the matrix
	 *-----------------------------------------------------------------*/
	template<class _PField>
	typename Toeplitz<typename _PField::CoeffField,_PField>::Element&
		Toeplitz<typename _PField::CoeffField,_PField>::det
		( Element& res ) const
	{
		return toeplitz_determinant( P, res, pdata, sysDim );
	}
	
	
	
	/*-----------------------------------------------------------------
	 *    Apply the matrix to a vector 
	 *----------------------------------------------------------------*/
	template <class _PField>
	template <class OutVector, class InVector>
	OutVector& Toeplitz<typename _PField::CoeffField,_PField>::apply( OutVector &v_out, 
									   const InVector& v_in) const
	{  
		
		if (v_out.size() != rowdim())
			std::cout << "\tToeplitz::apply()\t output vector not correct size, at "
					  << v_out.size() << ". System rowdim is" <<  rowdim() << std::endl;
		if ( v_in.size() != coldim() )
			std::cout << "\tToeplitz::apply()\t input vector not correct size at " 
					  << v_in.size() << ". System coldim is" <<  coldim() << std::endl;
		assert((v_out.size() == rowdim()) && 
			   (v_in.size() == coldim()))  ;
		
		Poly pOut, pIn;
		P.init( pIn, v_in );

#ifdef DBGMSGS
		std::cout << "\npX in is " << pxIn << std::endl;
		std::cout << "multiplied by " << pdata << std::endl;
#endif

		P.mul(pOut, pIn, pdata);

#ifdef DBGMSGS
		std::cout <<"pxOut is " << pxOut << std::endl;
#endif
		
		size_t N = rowdim();
		for( size_t i = 0; i < N; ++i )
			P.getCoeff(v_out[i], pOut, N-1+i);
		
		return v_out;
		
	}
	
	
	
	
	/*-----------------------------------------------------------------
	 *    Apply the transposed matrix to a vector
	 *----------------------------------------------------------------*/
	template <class _PField>
	template<class OutVector, class InVector>
	OutVector& Toeplitz<typename _PField::CoeffField,_PField>::applyTranspose( OutVector &v_out, 
												const InVector& v_in) const
	{  
		
		if (v_out.size() != coldim())
			std::cout << "\tToeplitz::apply()\t output vector not correct size, at "
					  << v_out.size() << ". System rowdim is" <<  coldim() << std::endl;
		if ( v_in.size() != rowdim() )
			std::cout << "\tToeplitz::apply()\t input vector not correct size at " 
					  << v_in.size() << ". System coldim is" <<  rowdim() << std::endl;
		assert((v_out.size() == coldim()) && 
			   (v_in.size() == rowdim()))  ;
		
		Poly pOut, pIn;
		P.init( pIn, v_in );

#ifdef DBGMSGS
		std::cout << "\npX in is " << pxIn << std::endl;
		std::cout << "multiplied by " << rpdata << std::endl;
#endif

		P.mul(pOut, pIn, rpdata);

#ifdef DBGMSGS
		std::cout <<"pxOut is " << pxOut << std::endl;
#endif
		
		size_t N = coldim();
		for( size_t i = 0; i < N; ++i )
			P.getCoeff(v_out[i], pOut, N-1+i);
		
		return v_out;
		
	}
	
} // namespace LinBox
