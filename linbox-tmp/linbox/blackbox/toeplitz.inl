/*-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
 *    toeplitz.inl     NTL_Toeplitz.cpp file
 *
 *    Copyright (C) 2002 Austin Lobo, B. David Saunders
 *    Author: Austin Lobo
 *    LinBox version 2001 and 2002
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========

 *    This file is included in the template description of ntl-Toeplitz.h
 *    it contains the implementations of templatized member functions in the
 *    partial template  specialization for toeplitz matrices that
 *    are manipulated in fields and rings according to the arithmetic
 *    in the ntl package from V. Shoup
 *-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+*/

#ifndef __LINBOX_bb_toeplitz_INL
#define __LINBOX_bb_toeplitz_INL

#include <iostream>
#include <fstream>
// #include <clinbox_check> // JGD 26.09.2003
#include "linbox/algorithms/toeplitz-det.h"

namespace LinBox
{
	/*-----------------------------------------------------------------
	 *----    Destructor
	 *----------------------------------------------------------------*/
	template <class _CField, class _PRing>
	inline ToeplitzBase<_CField,_PRing>::~ToeplitzBase()
	{
#ifdef DBGMSGS
		std::cout << "Toeplitz::~Toeplitz():\tDestroyed a " << rowdim() << "x"<< this->coldim()<<
		" Toeplitz matrix "<< std::endl;
#endif
	}//---- Destructor ---- [Tested 6/14/02 -- Works]



	/*-----------------------------------------------------------------
	 *----    Field-only Constructor
	 *----------------------------------------------------------------*/
	template <class _CField, class _PRing>
	ToeplitzBase<_CField, _PRing>::ToeplitzBase(const Field& F) :
		P(F), field_(&F)
	{
		sysDim =               // Default dimension is 0
		rowDim =               // Default row dim is 0
		this->colDim = 0;            // Default col dim is 0
		shape.shape(BlackboxSpecifier::TOEPLITZ);
#ifdef DBGMSGS
		std::cout << "Toeplitz::Toeplitz():\tCreated a " << rowDim << "x"<< this->colDim<<
		" Toeplitz matrix "<< std::endl;
#endif

	}//----- Field-only Constructor


#if 0
	/*-----------------------------------------------------------------
	 *----    Zero Parameter Constructor
	 *----------------------------------------------------------------*/
	template <class _CField, class _PRing>
	ToeplitzBase<_CField, _PRing>::ToeplitzBase() :
		P(0), field_(0)
	{
		sysDim =               // Default dimension is 0
		rowDim =               // Default row dim is 0
		this->colDim = 0;            // Default col dim is 0
		shape.shape(BlackboxSpecifier::TOEPLITZ);
#ifdef DBGMSGS
		std::cout << "Toeplitz::Toeplitz():\tCreated a " << rowDim << "x"<< this->colDim<<
		" Toeplitz matrix "<< std::endl;
#endif

	}//----- Zero Param Constructor ---- [Tested 6/14/02 -- Works]
#endif

	template <class _CField, class _PRing>
	ToeplitzBase<_CField, _PRing>::ToeplitzBase(const _PRing& PF) :
		P(PF), field_(&(PF.getCoeffField()))
	{
		sysDim = rowDim = this->colDim = 0;
		shape.shape(BlackboxSpecifier::TOEPLITZ);

	}//------ Polynomial Field constructor

	/*-----------------------------------------------------------------
	 *------ Polynomial constructor
	 *-----------------------------------------------------------------*/
	template< class _CField, class _PRing >
	ToeplitzBase<_CField,_PRing>::ToeplitzBase ( const PRing& PF
						     , const Poly& p
						     , size_t m
						     , size_t n ) :
		P(PF), field_(&(PF.getCoeffField())), rowDim(m), colDim(n), pdata(p)
	{
		shape.shape(BlackboxSpecifier::TOEPLITZ);
		if( n == 0 ) this->colDim = rowDim;
		if( rowDim >= this->colDim ) sysDim = rowDim;
		else sysDim = this->colDim;

		linbox_check( P.deg(p) <= rowDim + this->colDim - 2 );

		P.rev(rpdata, pdata);

		// Account for possible trailing zeroes
		if( P.deg(pdata) < rowDim + this->colDim - 2 ) {
			Poly x;
			P.assign(x,P.zero);
			P.setCoeff(x, (rowDim + this->colDim - 2 - P.deg(pdata)), field().one);
			P.mulin( rpdata, x );
		}


	}//------ Polynomial constructor

	/*-----------------------------------------------------------------
	 *----- Constructor With User-Supplied First Row And Column
	 *----------------------------------------------------------------*/
	template <class _PRing>
	void Toeplitz<typename _PRing::CoeffField, _PRing>::init_vector( const std::vector<Element>&v)
	{
		if ( (1 & v.size()) == 0)
		{
			std::cout << "There must be an ODD number of entries in the input vector " <<
			"The length given is " << v.size();
		}
		linbox_check( (1 & v.size()) == 1);

		this->P.init(this->pdata, v);
		this->P.rev(this->rpdata, this->pdata);

		// Account for possible trailing zeroes
		if( this->P.deg(this->pdata) < v.size() - 1 ) {
			Poly x;
			this->P.assign(x,this->P.zero);
			this->P.setCoeff(x, (v.size() - 1 - this->P.deg(this->pdata)), this->field().one);
			this->P.mulin( this->rpdata, x );
		}

		this->rowDim = this->colDim = this->sysDim = (v.size()+1)/2;

		//data = v;

#ifdef DBGMSGS
		std::cout << "Toeplitz::Toeplitz(F,V):\tCreated a " << rowDim << "x"<< this->colDim<<
		" Toeplitz matrix "<< std::endl;
#endif

	}//----- Constructor given a vector---- [Tested 6/14/02 -- Works]



	/*-----------------------------------------------------------------
	 *-----    Print The Matrix To Screen
	 *----------------------------------------------------------------*/
	template <class _PRing>
	std::ostream & Toeplitz<typename _PRing::CoeffField,_PRing>::write(std::ostream& os) const
	{

		int N;
		Element temp;

		os<< this->rowdim() << " " << this->coldim() << " " << this->shape.shape() << std::endl;
		N = (int) (this->rowdim() + this->coldim()) -1;

		if ( N < 20 )             // Print small matrices in dense format
		{
			int i;
			unsigned int j;
			for (i = (int)this->coldim()-1; i < N; i++)
			{
				for ( j = 0; j < this->coldim() ; j++)
					os << " " ;
				this->field().write(os,this->P.getCoeff(temp,this->pdata,static_cast<size_t>((unsigned int)i-j))) ;
				os << std::endl;
			}
		}
		else
		{                    // Print large matrices' first row and col
			os << "[";
			for (size_t ii = this->rowdim() + this->coldim() - 2; ii> 0;ii--)
				this->field().write(os, this->P.getCoeff(temp,this->pdata,ii) ) << " ";
			this->field().write(os,this->P.getCoeff(temp,this->pdata,0)) << "]\n";
			this->P.write(os, this->pdata) << std::endl;
		} //[v(2n-2),....,v(0)]; where v(0) is the top right entry of the matrix

		return os;
	} //---- write()----- [Tested 6/14/02 -- Works]


#if 0
	// 	/*-----------------------------------------------------------------
	// 	 *----    The infamous clone has been created here
	// 	 *----------------------------------------------------------------*/
	template <class Field, class Vector>
	BlackboxArchetype<Vector>* Toeplitz<Field, Vector>::clone() const
	{
		return new Toeplitz(*this);
	}// ------ This is not tested.
#endif

	/*-----------------------------------------------------------------
	 *----    Save To File, Given Destination Filename
	 *----------------------------------------------------------------*/
	template <class _PRing>
	void Toeplitz<typename _PRing::CoeffField, _PRing>::write( char *outFileName) const
	{
		Element temp;

		std::cout << "Printing toeplitz matrix to " << outFileName << std::endl;

		if ( outFileName == NULL )
			write();    // Print to stdout if no file is specified
		else
		{
			std::ofstream o_fp(outFileName, std::ios::out);
			o_fp << this->rowdim() << " " << this->coldim() << " " << this->shape.shape() << std::endl ;
			o_fp << "[";
			for (size_t i = this->rowdim() + this->coldim() - 1 ; i-- ; )
				this->field().write(o_fp,this->P.getCoeff(temp,this->pdata,i))
				<< " ";
			o_fp << "]\n";

			o_fp.close();
		}
		return;
	} // write(char *) [Tested 6/14/02 -- Works]


#if 0 //dated material with no known use -bds 2012Jul
	/*-----------------------------------------------------------------
	 *    Make the matrix upper triangular with determinant 1.
	 *    i.e. clear the last N-1 elements in the data vector
	 *----------------------------------------------------------------*/
	template <class _CField, class _PRing>
	void ToeplitzBase<_CField, _PRing>::setToUniModUT()
	{
		const PRing & PF = P.getCoeffField() ;

		for( size_t i = sysdim(); i <= P.deg(pdata); ++i )
			P.setCoeff(pdata,i,PF.zero);

		for( size_t i = 0; i < sysDim - 1; ++i )
			P.setCoeff(rpdata,i,PF.zero);

		P.setCoeff(pdata,sysDim-1,PF.one);
		P.setCoeff(rpdata,sysDim-1,PF.one);

		shape.shape(BlackboxSpecifier::UNIMOD_UT);
		return;
	}// [UNCOMMENTED PART Tested 6/14/02 -- Works]



	/*-----------------------------------------------------------------
	 *    Make matrix a unimodular Lower Triangular with det 1
	 *    i.e. clear the first N-1 elements in the data vector
	 *----------------------------------------------------------------*/
	template <class _CField, class _PRing>
	void ToeplitzBase<_CField, _PRing>::setToUniModLT()
	{
		const PRing & PF = P.getCoeffField() ;

		for( size_t i = sysDim; i <= P.deg(rpdata); ++i )
			P.setCoeff(rpdata,i,PF.zero);

		for( size_t i = 0; i < sysDim - 1; ++i )
			P.setCoeff(pdata,i,PF.zero);

		P.setCoeff(pdata,sysDim-1,PF.one);
		P.setCoeff(rpdata,sysDim-1,PF.one);

		shape.shape(BlackboxSpecifier::UNIMOD_LT);
		return;
	}// [UNCOMMENTED PART Tested 6/14/02 -- Works]

#endif

	/*-----------------------------------------------------------------
	 *     Compute the determinant of the matrix
	 *-----------------------------------------------------------------*/
	template<class _PRing>
	typename Toeplitz<typename _PRing::CoeffField,_PRing>::Element&
	Toeplitz<typename _PRing::CoeffField,_PRing>::det
	( Element& res ) const
	{
		return toeplitz_determinant( this->P, res, this->pdata, this->sysDim );
	}

	/*-----------------------------------------------------------------
	 *     Compute the trace of the matrix
	 *-----------------------------------------------------------------*/
	/*-----------------------------------------------------------------
	 *     Compute the trace of the matrix
	 *-----------------------------------------------------------------*/
	/*-----------------------------------------------------------------
	 *     Compute the trace of the matrix
	 *-----------------------------------------------------------------*/
	template<class _PRing>
	typename Toeplitz<typename _PRing::CoeffField,_PRing>::Element&
	Toeplitz<typename _PRing::CoeffField,_PRing>::trace
	( Element& res ) const
	{
		Element x; this->field().init(x, min(rowdim(), coldim()));
		this->P.getCoeff(res, pdata, rowdim()-1);
		return this->field().mulin(res, x);
	}



	/*-----------------------------------------------------------------
	 *    Apply the matrix to a vector
	 *----------------------------------------------------------------*/
	template <class _PRing>
	template <class OutVector, class InVector>
	OutVector& Toeplitz<typename _PRing::CoeffField,_PRing>::apply( OutVector &v_out,
									const InVector& v_in) const
	{

		linbox_check((v_out.size() == this->rowdim()) &&
			     (v_in.size() == this->coldim()))  ;

		Poly pOut, pIn;
		this->P.init( pIn, v_in );

#ifdef DBGMSGS
		std::cout << "\npX in is " << pxIn << std::endl;
		std::cout << "multiplied by " << this->pdata << std::endl;
#endif

		this->P.mul(pOut, pIn, this->pdata);

#ifdef DBGMSGS
		std::cout <<"pxOut is " << pxOut << std::endl;
#endif

		size_t N = this->rowdim();
		for( size_t i = 0; i < N; ++i )
			this->P.getCoeff(v_out[i], pOut, N-1+i);

		return v_out;

	}




	/*-----------------------------------------------------------------
	 *    Apply the transposed matrix to a vector
	 *----------------------------------------------------------------*/
	template <class _PRing>
	template<class OutVector, class InVector>
	OutVector& Toeplitz<typename _PRing::CoeffField,_PRing>::applyTranspose( OutVector &v_out,
										 const InVector& v_in) const
	{

		if (v_out.size() != this->coldim())
			std::cout << "\tToeplitz::apply()\t output vector not correct size, at "
			<< v_out.size() << ". System rowDim is" <<  this->coldim() << std::endl;
		if ( v_in.size() != this->rowdim() )
			std::cout << "\tToeplitz::apply()\t input vector not correct size at "
			<< v_in.size() << ". System colDim is" <<  this->rowdim() << std::endl;
		linbox_check((v_out.size() == this->coldim()) &&
			     (v_in.size() == this->rowdim()))  ;

		Poly pOut, pIn;
		this->P.init( pIn, v_in );

#ifdef DBGMSGS
		std::cout << "\npX in is " << pxIn << std::endl;
		std::cout << "multiplied by " << this->rpdata << std::endl;
#endif

		this->P.mul(pOut, pIn, this->rpdata);

#ifdef DBGMSGS
		std::cout <<"pxOut is " << pxOut << std::endl;
#endif

		size_t N = this->coldim();
		for( size_t i = 0; i < N; ++i )
			this->P.getCoeff(v_out[i], pOut, N-1+i);

		return v_out;

	}

} // namespace LinBox

#endif //__LINBOX_bb_toeplitz_INL

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
