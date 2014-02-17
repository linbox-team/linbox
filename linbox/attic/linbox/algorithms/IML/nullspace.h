
#ifndef __LINBOX_algorithms_iml_nullspace_H
#define __LINBOX_algorithms_iml_nullspace_H


namespace LinBox { namespace iml {


	//! @todo result can be references if _allocN, etc.. and ~ ok...
	template<class Ring, class Field=Ring>
	class Nullspace {
	private :
		BlasMatrix<Field>   & A ;
		BlasMatrix<Ring>    mp_N ;
		std::vector<size_t> rp ;
		std::vector<size_t> Pt ;
		Ring                R;
		typename Ring::Element mp_D ;
		size_t size ;
		size_t rank ;
		bool   komp ;
		unsigned ver1 ;
		typedef typename Field::Elemene Element ;

	private :
		void reconstruct();
		bool verify(unsigned verified);
		bool verify_compressed(unsigned verified);
		void justTry();
	public :
		Nullspace(Ring myR,BlasMatrix<Field> & myA) :
			A(myA),
			mp_N(R),
			rp(A.rowdim(),0),
			R(myR),
			mp_D(R.one),
			size(0),
			rank(0),
			komp(true),
			ver1(0)
		{}

		// mp_N, rp, mp_D
		size_t getNullspaceCompressed(int verified = 0) ;

		// 0 : no, 1 : proba, 2 : certif
		size_t getNullspace(int verified = 0) ;

		BlasMatrix<Ring> & refNullspace() { return mp_N ; }
		std::vector<size_t> & refP() { return rp ; }
		typename Ring::Element & refD () { return mp_D; }

		size_t getKernel(int verified = 1);


	};

#if 0
	template<class UnparaRing>
	size_t nullspace (
			  BlasMatrix<UnparaRing>        & N
			  ,const BlasMatrix<UnparaRing> & A
			  , const Tag::Side Side=Tag::Right
			 );
template<class ZZ>
size_t nullspace (
		BlasMatrix<ZZ>        & N
		,const BlasMatrix<ZZ> & A
		, const Tag::Side Side
		)
{
	if (Side == Tag::Right)
		return Protected::nullspaceRight(N,A);
	else {
		//! XXX throw error
		return -1 ;
	}
}
#endif



} // IML
} // LinBox


#include "nullspace.inl"

#endif // __LINBOX_algorithms_iml_nullspace_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

