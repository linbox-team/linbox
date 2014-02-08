#include "linbox-config.h"
#include "givaro/givpoly1.h"
#include "givaro/giv_randiter.h"

#ifndef __LINBOX_randiter_givaro_poly_H
#define __LINBOX_randiter_givaro_poly_H

namespace LinBox
{
	template<class Field>
	class GivaroPolyRandIter {
		typedef typename Field::Domain_t Domain;
		typedef typename Domain::Domain_t SubDomain;
		
		Field _pd;
		Givaro::GIV_randIter<SubDomain,integer> _randIter;
		
		integer _size;
		integer _seed;
	public:
		typedef typename Field::Element Element;
		
		GivaroPolyRandIter(Field pd, 
			const integer& size = 0,
			const integer& seed = 0) 
		{
			_pd = pd;
			_randIter = Givaro::GIV_randIter<SubDomain,integer>(pd.subdomain(), size, seed);
		}
		
		GivaroPolyRandIter(const GivaroPolyRandIter &R)
			: _pd(R._pd), _randIter(R._randIter), _size(R._size), _seed(R._seed) {}
		
		GivaroPolyRandIter &operator=(const GivaroPolyRandIter &R) {
			return *this;
		}
		
		Element &random(Element &a) const {
			return _pd.domain().random(_randIter, a);
		}
		
		Element &random(Element &a, Givaro::Degree d) const {
			return _pd.domain().random(_randIter, a, d);
		}
		
		Element &random(Element &a, long s) const {
			return _pd.domain().random(_randIter, a, s);
		}
	};
}

#endif // __LINBOX_randiter_givaro_poly_H
