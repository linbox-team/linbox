/** File linbox/randiter/ntl-ZZ.h
 *  Author: Zhendong Wan
 */

#ifndef ___RANDITER_NTL_ZZ_H__
#define ___RANDITER_NTL_ZZ_H__

#include <vector>

namespace LinBox {

	class NTL_ZZ;

	class NTL_ZZRandIter {
	public:
		typedef NTL::ZZ  Element;

		NTL_ZZRandIter (const NTL_ZZ& F, 
				const integer& size = 0, 
				const integer& seed = 0
				) {
			
			
			if (seed == integer(0)) NTL::SetSeed (NTL::to_ZZ(time(0)));

			else NTL::SetSeed(NTL::to_ZZ((unsigned int) seed));
		}
		
		Element& random (Element& x) const {
		
			NTL::RandomLen (x, 30);

			return x;
		} 
	}; 

} // namespace LinBox

#endif // 
