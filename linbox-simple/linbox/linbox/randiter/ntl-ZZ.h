/* File linbox/randiter/ntl-ZZ.h
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
					  
		  _size = NTL::to_ZZ(std::string(size).data());

		  if (seed == integer(0)) NTL::SetSeed (NTL::to_ZZ(time(NULL)));
		  
		  else NTL::SetSeed(NTL::to_ZZ(std::string(seed).data()));
		}
		
		Element& random (Element& x) const {
		
		  if (_size == 0)
		    NTL::RandomLen (x, 30);
		  else
		    NTL::RandomBnd (x, _size);

			return x;
		} 

	private:
		Element _size;
		
	}; 

} // namespace LinBox

#endif // 
