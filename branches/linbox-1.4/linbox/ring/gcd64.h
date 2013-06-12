/* Copyright (C) LinBox
 *
 *
 *
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 */

#ifndef __LINBOX_ring_gcd64_H
#define __LINBOX_ring_gcd64_H

template<typename Ints> Ints GCD2E64( const Ints p)
{
	if( p & 4294967295UL) {
		if( p & 65535UL) {
			if( p & 255UL) {
				if( p & 15UL) {
					if( p & 3UL) {
						if( p & 1UL) {
							return 1UL;
						}
						else {
							return 2UL;
						}
					}
					else {
						if( p & 7UL) {
							return 4UL;
						}
						else {
							return 8UL;
						}
					}
				}
				else {
					if( p & 63UL) {
						if( p & 31UL) {
							return 16UL;
						}
						else {
							return 32UL;
						}
					}
					else {
						if( p & 127UL) {
							return 64UL;
						}
						else {
							return 128UL;
						}
					}
				}
			}
			else {
				if( p & 4095UL) {
					if( p & 1023UL) {
						if( p & 511UL) {
							return 256UL;
						}
						else {
							return 512UL;
						}
					}
					else {
						if( p & 2047UL) {
							return 1024UL;
						}
						else {
							return 2048UL;
						}
					}
				}
				else {
					if( p & 16383UL) {
						if( p & 8191UL) {
							return 4096UL;
						}
						else {
							return 8192UL;
						}
					}
					else {
						if( p & 32767UL) {
							return 16384UL;
						}
						else {
							return 32768UL;
						}
					}
				}
			}
		}
		else {
			if( p & 16777215UL) {
				if( p & 1048575UL) {
					if( p & 262143UL) {
						if( p & 131071UL) {
							return 65536UL;
						}
						else {
							return 131072UL;
						}
					}
					else {
						if( p & 524287UL) {
							return 262144UL;
						}
						else {
							return 524288UL;
						}
					}
				}
				else {
					if( p & 4194303UL) {
						if( p & 2097151UL) {
							return 1048576UL;
						}
						else {
							return 2097152UL;
						}
					}
					else {
						if( p & 8388607UL) {
							return 4194304UL;
						}
						else {
							return 8388608UL;
						}
					}
				}
			}
			else {
				if( p & 268435455UL) {
					if( p & 67108863UL) {
						if( p & 33554431UL) {
							return 16777216UL;
						}
						else {
							return 33554432UL;
						}
					}
					else {
						if( p & 134217727UL) {
							return 67108864UL;
						}
						else {
							return 134217728UL;
						}
					}
				}
				else {
					if( p & 1073741823UL) {
						if( p & 536870911UL) {
							return 268435456UL;
						}
						else {
							return 536870912UL;
						}
					}
					else {
						if( p & 2147483647UL) {
							return 1073741824UL;
						}
						else {
							return 2147483648UL;
						}
					}
				}
			}
		}
	}
	else {
		if( p & 281474976710655ULL) {
			if( p & 1099511627775ULL) {
				if( p & 68719476735ULL) {
					if( p & 17179869183ULL) {
						if( p & 8589934591ULL) {
							return 4294967296ULL;
						}
						else {
							return 8589934592ULL;
						}
					}
					else {
						if( p & 34359738367ULL) {
							return 17179869184ULL;
						}
						else {
							return 34359738368ULL;
						}
					}
				}
				else {
					if( p & 274877906943ULL) {
						if( p & 137438953471ULL) {
							return 68719476736ULL;
						}
						else {
							return 137438953472ULL;
						}
					}
					else {
						if( p & 549755813887ULL) {
							return 274877906944ULL;
						}
						else {
							return 549755813888ULL;
						}
					}
				}
			}
			else {
				if( p & 17592186044415ULL) {
					if( p & 4398046511103ULL) {
						if( p & 2199023255551ULL) {
							return 1099511627776ULL;
						}
						else {
							return 2199023255552ULL;
						}
					}
					else {
						if( p & 8796093022207ULL) {
							return 4398046511104ULL;
						}
						else {
							return 8796093022208ULL;
						}
					}
				}
				else {
					if( p & 70368744177663ULL) {
						if( p & 35184372088831ULL) {
							return 17592186044416ULL;
						}
						else {
							return 35184372088832ULL;
						}
					}
					else {
						if( p & 140737488355327ULL) {
							return 70368744177664ULL;
						}
						else {
							return 140737488355328ULL;
						}
					}
				}
			}
		}
		else {
			if( p & 72057594037927935ULL) {
				if( p & 4503599627370495ULL) {
					if( p & 1125899906842623ULL) {
						if( p & 562949953421311ULL) {
							return 281474976710656ULL;
						}
						else {
							return 562949953421312ULL;
						}
					}
					else {
						if( p & 2251799813685247ULL) {
							return 1125899906842624ULL;
						}
						else {
							return 2251799813685248ULL;
						}
					}
				}
				else {
					if( p & 18014398509481983ULL) {
						if( p & 9007199254740991ULL) {
							return 4503599627370496ULL;
						}
						else {
							return 9007199254740992ULL;
						}
					}
					else {
						if( p & 36028797018963967ULL) {
							return 18014398509481984ULL;
						}
						else {
							return 36028797018963968ULL;
						}
					}
				}
			}
			else {
				if( p & 1152921504606846975ULL) {
					if( p & 288230376151711743ULL) {
						if( p & 144115188075855871ULL) {
							return 72057594037927936ULL;
						}
						else {
							return 144115188075855872ULL;
						}
					}
					else {
						if( p & 576460752303423487ULL) {
							return 288230376151711744ULL;
						}
						else {
							return 576460752303423488ULL;
						}
					}
				}
				else {
					if( p & 4611686018427387903ULL) {
						if( p & 2305843009213693951ULL) {
							return 1152921504606846976ULL;
						}
						else {
							return 2305843009213693952ULL;
						}
					}
					else {
						if( p & 9223372036854775807ULL) {
							return 4611686018427387904ULL;
						}
						else {
							return 9223372036854775808ULL;
						}
					}
				}
			}
		}
	}
}

#endif //__LINBOX_ring_gcd64_H

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

