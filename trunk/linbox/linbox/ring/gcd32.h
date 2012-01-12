/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
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

#ifndef __LINBOX_ring_gcd32_H
#define __LINBOX_ring_gcd32_H

template<typename Ints> Ints GCD2E32( const Ints p)
{
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

#endif //__LINBOX_ring_gcd32_H
