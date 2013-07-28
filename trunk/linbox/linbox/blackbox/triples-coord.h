/* Copyright (C) 2013 LinBox
 * Written by AJS <stachnik@udel.edu>
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

/*! @file   linbox/blackbox/triples-coord.h
 * @ingroup linbox/blackbox
 * @brief
 */

#ifndef __LINBOX_TRIPLES_COORD_H
#define __LINBOX_TRIPLES_COORD_H

#include <stdlib.h>
#include <fstream>
#include <limits.h>

namespace LinBox
{

typedef size_t Index;

union TriplesCoord {
        Index rowCol[2];
        Index blockIxArr[2]; // little endian

        TriplesCoord() {};

        TriplesCoord(Index row, Index col) {rowCol[0]=row;rowCol[1]=col;}

        inline TriplesCoord operator--() {
                linbox_check(!((blockIxArr[0]==0)&&(blockIxArr[1]==0)));
                if (blockIxArr[0] == 0) {
                        --blockIxArr[1];
                }
                --blockIxArr[0];
                return *this;
        }

        inline TriplesCoord operator++() {
                ++blockIxArr[0];
                if (blockIxArr[0] == 0) {
                        ++blockIxArr[1];
                }
                return *this;
        }
};

std::ostream& operator<<(std::ostream& out, const TriplesCoord& coord)
{
	out << coord.rowCol[0] << "," << coord.rowCol[1];
	return out;
}

inline TriplesCoord operator>>(const TriplesCoord& coord,unsigned int shift) {
	TriplesCoord retVal;
	retVal.blockIxArr[1]=coord.blockIxArr[1]>>shift;
        retVal.blockIxArr[0]=
                ((((1<<shift)-1)&coord.blockIxArr[1])<<(64-shift))|
                (coord.blockIxArr[0]>>shift);
	return retVal;
}

inline TriplesCoord operator+(const TriplesCoord& lhs,const TriplesCoord& rhs) {
	TriplesCoord retVal;
	retVal.blockIxArr[0]=lhs.blockIxArr[0]+rhs.blockIxArr[0];
	retVal.blockIxArr[1]=lhs.blockIxArr[1]+rhs.blockIxArr[1];
	if (retVal.blockIxArr[0]<lhs.blockIxArr[0]) {
		++(retVal.blockIxArr[1]);
	}
	return retVal;
}

inline bool operator==(const TriplesCoord& lhs,const TriplesCoord& rhs) {
        return (lhs.blockIxArr[1]==rhs.blockIxArr[1]) &&
                (lhs.blockIxArr[0]==rhs.blockIxArr[0]);

}

inline bool operator<(const TriplesCoord& lhs,const TriplesCoord& rhs) {
	return (lhs.blockIxArr[1]<rhs.blockIxArr[1])||
		((!(lhs.blockIxArr[1]>rhs.blockIxArr[1])) &&
		 (lhs.blockIxArr[0]<rhs.blockIxArr[0]));
}

inline TriplesCoord operator-(const TriplesCoord& lhs,const TriplesCoord& rhs) {
	TriplesCoord retVal;
	linbox_check(!(lhs<rhs));
	retVal.blockIxArr[0]=lhs.blockIxArr[0]-rhs.blockIxArr[0];
	retVal.blockIxArr[1]=lhs.blockIxArr[1]-rhs.blockIxArr[1];
	if (lhs.blockIxArr[0]<rhs.blockIxArr[0]) {
		--(retVal.blockIxArr[1]);
	}
	return retVal;
}

void coordFromBlock(TriplesCoord& coord)
{
	Index temp,final;
	Index localRow=0,localCol=0;
	if (coord.blockIxArr[1]!=0) {
		final=coord.blockIxArr[1];
		temp=(final^(final>>1))&0x2222222222222222;
		final^=temp^(temp<<1);
		temp=(final^(final>>2))&0x0C0C0C0C0C0C0C0C;
		final^=temp^(temp<<2);
		temp=(final^(final>>4))&0x00F000F000F000F0;
		final^=temp^(temp<<4);
		temp=(final^(final>>8))&0x0000FF000000FF00;
		final^=temp^(temp<<8);
		temp=(final^(final>>16))&0x00000000FFFF0000;
		final^=temp^(temp<<16);
		localRow=(Index)(final&(~0xFFFFFFFF));
		localCol=(Index)((final<<32)&(~0xFFFFFFFF));
	}
	final=coord.blockIxArr[0];
	temp=(final^(final>>1))&0x2222222222222222;
	final^=temp^(temp<<1);
	temp=(final^(final>>2))&0x0C0C0C0C0C0C0C0C;
	final^=temp^(temp<<2);
	temp=(final^(final>>4))&0x00F000F000F000F0;
	final^=temp^(temp<<4);
	temp=(final^(final>>8))&0x0000FF000000FF00;
	final^=temp^(temp<<8);
	temp=(final^(final>>16))&0x00000000FFFF0000;
	final^=temp^(temp<<16);
	localRow|=(Index)((final>>32)&0xFFFFFFFF);
	localCol|=(Index)(final&0xFFFFFFFF);
	coord.rowCol[0]=localRow;
	coord.rowCol[1]=localCol;
}

void coordToBlock(TriplesCoord& coord)
{
	Index temp,final;
	Index localRow=coord.rowCol[0],localCol=coord.rowCol[1];
	Index highOrderLocalRow=localRow&(~0xFFFFFFFF);
	Index highOrderLocalCol=localCol&(~0xFFFFFFFF);
	Index highOrderFinal;
	coord.blockIxArr[1]=0;
	if ((highOrderLocalRow != 0) || (highOrderLocalCol != 0)) {
		highOrderFinal=localRow|(localCol>>32);
		temp=(highOrderFinal^(highOrderFinal>>16))&0x00000000FFFF0000;
		highOrderFinal^=temp^(temp<<16);
		temp=(highOrderFinal^(highOrderFinal>>8))&0x0000FF000000FF00;
		highOrderFinal^=temp^(temp<<8);
		temp=(highOrderFinal^(highOrderFinal>>4))&0x00F000F000F000F0;
		highOrderFinal^=temp^(temp<<4);
		temp=(highOrderFinal^(highOrderFinal>>2))&0x0C0C0C0C0C0C0C0C;
		highOrderFinal^=temp^(temp<<2);
		temp=(highOrderFinal^(highOrderFinal>>1))&0x2222222222222222;
		highOrderFinal^=temp^(temp<<1);
		coord.blockIxArr[1]=highOrderFinal;
	}
	localCol&=0xFFFFFFFF;
	localRow&=0xFFFFFFFF;
	final=(localRow<<32)|localCol;
	temp=(final^(final>>16))&0x00000000FFFF0000;
	final^=temp^(temp<<16);
	temp=(final^(final>>8))&0x0000FF000000FF00;
	final^=temp^(temp<<8);
	temp=(final^(final>>4))&0x00F000F000F000F0;
	final^=temp^(temp<<4);
	temp=(final^(final>>2))&0x0C0C0C0C0C0C0C0C;
	final^=temp^(temp<<2);
	temp=(final^(final>>1))&0x2222222222222222;
	final^=temp^(temp<<1);
	coord.blockIxArr[0]=final;
}

}

#endif // __LINBOX_TRIPLES_COORD_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
