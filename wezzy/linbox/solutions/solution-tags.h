/* linbox/solutions/methods.h
 *
 * Written by -bds
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 *.
 */

/*! @file solutions/solution-tags.h
 * @ingroup solutions
 * @brief This allows files defining objects that have traits concerning several solutions
 to get the tags and traits with one inclusion.
 */

#ifndef __LINBOX_solution_tags_H
#define __LINBOX_solution_tags_H

namespace LinBox
{
	namespace SolutionTags {
		struct Generic{}; // use for the general case
		struct Local{}; // use if the object has a local function to perform the solution.
	}; // SolutionTags

	template<class BB> struct GetEntryCategory;
	template<class BB> struct TraceCategory;
	template<class BB> struct DetCategory;
	template<class BB> struct RankCategory;
}
#endif // __LINBOX_solution_tags_H


// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

