/* Copyright (C) 2013 LinBox
 * Written by Brice Boyer (briceboyer) <boyer.brice@gmail.com>
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

/*! @file   benchmarks/optimizer.h
 * @ingroup benchmarks
 *
 */

#ifndef __LINBOX_benchmarks_optimizer_H_
#define __LINBOX_benchmarks_optimizer_H_

namespace LinBox {

	/* optimiser from a graph. */

	/*
	 *   -----
	 *   if (a < threshold) then
	 *        toto()
	 *   else
	 *        titi()
	 *   ----
	 *
	 * OR (always better)
	 *
	 *   ----
	 *   toto();
	 *   ----
	 */

	// Optimizer opt(Data,Skeleton);
	// opt.run();
	// opt.report("file");

	/* optimiser from two timers */
	// Optimizer opt(class1,class2,Skeleton);
	// opt.run();
	// opt.report("file")
	//
	// opt.fit()
	//
	//
	class Optimizer {



		void report(const std::string & filename) ;

		void run() ;

		void fit();

	}
}

#endif // __LINBOX_benchmarks_optimizer_H_

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
