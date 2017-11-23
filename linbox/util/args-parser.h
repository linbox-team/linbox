/*
 * Copyright (C) 2014 the LinBox group
 *
 *  Written by Brice Boyer (briceboyer) <boyer.brice@gmail.com>
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
 *
 */

#ifndef __LINBOX_util_args_parser_H
#define __LINBOX_util_args_parser_H

#include <gmp++/gmp++.h>
#include <fflas-ffpack/utils/args-parser.h>

#include "linbox/util/commentator.h"
namespace LinBox
{
        void parseArguments (int argc, char **argv, Argument *args, bool printDefaults = true) {
        for (int i = 1; i < argc; ++i) {
            if (argv[i][0] == '-') {
                if (argv[i][1] == 0) {
                    LinBox::commentator().setReportStream (std::cout);
                    LinBox::commentator().setBriefReportStream (std::cout);
                } else {
                        // Skip the argument next to "-xxx"
                        // except if next argument is a switch
                    if ( ((i+1) < argc) &&
                         (argv[i+1][0] != '-') ) {
                        ++i;
                    }
                }
            } else {
                LinBox::commentator().setDefaultReportFile (argv[i]);
                LinBox::commentator().setBriefReportStream(std::cout);
            }
        }
        FFLAS::parseArguments(argc,argv,args,printDefaults);
    }
}
#endif // __LINBOX_util_args_parser_H
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
