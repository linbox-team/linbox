/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/tests/test-common.C
 * Copyright (C) 2001, 2002 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#include "linbox-config.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <cmath>

#include "linbox/util/commentator.h"
#include "linbox/field/archetype.h"

#include "test-common.h"

using namespace LinBox;

/* Display a help message on command usage */

void printHelpMessage (const char *program, Argument *args) 
{
	int i, l;

	// Skip past libtool prefix in program name
	if (!strncmp (program, "lt-", strlen ("lt-")))
		program += strlen ("lt-");

	std::cout << "Usage: " << program << " [options] [<report file>]" << std::endl;
	std::cout << std::endl;
	std::cout << "Where [options] are the following:" << std::endl;

	for (i = 0; args[i].c != '\0'; i++) {
		std::cout << "  " << args[i].example;
		l = 10 - strlen (args[i].example);
		do std::cout << ' '; while (--l > 0);
		std::cout << args[i].helpString << std::endl;
	}

	std::cout << "  -h or -?  Display this message" << std::endl;
	std::cout << std::endl;
	std::cout << "If <report file> is not given, then no detailed reporting is done. This is" << std::endl;
	std::cout << "suitable if you wish only to determine whether the tests succeeded." << std::endl;
	std::cout << std::endl;
	std::cout << "[1] N.B. This program does not verify the primality of Q, and does not use a" << std::endl;
	std::cout << "    field extension in the event that Q=p^n, n > 1" << std::endl;
	std::cout << std::endl;
}

/* Find an argument in the argument list for a character */

Argument *findArgument (Argument *args, char c) 
{
	int i;

	for (i = 0; args[i].c != '\0' && args[i].c != c; i++);

	if (args[i].c != '\0')
		return &(args[i]);
	else
		return (Argument *) 0;
}

/* Parse command line arguments */

void parseArguments (int argc, char **argv, Argument *args)
{
	int i;
	Argument *current;

	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			if (argv[i][1] == 'h' || argv[i][1] == '?') {
				printHelpMessage (argv[0], args);
				exit (1);
			}
			else if ((current = findArgument (args, argv[i][1])) != (Argument *) 0) {
				switch (current->type) {
				case TYPE_NONE:
					*(bool *) current->data = true;
					break;

				case TYPE_INT:
					*(int *) current->data = atoi (argv[i+1]);
					i++;
					break;

				case TYPE_INTEGER:
					*(integer *) current->data = atoi (argv[i+1]);
					i++;
					break;

				case TYPE_DOUBLE:
					*(double *) current->data = atof (argv[i+1]);
					i++;
					break;
				}
			} else {
				std::cerr << "ERROR: Bad argument " << argv[i] << std::endl;
				break;
			}
		} else {
			commentator.setDefaultReportFile (argv[i]);
			std::cout << "Writing report data to " << argv[i] << std::endl << std::endl;
			std::cout.flush ();
		}
	}
}

bool isPower (integer n, integer m)
{
	return (n == 1) || ((n % m) == 0) && isPower (n/m, m);
}

inline double incompleteGamma (double a, double x, double tol) 
{
	double xa_ex = pow (x, a) * exp (-x);
	double pi = 1.0;
	double xn = 1.0;
	double sigma = 0.0;
	double last_sigma;

	int n = 0;

	do {
		pi *= a + n;
		last_sigma = sigma;
		sigma += xn / pi;
		xn *= x;
		++n;
	} while (abs (sigma - last_sigma) >= tol);

	return sigma * xa_ex;
}

double chiSquaredCDF (double chi_sqr, double df)
{
	return incompleteGamma (df / 2.0, chi_sqr / 2.0, 1e-10) / exp (gamma (df / 2.0));
}

