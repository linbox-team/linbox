/* -*- mode: c; style: linux -*- */

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

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <cstring>

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

	cout << "Usage: " << program << " [options] [<report file>]" << endl;
	cout << endl;
	cout << "Where [options] are the following:" << endl;

	for (i = 0; args[i].c != '\0'; i++) {
		cout << "  " << args[i].example;
		l = 10 - strlen (args[i].example);
		do cout << ' '; while (--l > 0);
		cout << args[i].helpString << endl;
	}

	cout << "  -h or -?  Display this message" << endl;
	cout << endl;
	cout << "If <report file> is not given, then no detailed reporting is done. This is" << endl;
	cout << "suitable if you wish only to determine whether the tests succeeded." << endl;
	cout << endl;
	cout << "[1] N.B. This program does not verify the primality of Q, and does not use a" << endl;
	cout << "    field extension in the event that Q=p^n, n > 1" << endl;
	cout << endl;
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

void parseArguments (int argc, char **argv, ofstream &report, Argument *args)
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
				cerr << "ERROR: Bad argument " << argv[i] << endl;
				break;
			}
		} else {
			report.open (argv[i]);

			if (!report)
				cerr << "ERROR: Cannot open report file " << argv[i] << endl;
			else
				cout << "Writing report data to " << argv[i] << endl << endl;

			cout.flush ();
		}
	}
}

void test_header(char* T, ostream& report){
	cout << "Testing " <<  T <<  "...";
	cout.flush ();
	report << "Testing " << T << ":" << endl;
}

bool test_trailer(bool ret, ostream& report){
	if (ret) {
		cout << "passed" << endl;
		report << "Test passed" << endl << endl;
	} else {
		cout << "FAILED" << endl;
		report << "Test FAILED" << endl << endl;
	}
	cout.flush ();
}


