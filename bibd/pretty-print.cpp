/* pretty-print.cpp
 * Written by Bradford Hovinen <hovinen@ximian.com>
 *
 * Simple utility for pretty-printing a matrix given in a text file
 */

#include <iostream>
#include <fstream>
#include <stdlib.h>

#include "matrix.h"
#include "sparse-vector.h"

int main (int argc, char **argv) 
{
	ifstream input;
	unsigned int colwidth;
	matrix_t *matrix;

	if (argc < 2) {
		cerr << "Usage: pretty-print <filename> [<column width>]"
		     << endl;
		return -1;
	}

	if (argc >= 3)
		colwidth = atoi (argv[2]);
	else
		colwidth = 5;

	input.open (argv[1]);

	if (!input) {
		cerr << "Could not open file " << argv[1] << endl;
		return -1;
	}

	matrix = matrix_read (&sparse_vector_class, input);
	matrix_pretty_print (matrix, colwidth);
	matrix_destroy (matrix);

	input.close ();
	return 0;
}

