/* stage1.c
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * Program that given a singular matrix A with rank r, explicitly computes
 * permutation matrices P and Q so that PAQ has a nonsingular leading rxr
 * minor. Also applies the permutation on a vector b for convenience. Intended
 * for sparse computations. Outputs P, Q, PAQ, and Pb
 */

#include <stdlib.h>
#include <iostream>
#include <fstream>

#include "matrix.h"
#include "sparse-vector.h"

int main (int argc, char **argv) 
{
	ifstream input;
	ofstream output;
	matrix_t *A, *PAQ, *P, *Q;
	vector_t *b, *Pb;

	if (argc < 8) {
		cerr << "Usage: stage1 <matrix input> <vector input> "
			"<rank> <matrix output> <vector output> <p output> "
			"<q output>" << endl;
		return -1;
	}

	/* Read A */
	input.open (argv[1]);
	A = matrix_read (&sparse_vector_class, input);
	input.close ();

	/* Read b */
	input.open (argv[2]);
	b = vector_new (&sparse_vector_class, matrix_get_rows (A));
	vector_read (b, input);
	input.close ();

	/* Compute P, Q, PAQ, and Pb */
	PAQ = matrix_nonsingular_minor (A, atoi (argv[3]), &P, &Q);
	Pb = matrix_apply_to_vector (P, &sparse_vector_class, b);

	/* Debug */
#ifdef VERBOSE_DEBUG
	cout << "Orignal matrix:" << endl;
	matrix_pretty_print (A, 1);
	cout << "Final results:" << endl;
	matrix_pretty_print (PAQ, 1);
	matrix_pretty_print (P, 1);
	matrix_pretty_print (Q, 1);
#endif

	/* Output PAQ */
	output.open (argv[4]);
	matrix_write (PAQ, output);
	output.close ();

	/* Output Pb */
	output.open (argv[5]);
	vector_write (Pb, output);
	output.close ();

	/* Output P */
	output.open (argv[6]);
	matrix_write (P, output);
	output.close ();

	/* Output Q */
	output.open (argv[7]);
	matrix_write (Q, output);
	output.close ();

	return 0;
}
