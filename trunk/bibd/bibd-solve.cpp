/* bibd-solve.cpp
 *
 * Find a Moore-Penrose solution for the given BIBD problem
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
 */

#include <stdlib.h>
#include <iostream>
#include <fstream>

#include <string.h>

#include "bibd-matrix.h"
#include "moore-penrose.h"

static unsigned int extract_list (char *str, unsigned int *array, unsigned int maxlen) 
{
	unsigned int i = 0;

	while (*str && i < maxlen) {
		array[i++] = atoi (str);
		while (isdigit (*str)) str++;
		if (*str != '\0') str++;
	}

	return i;
}

typedef struct _report_state_t report_state_t;

struct _report_state_t
{
	unsigned int *subset_array;
	vector_t *soln;
	ostream *output;
};

static unsigned int *merge_lists (unsigned int *list1, unsigned int len1,
				  unsigned int *list2, unsigned int len2) 
{
	unsigned int *merged_list;
	unsigned int i = 0, j = 0, k = 0;

	merged_list = new unsigned int[len1 + len2];

	while (i < len1 && j < len2) {
		while (i < len1 && j < len2 && list1[i] <= list2[j])
			merged_list[k++] = list1[i++];
		while (j < len2 && i < len1 && list1[i] > list2[j])
			merged_list[k++] = list2[j++];
	}

	while (i < len1)
		merged_list[k++] = list1[i++];
	while (j < len2)
		merged_list[k++] = list2[j++];

	return merged_list;
}

static void modify_subset_array (unsigned int *subset_array,
				 unsigned int subset_array_len,
				 unsigned int *includes_list,
				 unsigned int includes_len,
				 unsigned int *excludes_list,
				 unsigned int excludes_len) 
{
	unsigned int *merged_list;
	int i;

	merged_list = merge_lists (includes_list, includes_len,
				   excludes_list, excludes_len);

	for (i = includes_len + excludes_len - 1; i >= 0; i--)
		memmove (subset_array + merged_list[i],
			 subset_array + merged_list[i] + 1,
			 (subset_array_len - merged_list[i] - 1) * sizeof (unsigned int));

	delete[] merged_list;
}

static void print_subset (ostream &out, unsigned int bitmap) 
{
	unsigned int i = 1;

	out << "(";
	while (bitmap != 0) {
		if (bitmap & 1) {
			out << i;
			if ((bitmap >> 1) != 0) out << " ";
		}

		bitmap >>= 1;
		i++;
	}
	out << ")";
}

static int report_cb (vector_t *vector, unsigned int column, Element &entry,
		      report_state_t *state) 
{
	state->output->width (3);
	(*state->output) << column << " ";
	print_subset (*state->output, state->subset_array[column]);
	(*state->output) << " ";
	field.write (*state->output, entry);
	(*state->output) << endl;
	return 0;
}

int main (int argc, char **argv) 
{
	int v, k, lambda;
	matrix_t *bibd_matrix, *bibd_matrix_1;
	vector_t *bibd_vector, *bibd_vector_1;
	vector_t *soln;
	unsigned int *subset_array, array_len;
	ofstream out;
	report_state_t state;

	unsigned int includes_list[256];
	unsigned int excludes_list[256];
	unsigned int includes_len = 0, excludes_len = 0;

	if (argc < 5) {
		cerr << "Usage: bibd-solve v k lambda output [includes list] [excludes list]" << endl;
		return -1;
	}

	v = atoi (argv[1]);
	k = atoi (argv[2]);
	lambda = atoi (argv[3]);

	if (argc >= 6) includes_len = extract_list (argv[5], includes_list, 256);
	if (argc >= 7) excludes_len = extract_list (argv[6], excludes_list, 256);

	bibd_matrix = make_bibd_matrix (v, k);
	bibd_vector = make_bibd_vector (v, k, lambda);
	bibd_matrix_1 = refine_bibd_matrix (bibd_matrix, 
					    includes_list, includes_len,
					    excludes_list, excludes_len);
	bibd_vector_1 = refine_bibd_vector (bibd_vector, includes_list, includes_len, bibd_matrix);
	matrix_destroy (bibd_matrix);
	vector_unref (bibd_vector);
	soln = full_moore_penrose (bibd_matrix_1, bibd_vector_1, matrix_get_rows (bibd_matrix_1));

	subset_array = get_subset_array (v, k, &array_len);

	out.open (argv[4]);
	state.subset_array = subset_array;
	state.soln = soln;
	state.output = &out;
	modify_subset_array (subset_array, array_len,
			     includes_list, includes_len,
			     excludes_list, excludes_len);
	vector_foreach_entry (soln, (entry_cb_t) report_cb, 0, &state);
	out.close ();

	return 0;
}
