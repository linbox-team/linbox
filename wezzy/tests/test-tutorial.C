#include <cstdlib>
#include <stdio.h>

int main(){
	bool pass = true;
	int d;
	FILE * infile = fopen("tutorial-1.in", "w");
	fprintf(infile, "2 2 d\n");
	fprintf(infile, "3 1\n");
	fprintf(infile, "1 2\n");
	fclose(infile);

	system("./tutorial-1 <tutorial-1.in >tutorial-1.out");

	FILE * outfile = fopen("tutorial-1.out", "r");
	fscanf(outfile, "the determinant is %d", &d);
	fclose(outfile);

	system("rm tutorial-1.in tutorial-1.out");

	pass = pass and d == 5;
	return pass? 0 : -1;
}
