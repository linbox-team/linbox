%module linboxfunctions
%{
	#include <iostream>

extern bool det(std::istream& matrix_in, std::ostream& det_out);
extern long long detFiles(char *matfile);
extern bool rank(std::istream& matrix_in, std::ostream& rank_out);
extern int rankFiles(char *matfile);
extern bool val(std::istream& matrix_in, std::ostream& val_out);
extern int valFiles(char *matfile);
extern bool trace(std::istream& matrix_in, std::ostream& trace_out);
extern int traceFiles(char *matfile);
%}

extern bool det(std::istream& matrix_in, std::ostream& det_out);
extern long long detFiles(char *matfile);
extern bool rank(std::istream& matrix_in, std::ostream& rank_out);
extern int rankFiles(char *matfile);
extern bool val(std::istream& matrix_in, std::ostream& trace_out);
extern int valFiles(char *matfile);
extern bool trace(std::istream& matrix_in, std::ostream& trace_out);
extern int traceFiles(char *matfile);
