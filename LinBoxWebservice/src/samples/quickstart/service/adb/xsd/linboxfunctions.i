%module linboxfunctions
%{
	#include <iostream>
	#include <fstream>

extern bool det(std::istream& matrix_in, std::ostream& det_out);
extern char* detFiles(char *mat);
extern bool rank(std::istream& matrix_in, std::ostream& rank_out);
extern const char* rankFiles(char *mat);
extern int estimateRankTime(char *mat);
extern bool val(std::istream& matrix_in, std::ostream& val_out);
extern const char* valFiles(char *mat);
extern bool trace(std::istream& matrix_in, std::ostream& trace_out);
extern const char* traceFiles(char *mat);
extern bool smithNormalForm(std::istream& matrix_in, std::ostream& snf_out);
extern const char* smithNormalFormFiles(char *mat);
%}

extern bool det(std::istream& matrix_in, std::ostream& det_out);
extern char* detFiles(char *mat);
extern bool rank(std::istream& matrix_in, std::ostream& rank_out);
extern const char* rankFiles(char *mat);
extern int estimateRankTime(char *matfile);
extern bool val(std::istream& matrix_in, std::ostream& trace_out);
extern const char* valFiles(char *mat);
extern bool trace(std::istream& matrix_in, std::ostream& trace_out);
extern const char* traceFiles(char *mat);
extern bool smithNormalForm(std::istream& matrix_in, std::ostream& snf_out);
extern const char* smithNormalFormFiles(char *mat);
