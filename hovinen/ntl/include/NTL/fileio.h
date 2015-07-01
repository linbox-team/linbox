
#ifndef NTL_fileio__H
#define NTL_fileio__H

#include <NTL/tools.h>

#include <fstream.h>


void OpenWrite(ofstream& s, const char *name);

// opens file for writing...aborts if fails

void OpenRead(ifstream& s, const char *name);

// opens file for reading

char *FileName(const char* stem, const char *ext);

// builds the name "stem.ext"

char *FileName(const char* stem, const char *ext, long d);

// builds the name stem.ext.d


#endif


