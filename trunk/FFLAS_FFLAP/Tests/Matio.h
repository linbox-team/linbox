#include <stdio.h>
#include <stdlib.h>

// Reading and writing matrices over double

// Reading a matrice from a (eventually zipped) file
double * read_dbl(char * mat_file,int* tni,int* tnj) {
  char *UT, *File_Name;
  int is_gzipped = 0;
  size_t s = strlen(mat_file);
  double* X;
  if ((mat_file[--s] == 'z') && 
      (mat_file[--s] == 'g') && 
      (mat_file[--s] == '.')) {
    is_gzipped = 1;
    File_Name = "/tmp/bbXXXXXX_";
    mkstemp(File_Name);
    UT = new char[s+34+strlen(File_Name)];
    sprintf(UT,"gunzip -c %s > %s", mat_file, File_Name);
    system(UT);
    sprintf(UT,"\\rm %s", File_Name);
  } else
    File_Name = mat_file;
  
  FILE* FileDes = fopen(File_Name, "r");
  if (FileDes != NULL) {
    char * tmp = new char[200];// usigned long tni, tnj;
    fscanf(FileDes,"%d %d %s\n",tni, tnj, &tmp) ;
    int n=*tni;
    int p=*tnj;
    X = new double[n*p];
    for (int i=0;i<n*p;++i)
      X[i] = (double) 0;
    long i,j; long val;
    fscanf(FileDes,"%ld %ld %ld\n",&i, &j, &val) ;
    while(i && j) {
      X[p*(i-1)+j-1] = (double) val;
      fscanf(FileDes,"%ld %ld %ld\n",&i, &j, &val) ;
    }
  }
  
  fclose(FileDes);
  if (is_gzipped) system(UT);        
  return X;
}

// Displays a matrix 
std::ostream& write_dbl(std::ostream& c, 
		   double* E,
		   int n, int m, int id){
  
  for (int i = 0; i<n;++i){
    for (int j=0; j<m;++j)
      c << *(E+j+id*i) << " ";
    c << std::endl;
  }
  c << std::endl;
}

// Reading and writing matrices over field

// Reading a matrice from a (eventually zipped) file
template<class Field>
typename Field::element * read_field(const Field& F,char * mat_file,int* tni,int* tnj)
{
  char *UT, *File_Name;
  int is_gzipped = 0;
  size_t s = strlen(mat_file);
 
  typename Field::element * X;
  if ((mat_file[--s] == 'z') && 
      (mat_file[--s] == 'g') && 
      (mat_file[--s] == '.')) {
    is_gzipped = 1;
    File_Name = "/tmp/bbXXXXXX_";
    mkstemp(File_Name);
    UT = new char[s+34+strlen(File_Name)];
    sprintf(UT,"gunzip -c %s > %s", mat_file, File_Name);
    system(UT);
    sprintf(UT,"\\rm %s", File_Name);
  } else
    File_Name = mat_file;
  FILE* FileDes = fopen(File_Name, "r");
  if (FileDes != NULL) {
    char * tmp = new char[200];// usigned long tni, tnj;
    fscanf(FileDes,"%d %d %s\n",tni, tnj, &tmp) ;
    int n=*tni;
    int p=*tnj;
    X = new typename Field::element[n*p];
    for (int i=0;i<n*p;++i)
      X[i] = F.zero;
    long i,j; long val;
    fscanf(FileDes,"%ld %ld %ld\n",&i, &j, &val) ;
    while(i && j) {
      F.init(X[p*(i-1)+j-1],val);
      fscanf(FileDes,"%ld %ld %ld\n",&i, &j, &val) ;
    }
  }
  
  fclose(FileDes);
  if (is_gzipped) system(UT);        
  return X;
}

// Displays a matrix
template<class Field>
std::ostream& write_field(const Field& F,std::ostream& c, 
		     const typename Field::element* E,
		     int n, int m, int id){
  
  long tmp;
  for (int i = 0; i<n;++i){
    for (int j=0; j<m;++j){
      F.convert(tmp,*(E+j+id*i));
      c << tmp << " ";}
    c << std::endl;
  }
  c << std::endl;
}
