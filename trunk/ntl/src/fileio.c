
#include <NTL/fileio.h>

#include <string.h>
#include <stdlib.h>


void OpenWrite(ofstream& s, const char *name)
{
   s.open(name, ios::out);

   if (!s) {
      cerr << "open write error: " << name << "\n";
      abort();
   }
}


void OpenRead(ifstream& s, const char *name)
{
   s.open(name, ios::in);
   if (!s) {
      cerr << "open read error: " << name << "\n";
      abort();
   }
}

static char sbuf[256];

char *FileName(const char* stem, const char *ext)
{
   strcpy(sbuf, stem);
   strcat(sbuf, "-");
   strcat(sbuf, ext);

   return sbuf;
}

char *FileName(const char* stem, const char *ext, long d)
{
   strcpy(sbuf, stem);
   strcat(sbuf, "-");
   strcat(sbuf, ext);
   strcat(sbuf, "-");

   char dbuf[6];
   dbuf[5] = '\0';
   long i, dd;
   dd = d;
   for (i = 4; i >= 0; i--) {
      dbuf[i] = (char) ((dd % 10) + '0');
      dd = dd / 10;
   }

   strcat(sbuf, dbuf);

   return sbuf;
}
