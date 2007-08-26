#include <iostream>
#include <linbox/util/commentator.h>
#include "test-common.h"
#include <linbox/vector/subiterator.h>
using namespace LinBox;
bool test() 
;
int main(int argc, char** argv)
{	
   static Argument args[] = {
   };
   parseArguments (argc, argv, args);
	commentator.start("Subiterator test suite", "Subiterator");

   bool pass = test();
	commentator.stop("Subiterator test suite");
   return pass? 0 : -1;
}

bool test()
{
	char* title = "Subiterator test";
	commentator.start(title, title, 1);
	ostream &report = commentator.report 
		(Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	std::vector<int> v; 
	for (int i = 1; i <= 10; ++i) v.push_back(i);

	Subiterator<std::vector<int>::iterator> s(v.begin()+2, 2);
	Subiterator<std::vector<int>::iterator> t(v.end()+2, 2);
	bool res = true;
	if (*s != 3) {report << 1 << endl; res=false;}
	if (*(++s) != 5) {report << 2 << endl; res=false;}
	if (*(s++) != 5) {report << 3 << endl; res=false;}
	if (*(s) != 7) {report << 4 << endl; res=false;}
	if (*(--s) != 5) {report << 5 << endl; res=false;}
	if (*(s--) != 5) {report << 6 << endl; res=false;}
	if (*(s) != 3) {report << 7 << endl; res=false;}
	if (*(s+1) != 5) {report << 8 << endl; res=false;}
	if (*(s-1) != 1) {report << 9 << endl; res=false;}
	if ((t-s) != 5) {report << 10 << endl; res=false;}
	//	report << "10 " << t-s << endl; 
	if ((s[1]) != 5) {report << 11 << endl; res=false;}

	if (s == t) {report << 12 << endl; res=false;}
	if (!(s != t)) {report << 13 << endl; res=false;}
	if (!(s <= t)) {report << 14 << endl; res=false;}
	if (s >= t) {report << 15 << endl; res=false;}
	if (!(s < t)) {report << 16 << endl; res=false;}
	if (s > t) {report << 17 << endl; res=false;}

	*s = 0;
	if (v[2] != 0) {report << 18 << endl; res=false;}
	s[3]=0;
	if (v[8] != 0) {report << 19 << endl; res=false;}
	std::swap(s, t);
	if (s < t) {report << 20 << endl; res=false;}
	s.swap(t);
	if (s > t) {report << 21 << endl; res=false;}


	commentator.stop (MSG_STATUS (res), (const char *) 0, title);
	return res;
}

