/* File: test-bitonic-sort.C
 *  Author: Zhendong Wan
 */

#include <algorithm>
#include <linbox/algorithms/bitonic-sort.h>
#include <linbox/util/commentator.h>
#include <linbox/vector/stream.h>
#include "test-common.h"

using namespace LinBox;
 
class Comparator {

	public:

	void operator()(int& min, int& max) const {

		if (min > max)

			std::swap (min, max);
	}

};

bool testRandom (int s, int iteration) {
 
	using namespace std;
	
        commentator.start ("test bitonic sort", "test bitonic sort", iteration);

	ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);  

	bool iter_passed = true;
	bool ret = true;

	Comparator comp;

	for (int i = 0; i < iteration; ++ i) {

		commentator.startIteration (i);
	
		std::vector<int> v(s), d(s);

	
		std::vector<int>::iterator p1, p2;

		
		for (p1 = v. begin(), p2 = d.begin(); p1 != v.end(); ++ p1, ++p2) 
			
			*p1 = *p2 = rand();
		
		report << "Input vector:  ";
	

		for (p1 = v.begin(); p1 != v.end(); ++ p1)

			report << *p1 << " ";
		
		report << endl;


		bitonicSort(v.begin(), v.end(), comp);

		

		report << "Computed sequence by bitonic sorting network:\n";
		
		
		for (p1 = v.begin(); p1 != v.end(); ++ p1)

			report << *p1 << " ";
		
		report << '\n';

		stable_sort(d.begin(), d.end());
	
		report << "Expected sequence after sorting:\n";
		
		
		for (p1 = d.begin(); p1 != d.end(); ++ p1)

			report << *p1 << " ";
		
		report << '\n';
	
	

		for (p1 = v.begin(), p2 = d. begin(); p1 != v.end(); ++ p1, ++ p2) 

			if (*p1 != *p2) {
		
				ret = iter_passed = false;

				break;
			}
	
		if (!iter_passed) 
		
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: Computed Smith form is incorrect" << endl;
		
	
		
		commentator.stop ("done");
		
		commentator.progress ();
	
	}
	 
	
	  commentator.stop (MSG_STATUS (ret), (const char *) 0, "testRandom");
                                                                                                        
	  return ret;

}

int main(int argc, char** argv) {
                                                                                                        
        using namespace LinBox;
                                                                                                        
        bool pass = true;
                                                                                                        
        static size_t n = 512;
                                                                                                        
        static int iterations = 10;
                                                                                                        
        static Argument args[] = {
                { 'n', "-n N", "Set size of sequnce to N (default 512, must be power of 2)",  TYPE_INT,     &n },
                { 'i', "-i I", "Perform each test for I iterations (default 10)"
,           TYPE_INT,     &iterations },
        };
                                                                                                        
                                                                                                        
        parseArguments (argc, argv, args);
                                                                                                        
	std::cout << std::endl << "Sort network test suite:\n";

        commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (5);

	if (!testRandom(n, iterations)) pass = false;

        return pass ? 0 : -1;
                                                                                                        
}
