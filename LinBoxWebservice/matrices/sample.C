#include <iostream>
using namespace std;
int main(int argc, char* argv[])
{
	int n = 10;
	if (argc == 1) cerr << "Usage: " << argv[0] << " n" << endl << "to generate a (0..9)-random n by n dense matrix." << endl;
	if (argc > 1) n = atoi(argv[1]); 
	cout << n << " " << n << " D" << endl;
	for(int i = 0; i < n; ++i) { 
		for(int j = 0; j < n; ++j) cout << rand()%10 << " ";
		cout << endl;
	}
}
