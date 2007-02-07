#include <iostream>
using namespace std;


// Makes an n by n matrix whose determinant is considerably less than the 
// Hadamard bound.

int main(int argc, char* argv[])
{
    if (argc != 2 ) {
        cerr << "Usage: bigmat <n>, where <n> is the size you like." << endl;
		        return -1;
				    }

	int n = atoi(argv[1]);
	cout << n << " " << n << " M" << endl;
	cout << "1 1 3" << endl;
	cout << "1 " << n << " 2" << endl;
	for (int i = 2; i <=n; ++i) 
	{
		cout << i << " " << i-1 << " " << 2 << endl;
		cout << i << " " << i << " " << 3 << endl;
	}
	cout << "0 0 0" << endl;
}
