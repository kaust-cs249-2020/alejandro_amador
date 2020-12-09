#include"parsiTree.h"


int main()
{
	int n = 20;
	vector<vector<float>> lol = readMatrix("input.txt", n, n, 7.f);
	vector<vector<int>> test = UPGMAcluster(lol);
	vector<vector<int>> finally = transformUPGMA(test, n);
	for (auto& x : finally)
	{
		printVector(x);
	}
	return 0;
}