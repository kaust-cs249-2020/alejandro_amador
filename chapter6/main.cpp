#include "genome.h"


int main()
{
	genome t;
	string g1 = t.loadString("data");
	string g2 = t.loadString("data2");
	int k = 30;
	int la = t.numberSharedKmers(g1, g2, k);
	cout << la << endl;
	return 0;
}