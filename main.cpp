#include "amino.h"

int main()
{
	amino test("ex", totalAlphabet);
	vector<string> lala = test.convolutionCycloPeptideSequencing(20, 1000);
	vector<vector<string>> lo = { lala };
	vector<vector<int>> ka = test.peptidesToMass(lo);
	printPeptideSequence(ka);
	return 0;
}