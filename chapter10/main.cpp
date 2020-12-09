#include "profileMarkov.h"

int main()
{
	markov x("input");
	string msg = "HHTT";
	cout << x.viterbiAlgorithm(msg) <<endl;
	return 0;
}