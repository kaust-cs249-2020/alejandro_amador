#include "align.h"

int main()
{
	alignment test("input", coin);
	string a = "ACACA";
	string b = "CACAC";
	int score;

	//Global
	vector<vector<steps>> back1 = test.generalBackTrack(a, b, 10, score);
	vector<string> re1 = test.generalLCS(back1, a, b, a.size(), b.size());
	cout << score << endl;
	printVectorRows(re1);
	exit(-3);

	//Local
	vector<vector<string>> back2 = test.taxiBackTrack(a, b, 1);
	vector<string> re2 = test.taxiLCS(back2, a, b, a.size(), b.size());
	printVectorRows(re2);

	//Fit
	//vector<vector<string>> back3 = test.fittingBackTrack(a, b, 5, score);
	//vector<string> re3 = test.fittingString(back3, a, b, a.size(), b.size());
	//cout << score << endl;
	//printVectorRows(re3);

	int level;
	for (int i = -10; i < 11; i++)
	{
		vector<vector<vector<string>>> back4 = test.gapBackTrack(a, b, 1, i, score, level);
		cout << i << "," << score << endl;
	}
//	vector<string> lala = test.gapAllignment(back4, a, b, a.size(), b.size(), level);
//	cout << score << endl;
//	printVectorRows(lala);
	return 0;
}