#pragma once
#include "directedGraph.h"
#include "weightedGraph.h"

enum problemType
{
	coin,
	other
};

enum steps
{
	down,
	rights,
	diagonal,
	jumpStart,
	jumpEnd
};

enum multiSteps
{
	moveX,
	moveY,
	moveZ,
	diagX,
	diagY,
	diagZ,
	superDiag,
};

class alignment
{
private:
	vector<int> m_data;
	vector<vector<int>> m_scoreMatrix;
public:
	alignment() { ; };
	alignment(string name, problemType type);
	int matchScore(char x);
	int mismatchScore(char x, char y);
	int DPChange(int money);
	int ManhattanTourist(int n, int m);
	vector<vector<steps>> LCSBackTrack(string v, string w);
	string outputLCS(vector<vector<steps>> backtrack, string v, int i, int j);
	vector<vector<steps>> generalBackTrack(string v, string w, int sigma, int& score);
	vector<string> generalLCS(vector<vector<steps>> backtrack, string v, string w, int i, int j);
	vector<vector<vector<multiSteps>>> multiDimBackTrack(string u, string v, string w, int& score);
	vector<string> multiDimLCS(vector<vector<vector<multiSteps>>> back, string u, string v, string w, int i, int j, int k);
	vector<vector<string>> taxiBackTrack(string v, string w, int sigma);
	vector<string> taxiLCS(vector<vector<string>> backtrack, string v, string w, int i, int j);
	vector<vector<string>> fittingBackTrack(string v, string w, int sigma, int& score);
	vector<string> fittingString(vector<vector<string>> backtrack, string v, string w, int i, int j);
	vector<vector<string>> overlapBackTrack(string v, string w, int sigma, int& score);
	vector<string> overlapAllignment(vector<vector<string>> backtrack, string v, string w, int i, int j);
	vector<vector<vector<string>>> gapBackTrack(string v, string w, int sigma, int epsi, int& score, int& level);
	vector<string> gapAllignment(vector<vector<vector<string>>> backtrack, string v, string w, int i, int j, int k);
	steps middleEdge(string realV, string realW, int top, int bottom, int left, int right, int sigma, int& point);
	void linearSpaceAlignment(string v, string w, int top, int bottom, int left, int right, int sigma, vector<steps>& res);
	vector<string> pathTostrings(string str1, string str2, vector<steps> path, int sigma, int& score);
	int editDistance(string str1, string str2);
	void cutHairString(vector<string>& match, int& score);
};
