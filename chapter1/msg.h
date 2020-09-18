#pragma once

#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<algorithm>
#include<map>

using namespace std;


class msg
{
private:
	string m_inputString;
public:
	msg(string filename);
	void print();
	int findFrequency(string pattern);
	int findMismatchFrequency(string pattern, int distance);
	vector<string> frequentWords(int k);
	vector<string> frequentWordsWithMismatch(int k, int d);
	vector<string> betterFrequentWords(int k);
	vector<string> frequentWordsMRC(int k, int d);
	vector<int> patternMatching(string pattern);
	vector<int> approxPatternMatching(string pattern, int distance);
	vector<string> findClumps(int k, int L, int t);
	vector<int> frequencyPlot(string pattern, int d);
	vector<int> skewDiagram();
	vector<int> skewMinPos();
	~msg();
};

