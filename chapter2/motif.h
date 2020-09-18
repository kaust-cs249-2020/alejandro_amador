#pragma once
#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<algorithm>
#include<map>
#include<numeric>
#include<ctime>
#include <cstdlib>
#include "agrozUtil.h"

enum algorithm
{
	MEDIAN_STRING,
	RANDOMIZED_MOTIF_SEARCH,
	GIBBS_SAMPLING
};

using namespace std;

class motif
{
private:
	vector<string> m_dna;
public:
	motif() { ; };
	motif(string filename);
	string consensus(int k, algorithm option, int runs = 1);
	vector<string> greedyMotifSearch(int k);
	vector<string> initialBestMotifs(int k);
	vector<string> randomizedMotifSearch(int k, int runs);
	vector<string> gibbsSampler(int k, int pop, int runs);
	vector<string> initialRandomMotifs(int k);
	vector<string> motifsMostProbable(vector<vector<float>> profile);
	string profileMostProbable(vector<vector<float>> profile, int position);
	string profileRandomly(vector<vector<float>> profile, int position);
	string orderLine();
	string medianString(int k);
	vector<string> motifEnumeration(int k, int d);
	void print();
	bool appearsMismatch(string pattern, int distance);
	int numberAppearsMismatch(string pattern, int distance);
	int DistanceBetweenPatternAndStrings(string pattern);
};

