#pragma once
#include "suffixArray.h"


class BWT
{
private:
	string sentence;
	vector<pair<char, int>> first;
	vector<pair<char, int>> last;
	map<int, int> lastToFirst;
	vector<map<char, int>> count;
	map<char, int> firstOcurrence;
	suffixArray sufijo;
public:
	//Default constructor
	BWT() { ; };
	//Constructor given text input
	BWT(string text);
	//Matching of a pattern to the string, counts total number of matches
	int BWMatching(string pattern);
	//Given a vector of patterns, gives vector of number of matchings for them
	vector<int> multipleBWMatching(vector<string>& patterns);
	//Better version of BMW matching
	int betterBWMatching(string pattern);
	//Multiple version of better BMW
	vector<int> multipleBetterBWM(vector<string>& patterns);
	//Better BMW but tells the starting positions
	vector<int> betterBWMPositions(string& pattern);
	//Multiple BMW for positions
	vector<int> multipleBWMPositions(vector<string>& patterns);
	//Naive text find with mismatches
	vector<int> naiveMismatch(string& pattern, int d);
	//Multiple naive search with mismatches
	vector<int> multepleNaiveMismatch(vector<string>& patterns, int d);

	//Misc Functions

	//Checks from top to bottom in lastcolumn to find a symbol
	bool checkSymbolHere(int top, int bot, int& topIndex, int& botIndex, char& symbol);
	//Checks from top to bottom in lastColumn but only checks 
	bool checkSymbolHereNoPos(int top, int bot, char& symbol);
};

//Gives the BWT
string transformBW(string& text);
//Given BWT it applies inverse transform
string inverseBW(string& text);
