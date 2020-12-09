#pragma once
#include "suffixTree.h"

enum nodeColor
{
	gray,
	blue,
	red,
	purple
};

class coloredSuffix : public suffixTree
{
private:
	map<int, nodeColor> nodeCanvas;
	string sentence1;
	string sentence2;
public:

	//Default constructor
	coloredSuffix() { ; };
	//Constructor just calls suffixTree constructor on both words
	coloredSuffix(string text1, string text2);
	//Constructor given book's example
	coloredSuffix(string name);
	//Print the colors
	void printColors();
	//Coloring the tree
	void giveMeColor();
	//Longest shared substring problem
	string longestSharedSubstringProblem();
	//Shortest non-shared substring problem
	string shortestNonSharedSubstringProblem();
	//DFS for strings from root but also saves color of final node
	void stringsColorDFS(int x, map<int, bool>& visited, vector<string>& list, vector<nodeColor>& colList, string& soFar);
	//DFS for strings from root with color but ends in leafes (only blue ones)
	void stringsBlueDFS(int x, map<int, bool>& visited, vector<string>& list, string& soFar);

	//Misc Functions
	//Given a sentence it finds whether it is blue or red
	nodeColor findMyColor(string x);
	//Finds ripe nodes (gray with no gray children)
	vector<int> findRipeNodes();
};