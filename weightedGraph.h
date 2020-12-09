#pragma once
#include "agrozUtil.h"

class wegraph
{
private:
	vector<vector<vector<int>>> m_adjList;
public:
	wegraph() { ; };
	wegraph(string name, int n);
	vector<int> longestPath(int start, int end);
	vector<vector<int>> antecesores(int node);
	void print();
};