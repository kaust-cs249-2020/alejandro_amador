#pragma once
#include"agrozUtil.h"


class digraph
{
private:
	vector<vector<int>> m_adjList;
public:
	digraph() { ; };
	digraph(string name, int n);
	vector<int> topologicalOrdering();
	void print();
};
