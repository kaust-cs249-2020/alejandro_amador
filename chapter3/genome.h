#pragma once

#include "agrozUtil.h"

enum classType
{
	patterns,
	graph,
	KDgraph,
	temp
};

class genome
{
private:
	vector<string> m_dna;
	vector<vector<string>> m_graph;
public:
	genome() { ; };
	genome(int k);
	genome(string filename, classType type);
	string universalString();
	vector<string> composition(int pos, int k);
	vector<vector<string>> overlap();
	vector<vector<string>> pathGraph(int k);
	vector<vector<string>> deBrujin(int k);
	vector<vector<string>> pathFromKmers();
	vector<vector<string>> deBrujinFromKmers();
	vector<string> contigsFromKmers();
	vector<string> eulerianCycle();
	vector<string> generateNodes();
	vector<string> eulerianPath();
	vector<vector<string>> eulerianPathKD();
	string stringReconstruction();
	string stringFromKDmers(int k, int d);
	vector<vector<string>> generateKDmers(int k, int d);
};