#pragma once
#include "agrozUtil.h"

class genome
{
private:
public:
	//Empty constructor
	genome() { ; };
	//Given a permutation it computes reverse distance
	int greedySorting(vector<int>& p);
	//Given a permutation it computes number of breaking points
	int breakingPoints(vector<int> p);
	//Transforms chromosome to head/tail nodes for graph
	vector<int> chromosomeToCycle(vector<int> chrome);
	//Transform complete cycle to chromosome
	vector<int> cycleToChromosome(vector<int> nodes);
	//Transform red cycle to chromosome
	vector<int> redCycleToChromosome(vector<int> nodesRed);
	//Given a genome it returns its colored edges
	vector<vector<int>> coloredEdges(vector<vector<int>> p);
	//Given a genome it returns its black edges
	vector<vector<int>> blackEdges(vector<vector<int>> p);
	//Transforms graph black/color back to genome
	vector<vector<int>> graphtoGenome(vector<vector<int>> graph, vector<vector<int>> reds);
	//Transforms red graph back to genome
	vector<vector<int>> redGraphToGenome(vector<vector<int>> graph);
	//Given two genomes it returns their 2-Break distance
	int twoBreakDistance(vector<vector<int>> P, vector<vector<int>> Q);
	//Given two genomes it prints the 2-break from one to another
	void shortestRearrangementScenario(vector<vector<int>> P, vector<vector<int>> Q);
	//given two strings finds all shared k-mers
	vector<vector<int>> sharedKmers(string g1, string g2, int k);
	//given two strings finds number of shared k-mers
	int numberSharedKmers(string g1, string g2, int k);


	//Misc functions.
	void kSortingReversal(vector<int>& p, int k);
	//This for reading permutations/cycles
	vector<int> readPermutation(string name);
	//This for reading genomes as colection of permutation
	vector<vector<int>> readGenome(string name);
	//This for reading a genomegraph as collection of nodes
	vector<vector<int>> readGraph(string name);
	//This for reading multiple genomes as collections of permutations
	vector<vector<vector<int>>> readMultipleGenomes(string name);
	//Performs a 2break on a graph
	vector<vector<int>> twoBreakGraph(vector<vector<int>> g, int a, int b, int c, int d);
	//Given a genome it applies a 2-break on it
	vector<vector<int>> twoBreakGenome(vector<vector<int>> p, int a, int b, int c, int d);
	//Loads a string from a file
	string loadString(string name);
};
