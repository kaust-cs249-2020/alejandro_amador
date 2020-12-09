#pragma once
#include"agrozUtil.h"
#include<queue>

class weightedGraph
{
protected:
	vector<vector<float>> adjList;
public:
	weightedGraph() { ; };
	//Constructor loads weighted graph as book's format
	weightedGraph(string name);
	//Constructor given the adjacency list
	weightedGraph(vector<vector<float>> mAdj) { adjList = mAdj; };
	//Constructor for an isolated tree of n nodes
	weightedGraph(int n);
	//Prints adj matrix
	void print();
	//Computes path dstance between every pair of nodes (tree)
	vector<vector<float>> distMat();
	//Prints a weigthed graph as books's format (its adjList)
	void printBookFormat();

	//Mics functions
	//Appends back a leaf to a potentially new node
	void limbAppender(int v, int n, int A, int B, int DA, int DB, int limb);
	//Given a node and a ditance, it looks for a node at that distance from the source
	//If the node is not found then creates a new one pop+1
	int findNodeAtDistance(int i, int k, int d, int& pop, int& limA, int& limB, int& limDA, int& limDB);
	//Given two nodes it returns a path between them (tree)
	vector<int> pathNodesTree(int a, int b);
	//Gives all adjacent nodes of a given node
	vector<int> adjacentsToNode(int node);
	//Given a path it returns its total weight
	int pathWeight(vector<int> path);
	//Given a node it returns the list of its in-neighbors
	vector<int> inNeighborsNode(int node);
	//Given a node it returns the list of its out-neighbors
	vector<int> outNeighborsNode(int node);
	//It returns true if a node is a leaf (degree 1)
	bool isNodeLeaf(int node);
	//Given two nodes it returns its path distance
	int pathDistanceNodes(int a, int b);
	//Returns all leaves in the graph
	vector<int> leavesFinder();
	//Pushes a new node to the graph (isolated node at the end of the adj list, returns its index)
	int pushNode();
	//Sets the distance between two nodes (it connects/disconnects nodes)
	void connectNodes(int u, int v, float dist);
	//Connects a node to a set of nodes (cluster) with "distance" -2 for reference
	void connectNodeToCluster(int u, vector<int>& cluster, map<int, int>& mapi);
	//Helper function, given a vector of ages it sets the correspondent adjacencies
	void setEdgesFromAge(vector<float>& ages);
	//Gets number of nodes
	int numberNodes() { return adjList.size(); };
};

//Functions outside class 

//Loads distance matrix of size nXn
vector<vector<float>> loadDistMatrix(string name, int n);
//Solves the limb length problem
int limbLength(vector<vector<float>>& distMatrix, int j);
//Algorithm for phylogeny reconstruction of tree given distance matrix
weightedGraph additivePhylogeny(vector<vector<float>>& D, int& pop);
//UPGMA algorithm for creating tree given a distance matrix
weightedGraph UPGMA(vector<vector<float>>& D);
//Modified UPGMA algorithm for clustering
vector<vector<int>> UPGMAcluster(vector<vector<float>>& D);
//Joining algorithm for creating tree given a distance matrix
weightedGraph neighborJoining(vector<vector<float>>& D, int pop = -1, map<int,int> dictio = {});

//Transformation UPGMA into readable book
vector<vector<int>> transformUPGMA(vector<vector<int>>& backs, int n);