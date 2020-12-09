#pragma once
#include"weightedGraph.h"

enum treeType
{
	rooted,
	unrooted
};

class parsiTree : public weightedGraph
{
private:
	map<int,string> alphaTrans;
	vector<char> alphabet;
	vector<vector<float>> directed;
public:
	parsiTree() { ; };
	//Constructor loads adj list of the tree given the list and number of leaves
	parsiTree(string name, treeType tipo);
	//Constructor takes a forest and transforms back into tree
	parsiTree(vector<parsiTree> bosque);
	//Inits alphabet
	void initAlpha();
	//Prints directed info
	void printDirected();
	//Prints in the maps format
	void printMap();
	//Prints the actual map
	void printTheMap();
	//Prints the tree without weight, only connections
	void printNoWeight();
	//Gets alphabet
	vector<char> getAlphabet() { return alphabet; };
	//Gets char of a node (tree with map from int to char space)
	char getChar(int pos);
	//Push char of a node (tree with map from int to char space)
	void pushMap(int pos, string x);
	//Given a node it gives all of its neighbors
	vector<int> getNeighborsNode(int pos);
	//Gives the children of a node
	vector<int> getChildrenNode(int pos);
	//Gives the father of a node
	int getFatherNode(int pos);

	//Helper functions
	//Gets the chars associated to all leaves
	vector<char> charLeaves();
	//Gets list of ripes on T
	vector<int> getRipes(vector<int>& tags);
	//Returns size of string in the map
	int sizeMapString();
	//Creates subTree Ti with the map to the ith char
	parsiTree subCharTree(int pos);
	//Creates vector of possibleSubtrees
	vector<parsiTree> subCharForest();
	//Adds a root to the tree (pushes a new node at the top)
	void addRoot();
	//Removes a root to the tree (on the last pos)
	void removeRoot();
	//Given and edge e (from a to b) it returns the two nearest neighbors
	vector<parsiTree> nearestNeighborsProblem(int a, int b);
	//Gets internal nodes of the tree
	vector<int> internalNodes();
	//Gets the internal edges of the tree as pairs of nodes
	vector<vector<int>> internalEdges();
	//Disconect pair of nodes (directed edge)
	void disconectDirected(int a, int b);
	//Conects pair of nodes (drected edge, in given order)
	void connectDirected(int a, int b);
	//Fixes parent-child relations
	void fixChildren();
	//Gives children when the father is known
	vector<int> getChildrenKnownFather(int x);
};

//Other functions
//SmallParsimony for unique char tree
int smallParsimony(parsiTree T, parsiTree& Tfilled);
//SmallParsimony for string tree
int smallParsimonyString(parsiTree T, parsiTree& Tfilled);
//SmallParsimony for unrooted tree
int smallParsimonyUnrooted(parsiTree T, parsiTree& Tfilled);
//Heuristic local-search for parsimony score
parsiTree nearestNeighborInterchange(parsiTree T);
