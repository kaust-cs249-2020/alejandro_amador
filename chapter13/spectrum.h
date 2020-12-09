#pragma once
#include "mathStuff.h"

class spectrum
{
private: 
	vector<vector<pair<int, int>>> adjList;
	map<string, int> massTable;
	map<int, string> nameTable;
	map<int, int> nodeTranslator;
	vector<int> linearSpectrum;
public:
	spectrum() { loadMassTable(); };
	spectrum(string name);
	void loadMassTable();
	void printBookFormat();
	//Decodes the spectrum 
	string decodingIdealSpectrum();
	//Peptide Sequencing Problem
	string peptideSequencingProblem();
	//Finds all paths from u to v
	vector<vector<int>> allPaths(int u, int v);
	//Helper function all paths
	void DFS(int& u, int& v, vector<int> currentPath, vector<vector<int>>& simplePaths, vector<bool> visited);
	//Misc Functions
	//Returns out neighbors
	vector<int> getOutNeighbors(int u);
	//Transforms path to peptide
	string pathToPeptide(vector<int>& path);
	//Given a peptide it returns its spectrum
	vector<int> peptideToSpectrum(string peptide);
	//Given a peptide it returns its total mass
	int massFragment(string fragment);
	//Transforms peptide to peptideVector
	vector<int> peptideIntoVector(string peptide);
	//Transforms peptideVector into peptide
	string vectorIntoPeptide(vector<int>& peptiVector);
	//Given a path it returns its weight for final problem
	int pathWeight(vector<int>& path);
};