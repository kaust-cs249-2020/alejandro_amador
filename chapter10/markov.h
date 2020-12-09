#pragma once
#include"mathStuff.h"

class markov
{
protected:
	vector<vector<double>> tranMat;
	vector<vector<double>> emmiMat;
	map<string, int>  statePos;
	map<char, int> charPos;
	vector<string> stateList;
	vector<char> charList;
public:
	//Default constructor
	markov() { ; };
	//Load transition matrix from file
	markov(string name);
	//Given a read it computes its profile HMM (takes pre/post treshold)
	markov(vector<char>& pre, vector<char>& post, vector<bool>& deleted, map<char,int>& alphabet);
	//Given a path it computes its prob
	double pathProb(string path);
	//Given an emmiter string and a path, gives the prob of having being emmited by that path
	double emmisionPathProb(string msg, string path);
	//Viterbi algorithm for solving decoding problem 
	//(given msg which path of states was more likely to have produced it)
	string viterbiAlgorithm(string msg);
	//Outcome likelihood problem
	//Given msg what is the prob that HMM emits it
	double outcomeLikelihoodProblem(string msg);
	//Given a pseucocount value sigma it normalices
	void pseudoCountNormalize(double sigma);
	//Soft decoding algorithm
	vector<vector<double>> softDecodingProblem(string msg);
	//Soft decoding but gives you Gamma** for the Baum-Welch algorithm
	void softDecodingGammas(string msg, vector<vector<double>>& gamma, vector<vector<double>>& gammaStar);
	//M-step for baum-welch
	void mStepBaum(vector<vector<double>>& gamma, vector<vector<double>>& gammaStar, string msg);
	//Baum-Welch learning
	void baumWelchLearning(int numIter, string msg);

	//Misc Functions
	//Prints the tran mat
	void printTran();
	//Prints the states with positions
	void printStates();
	//Prints the emmision matrix
	void printEmmision();
	//Prints the alphabet with positions
	void printChars();
	//Prints only connected states
	void printConnected();
	//Get state list
	vector<string> getStateList() { return stateList; };
	//Getter state pos
	map<string, int> getStatePos() { return statePos; };
	//Getter tran mat entry
	double getTranAt(int i, int j) { return tranMat[i][j]; };
	//Getter emmision mat entry
	double getEmitAt(int i, int j) { return emmiMat[i][j]; };
	//Getter emmision mat row
	vector<double> getEmitRow(int i) { return emmiMat[i]; };
	//Getter transmtion mat row
	vector<double> getTranRow(int i) { return tranMat[i]; };
};