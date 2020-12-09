#pragma once
#include "markov.h"

class profileMarkov : public markov
{
private:
	double theta;
	vector<vector<char>> alignment;
	vector<vector<char>> alignmentStar;
	vector<bool> deletedCols;
public:
	//Default constructor
	profileMarkov() { ; };
	//Constructor given file
	profileMarkov(string name);
	//Sequence alignment with profile (viterbi)
	vector<string> sequenceAlignmentViterbi(string msg);
	//Estimation of parameters (transition/emition matrices) given msg and path
	void parameterEstimationProblem(string name);
	//Estimation of parameters but given only the path and the msg (states and chars already saved)
	void parameterEstimationProblemForViterbi(string msg, string path);
	//Viterbi learning algorithm
	void viterbiLearning(string name);

	//Misc functions
	void readFileProfile(string name);
};