#pragma once
#include "agrozUtil.h"

enum dataType
{
	dna,
	spectrum,
	totalAlphabet
};

class amino
{
private:
	vector<string> m_rna;
	vector<vector<string>> m_trans;
	vector<int> m_spectrum;
	map<string, int> m_masses;
public:
	amino() { ; };
	amino(string filename, dataType type);
	void loadTranslation();
	void loadMasses();
	void loadAllMass();
	string codonToAmino(string codon);
	string rnaToAmino(string rna);
	void multiplicityCount();
	vector<string> substringsEncoding(string aminos);
	vector<int> linearSpectrum(vector<string> peptide);
	vector<int> cyclicSpectrum(vector<string> peptide);
	bool isConsistent(vector<int> subset);
	long long int numberOfPeptides(int m);
	int peptideMass(vector<string> peptide);
	void expand(vector<vector<string>>& peptides);
	vector<vector<string>> cycloPeptideSequencing();
	vector<string> convolutionCycloPeptideSequencing(int M, int N);
	vector<string> leaderBoardCycloPeptideSequencing(int N);
	vector<vector<int>> peptidesToMass(vector<vector<string>> peptides);
	vector<int> peptideToMass(vector<string> peptide);
	int score(vector<string> peptide);
	int linearScore(vector<string> peptide);
	void trim(vector<vector<string>>& leaderboard, int N);
	vector<int> convolution();
	void extendedAlphabet(int M);
	void pushScore(vector<string> p, vector<vector<string>>& kings);

};

