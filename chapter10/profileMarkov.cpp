#include "profileMarkov.h"

//Given two vectors it sums their entries on the first one
void vectorSum(vector<double>& x, vector<double> y);
//Normalices a vector
void normaliceVector(vector<double>& x);
//Very specific function for the new viterbi alignment algorithm, only checks last three entries
void maxValueVectorNew(double& val, int& pos, vector<double>& x);
//Given a path, computes number of transitions from l to k on it
int numberTransitionsInPath(int i, int j, string& path, map<string, int>& staPos);
//Given a path, computes sum of transitions from i to all states
int numberOutPath(int i, string& path, map<string, int>& staPos);
//Given a path, a state k and a symbol sym, it computes the number of time that symbol
//is emmited when the path is on state k
int numberEmmisionsInState(string k, char sym, string& path, string& msg);
//Given a path and a state k, gives sum over all emmited symbols
int numberTotalEmmisionsInState(string k, string& path);

profileMarkov::profileMarkov(string name)
{
	//First reads everything from the file (alphabet, alignments and theta)
	readFileProfile(name);

	//Now assembles HMM from data

	//First creates alignment*
	deletedCols = vector<bool>(alignment[0].size(), false);
	alignmentStar = alignment;

	int numberAligns = alignment.size();
	for (int i = alignment[0].size() - 1; i >= 0; i--)
	{
		vector<char> column = getColumn(alignment, i);
		int counter = count(column.begin(), column.end(), '-');
		double score = double(counter) / numberAligns;
		if (score >= theta)
		{
			eraseColumn(alignmentStar, i);
			deletedCols[i] = true;
		}
	}

	//Creates profile HMM for every read
	//This seems very dumb, maybe just do it on the fly
	vector<markov> readsHMM;
	for (int i = 0; i < alignment.size(); i++)
	{
		vector<char> preTheta = alignment[i];
		vector<char> posTheta = alignmentStar[i];
		markov temp(preTheta, posTheta, deletedCols, charPos);
		readsHMM.push_back(temp);
	}

	//Computes tran matrix from all the hmms
	stateList = readsHMM[0].getStateList();
	statePos = readsHMM[0].getStatePos();
	int numStates = stateList.size();
	int numChars = charList.size();
	tranMat = vector<vector<double>>(numStates, vector<double>(numStates, 0));
	emmiMat = vector<vector<double>>(numStates, vector<double>(numChars, 0));

	for (int i = 0; i < numStates; i++)
	{
		vector<double> stateTran(numStates, 0);
		for (auto& m : readsHMM)
		{
			vectorSum(stateTran, m.getTranRow(i));
		}

		normaliceVector(stateTran);
		tranMat[i] = stateTran;
	}


	//Now the emmiting matrix
	for (int i = 1; i < numStates - 1; i++)
	{
		vector<double> stateEmit(numChars, 0);
		for (auto& m : readsHMM)
		{
			vectorSum(stateEmit, m.getEmitRow(i));
		}

		normaliceVector(stateEmit);
		emmiMat[i] = stateEmit;
	}

}


//By now we dont take care of beginning with a deletion state just for mind sanity
vector<string> profileMarkov::sequenceAlignmentViterbi(string msg)
{
	int numStates = tranMat.size();

	//First transition
	char emmited = msg.at(0);
	int emmitedPos = charPos[emmited];
	vector<double> lastS(numStates, 0);
	vector<vector<int>> backTrack(numStates, vector<int>());
	for (int i = 0; i < numStates; i++)
	{
		backTrack[i] = { i };
		double localS;
		char stateType = stateList[i].at(0);
		if (stateType == 'D')
		{
			localS = tranMat[0][i];
		}
		else
		{
			localS = tranMat[0][i] * emmiMat[i][emmitedPos];
		}

		lastS[i] = localS;
	}

	//Now fills manually D1. It either comes from source or from I0
	if (lastS[1] * tranMat[1][3] > lastS[3])
	{
		lastS[3] = lastS[1] * tranMat[1][3];
		backTrack[3] = { 1,3 };
	}

	//Next, in the first column,  Di can come from D_{i-1} or M_{i-1}
	for (int i = 6; i < numStates; i++)
	{
		if (stateList[i].at(0) == 'D')
		{
			double lastM = lastS[i - 4] * tranMat[i - 4][i];
			double lastD = lastS[i - 3] * tranMat[i - 3][i];
			if (lastM > lastD)
			{
				vector<int> currentPath = backTrack[i - 4];
				currentPath.push_back(i);
				backTrack[i] = currentPath;
				lastS[i] = lastM;
			}
			else
			{
				vector<int> currentPath = backTrack[i - 3];
				currentPath.push_back(i);
				backTrack[i] = currentPath;
				lastS[i] = lastD;
			}
		}
	}

	//Other transitions
	//First fills S_{ik} for states other than deletions, luckily deletions can only travel to other
	//deletions within the same column, so we can fill all other states first
	for (int i = 1; i < msg.size(); i++)
	{
		emmited = msg.at(i);
		emmitedPos = charPos[emmited];
		vector<vector<int>> tempBacktrack(backTrack);
		vector<double> tempS(lastS);
		for (int j = 1; j < numStates - 1; j++)
		{
			if (stateList[j].at(0) != 'D')
			{
				double maxS = 0;
				int maxState;
				for (int k = 1; k < numStates - 1; k++)
				{
					double localS = lastS[k] * tranMat[k][j] * emmiMat[j][emmitedPos];
					if (localS > maxS)
					{
						maxS = localS;
						maxState = k;
					}
				}

				tempS[j] = maxS;
				vector<int> currentPath = backTrack[maxState];
				currentPath.push_back(j);
				tempBacktrack[j] = currentPath;
			}
		}

		backTrack = tempBacktrack;
		lastS = tempS;

		//After filling all other states, we can start filling deletions on the same column
		//First, D1 can only come from I0
		lastS[3] = lastS[1] * tranMat[1][3];
		vector<int> localPath = backTrack[1];
		localPath.push_back(3);
		backTrack[3] = localPath;

		//The next Di can come from previous I,D,M on the column
		for (int j = 5; j < numStates - 1; j++)
		{
			if (stateList[j].at(0) == 'D')
			{
				double maxS = 0;
				int maxState;
				for (int k = 1; k <= j; k++)
				{
					double localS = lastS[k] * tranMat[k][j];
					if (localS > maxS)
					{
						maxS = localS;
						maxState = k;
					}
				}

				lastS[j] = maxS;
				vector<int> currentPath = backTrack[maxState];
				currentPath.push_back(j);
				backTrack[j] = currentPath;
			}
		}
	}

	//Gets maximum path to the sink (only three states can go to the sink)
	vector<string> path;
	int pathPos;
	double pathMax;
	maxValueVectorNew(pathMax, pathPos, lastS);
	vector<int> pathStates = backTrack[pathPos];
	for (auto& x : pathStates)
	{
		path.push_back(stateList[x]);
	}

	return path;
}

void profileMarkov::parameterEstimationProblem(string name)
{
	//FILE READING
	int numChars = 0;
	int numStates = 0;
	string msg;
	string path;
	ifstream file;
	file.open(name + ".txt");
	string temp;

	getline(file, msg);
	getline(file, temp);
	getline(file, temp);
	strstream s;
	s << temp;
	while (true)
	{
		string newChar;
		s >> newChar;
		if (newChar.empty())
		{
			break;
		}
		else
		{
			char realNew = newChar.at(0);
			charList.push_back(realNew);
			charPos[realNew] = numChars;
			numChars++;
		}
	}

	getline(file, temp);
	getline(file, path);
	getline(file, temp);
	getline(file, temp);
	s.clear();
	s << temp;
	while (true)
	{
		string newState;
		s >> newState;
		if (newState.empty())
		{
			break;
		}
		else
		{
			stateList.push_back(newState);
			statePos[newState] = numStates;
			numStates++;
		}
	}

	file.close();
	tranMat = vector<vector<double>>(numStates, vector<double>(numStates, 0));
	emmiMat = vector<vector<double>>(numChars, vector<double>(numChars, 0));

	//ACTUAL ALGORITHM
	for (int i = 0; i < numStates; i++)
	{
		double normalization = 1.f / numberOutPath(i, path, statePos);
		if (normalization > 0)
		{
			for (int j = 0; j < numStates; j++)
			{
				tranMat[i][j] = normalization * numberTransitionsInPath(i, j, path, statePos);
			}
		}
		else
		{
			normalization = 1.f / numStates;
			for (int j = 0; j < numStates; j++)
			{
				tranMat[i][j] = normalization;
			}
		}
	}

	for (int i = 0; i < numStates; i++)
	{
		string state = stateList[i];
		double normalization = 1.f / numberTotalEmmisionsInState(state, path);
		if (normalization > 0)
		{
			for (int j = 0; j < numChars; j++)
			{
				emmiMat[i][j] = normalization * numberEmmisionsInState(state, charList[j], path, msg);
			}
		}
		else
		{
			normalization = 1.f / numChars;
			for (int j = 0; j < numChars; j++)
			{
				emmiMat[i][j] = normalization;
			}
		}
	}
}

void profileMarkov::parameterEstimationProblemForViterbi(string msg, string path)
{
	int numStates = stateList.size();
	int numChars = charList.size();

	//ACTUAL ALGORITHM
	for (int i = 0; i < numStates; i++)
	{
		double normalization = 1.f / numberOutPath(i, path, statePos);
		if (normalization > 0)
		{
			for (int j = 0; j < numStates; j++)
			{
				tranMat[i][j] = normalization * numberTransitionsInPath(i, j, path, statePos);
			}
		}
		else
		{
			normalization = 1.f / numStates;
			for (int j = 0; j < numStates; j++)
			{
				tranMat[i][j] = normalization;
			}
		}
	}

	for (int i = 0; i < numStates; i++)
	{
		string state = stateList[i];
		double normalization = 1.f / numberTotalEmmisionsInState(state, path);
		if (normalization > 0)
		{
			for (int j = 0; j < numChars; j++)
			{
				emmiMat[i][j] = normalization * numberEmmisionsInState(state, charList[j], path, msg);
			}
		}
		else
		{
			normalization = 1.f / numChars;
			for (int j = 0; j < numChars; j++)
			{
				emmiMat[i][j] = normalization;
			}
		}
	}
}

void profileMarkov::viterbiLearning(string name)
{
	//FILE READING
	int numChars = 0;
	int numStates = 0;
	string msg;
	ifstream file;
	file.open(name + ".txt");
	string temp;

	getline(file, temp);
	int numIter = stoi(temp);
	getline(file, temp);
	getline(file, msg);
	getline(file, temp);
	getline(file, temp);
	strstream s;
	s << temp;
	while (true)
	{
		string temp;
		s >> temp;
		if (temp.empty())
		{
			break;
		}
		else
		{
			char newChar = temp.at(0);
			charList.push_back(newChar);
			charPos[newChar] = numChars;
			numChars++;
		}
	}

	getline(file, temp);
	getline(file, temp);
	s.clear();
	s << temp;
	while (true)
	{
		string temp;
		s >> temp;
		if (temp.empty())
		{
			break;
		}
		else
		{
			stateList.push_back(temp);
			statePos[temp] = numStates;
			numStates++;
		}
	}

	tranMat = vector<vector<double>>(numStates, vector<double>(numStates, 0));
	emmiMat = vector<vector<double>>(numStates, vector<double>(numChars, 0));

	getline(file, temp);
	getline(file, temp);
	for (int i = 0; i < numStates; i++)
	{
		s.clear();
		getline(file, temp);
		string tempState;
		s << temp;
		s >> tempState;
		for (int j = 0; j < numStates; j++)
		{
			s >> tranMat[i][j];
		}
	}

	getline(file, temp);
	getline(file, temp);
	for (int i = 0; i < numStates; i++)
	{
		s.clear();
		getline(file, temp);
		string tempState;
		s << temp;
		s >> tempState;
		for (int j = 0; j < numChars; j++)
		{
			s >> emmiMat[i][j];
		}
	}

	//Actual algorithm
	for (int i = 0; i < numIter; i++)
	{
		string path = viterbiAlgorithm(msg);
		parameterEstimationProblemForViterbi(msg, path);
	}
}

void profileMarkov::readFileProfile(string name)
{
	//First reads everything from the file
	int numChars = 0;
	ifstream file;
	file.open(name + ".txt");
	string temp;

	//First gets treshold
	getline(file, temp);
	strstream s;
	s << temp;
	s >> theta;

	//Then gets alphabet
	getline(file, temp);
	getline(file, temp);
	s.clear();
	s << temp;
	while (true)
	{
		string local;
		s >> local;
		if (local.empty())
		{
			break;
		}
		else
		{
			charPos[local.at(0)] = numChars;
			charList.push_back(local.at(0));
			numChars++;
		}
	}

	getline(file, temp);
	//Finally it gets the alignment
	while (getline(file, temp))
	{
		vector<char> reading;
		while (!temp.empty())
		{
			reading.push_back(temp.at(0));
			temp.erase(temp.begin());
		}

		alignment.push_back(reading);
	}

	file.close();
}

void vectorSum(vector<double>& x, vector<double> y)
{
	for (int i = 0; i < x.size(); i++)
	{
		x[i] += y[i];
	}
}

void normaliceVector(vector<double>& x)
{
	double alpha = 0;
	for (auto& r : x)
	{
		alpha += r;
	}

	if (alpha != 0)
	{
		for (int i = 0; i < x.size(); i++)
		{
			x[i] /= alpha;
		}
	}
}

void maxValueVectorNew(double & val, int & pos, vector<double>& x)
{
	int totalSize = x.size();
	val = x[totalSize - 4];
	pos = totalSize - 4;
	for (int i = totalSize - 3; i < totalSize - 1; i++)
	{
		if (x[i] > val)
		{
			val = x[i];
			pos = i;
		}
	}
}

int numberTransitionsInPath(int i, int j, string & path, map<string, int>& staPos)
{
	int counter = 0;
	for (int k = 0; k < path.size() - 1; k++)
	{
		int stateFrom = staPos[path.substr(k, 1)];
		int stateTo = staPos[path.substr(k + 1, 1)];
		if (stateFrom == i && stateTo == j)
		{
			counter++;
		}
	}

	return counter;
}

int numberOutPath(int i, string & path, map<string, int>& staPos)
{
	int counter = 0;
	for (int k = 0; k < path.size() - 1; k++)
	{
		int stateFrom = staPos[path.substr(k, 1)];
		if (stateFrom == i)
		{
			counter++;
		}
	}

	if (counter == 0)
	{
		return -1;
	}
	else
	{
		return counter;
	}
}

int numberEmmisionsInState(string k, char sym, string & path, string & msg)
{
	int counter = 0;
	for (int i = 0; i < path.size(); i++)
	{
		string tempState = path.substr(i, 1);
		char tempEmmi = msg.at(i);
		if (tempState == k && tempEmmi == sym)
		{
			counter++;
		}
	}
	return counter;
}

int numberTotalEmmisionsInState(string k, string & path)
{
	int counter = 0;
	for (int i = 0; i < path.size(); i++)
	{
		string tempState = path.substr(i, 1);
		if (tempState == k)
		{
			counter++;
		}
	}

	if (counter == 0)
	{
		return -1;
	}
	else
	{
		return counter;
	}
}
