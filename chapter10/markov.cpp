#include "markov.h"

//Normalices a vector from position a to b
void normaliceVectorRange(vector<double>& x, int a, int b, double theta);
//Returns max value in a double vector and also the position
void maxValueVector(double& val, int& pos, vector<double>& x);
//Gives T_{l,k} for baumm learning
double computeT(int l, int k, vector<vector<double>>& gamma);
//Gives E_{k}(b) for baumm learning
double computeE(int k, char b, string msg, vector<vector<double>>& gamma);
//Gives normalization for E_{k}
double normalE(vector<char>& alphabet, int k, string msg, vector<vector<double>>& gamma);
//Gives normalization for T_{l}
double normalT(int l, vector<vector<double>>& gammaStar);

markov::markov(string name)
{
	//READS BOTH TRANSITION AND EMMISION MATRIX
	int numStates = 0;
	int numChars = 0;
	ifstream file;
	file.open(name + ".txt");
	string temp;

	//First gets alphabet
	getline(file, temp);
	getline(file, temp);
	strstream s;
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

	//Now gets states
	s.clear();
	getline(file, temp);
	getline(file, temp);
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
			statePos[local] = numStates;
			stateList.push_back(local);
			numStates++;
		}
	}

	//Finally it gets the matrices
	getline(file, temp);
	getline(file, temp);
	tranMat = vector<vector<double>>(numStates, vector<double>(numStates, 0));
	emmiMat = vector<vector<double>>(numStates, vector<double>(numChars, 0));

	for (int i = 0; i < numStates; i++)
	{
		s.clear();
		getline(file, temp);
		s << temp;
		string state;
		s >> state;
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
		s << temp;
		string state;
		s >> state;
		for (int j = 0; j < numChars; j++)
		{
			s >> emmiMat[i][j];
		}
	}


	file.close();


	////READS ONLY TRANSITION MATRIX
	//int size = 0;
	//ifstream file;
	//file.open(name + ".txt");
	//string temp;
	//getline(file, temp);
	//getline(file, temp);

	////Gets states
	//strstream s;
	//s << temp;
	//while (true)
	//{
	//	string local;
	//	s >> local;
	//	if (local.empty())
	//	{
	//		break;
	//	}
	//	else
	//	{
	//		statePos[local] = size;
	//		size++;
	//	}
	//}

	//tranMat = vector<vector<double>>(size, vector<double>(size, 0));
	////Gets transition matrix, first two trash lines
	//getline(file, temp);
	//getline(file, temp);
	//for (int i = 0; i < size; i++)
	//{
	//	s.clear();
	//	getline(file, temp);
	//	s << temp;
	//	string state;
	//	s >> state;
	//	for (int j = 0; j < size; j++)
	//	{
	//		s >> tranMat[i][j];
	//	}
	//}

	////READS ONLY EMISSION MATRIX
	//int numStates = 0;
	//int numChars = 0;
	//ifstream file;
	//file.open(name + ".txt");
	//string temp;

	////First gets alphabet
	//getline(file, temp);
	//getline(file, temp);
	//strstream s;
	//s << temp;
	//while (true)
	//{
	//	string local;
	//	s >> local;
	//	if (local.empty())
	//	{
	//		break;
	//	}
	//	else
	//	{
	//		charPos[local.at(0)] = numChars;
	//		numChars++;
	//	}
	//}

	////Now gets states
	//s.clear();
	//getline(file, temp);
	//getline(file, temp);
	//getline(file, temp);
	//getline(file, temp);
	//s << temp;
	//while (true)
	//{
	//	string local;
	//	s >> local;
	//	if (local.empty())
	//	{
	//		break;
	//	}
	//	else
	//	{
	//		statePos[local] = numStates;
	//		numStates++;
	//	}
	//}

	////Finally gets emmision matrix
	//getline(file, temp);
	//getline(file, temp);
	//emmiMat = vector<vector<double>>(numStates, vector<double>(numChars, 0));
	//for (int i = 0; i < numStates; i++)
	//{
	//	s.clear();
	//	getline(file, temp);
	//	s << temp;
	//	string state;
	//	s >> state;
	//	for (int j = 0; j < numChars; j++)
	//	{
	//		s >> emmiMat[i][j];
	//	}
	//}
}

markov::markov(vector<char>& pre, vector<char>& post, vector<bool>& deleted, map<char,int> &alphabet)
{
	//Total nodes is internal plus source, sink and Io
	int numChars = alphabet.size();
	int numberInternal = post.size();
	int size = 3 * numberInternal + 3;
	tranMat = vector<vector<double>>(size, vector<double>(size, 0));
	emmiMat = vector<vector<double>>(size, vector<double>(numChars, 0));

	//Gets states
	statePos["S"] = 0;
	statePos["I0"] = 1;
	stateList.push_back("S");
	stateList.push_back("I0");
	for (int i = 0; i < numberInternal; i++)
	{
		string inser = "I" + to_string(i + 1);
		string delet = "D" + to_string(i + 1);
		string match = "M" + to_string(i + 1);
		statePos[match] = 1 + 3 * i + 1;
		statePos[delet] = 1 + 3 * i + 2;
		statePos[inser] = 1 + 3 * i + 3;
		stateList.push_back(match);
		stateList.push_back(delet);
		stateList.push_back(inser);
	}

	statePos["E"] = size - 1;
	stateList.push_back("E");

	//Fills transition matrix (1 if visited, 0 otherwise)
	int lastPos = 0;
	int prePosDiff = 0;
	for (int i = 0; i < pre.size(); i++)
	{
		string goesTo;
		int emmits; 
		if (!deleted[i])
		{
			if (pre.at(i) == '-')
			{
				goesTo = "D" + to_string(i + 1 - prePosDiff);
			}
			else
			{
				emmits = alphabet[pre.at(i)];
				goesTo = "M" + to_string(i + 1 - prePosDiff);
			}

			int goesPos = statePos[goesTo];
			tranMat[lastPos][goesPos] += 1;
			if (pre.at(i) != '-')
			{
				emmits = alphabet[pre.at(i)];
				emmiMat[goesPos][emmits] += 1;
			}
			lastPos = goesPos;
		}
		else
		{
			if (pre.at(i) != '-')
			{
				goesTo = "I" + to_string(i - prePosDiff);
				emmits = alphabet[pre.at(i)];
				int goesPos = statePos[goesTo];
				tranMat[lastPos][goesPos] += 1;
				emmiMat[goesPos][emmits] += 1;
				lastPos = goesPos;
			}

			prePosDiff++;
		}
	}

	tranMat[lastPos][size-1] = 1;


}

double markov::pathProb(string path)
{
	double proba = 1;
	string lastState = path.substr(0,1);
	int lastPos = statePos[lastState];
	string nextState;
	int nextPos;
	for (int i = 1; i < path.size(); i++)
	{
		nextState = path.substr(i, 1);
		int nextPos = statePos[nextState];
		proba *= tranMat[lastPos][nextPos];

		lastPos = nextPos;
	}

	proba /= statePos.size();
	return proba;
}

double markov::emmisionPathProb(string msg, string path)
{
	double proba = 1;
	for (int i = 0; i < msg.size(); i++)
	{
		int state = statePos[path.substr(i, 1)];
		int x = charPos[msg.at(i)];
		proba *= emmiMat[state][x];
	}

	return proba;
}

string markov::viterbiAlgorithm(string msg)
{
	int numStates = tranMat.size();

	//First transition
	char emmited = msg.at(0);
	int emmitedPos = charPos[emmited];
	double firstProb = double(1.0) / numStates;
	vector<double> lastS(numStates, 0);
	vector<vector<int>> backTrack(numStates, vector<int>());
	for (int i = 0; i < numStates; i++)
	{
		backTrack[i].push_back(i);
		double localS = firstProb * emmiMat[i][emmitedPos];
		lastS[i] = localS;
	}

	//Other transitions
	for (int i = 1; i < msg.size(); i++)
	{
		emmited = msg.at(i);
		emmitedPos = charPos[emmited];
		vector<vector<int>> tempBacktrack(backTrack);
		vector<double> tempS(lastS);
		for (int j = 0; j < numStates; j++)
		{
			double maxS = 0;
			int maxState;
			for (int k = 0; k < numStates; k++)
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

		backTrack = tempBacktrack;
		lastS = tempS;
	}

	//Gets maximum path
	string path;
	int pathPos;
	double pathMax;
	maxValueVector(pathMax, pathPos, lastS);
	vector<int> pathStates = backTrack[pathPos];
	for (auto& x : pathStates)
	{
		path += stateList[x];
	}

	return path;
}

double markov::outcomeLikelihoodProblem(string msg)
{
	int numStates = tranMat.size();

	//First transition
	char emmited = msg.at(0);
	int emmitedPos = charPos[emmited];
	double firstProb = double(1.0) / numStates;
	vector<double> lastS(numStates, 0);
	for (int i = 0; i < numStates; i++)
	{
		double localS = firstProb * emmiMat[i][emmitedPos];
		lastS[i] = localS;
	}

	//Other transitions
	for (int i = 1; i < msg.size(); i++)
	{
		emmited = msg.at(i);
		emmitedPos = charPos[emmited];
		vector<double> tempS(lastS);
		for (int j = 0; j < numStates; j++)
		{
			double maxS = 0;
			int maxState;
			for (int k = 0; k < numStates; k++)
			{
				double localS = lastS[k] * tranMat[k][j] * emmiMat[j][emmitedPos];
				maxS += localS;
			}

			tempS[j] = maxS;
		}

		lastS = tempS;
	}

	//Gets maximum path
	double pathMax = 0;
	for (auto& x : lastS)
	{
		pathMax += x;
	}

	return pathMax;
}

void markov::pseudoCountNormalize(double sigma)
{
	int numberStates = stateList.size();

	//First normalices transition matrix

	//Outermost entries
	normaliceVectorRange(tranMat[0], 1, 3, sigma);
	normaliceVectorRange(tranMat[1], 1, 3, sigma);
	normaliceVectorRange(tranMat[numberStates - 4], numberStates - 2, numberStates - 1, sigma);
	normaliceVectorRange(tranMat[numberStates - 3], numberStates - 2, numberStates - 1, sigma);
	normaliceVectorRange(tranMat[numberStates - 2], numberStates - 2, numberStates - 1, sigma);


	//Inner states
	int innerHolder = 1;
	int innerPlus = 0;
	for (int i = 2; i < numberStates - 4; i++)
	{
		int a = 3 * innerHolder + 1;
		int b = a + 2;
		normaliceVectorRange(tranMat[i], a, b, sigma);
		innerPlus++;
		if (innerPlus == 3)
		{
			innerHolder++;
			innerPlus = 0;
		}

	}

	//Now normalices emmision matrix

	int numChars = charList.size();

	for (int i = 1; i < numberStates - 1; i++)
	{
		string tempState = stateList[i];
		if (tempState.at(0) != 'D')
		{
			normaliceVectorRange(emmiMat[i], 0, numChars - 1, sigma);
		}
	}

}

vector<vector<double>> markov::softDecodingProblem(string msg)
{
	//First transition
	int numStates = tranMat.size();
	int msgSize = msg.size();
	char emmited = msg.at(0);
	int emmitedPos = charPos[emmited];
	double firstProb = double(1.0) / numStates;
	vector<vector<double>> forward = vector<vector<double>>(msgSize, vector<double>(numStates, 0));
	for (int i = 0; i < numStates; i++)
	{
		double localS = firstProb * emmiMat[i][emmitedPos];
		forward[0][i] = localS;
	}

	//Other transitions
	for (int i = 1; i < msg.size(); i++)
	{
		emmited = msg.at(i);
		emmitedPos = charPos[emmited];
		vector<double> tempS(forward[i-1]);
		for (int j = 0; j < numStates; j++)
		{
			double maxS = 0;
			int maxState;
			for (int k = 0; k < numStates; k++)
			{
				double localS = forward[i-1][k] * tranMat[k][j] * emmiMat[j][emmitedPos];
				maxS += localS;
			}

			tempS[j] = maxS;
		}

		forward[i] = tempS;
	}

	//Now does backward, to get the back graph we only reverse the msg
	//First transition
	vector<vector<double>> backward = vector<vector<double>>(msgSize, vector<double>(numStates, 0));
	for (int i = 0; i < numStates; i++)
	{
		double localS = firstProb;
		backward[msgSize - 1][i] = 1;
	}

	//Other transitions
	for (int i = msgSize - 2; i >= 0; i--)
	{
		emmited = msg.at(i + 1);
		emmitedPos = charPos[emmited];
		vector<double> tempS(backward[i + 1]);
		for (int j = 0; j < numStates; j++)
		{
			double maxS = 0;
			int maxState;
			for (int k = 0; k < numStates; k++)
			{
				double localS = backward[i + 1][k] * tranMat[j][k] * emmiMat[k][emmitedPos];
				maxS += localS;
			}

			tempS[j] = maxS;
		}

		backward[i] = tempS;
	}


	vector<vector<double>> soft = vector<vector<double>>(msgSize, vector<double>(numStates, 0));
	double sink = 0;
	for (int i = 0; i < numStates; i++)
	{
		sink += forward.back()[i];
	}

	//Gets all conditional probs
	for (int i = 0; i < msgSize; i++)
	{
		for (int j = 0; j < numStates; j++)
		{
			soft[i][j] = forward[i][j] * backward[i][j] / sink;
		}
	}

	return soft;
}

void markov::softDecodingGammas(string msg, vector<vector<double>>& gamma, vector<vector<double>>& gammaStar)
{
	//First transition
	int numStates = tranMat.size();
	int msgSize = msg.size();
	char emmited = msg.at(0);
	int emmitedPos = charPos[emmited];
	double firstProb = double(1.0) / numStates;
	vector<vector<double>> forward = vector<vector<double>>(msgSize, vector<double>(numStates, 0));
	for (int i = 0; i < numStates; i++)
	{
		double localS = firstProb * emmiMat[i][emmitedPos];
		forward[0][i] = localS;
	}

	//Other transitions
	for (int i = 1; i < msg.size(); i++)
	{
		emmited = msg.at(i);
		emmitedPos = charPos[emmited];
		vector<double> tempS(forward[i - 1]);
		for (int j = 0; j < numStates; j++)
		{
			double maxS = 0;
			int maxState;
			for (int k = 0; k < numStates; k++)
			{
				double localS = forward[i - 1][k] * tranMat[k][j] * emmiMat[j][emmitedPos];
				maxS += localS;
			}

			tempS[j] = maxS;
		}

		forward[i] = tempS;
	}

	//Now does backward, to get the back graph we only reverse the msg
	//First transition
	vector<vector<double>> backward = vector<vector<double>>(msgSize, vector<double>(numStates, 0));
	for (int i = 0; i < numStates; i++)
	{
		double localS = firstProb;
		backward[msgSize - 1][i] = 1;
	}

	//Other transitions
	for (int i = msgSize - 2; i >= 0; i--)
	{
		emmited = msg.at(i + 1);
		emmitedPos = charPos[emmited];
		vector<double> tempS(backward[i + 1]);
		for (int j = 0; j < numStates; j++)
		{
			double maxS = 0;
			int maxState;
			for (int k = 0; k < numStates; k++)
			{
				double localS = backward[i + 1][k] * tranMat[j][k] * emmiMat[k][emmitedPos];
				maxS += localS;
			}

			tempS[j] = maxS;
		}

		backward[i] = tempS;
	}


	vector<vector<double>> soft = vector<vector<double>>(numStates*numStates, vector<double>(msgSize - 1, 0));
	double sink = 0;
	for (int i = 0; i < numStates; i++)
	{
		sink += forward.back()[i];
	}

	//Fills Gamma**
	for (int i = 0; i < numStates; i++)
	{
		for (int j = 0; j < numStates; j++)
		{
			for (int k = 0; k < msgSize - 1; k++)
			{
				emmited = msg.at(k + 1);
				emmitedPos = charPos[emmited];
				double tempWeight = tranMat[i][j] * emmiMat[j][emmitedPos];
				soft[numStates*i + j][k] = forward[k][i] * backward[k + 1][j] * tempWeight / sink;
			}
		}
	}

	vector<vector<double>> saft = vector<vector<double>>(numStates, vector<double>(msgSize, 0));

	//Fills gamma
	for (int i = 0; i < numStates; i++)
	{
		for (int j = 0; j < msgSize; j++)
		{
			saft[i][j] = forward[j][i] * backward[j][i] / sink;
		}
	}

	gammaStar = soft;
	gamma = saft;
}

void markov::mStepBaum(vector<vector<double>>& gamma, vector<vector<double>>& gammaStar, string msg)
{
	int numStates = stateList.size();
	int numChars = charList.size();
	for (int i = 0; i < numStates; i++)
	{
		double normal = normalT(i, gammaStar);
		for (int j = 0; j < numStates; j++)
		{
			tranMat[i][j] = normal * computeT(i, j, gammaStar);
		}
	}

	for (int i = 0; i < numStates; i++)
	{
		double normal = normalE(charList, i, msg, gamma);
		for (int j = 0; j < numChars; j++)
		{
			emmiMat[i][j] = normal * computeE(i, charList[j], msg, gamma);
		}
	}
}

void markov::baumWelchLearning(int numIter, string msg)
{
	vector<vector<double>> gamma;
	vector<vector<double>> gammaStar;
	for (int i = 0; i < numIter; i++)
	{
		softDecodingGammas(msg, gamma, gammaStar);
		mStepBaum(gamma, gammaStar, msg);
	}
}

void markov::printTran()
{

	cout << setprecision(3) << fixed;
	int numStates = tranMat.size();
	int numChars = charList.size();
	cout << " ";
	for (int i = 0; i < numStates; i++)
	{
		cout << stateList[i] << " ";
	}
	cout << endl;
	for (int i = 0; i < numStates; i++)
	{
		cout << stateList[i] << " " << left;
		for (int j = 0; j < numStates; j++)
		{
			cout << tranMat[i][j] << " ";
		}
		cout << endl;
	}
}

void markov::printStates()
{
	for (int i=0; i<stateList.size(); i++)
	{
		string namae = stateList[i];
		cout << namae << "," << statePos[namae] << endl;
	}
}

void markov::printEmmision()
{
	int numStates = tranMat.size();
	int numChars = charList.size();
	cout <<  " ";
	for (int i = 0; i < numChars; i++)
	{
		cout << charList[i] << " ";
	}
	cout << endl;
	for (int i = 0; i < numStates; i++)
	{
		cout << stateList[i] << " ";
		for (int j = 0; j < numChars; j++)
		{
			cout << emmiMat[i][j] << " ";
		}
		cout << endl;
	}
}

void markov::printChars()
{
	printMap(charPos);
}

void markov::printConnected()
{
	int size = tranMat.size();
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			if (tranMat[i][j] != 0)
			{
				cout << stateList[i] << "->" << stateList[j] << endl;
			}
		}
	}
}

void normaliceVectorRange(vector<double>& x, int a, int b, double theta)
{
	double normal = 0;
	for (int i = a; i < b + 1; i++)
	{
		x[i] += theta;
		normal += x[i];
	}

	normal = 1.f / normal;
	for (int i = a; i < b + 1; i++)
	{
		x[i] *= normal;
	}
}

void maxValueVector(double & val, int & pos, vector<double>& x)
{
	val = x[0];
	pos = 0;
	for (int i = 1; i < x.size(); i++)
	{
		if (x[i] > val)
		{
			val = x[i];
			pos = i;
		}
	}
}

double computeT(int l, int k, vector<vector<double>>& gamma)
{
	double result = 0;
	int size = gamma[0].size();
	int numStates = sqrt(gamma.size());
	for (int i = 0; i < size; i++)
	{
		result += gamma[numStates*l + k][i];
	}

	return result;
}

double computeE(int k, char b, string msg, vector<vector<double>>& gamma)
{
	double result = 0;
	for (int i = 0; i < msg.size(); i++)
	{
		if (msg.at(i) == b)
		{
			result += gamma[k][i];
		}
	}
	return result;
}

double normalE(vector<char>& alphabet, int k, string msg, vector<vector<double>>& gamma)
{
	double normal = 0;
	for (auto& x : alphabet)
	{
		normal += computeE(k, x, msg, gamma);
	}

	normal = 1.l / normal;

	return normal;
}

double normalT(int l, vector<vector<double>>& gammaStar)
{
	double normal = 0;
	int numStates = sqrt(gammaStar.size());
	for (int i = 0; i < numStates; i++)
	{
		normal += computeT(l, i, gammaStar);
	}
	normal = 1.l / normal;

	return normal;
}
