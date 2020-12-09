#include "parsiTree.h"

//Helper function to read node adjacency
void readNodeConection(string line, int& tempA, int& tempB, string& namae, int& leaves);
//Helper function read unrooted tree
void readNodeUnrooted(string line, int& tempA, int& tempB, string& namae, int& leaves, bool& ignore);
//Finds min value of 4 ints
int min4(int a, int b, int c, int d);
//Compare chars
int alphaChar(char a, char b);
//Gets the "parsimony score" of a node
int parsiScoreNode(parsiTree& T, vector<vector<int>>& parsiNodes, vector<char>& alphabet, char k, int v);
//Given the parsy scores of a node it gives the min
int minParsiOption(vector<vector<int>>& parsiNodes, int pos);
//Given the parsy scores of a node it gives the possible min
vector<int> minParsiOptions(vector<vector<int>>& parsiNodes, int pos);
//First try for backtrack
void tryBack(vector<vector<int>>& parsiNodes, parsiTree& T, vector<char>& alphabet);
//computes hamming ditance
int compareStrings(string s1, string s2);
//Deletes element from vector
void deleteElement(vector<int>& x, int val);

parsiTree::parsiTree(string name, treeType tipo)
{
	switch (tipo)
	{
	case rooted:
	{
		vector<vector<int>> nodesHolder;
		int maxNode = 0;
		int leaves = 0;

		ifstream file(name + ".txt");
		string temp;
		while (getline(file, temp))
		{
			int tempA, tempB;
			string tempName = "";
			readNodeConection(temp, tempA, tempB, tempName, leaves);
			vector<int> hold = { tempA, tempB };
			nodesHolder.push_back(hold);
			maxNode = max(maxNode, max(tempA, tempB));

			//Updates map
			if (tempName.size() != 0)
			{
				alphaTrans[tempB] = tempName;
			}
		}

		adjList = vector<vector<float>>(maxNode + 1, vector<float>(maxNode + 1, -1));
		directed = vector<vector<float>>(maxNode + 1, vector<float>(maxNode + 1, -1));
		for (auto& x : nodesHolder)
		{
			int tempA = x[0];
			int tempB = x[1];
			int tempWeight = -2;
			adjList[tempA][tempB] = tempWeight;
			adjList[tempB][tempA] = tempWeight;


			//Directed info
			directed[tempA][tempB] = tempWeight;
		}

		//Updates map
		for (int i = leaves; i < adjList.size(); i++)
		{
			alphaTrans[i] = "N";
		}
		file.close();

		break;
	}
	case unrooted:
	{
		vector<vector<int>> nodesHolder;
		int maxNode = 0;
		int leaves = 0;

		ifstream file(name + ".txt");
		string temp;
		bool ignore;
		while (getline(file, temp))
		{
			int tempA, tempB;
			string tempName = "";
			readNodeUnrooted(temp, tempA, tempB, tempName, leaves, ignore);
			if (!ignore)
			{
				vector<int> hold = { tempA, tempB };
				nodesHolder.push_back(hold);
				maxNode = max(maxNode, max(tempA, tempB));

				//Updates map
				if (tempName.size() != 0)
				{
					alphaTrans[tempB] = tempName;
				}
			}
		}

		adjList = vector<vector<float>>(maxNode + 1, vector<float>(maxNode + 1, -1));
		directed = vector<vector<float>>(maxNode + 1, vector<float>(maxNode + 1, -1));
		for (auto& x : nodesHolder)
		{
			int tempA = x[0];
			int tempB = x[1];
			int tempWeight = -2;
			adjList[tempA][tempB] = tempWeight;
			adjList[tempB][tempA] = tempWeight;
		}

		//Updates map
		for (int i = leaves; i < adjList.size(); i++)
		{
			alphaTrans[i] = "N";
		}
		file.close();
		break;
	}
	default:
		break;
	}

	initAlpha();
}

parsiTree::parsiTree(vector<parsiTree> bosque)
{
	*this = bosque[0];
	int mapSize = alphaTrans.size();
	for (int i = 1; i < bosque.size(); i++)
	{
		//Updates adjacency
		for (int l = 0; l < mapSize; l++)
		{
			for (int m = 0; m < mapSize; m++)
			{
				int temp = bosque[i].adjList[l][m];
				if (temp >= 0)
				{
					adjList[l][m] += temp;
				}
			}
		}

		//Updates map
		for (int j = 0; j < mapSize; j++)
		{
			string temp = bosque[i].alphaTrans[j];
			alphaTrans[j] += temp;
		}
	}


}

void parsiTree::initAlpha()
{
	alphabet = { 'A', 'C', 'G', 'T' };
}

void parsiTree::printDirected()
{
	int size = adjList.size();
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			if (directed[i][j] != -1)
			{
				cout << i << "->" << j << endl;
			}
		}
	}
}

void parsiTree::printMap()
{
	int size = adjList.size();
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			if (adjList[i][j] != -1)
			{
				cout << alphaTrans[i] << "->" << alphaTrans[j] << ":" << adjList[i][j] << endl;
			}
		}
	}
}

void parsiTree::printTheMap()
{
	for (auto& x : alphaTrans)
	{
		cout << x.first << " " << x.second << endl;
	}
}

void parsiTree::printNoWeight()
{
	int size = adjList.size();
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			if (adjList[i][j] != -1)
			{
				cout << i << "->" << j << endl;
			}
		}
	}
}

char parsiTree::getChar(int pos)
{
	char result = alphaTrans[pos].at(0);
	return result;
}

void parsiTree::pushMap(int pos, string x)
{
	alphaTrans[pos] = x;
}

vector<int> parsiTree::getNeighborsNode(int pos)
{
	vector<int> neighbors;
	for (int i = 0; i < adjList.size(); i++)
	{
		if (adjList[pos][i] != -1)
		{
			neighbors.push_back(i);
		}
	}

	return neighbors;
}

vector<int> parsiTree::getChildrenNode(int pos)
{
	vector<int> neighbors;
	int size = directed.size();
	for (int i = 0; i < size; i++)
	{
		if (directed[pos][i] != -1)
		{
			neighbors.push_back(i);
		}
	}

	return neighbors;
}

int parsiTree::getFatherNode(int pos)
{
	int neighbors = -1;
	for (int i = 0; i < directed.size(); i++)
	{
		if (directed[i][pos] != -1)
		{
			neighbors = i;
			break;
		}
	}

	return neighbors;
}

vector<char> parsiTree::charLeaves()
{
	vector<char> result;
	for (auto& x : alphaTrans)
	{
		char temp = x.second.at(0);
		if (temp != 'N')
		{
			result.push_back(temp);
		}
		else
		{
			break;
		}
	}
	return result;
}

vector<int> parsiTree::getRipes(vector<int>& tags)
{
	vector<int> result;
	for (int i = 0; i < tags.size(); i++)
	{
		if (tags[i] == 0)
		{
			vector<int> neighbors = getChildrenNode(i);
			bool ripe = true;
			for (auto& x : neighbors)
			{
				if (tags[x] != 1)
				{
					ripe = false;
					break;
				}
			}

			if (ripe)
			{
				result.push_back(i);
			}
		}
	}


	return result;
}

int parsiTree::sizeMapString()
{
	return alphaTrans[0].size();
}

parsiTree parsiTree::subCharTree(int pos)
{
	parsiTree subTree = *this;
	for (auto& x : subTree.alphaTrans)
	{
		if (x.second != "N")
		{
			string cut = x.second.substr(pos, 1);
			x.second = cut;
		}
	}

	return subTree;
}

vector<parsiTree> parsiTree::subCharForest()
{
	int size = alphaTrans[0].size();
	vector<parsiTree> bosque;
	for (int i = 0; i < size; i++)
	{
		parsiTree temp = subCharTree(i);
		bosque.push_back(temp);
	}

	return bosque;
}

void parsiTree::addRoot()
{
	int size = adjList.size();
	//Updates map
	alphaTrans[size] = "N";
	//Updates adj
	pushNode();
	//Connects
	int child = getNeighborsNode(size - 1).back();
	connectNodes(size-1, child, -1);
	connectNodes(size, child, -2);
	connectNodes(size, size - 1, -2);

	//updates directed
	for (int i = 0; i < size; i++)
	{
		directed[i].push_back(-1);
	}
	vector<float> newNode(size + 1, -1);
	directed.push_back(newNode);
	fixChildren();
}

void parsiTree::removeRoot()
{
	int size = adjList.size();
	//Temp new conection
	vector<int> brothers = getChildrenNode(size-1);

	//Updates map
	alphaTrans.erase(size-1);
	//Updates adj
	for (int i = 0; i < size; i++)
	{
		adjList[i].pop_back();
	}
	adjList.pop_back();

	//Connects
	int b0 = brothers[0];
	int b1 = brothers[1];
	string bs1 = alphaTrans[b0];
	string bs2 = alphaTrans[b1];
	int weight = compareStrings(bs1, bs2);
	connectNodes(b0, b1, weight);

	//updates directed
	directed = vector<vector<float>>(size-1, vector<float>(size-1, -1));


}

vector<parsiTree> parsiTree::nearestNeighborsProblem(int a, int b)
{
	vector<int> neighA = getNeighborsNode(a);
	vector<int> neighB = getNeighborsNode(b);
	deleteElement(neighA, b);
	deleteElement(neighB, a);

	int w = neighA[0];
	int x = neighA[1];
	int y = neighB[0];
	int z = neighB[1];

	parsiTree N1 = *this;
	N1.connectNodes(a, x, -1);
	N1.connectNodes(a, y, -2);
	N1.connectNodes(b, y, -1);
	N1.connectNodes(b, x, -2);

	parsiTree N2 = *this;
	N2.connectNodes(a, x, -1);
	N2.connectNodes(a, z, -2);
	N2.connectNodes(b, z, -1);
	N2.connectNodes(b, x, -2);

	vector<parsiTree> result = { N1, N2 };
	return result;
}

vector<int> parsiTree::internalNodes()
{
	int size = adjList.size();
	vector<int> intern;
	for (int i = 0; i < size; i++)
	{
		int neighbors = 0;
		for (int j = 0; j < size; j++)
		{
			if (adjList[i][j] != -1)
			{
				neighbors++;
			}
		}

		if (neighbors > 1)
		{
			intern.push_back(i);
		}
	}

	return intern;
}

vector<vector<int>> parsiTree::internalEdges()
{
	vector<int> interns = internalNodes();
	vector<vector<int>> edges;
	int size = interns.size();
	for (int i = 0; i < size-1; i++)
	{
		int a = interns[i];
		for (int j = i + 1; j < size; j++)
		{
			int b = interns[j];
			if (adjList[a][b] != -1)
			{
				vector<int> couple = { a,b };
				edges.push_back(couple);
			}
		}
	}
	return edges;
}

void parsiTree::disconectDirected(int a, int b)
{
	int temp = directed[a][b];
	if (temp != -1)
	{
		directed[a][b] = -1;
	}
	else
	{
		directed[b][a] = -1;
	}
}

void parsiTree::connectDirected(int a, int b)
{
	directed[a][b] = -2;
}

void parsiTree::fixChildren()
{
	int father = directed.size() - 1;
	vector<int> parents;
	vector<int> newParents = { father };
	while (!newParents.empty())
	{
		parents = newParents;
		newParents.clear();
		//Properly connects father to son
		for (auto& x : parents)
		{
			vector<int> children = getChildrenKnownFather(x);
			for (auto& c : children)
			{
				directed[x][c] = -2;
				directed[c][x] = -1;
				newParents.push_back(c);
			}
		}
	}
}

vector<int> parsiTree::getChildrenKnownFather(int x)
{
	vector<int> neighbors = getNeighborsNode(x);
	int father = getFatherNode(x);
	deleteElement(neighbors, father);
	return neighbors;
}

void readNodeConection(string line, int & tempA, int & tempB, string & namae, int& leaves)
{
	int lim1 = line.find('-');
	int t1 = line.find('A');
	int t2 = line.find('C');
	int t3 = line.find('G');
	int t4 = line.find('T');

	tempA = stoi(line.substr(0, lim1));
	if (t1 != string::npos || t2 != string::npos || t3 != string::npos || t4 != string::npos)
	{
		tempB = leaves;
		namae = line.substr(lim1 + 2);
		leaves++;
	}
	else
	{
		tempB = stoi(line.substr(lim1 + 2));
	}
}

void readNodeUnrooted(string line, int & tempA, int & tempB, string & namae, int & leaves, bool& ignore)
{
	//Finds if it has to ignore it or not
	ignore = false;

	//It ignores whenever we have i->j for i<j
	int tem1 = line.find('A');
	int tem2 = line.find('C');
	int tem3 = line.find('G');
	int tem4 = line.find('T');
	if (tem1 == 0 || tem2 == 0 || tem3 == 0 || tem4 == 0)
	{
		ignore = true;
		return;
	}


	//If no ignore then continues
	int lim1 = line.find('-');
	int t1 = line.find('A');
	int t2 = line.find('C');
	int t3 = line.find('G');
	int t4 = line.find('T');

	tempA = stoi(line.substr(0, lim1));
	if (t1 != string::npos || t2 != string::npos || t3 != string::npos || t4 != string::npos)
	{
		tempB = leaves;
		namae = line.substr(lim1 + 2);
		leaves++;
	}
	else
	{
		tempB = stoi(line.substr(lim1 + 2));
		if (tempA < tempB)
		{
			ignore = true;
		}
	}
}

int min4(int a, int b, int c, int d)
{
	return min(a,min(b,min(c,d)));
}

int alphaChar(char a, char b)
{
	int result;
	if (a == b)
	{
		result = 0;
	}
	else
	{
		result = 1;
	}

	return result;
}

int parsiScoreNode(parsiTree & T, vector<vector<int>>& parsiNodes, vector<char>& alphabet, char k, int v)
{
	int min1 = 9999;
	int min2 = 9999;
	vector<int> children = T.getChildrenNode(v);
	int a = children[0];
	int b = children[1];
	for (int i = 0; i < alphabet.size(); i++)
	{
		char symbol = alphabet[i];
		int temp1 = parsiNodes[a][i] + alphaChar(k, symbol);
		int temp2 = parsiNodes[b][i] + alphaChar(k, symbol);
		if (temp1 < min1)
		{
			min1 = temp1;
		}
		if (temp2 < min2)
		{
			min2 = temp2;
		}
	}

	return min1 + min2;
}

int minParsiOption(vector<vector<int>>& parsiNodes, int pos)
{
	int result = parsiNodes[pos][0];
	for (int i = 1; i < parsiNodes[pos].size(); i++)
	{
		int temp = parsiNodes[pos][i];
		if (temp < result)
		{
			result = temp;
		}
	}

	return result;
}

vector<int> minParsiOptions(vector<vector<int>>& parsiNodes, int pos)
{
	int minPos = minParsiOption(parsiNodes, pos);
	vector<int> result;
	for (int i = 0; i < parsiNodes[pos].size(); i++)
	{
		int temp = parsiNodes[pos][i];
		if (temp == minPos)
		{
			result.push_back(i);
		}
	}

	return result;
}

void tryBack(vector<vector<int>>& parsiNodes, parsiTree& T, vector<char>& alphabet)
{
	int size = T.numberNodes();
	vector<int> best = minParsiOptions(parsiNodes, size - 1);
	char bestChar = alphabet[best[0]];
	string root(1, bestChar);
	T.pushMap(size - 1, root);

	for (int i = size - 2; i >= 0; i--)
	{
		int father = T.getFatherNode(i);
		char fatherChar = T.getChar(father);
		vector<int> maybe = minParsiOptions(parsiNodes, i);
		if (maybe.size() == 1)
		{
			char sonChar = alphabet[maybe[0]];
			string mapo(1, sonChar);
			T.pushMap(i, mapo);
			T.connectNodes(father, i, alphaChar(fatherChar, sonChar));
		}
		else
		{
			bool tie = true;
			for (auto& x : maybe)
			{
				char maybeChar = alphabet[x];
				if (maybeChar == fatherChar)
				{
					string mapo(1, fatherChar);
					T.pushMap(i, mapo);
					T.connectNodes(father, i, 0);
					tie = false;
					break;
				}
			}
			if (tie)
			{
				char sonChar = alphabet[maybe[0]];
				string mapo(1, sonChar);
				T.pushMap(i, mapo);
				T.connectNodes(father, i, alphaChar(fatherChar, sonChar));
			}
		}
	}
}

int compareStrings(string s1, string s2)
{
	int result = 0;
	for (int i = 0; i < s1.size(); i++)
	{
		if (s1.at(i) != s2.at(i))
		{
			result++;
		}
	}
	return result;
}

void deleteElement(vector<int>& x, int val)
{
	vector<int>::iterator pos = find(x.begin(), x.end(), val);
	if (pos != x.end())
	{
		x.erase(pos);
	}
}


int smallParsimony(parsiTree T, parsiTree& TFilled)
{
	//Tree for storing the result
	TFilled = T;

	//Gets all data from the tree
	vector<char> leaves = T.charLeaves();
	vector<char> alphabet = T.getAlphabet();
	int alphaSize = alphabet.size();
	int size = T.numberNodes();
	int numLeaves = leaves.size();

	//Algorithm
	vector<vector<int>> nodeParsimony(size, vector<int>(alphaSize, 999999));
	vector<int> tag(size, 0);

	//Init tags
	for (int i = 0; i < numLeaves; i++)
	{
		tag[i] = 1;
		for (int j=0; j<alphaSize; j++)
		{
			if (T.getChar(i) == alphabet[j])
			{
				nodeParsimony[i][j] = 0;
			}
		}
	}

	vector<int> ripes = T.getRipes(tag);
	while (!ripes.empty())
	{
		int v = ripes[0];
		tag[v] = 1;
		for (int i = 0; i < alphaSize; i++)
		{
			char k = alphabet[i];
			int s_k_v = parsiScoreNode(T, nodeParsimony, alphabet, k, v);
			nodeParsimony[v][i] = s_k_v;
		}

		ripes = T.getRipes(tag);
	}

	int result = minParsiOption(nodeParsimony, size - 1);

	//Backtrack
	tryBack(nodeParsimony, TFilled, alphabet);

	return result;
}

int smallParsimonyString(parsiTree T, parsiTree & Tfilled)
{
	vector<parsiTree> bosque = T.subCharForest();
	vector<parsiTree> bosqueFilled;
	int totalScore = 0;
	for (auto& x : bosque)
	{
		parsiTree temp;
		totalScore += smallParsimony(x, temp);
		bosqueFilled.push_back(temp);
	}

	Tfilled = parsiTree(bosqueFilled);
	return totalScore;
}

int smallParsimonyUnrooted(parsiTree T, parsiTree & Tfilled)
{
	parsiTree tempo = T;
	tempo.addRoot();
	int result = smallParsimonyString(tempo, Tfilled);
	Tfilled.removeRoot();
	return result;
}

parsiTree nearestNeighborInterchange(parsiTree T)
{
	int score = 999999;
	parsiTree tree = T;
	parsiTree R;
	parsiTree solved;
	int newScore = smallParsimonyUnrooted(tree, R);
	solved = R;
	parsiTree newT = tree;
	while(newScore < score)
	{
		cout << newScore << endl;
		solved.printMap();
		cout << endl;
		score = newScore;
		tree = newT;
		vector<vector<int>> internal = tree.internalEdges();
		for (auto& inter : internal)
		{
			vector<parsiTree> neighbors = tree.nearestNeighborsProblem(inter[0], inter[1]);
			for (auto& arbol : neighbors)
			{
				int tempScore = smallParsimonyUnrooted(arbol, R);
				if(tempScore < newScore) 
				{
					newScore = tempScore;
					newT = arbol;
					solved = R;
				}
			}
		}
	}

	return R;
}
