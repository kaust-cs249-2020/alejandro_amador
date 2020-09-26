#include "genome.h"

bool edgeCheck(string u, string v);
vector<string> glueStrings(vector<string> a, vector<string> b);
void glueVertices(int a, int b, vector<vector<string>>& adjacency);
void randomCycle(vector<vector<string>>& unused, vector<string>& cycle, vector<string>& nodes, vector<string>& nodesWithEdges);
void randomCycleNew(vector<vector<string>>& unused, vector<string>& oldCycle, vector<string>& nodes, vector<string>& nodesWithEdges);
void initialPath(vector<vector<string>>& unused, vector<string>& cycle, vector<string>& nodes, vector<string>& nodesWithEdges);
int findElementIndex(vector<string>& vec, string& x);
void mergeCycles(vector<string>& cycleA, vector<string>& cycleB, int pos);
int getUnbalanced(vector<vector<string>>& unused, vector<string>& nodes);
void getDegrees(vector<vector<string>>& unused, string& node, int& in, int& out);
vector<string> notAlone(vector<vector<string>>& unused, vector<string>& nodes);
string pathToGenome(vector<string> path);
string stringSpelledByGappedPatterns(vector<vector<string>> matrix, int k, int d);
string floatToBinary(int k, int n);
void lexiOrder(vector<vector<string>>& matrix);
vector<vector<string>> pairedCompositionGraph(vector<vector<string>> matrix);
vector<vector<string>> deBrujinFromKDmers(string filename);
vector<vector<string>> maximalNonBranchingPaths(vector<vector<string>> grafo);
vector<string> maximalToContigs(vector<vector<string>> paths);

genome::genome(int k)
{
	int size = pow(2, k);
	for (int i = 0; i < size; i++)
	{
		vector<string> temp;
		string node = floatToBinary(k, i);
		temp.push_back(node);
		string neighborOne = node.substr(1) + "0";
		string neighborTwo = node.substr(1) + "1";
		temp.push_back(neighborOne);
		temp.push_back(neighborTwo);
		m_graph.push_back(temp);
	}
}

genome::genome(string filename, classType type)
{
	if (type == patterns)
	{
		ifstream file(filename);
		string temp;
		while (getline(file, temp))
		{
			m_dna.push_back(temp);
		}
		file.close();
	}

	else if(type == graph)
	{
		m_graph = readAdjacency(filename);
	}
	else if (type == KDgraph)
	{
		m_graph = deBrujinFromKDmers(filename);
		int error = 0;
		for (int i = 0; i < m_graph.size(); i++)
		{
			int in, out;
			string node = m_graph[i][0];
			getDegrees(m_graph, node, in, out);
			if (in!=out)
			{
				error++;
			}
		}

		cout << error << endl;
	}
	else if (type == temp)
	{
		m_graph = readAdjacency(filename);
		vector<vector<string>> jaja = maximalNonBranchingPaths(m_graph);
		for (int i = 0; i < jaja.size(); i++)
		{
			printEuler(jaja[i]);
		}
	}
}

string genome::universalString()
{
	string result;
	vector<string> eulerian = eulerianCycle();
	result = pathToGenome(eulerian);
	result = result.substr(0,result.size() - m_graph[0][0].size());
	return result;
}

vector<string> genome::composition(int pos, int k)
{
	vector<string> result;
	int size = m_dna[pos].size() - k;
	for (int i = 0; i <= size; i++)
	{
		string temp = m_dna[pos].substr(i, k);
		result.push_back(temp);
	}
	return result;
}

string pathToGenome(vector<string> path)
{
	int numberKmers = path.size();
	int k = path[0].size() - 1;
	string result = path[0];
	for (int i = 1; i < numberKmers; i++)
	{
		result += path[i].at(k);
	}
	return result;
}

vector<vector<string>> genome::overlap()
{
	vector<vector<string>> result;
	int size = m_dna.size();
	for (int i = 0; i < size; i++)
	{
		vector<string> vertex = { m_dna[i] };
		for (int j = 0; j < size; j++)
		{
			if (j != i && edgeCheck(m_dna[i], m_dna[j]))
			{
				vertex.push_back(m_dna[j]);
			}
		}

		result.push_back(vertex);
	}
	return result;
}

vector<vector<string>> genome::pathGraph(int k)
{
	vector<vector<string>> result;
	int size = m_dna.size();
	int sizeDna = m_dna[0].size();
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j <= sizeDna - k; j++)
		{
			string vertex = m_dna[i].substr(j, k - 1);
			string outVertex = m_dna[i].substr(j + 1, k - 1);
			vector<string> row = { vertex, outVertex };
			result.push_back(row);
		}
	}
	return result;
}

vector<vector<string>> genome::deBrujin(int k)
{
	vector<vector<string>> result = pathGraph(k);
	int check = 0;
	while (check < result.size())
	{
		string vertexA = result[check][0];
		int check2 = check + 1;
		while (check2 < result.size())
		{
			string vertexB = result[check2][0];
			if (vertexA == vertexB)
			{
				glueVertices(check, check2, result);
			}
			else
			{
				check2++;
			}
		}
		check++;
	}
	return result;
}

vector<vector<string>> genome::pathFromKmers()
{
	vector<vector<string>> result;
	int size = m_dna.size();
	int k = m_dna[0].size();
	for (int i = 0; i < size; i++)
	{
		string vertexA = m_dna[i].substr(0, k - 1);
		string vertexB = m_dna[i].substr(1, k - 1);
		vector<string> edge = { vertexA, vertexB };
		result.push_back(edge);
	}
	return result;
}

vector<vector<string>> genome::deBrujinFromKmers()
{
	vector<vector<string>> result = pathFromKmers();
	int check = 0;
	while (check < result.size())
	{
		string vertexA = result[check][0];
		int check2 = check + 1;
		while (check2 < result.size())
		{
			string vertexB = result[check2][0];
			if (vertexA == vertexB)
			{
				glueVertices(check, check2, result);
			}
			else
			{
				check2++;
			}
		}
		check++;
	}

	return result;
}

vector<string> genome::contigsFromKmers()
{
	vector<string> result;
	vector<vector<string>> deBrujin = deBrujinFromKmers();
	vector<vector<string>> maximalPaths = maximalNonBranchingPaths(deBrujin);
	result = maximalToContigs(maximalPaths);
	return result;
}

vector<string> genome::eulerianCycle()
{
	vector<string> result;
	vector<string> nodes = generateNodes();
	vector<vector<string>> unused = m_graph;
	vector<string> nodesWithEdges = nodes;
	randomCycle(unused, result, nodes, nodesWithEdges);
	while (nodesWithEdges.size() > 0)
	{
		randomCycleNew(unused, result, nodes, nodesWithEdges);
	}
	return result;
}

vector<string> genome::generateNodes()
{
	vector<string> result;
	int size = m_graph.size();
	for (int i = 0; i < size; i++)
	{
		string node = m_graph[i][0];
		result.push_back(node);
	}
	return result;
}

vector<string> genome::eulerianPath()
{
	vector<string> result;
	vector<string> nodes = generateNodes();
	vector<vector<string>> unused = m_graph;
	vector<string> nodesWithEdges = nodes;
	initialPath(unused, result, nodes, nodesWithEdges);
	while (nodesWithEdges.size() > 0)
	{
		randomCycleNew(unused, result, nodes, nodesWithEdges);
	}
	return result;
}

vector<vector<string>> genome::eulerianPathKD()
{
	vector<vector<string>> result;
	vector<string> usualPath;
	vector<string> nodes = generateNodes();
	vector<vector<string>> unused = m_graph;
	vector<string> nodesWithEdges = nodes;
	initialPath(unused, usualPath, nodes, nodesWithEdges);
	while (nodesWithEdges.size() > 0)
	{
		randomCycleNew(unused, usualPath, nodes, nodesWithEdges);
	}

	int k = usualPath[0].size() / 2;
	for (int i = 0; i < usualPath.size(); i++)
	{
		string edgeA = usualPath[i].substr(0, k);
		string edgeB = usualPath[i].substr(k);
		vector<string> row = { edgeA, edgeB };
		result.push_back(row);
	}
	return result;
}

string genome::stringReconstruction()
{
	m_graph = deBrujinFromKmers();
	vector<string> eulerian = eulerianPath();
	string result = pathToGenome(eulerian);
	return result;
}

string genome::stringFromKDmers(int k, int d)
{
	vector<vector<string>> path = eulerianPathKD();
	string result = stringSpelledByGappedPatterns(path, k, d);
	return result;
}

vector<vector<string>> genome::generateKDmers(int k, int d)
{
	vector<vector<string>> result;
	int size = m_dna[0].size() - 2*k - d;
	for (int i = 0; i <= size; i++)
	{
		vector<string> temp = {m_dna[0].substr(i,k), m_dna[0].substr(i + k + d, k)};
		result.push_back(temp);
	}
	lexiOrder(result);
	return result;
}

bool edgeCheck(string u, string v)
{
	int size = u.size();
	bool result = (u.substr(1) == v.substr(0, size - 1));
	return result;
}

void glueVertices(int a, int b, vector<vector<string>>& adjacency)
{
	vector<string> gluedVertex = glueStrings(adjacency[a], adjacency[b]);
	adjacency[a] = gluedVertex;
	adjacency.erase(adjacency.begin() + b);
}

vector<string> glueStrings(vector<string> a, vector<string> b)
{
	vector<string> glued = a;
	glued.insert(glued.end(), b.begin() + 1, b.end());
	return glued;
}

void randomCycle(vector<vector<string>>& unused, vector<string>& cycle, vector<string>& nodes, vector<string>& nodesWithEdges)
{
	srand(time(0));
	int numberNodes = unused.size();
	int startPos = rand() % numberNodes;
	string startNode = nodes[startPos];
	cycle.push_back(startNode);
	int lastNode = startPos;
	while (true)
	{
		int numberEdges = unused[lastNode].size() - 1;
		if (numberEdges == 1)
		{
			int remove = findElementIndex(nodesWithEdges, unused[lastNode][0]);
			nodesWithEdges.erase(nodesWithEdges.begin() + remove);
		}
		int nextPos = rand() % numberEdges + 1;
		string nextNode = unused[lastNode][nextPos];
		cycle.push_back(nextNode);
		unused[lastNode].erase(unused[lastNode].begin() + nextPos);
		lastNode = findElementIndex(nodes, nextNode);
		if (lastNode == startPos)
		{
			return;
		}
	}
}

void randomCycleNew(vector<vector<string>>& unused, vector<string>& oldCycle, vector<string>& nodes, vector<string>& nodesWithEdges)
{
	srand(time(0));
	vector<string> cycle;
	int startPos;
	int connect;
	for (int i = 0; i < oldCycle.size(); i++)
	{
		int tempIndex = findElementIndex(nodes, oldCycle[i]);
		if (unused[tempIndex].size() > 1)
		{
			startPos = tempIndex;
			connect = i;
			break;
		}
	}
	int numberNodes = unused.size();
	string startNode = nodes[startPos];
	cycle.push_back(startNode);
	int lastNode = startPos;
	while (true)
	{
		int numberEdges = unused[lastNode].size() - 1;
		if (numberEdges == 1)
		{
			int remove = findElementIndex(nodesWithEdges, unused[lastNode][0]);
			nodesWithEdges.erase(nodesWithEdges.begin() + remove);
		}
		int nextPos = rand() % numberEdges + 1;
		string nextNode = unused[lastNode][nextPos];
		cycle.push_back(nextNode);
		unused[lastNode].erase(unused[lastNode].begin() + nextPos);
		lastNode = findElementIndex(nodes, nextNode);
		if (lastNode == startPos)
		{
			mergeCycles(oldCycle, cycle, connect);
			return;
		}
	}
}

int findElementIndex(vector<string>& vec, string& x)
{
	int index;
	vector<string>::iterator it = find(vec.begin(), vec.end(), x);
	if (it == vec.end())
	{
		index = -1;
	}
	else
	{
		index = distance(vec.begin(), it);
	}
	return index;
}

void mergeCycles(vector<string>& cycleA, vector<string>& cycleB, int pos)
{
	cycleA.erase(cycleA.begin() + pos);
	cycleA.insert(cycleA.begin() + pos, cycleB.begin(), cycleB.end());
}

void initialPath(vector<vector<string>>& unused, vector<string>& cycle, vector<string>& nodes, vector<string>& nodesWithEdges)
{
	srand(time(0));
	nodesWithEdges = notAlone(unused, nodes);
	int startPos = getUnbalanced(unused, nodes);
	string startNode = nodes[startPos];
	cycle.push_back(startNode);
	int lastNode = startPos;
	while (true)
	{
		int numberEdges = unused[lastNode].size() - 1;
		if (numberEdges == 1)
		{
			int remove = findElementIndex(nodesWithEdges, unused[lastNode][0]);
			nodesWithEdges.erase(nodesWithEdges.begin() + remove);
		}
		int nextPos = rand() % numberEdges + 1;
		string nextNode = unused[lastNode][nextPos];
		cycle.push_back(nextNode);
		unused[lastNode].erase(unused[lastNode].begin() + nextPos);
		lastNode = findElementIndex(nodes, nextNode);
		int endLine = findElementIndex(nodesWithEdges, nextNode);
		if (endLine < 0)
		{
			return;
		}
	}
}

int getUnbalanced(vector<vector<string>>& unused, vector<string>& nodes)
{
	for (int i = 0; i < nodes.size(); i++)
	{
		string node = nodes[i];
		int outDegree = 0;
		int inDegree = 0;
		getDegrees(unused, node, inDegree, outDegree);
		if (inDegree < outDegree)
		{
			return i;
		}
	}
}

void getDegrees(vector<vector<string>>& unused, string& node, int& in, int& out)
{
	in = 0;
	out = 0;
	for (int j = 0; j < unused.size(); j++)
	{
		if (unused[j][0] == node)
		{
			out += (unused[j].size() - 1);
		}

		for (auto s : unused[j])
		{
			if (s == node)
			{
				in++;
			}
		}
	}

	in--;
}

vector<string> notAlone(vector<vector<string>>& unused, vector<string>& nodes)
{
	vector<string> result;
	for (int i = 0; i < nodes.size(); i++)
	{
		string node = nodes[i];
		int outDegree = 0;
		int inDegree = 0;
		getDegrees(unused, node, inDegree, outDegree);
		if (outDegree > 0)
		{
			result.push_back(node);
		}
	}

	return result;
}

string floatToBinary(int k, int n)
{
	string binaryString;
	string conversor = "01";
	int temp = n/2;
	int digit = n%2;
	binaryString.push_back(conversor[digit]);
	while (temp != 0)
	{
		digit = temp % 2;
		binaryString = conversor[digit] + binaryString;
		temp = temp / 2;
	}
	int size = k - binaryString.size();
	for (int i = 0; i < size; i++)
	{
		binaryString = "0" + binaryString;
	}

	return binaryString;
}

void lexiOrder(vector<vector<string>>& matrix)
{
	int size = matrix.size();
	for (int i = 0; i < (size - 1); i++)
	{
		for (int j = i+1; j < size; j++)
		{
			if (matrix[i][0] > matrix[j][0])
			{
				vector<string> temp = matrix[i];
				matrix[i] = matrix[j];
				matrix[j] = temp;
			}
		}
	}
}

vector<vector<string>> pairedCompositionGraph(vector<vector<string>> matrix)
{
	int k = matrix[0][0].size();
	vector<vector<string>> result;
	int size = matrix.size();
	for (int i = 0; i < size; i++)
	{
		vector<string> row;
		string nodeA = matrix[i][0].substr(0,k-1) + matrix[i][1].substr(0,k-1);
		string nodeB = matrix[i][0].substr(1) + matrix[i][1].substr(1);
		row.push_back(nodeA);
		row.push_back(nodeB);
		result.push_back(row);
	}
	return result;
}

string stringSpelledByGappedPatterns(vector<vector<string>> matrix, int k, int d)
{
	vector<string> firstPatterns;
	vector<string> secondPatterns;
	for (int i = 0; i < matrix.size(); i++)
	{
		firstPatterns.push_back(matrix[i][0]);
		secondPatterns.push_back(matrix[i][1]);
	}
	string prefixString = pathToGenome(firstPatterns);
	string suffixString = pathToGenome(secondPatterns);
	for (int i = k + d ; i < prefixString.size(); i++)
	{
		if (prefixString[i] != suffixString[i - k - d])
		{
			cout << "There is no string spelled by the gapped patterns." << endl;
			exit(-10);
		}
	}
	prefixString += suffixString.substr(suffixString.size() - k - d);
	return prefixString;
}

vector<vector<string>> deBrujinFromKDmers(string filename)
{
	vector<vector<string>> result = readKDmer(filename);
	result = pairedCompositionGraph(result);
	int check = 0;
	while (check < result.size())
	{
		string vertexA = result[check][0];
		int check2 = check + 1;
		while (check2 < result.size())
		{
			string vertexB = result[check2][0];
			if (vertexA == vertexB)
			{
				glueVertices(check, check2, result);
			}
			else
			{
				check2++;
			}
		}
		check++;
	}
	return result;
}

vector<vector<string>> maximalNonBranchingPaths(vector<vector<string>> grafo)
{
	vector<vector<string>> graph = grafo;
	vector<vector<string>> unused = grafo;
	vector<string> nodes;
	for (int i = 0; i < graph.size(); i++)
	{
		nodes.push_back(graph[i][0]);
	}
	vector<vector<string>> result;
	for (int i = 0; i < graph.size(); i++)
	{
		string v = graph[i][0];
		int inD, outD;
		getDegrees(graph, v, inD, outD);
		if (inD != 1 || outD != 1)
		{
			if (outD > 0)
			{
				for (int j = 1; j < graph[i].size() ; j++)
				{
					string w = graph[i][j];
					int indox = findElementIndex(unused[i], w);
					unused[i].erase(unused[i].begin() + indox);
					int inW, outW;
					vector<string> nonBranchingPath = { v, w };
					getDegrees(graph, w, inW, outW);
					while (inW == 1 && outW == 1)
					{
						int index = findElementIndex(nodes, w);
						string u = graph[index][1];
						nonBranchingPath.push_back(u);
						getDegrees(graph, u, inW, outW);
						w = u;
						unused[index].pop_back();
					}
					result.push_back(nonBranchingPath);
				}
			}
		}
	}


	for (int j = 0; j < graph.size(); j++)
	{
		string v = unused[j][0];
		int inD, outD;
		getDegrees(unused, v, inD, outD);
		if (outD == 1)
		{
			string w = unused[j][1];
			int inW, outW;
			vector<string> cycle = { v,w };
			unused[j].pop_back();
			getDegrees(unused, w, inW, outW);
			while (outW == 1)
			{
				int index = findElementIndex(nodes, w);
				string u = unused[index][1];
				cycle.push_back(u);
				unused[index].pop_back();
				getDegrees(unused, u, inW, outW);
				w = u;
			}
			result.push_back(cycle);
		}
	}
	return result;
}

vector<string> maximalToContigs(vector<vector<string>> paths)
{
	vector<string> result;
	for (int i = 0; i < paths.size(); i++)
	{
		string temp;
		temp = pathToGenome(paths[i]);
		result.push_back(temp);
	}
	return result;
}