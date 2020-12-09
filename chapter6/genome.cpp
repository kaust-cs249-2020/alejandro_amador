#include "genome.h"

//Identifies where the next cycle ends
int cycleNextPos(vector<vector<int>>& graph, int x);
//Given a red graph it returns all of its disconnected cycles
vector<vector<int>> redGraphCycles(vector<vector<int>> graph);
//Given a graph it returns all of its disconnected cycles
vector<vector<int>> graphCycles(vector<vector<int>> graph, vector<vector<int>> red);
//Transforms red edges to black edges
void transformColors(vector<int>& color);
//Given two graphs finds number of alternating cycles between them
int numberCyclesBi(vector<vector<int>>& g1, vector<vector<int>>& g2);
//Given an edge finds possible new edge in a graph for continuing cycle
vector<int> nextEdge(vector<vector<int>>& g, vector<int> edge);
//Removes an edge from a graph
int edgeRemove(vector<vector<int>>& g, vector<int> edge);
//Pushes an edge to a graph
void edgePush(vector<vector<int>>& g, vector<int> edge, int pos);
//Changes cycles to only black edges, also takes into account first one cycle starting with black (fixes it) and the 
//second starting with red.
void fixCycles(vector<vector<int>>& c);
//Function I do not understand at all, some weird fixes
void fixGenome(vector<vector<int>>& g);
//Finds nontrivial cycle
int findNonTrivial(vector<vector<int>> cycles);
//Given a nontrivial cycle it finds a path starting with red
vector<int> findRedPatho(vector<int> nonTrivial, vector<vector<int>> redEdges);
//Finds if an edge is in the graph and gives its position
int findEdge(vector<vector<int>> g, vector<int> edge);
//Given a kmer returns its reverse complement
string reverseComp(string g);
//Given a kmer and a string, returns all its positions if found
vector<int> findKmer(string g, string kmer);
//Given a string creates a map with all of its kmers
map<string, int> findKmers(string g, int k);
//Given a kmerTable and a new kmer it gives how many times it appears
void findKmerAppereances(map<string, int>& table, string kmer, int& counter);


//Given a permutation it computes reverse distance
int genome::greedySorting(vector<int>& p)
{
	int reversalDistance = 0;
	int size = p.size();
	for (int i = 0; i < size; i++)
	{
		if (p[i] != (i + 1))
		{
			kSortingReversal(p, i);
			reversalDistance++;
			printVector(p);
		}

		if (p[i] < 0)
		{
			p[i] *= -1;
			reversalDistance++;
		    printVector(p);
		}
	}

	return reversalDistance;
}

int genome::breakingPoints(vector<int> p)
{
	int size = p.size();
	int score = 0;
	if (p[0] != 1)
	{
		score++;
	}
	for (int i = 1; i < size; i++)
	{
		if (p[i] - p[i - 1] != 1)
		{
			score++;
		}
	}
	if (size + 1 - p[size - 1]  != 1)
	{
		score++;
	}
	return score;
}

vector<int> genome::chromosomeToCycle(vector<int> chrome)
{
	int size = chrome.size();
	vector<int> result(2 * size, 0);
	for (int j = 0; j < size; j++)
	{
		int i = chrome[j];
		if (i > 0)
		{
			result[2 * j + 1] = 2 * i;
			result[2 * j] = 2 * i - 1;
		}
		else
		{
			result[2 * j + 1] = -2 * i - 1;
			result[2 * j] = -2 * i;
		}
	}
	return result;
}

vector<int> genome::cycleToChromosome(vector<int> nodes)
{
	int size = nodes.size();
	vector<int> chrome(size / 2, 0);
	for (int i = 0; i < size/2; i++)
	{
		if (nodes[2 * i] < nodes[2 * i + 1])
		{
			chrome[i] = nodes[2 * i + 1] / 2;
		}
		else
		{
			chrome[i] = -nodes[2 * i] / 2;
		}
	}

	return chrome;
}

vector<int> genome::redCycleToChromosome(vector<int> nodesRed)
{

	vector<int> nodes = nodesRed;
	transformColors(nodes);
	int size = nodes.size();
	vector<int> chrome(size / 2, 0);
	for (int i = 0; i < size / 2; i++)
	{
		if (nodes[2 * i] < nodes[2 * i + 1])
		{
			chrome[i] = nodes[2 * i + 1] / 2;
		}
		else
		{
			chrome[i] = -nodes[2 * i] / 2;
		}
	}

	return chrome;
}

vector<vector<int>> genome::coloredEdges(vector<vector<int>> p)
{
	vector<vector<int>> edges;
	for (auto &c : p)
	{
		vector<int> nodes = chromosomeToCycle(c);
		nodes.push_back(nodes[0]);
		for (int j = 1; j <= c.size(); j++)
		{
			vector<int> tempEdge = { nodes[2 * j - 1], nodes[2 * j] };
			edges.push_back(tempEdge);
		}
	}

	return edges;
}

vector<vector<int>> genome::blackEdges(vector<vector<int>> p)
{
	vector<vector<int>> edges;
	for (auto &c : p)
	{
		vector<int> nodes = chromosomeToCycle(c);
		for (int i = 0; i < nodes.size()/2; i++)
		{
			vector<int> temp = { nodes[2 * i], nodes[2 * i + 1] };
			edges.push_back(temp);
		}
	}
	return edges;
}

vector<vector<int>> genome::graphtoGenome(vector<vector<int>> graph, vector<vector<int>> reds)
{
	vector<vector<int>> P;
	vector<vector<int>> cycles = graphCycles(graph, reds);
	fixCycles(cycles);
	for (auto &c : cycles)
	{
		vector<int> chrome = cycleToChromosome(c);
		P.push_back(chrome);
	}
	return P;
}

vector<vector<int>> genome::redGraphToGenome(vector<vector<int>> graph)
{
	vector<vector<int>> P;
	vector<vector<int>> cycles = redGraphCycles(graph);
	for (auto &c : cycles)
	{
		vector<int> chrome = redCycleToChromosome(c);
		P.push_back(chrome);
	}
	return P;
}

int genome::twoBreakDistance(vector<vector<int>> P, vector<vector<int>> Q)
{
	vector<vector<int>> colorP = coloredEdges(P);
	vector<vector<int>> colorB = coloredEdges(Q);

	int block = colorP.size();
	int cycles = numberCyclesBi(colorP, colorB);
	int dist = block - cycles;

	return dist;
}

void genome::shortestRearrangementScenario(vector<vector<int>> P, vector<vector<int>> Q)
{
	vector<vector<int>> tempGen = P;
	printGenome(tempGen);
	vector<vector<int>> redEdges = coloredEdges(P);
	vector<vector<int>> blueEdges = coloredEdges(Q);
	vector<vector<int>> breakGraph = redEdges;
	breakGraph.insert(breakGraph.begin(), blueEdges.begin(), blueEdges.end());

	vector<vector<int>> cycles = graphCycles(breakGraph, redEdges);
	int cyclePos = findNonTrivial(cycles);
	while (cyclePos != -1)
	{
		vector<int> nonTrivial = cycles[cyclePos];
		vector<int> redPatho = findRedPatho(nonTrivial, redEdges);
		vector<int> e12 = { redPatho[0], redPatho[1] };
		vector<int> e34 = { redPatho[2], redPatho[3] };
		vector<int> e14 = { redPatho[0], redPatho[3] };
		vector<int> e23 = { redPatho[1], redPatho[2] };


		int p1 = edgeRemove(redEdges, e12);
		edgePush(redEdges, e14, p1);
		int p2 = edgeRemove(redEdges, e34);
		edgePush(redEdges, e23, p2);


		breakGraph = redEdges;
		breakGraph.insert(breakGraph.begin(), blueEdges.begin(), blueEdges.end());
		cycles = graphCycles(breakGraph, redEdges);
		cyclePos = findNonTrivial(cycles);

		tempGen = twoBreakGenome(tempGen, redPatho[0], redPatho[1], redPatho[3], redPatho[2]);
		printGenome(tempGen);
	}
}

vector<vector<int>> genome::sharedKmers(string g1, string g2, int k)
{
	vector<vector<int>> positions;
	for (int i = 0; i <= g1.size() - k; i++)
	{
		string temp = g1.substr(i, k);
		vector<int> pos = findKmer(g2, temp);

		for (int j = 0; j < pos.size(); j++)
		{
			vector<int> found = { i, pos[j] };
			positions.push_back(found);
		}
	}
	return positions;
}

int genome::numberSharedKmers(string g1, string g2, int k)
{
	map<string, int> kmerTable = findKmers(g2, k);
	int count = 0;
	for (int i = 0; i <= g1.size() - k; i++)
	{
		string temp = g1.substr(i, k);
		findKmerAppereances(kmerTable, temp, count);
	}
	return count;
}

void genome::kSortingReversal(vector<int>& p, int k)
{
	vector<int>::iterator it1 = find(p.begin() + k, p.end(), k + 1);
	vector<int>::iterator it2 = find(p.begin() + k, p.end(), -k - 1);
	vector<int>::iterator it = p.begin() + k;
	if (it1 != p.end())
	{
		reverse(p.begin() + k, it1 + 1);
		for (it; it != it1 + 1; it++)
		{
			*it *= -1;
		}
	}
	else if (it2 != p.end())
	{
		reverse(p.begin() + k, it2 + 1);
		for (it; it != it2 + 1; it++)
		{
			*it *= -1;
		}
	}
	else
	{
		reverse(p.begin() + k, it1);
		for (it; it != it1; it++)
		{
			*it *= -1;
		}
	}

}

vector<int> genome::readPermutation(string name)
{
	ifstream file(name+".txt");
	vector<int> result;
	int temp;
	while (file >> temp)
	{
		result.push_back(temp);
	}
	file.close();
	return result;
}

vector<vector<int>> genome::readGenome(string name)
{
	ifstream file(name + ".txt");
	vector<vector<int>> result;
	vector<int> line;
	string tempStr;
	while (file >> tempStr)
	{
		int l1 = tempStr.find('(');
		int l2 = tempStr.find(')');
		if (l1 != string::npos && l2== string::npos)
		{
			line.push_back(stoi(tempStr.substr(1)));
		}
		else if (l1 == string::npos && l2 == string::npos)
		{
			line.push_back(stoi(tempStr));
		}
		else if (l1 != string::npos && l2 != string::npos)
		{
			int end = stoi(tempStr.substr(0, l2));
			line.push_back(end);
			result.push_back(line);

			line = {};
			line.push_back(stoi(tempStr.substr(l1 + 1)));
		}
		else if (l1 == string::npos && l2 != string::npos)
		{
			int end = stoi(tempStr.substr(0, l2));
			line.push_back(end);
			result.push_back(line);
		}
	}

	file.close();
	return result;
}

vector<vector<int>> genome::readGraph(string name)
{
	ifstream file(name + ".txt");
	vector<vector<int>> result;
	vector<int> line;
	string tempStr;
	while (file >> tempStr)
	{
		int l1 = tempStr.find('(');
		int l2 = tempStr.find(')');
		if (l1 != string::npos && l2 == string::npos)
		{
			line.push_back(stoi(tempStr.substr(1, tempStr.size() - 2)));
		}
		else if (l1 == string::npos && l2 != string::npos)
		{
			int size;
			if (tempStr.find(",") == string::npos)
			{
				size = tempStr.size() - 1;
			}
			else
			{
				size = tempStr.size() - 2;
			}
			int end = stoi(tempStr.substr(0, size));
			line.push_back(end);
			result.push_back(line);

			line = {};
		}
	}

	file.close();
	return result;
}

vector<vector<vector<int>>> genome::readMultipleGenomes(string name)
{
	ifstream file(name + ".txt");
	vector<vector<vector<int>>> result;
	string row;

	while (getline(file, row))
	{
		vector<vector<int>> tGen;
		vector<int> line;
		string tempStr;
		istringstream iss(row);
		while (iss >> tempStr)
		{
			int l1 = tempStr.find('(');
			int l2 = tempStr.find(')');
			if (l1 != string::npos && l2 == string::npos)
			{
				line.push_back(stoi(tempStr.substr(1)));
			}
			else if (l1 == string::npos && l2 == string::npos)
			{
				line.push_back(stoi(tempStr));
			}
			else if (l1 != string::npos && l2 != string::npos)
			{
				int end = stoi(tempStr.substr(0, l2));
				line.push_back(end);
				tGen.push_back(line);

				line = {};
				line.push_back(stoi(tempStr.substr(l1 + 1)));
			}
			else if (l1 == string::npos && l2 != string::npos)
			{
				int end = stoi(tempStr.substr(0, l2));
				line.push_back(end);
				tGen.push_back(line);
			}
		}

		result.push_back(tGen);
	}


	file.close();
	return result;
}

vector<vector<int>> genome::twoBreakGraph(vector<vector<int>> g, int a, int b, int c, int d)
{
	vector<vector<int>> result = g;
	vector<int> edge1 = { a,b };
	vector<int> edge2 = { c,d };
	vector<int> edge3 = { c,a };
	vector<int> edge4 = { d,b };
	int pos1 = edgeRemove(result, edge1);
	edgePush(result, edge4, pos1);
	int pos2 = edgeRemove(result, edge2);
	edgePush(result, edge3, pos2);
	return result;
}

vector<vector<int>> genome::twoBreakGenome(vector<vector<int>> P, int a, int b, int c, int d)
{
	vector<vector<int>> result;
	vector<vector<int>> genomeB = blackEdges(P);
	vector<vector<int>> genomeR = coloredEdges(P);
	vector<vector<int>> genomeGraph = genomeR;
	genomeGraph.insert(genomeGraph.end(), genomeB.begin(), genomeB.end());
	genomeGraph = twoBreakGraph(genomeGraph, a, b, c, d);
	result = graphtoGenome(genomeGraph, genomeR);

	return result;
}

string genome::loadString(string name)
{
	ifstream file(name + ".txt");
	string s;
	getline(file, s);
	file.close();
	return s;
}

int cycleNextPos(vector<vector<int>>& graph, int x)
{
	int pos = -99999;
	for (int i = 1; i < graph.size(); i++)
	{
		if (graph[i][0] == x)
		{
			pos = i;
			break;
		}
		else if (graph[i][1] == x)
		{
			pos = -i;
			break;
		}
	}

	return pos;
}

vector<vector<int>> redGraphCycles(vector<vector<int>> graph)
{
	vector<vector<int>> g = graph;
	vector<vector<int>> result;
	vector<int> temp;
	while (g.size() != 0)
	{
		temp = {};
		temp.push_back(g[0][0]);
		temp.push_back(g[0][1]);
		g.erase(g.begin());
		while (true)
		{
			int m1 = g[0][0];
			int m2 = g[0][1];
			if (m1 > m2)
			{
				temp.push_back(m1);
				temp.push_back(m2);
				g.erase(g.begin());
				break;
			}
			else
			{
				temp.push_back(m1);
				temp.push_back(m2);
				g.erase(g.begin());
			}
		}
		result.push_back(temp);
	}
	return result;
}

vector<vector<int>> graphCycles(vector<vector<int>> graph, vector<vector<int>> red)
{
	vector<vector<int>> r = red;
	vector<vector<int>> g = graph;
	vector<vector<int>> result;
	vector<int> temp;
	while (g.size()!=0)
	{
		temp = {};
		temp.push_back(g[0][0]);
		temp.push_back(g[0][1]);
		int nexPos = cycleNextPos(g, temp.back());
		while (nexPos!=-99999)
		{
			if (nexPos > 0)
			{
				temp.push_back(g[nexPos][0]);
				temp.push_back(g[nexPos][1]);
				g.erase(g.begin() + nexPos);
				nexPos = cycleNextPos(g, temp.back());
			}
			else if (nexPos < 0)
			{
				temp.push_back(g[-nexPos][1]);
				temp.push_back(g[-nexPos][0]);
				g.erase(g.begin() - nexPos);
				nexPos = cycleNextPos(g, temp.back());
			}
		}
		g.erase(g.begin());
		result.push_back(temp);
	}

	return result;
}

void transformColors(vector<int>& color)
{
	color.insert(color.begin(), color.back());
	color.pop_back();
}

int numberCyclesBi(vector<vector<int>>& g1, vector<vector<int>>& g2)
{
	int score = 0;
	while (g1.size() != 0)
	{
		vector<vector<int>> cycle;
		vector<int> current = g1[0];
		vector<int> next;
		cycle.push_back(current);
		g1.erase(g1.begin());
		int inter = 1;
		while (true)
		{
			if (inter % 2 == 0)
			{
				next = nextEdge(g1, current);
				cycle.push_back(next);
				current = next;
				inter++;
			}
			else
			{
				next = nextEdge(g2, current);
				cycle.push_back(next);
				current = next;
				inter++;
			}

			if (current[1] == cycle[0][0])
			{
				score++;
				break;
			}
		}
	}
	return score;
}

vector<int> nextEdge(vector<vector<int>>& g, vector<int> edge)
{
	vector<int> next;
	for (int i = 0; i < g.size(); i++)
	{
		if (edge[1] == g[i][0])
		{
			next = g[i];
			g.erase(g.begin() + i);
			break;
		}
		else if (edge[1] == g[i][1])
		{
			next = { g[i][1], g[i][0] };
			g.erase(g.begin() + i);
			break;
		}
	}
	return next;
}

int edgeRemove(vector<vector<int>>& g, vector<int> edge)
{
	for (int i = 0; i < g.size(); i++)
	{
		vector<int> temp = g[i];
		vector<int> edgeVar = { edge[1], edge[0] };
		if (temp == edge || temp == edgeVar)
		{
			g.erase(g.begin() + i);
			return i;
		}
	}
}

void edgePush(vector<vector<int>>& g, vector<int> edge, int pos)
{
	g.insert(g.begin() + pos, edge);
}

void fixCycles(vector<vector<int>>& c)
{
	for (int i = 0; i < c.size(); i++)
	{
		//First fixes position of the first cycle so that it starts with black
		c[i].insert(c[i].begin(), c[i].back());
		c[i].pop_back();
		c[i].insert(c[i].begin(), c[i].back());
		c[i].pop_back();

		//Now transforms them to only black
		int size1 = c[i].size();
		for (int j = size1 - 1; j >= 0; j -= 4)
		{
			c[i].erase(c[i].begin() + j);
			c[i].erase(c[i].begin() + j - 1);
		}
	}
}

void fixGenome(vector<vector<int>>& g)
{
	vector<vector<int>> temp = g;

	g[0].insert(g[0].begin(), g[0].back());
	g[0].pop_back();

	for (int i = 0; i < g[1].size(); i++)
	{
		g[1][i] = -temp[1][g[1].size() - 1 - i];
	}

	g[1].insert(g[1].begin(), g[1].back());
	g[1].pop_back();
}

int findNonTrivial(vector<vector<int>> cycles)
{
	int pos = -1;
	for (int i = 0; i < cycles.size(); i++)
	{
		int temp = cycles[i].size() - 2;
		if (temp > 2)
		{
			pos = i;
			break;
		}
	}

	return pos;
}

vector<int> findRedPatho(vector<int> nonTrivial, vector<vector<int>> redEdges)
{
	vector<int> path;
	for (int i = 0; i < 2; i++)
	{
		vector<int> check = { nonTrivial[2 * i], nonTrivial[2 * i + 1] };
		int pos = findEdge(redEdges, check);
		if (pos != -1)
		{
			path = { nonTrivial[2 * i], nonTrivial[2 * i + 1], nonTrivial[2 * i + 3], nonTrivial[2 * i + 5] };
			break;
		}
	}

	return path;
}

int findEdge(vector<vector<int>> g, vector<int> edge)
{
	int pos = -1;
	int size = g.size();
	for (int i = 0; i < size; i++)
	{
		vector<int> check1 = g[i];
		vector<int> check2 = { g[i][1], g[i][0] };
		if (edge == check1 || edge == check2)
		{
			pos = i;
			break;
		}
		
	}
	return pos;
}

string reverseComp(string g)
{
	string r = g;
	reverse(r.begin(), r.end());
	for (auto& c : r)
	{
		switch (c)
		{
		case('A'):
		{
			c = 'T';
			break;
		}
		case('T'):
		{
			c = 'A';
			break;
		}
		case('C'):
		{
			c = 'G';
			break;
		}
		case('G'):
		{
			c = 'C';
			break;
		}
		}
	}

	return r;
}

vector<int> findKmer(string g, string kmer)
{
	vector<int> pos;
	string r = reverseComp(kmer);
	int total = g.size();
	int local = kmer.size();
	for (int i = 0; i <= total - local; i++)
	{
		string temp = g.substr(i, local);
		if (temp == kmer || temp == r)
		{
			pos.push_back(i);
		}
	}
	return pos;
}

map<string, int> findKmers(string g, int k)
{
	map<string, int> result;
	int size = g.size();
	for (int i = 0; i <= size - k; i++)
	{
		string temp = g.substr(i, k);
		if (result.find(temp) == result.end())
		{
			result.insert({ temp, 1 });
		}
		else
		{
			result[temp]++;
		}
	}
	return result;
}

void findKmerAppereances(map<string, int>& table, string kmer, int & counter)
{
	string r = reverseComp(kmer);
	counter += table[kmer];
	if (r!=kmer)
	{
		counter += table[r];
	}


}
