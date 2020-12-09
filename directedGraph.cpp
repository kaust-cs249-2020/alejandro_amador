#include "directedGraph.h"

vector<int> nodesWithNoOut(vector<vector<int>>& adj);
vector<int> nodesWithNoIn(vector<vector<int>>& adj);
void eraseEdge(int node, int edge, vector<vector<int>>& adj);
int inDegree(int node, vector<vector<int>>& adj);
int numberEdges(vector<vector<int>>& adj);

digraph::digraph(string name, int n)
{
	vector<vector<string>> temp = readAdjacency(name+".txt");
	for (int i = 0; i < n; i++)
	{
		vector<int> tempRow = { i };
		for (int j = 0; j < temp.size(); j++)
		{
			if (stoi(temp[j][0]) == i)
			{
				for (int k = 1; k < temp[j].size(); k++)
				{
					tempRow.push_back(stoi(temp[j][k]));
				}
				break;
			}
		}
		m_adjList.push_back(tempRow);
	}
}

vector<int> digraph::topologicalOrdering()
{
	vector<vector<int>> adj = m_adjList;
	vector<int> list;
	vector<int> candidates = nodesWithNoIn(adj);
	while (candidates.size() != 0)
	{
		int a = candidates[0];
		list.push_back(a);
		candidates.erase(candidates.begin());
		vector<int> outNodes = adj[a];
		for (int i = 1; i < outNodes.size(); i++)
		{
			int b = outNodes[i];
			eraseEdge(a, b, adj);
			if (inDegree(b, adj)==0)
			{
				candidates.push_back(b);
			}
		}
	}

	if (numberEdges(adj) != 0)
	{
		cout << "the input graph is not a DAG" << endl;
		printAdjacent(adj);
		exit(-2);
	}

	return list;
}


void digraph::print()
{
	for (int i = 0; i < m_adjList.size(); i++)
	{
		printVector(m_adjList[i]);
	}
}

vector<int> nodesWithNoOut(vector<vector<int>>& adj)
{
	vector<int> result;
	for (int i = 0; i < adj.size(); i++)
	{
		if (adj[i].size() == 1)
		{
			result.push_back(i);
		}
	}
	return result;
}

void eraseEdge(int node, int edge, vector<vector<int>>& adj)
{
	vector<int>::iterator it = find(adj[node].begin(), adj[node].end(), edge);
	adj[node].erase(it);
}

int inDegree(int node, vector<vector<int>>& adj)
{
	int count = 0;
	for (int i = 0; i < adj.size(); i++)
	{
		vector<int> tempRow = adj[i];
		for (int j = 1; j < tempRow.size(); j++)
		{
			if (tempRow[j] == node)
			{
				count++;
				break;
			}
		}
	}

	return count;
}

int numberEdges(vector<vector<int>>& adj)
{
	int count = 0;
	for (int i = 0; i < adj.size(); i++)
	{
		count += adj[i].size() - 1;
	}

	return count;
}

vector<int> nodesWithNoIn(vector<vector<int>>& adj)
{
	vector<int> result;
	for (int i = 0; i < adj.size(); i++)
	{
		if (inDegree(i, adj) == 0)
		{
			result.push_back(i);
		}
	}

	return result;
}