#include "weightedGraph.h"

int findUsed(vector<vector<int>> used, int node);

wegraph::wegraph(string name, int n)
{
	vector<vector<vector<int>>> result(n, vector<vector<int>>());
	ifstream file(name+".txt");
	string temp;
	while (getline(file, temp))
	{
		int l1 = temp.find(">");
		int l2 = temp.find(":", l1);
		int node = stoi(temp.substr(0, l1 - 1));
		int next = stoi(temp.substr(l1 + 1, l2 - l1));
		int weight = stoi(temp.substr(l2 + 1));
		vector<int> info = { next, weight };
		result[node].push_back(info);
	}
	file.close();
	m_adjList = result;
}

vector<int> wegraph::longestPath(int start, int end)
{
	vector<int> path = {end};
	vector<vector<int>> usedOnes;
	vector<int> first = { start, start };
	usedOnes.push_back(first);

	//Creates S and keeps track of used.
	int size = m_adjList.size();
	vector<int> s(size, 0);
	int node = start + 1;
	while (node <= end)
	{
		vector<int> tempUsed = { 0,node };
		vector<vector<int>> candidates = antecesores(node);
		int localS = 0;
		for (int i = 0; i < candidates.size(); i++)
		{
			int tempS = s[candidates[i][0]] + candidates[i][1];
			if (tempS > localS)
			{
				tempUsed[0] = candidates[i][0];
				localS = tempS;
			}
		}
		s[node] = localS;
		node++;
		usedOnes.push_back(tempUsed);
	}


	cout << s.back() << endl;
	//Backtracking
	int step = end;
	while (step > start + 1)
	{
		int temp = findUsed(usedOnes, step);
		path.insert(path.begin(), temp);
		step = temp;
	}

	return path;
}

vector<vector<int>> wegraph::antecesores(int node)
{
	vector<vector<int>> result;
	for (int i = 0; i < m_adjList.size(); i++)
	{
		vector<vector<int>> temp = m_adjList[i];
		for (int j = 0; j < temp.size(); j++)
		{
			if (temp[j][0] == node)
			{
				vector<int> info = { i, temp[j][1] };
				result.push_back(info);
				break;
			}
		}
	}
	return result;
}


void wegraph::print()
{
	for (int i = 0; i < m_adjList.size(); i++)
	{
		vector<vector<int>> temp = m_adjList[i];
		if (temp.size() == 0)
		{
			cout << i << endl;
		}
		else
		{
			cout << i << " ";
			for (int j = 0; j < temp.size(); j++)
			{
				cout << temp[j][0] << ":" << temp[j][1] << " ";
			}
			cout << endl;
		}
	}
}


int findUsed(vector<vector<int>> used, int node)
{
	int result = 0;
	for (int i = 0; i < used.size(); i++)
	{
		if (used[i][1] == node)
		{
			result = used[i][0];
			break;
		}
	}
	return result;
}