#include "spectrum.h"

spectrum::spectrum(string name)
{
	loadMassTable();
	int size = 1;
	nodeTranslator[0] = 0;
	ifstream file;
	file.open(name + ".txt");
	strstream s;
	string temp;
	getline(file, temp);
	s << temp;
	while (true)
	{
		string lala;
		s >> lala;
		if (lala.empty())
		{
			break;
		}
		else
		{
			int pept = stoi(lala);
			nodeTranslator[size] = pept;
			size++;
		}
	}

	adjList = vector<vector<pair<int, int>>>(size);
	//Makes possible connections
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			int weight = j - i;
			if (nameTable.find(weight) != nameTable.end())
			{
				pair<int, int> edge = { j,weight };
				adjList[i].push_back(edge);
			}
		}
	}
}

void spectrum::loadMassTable()
{
	ifstream file;
	file.open("massTable.txt");
	string temp;
	strstream s;
	while (getline(file, temp))
	{
		s.clear();
		s << temp;
		string pepto;
		int peptoMass;
		s >> pepto;
		s >> peptoMass;
		massTable[pepto] = peptoMass;
		nameTable[peptoMass] = pepto;
	}
	file.close();

}

void spectrum::printBookFormat()
{
	for (int i = 0; i < adjList.size(); i++)
	{
		for (auto& p : adjList[i])
		{
			cout << nodeTranslator[i] << "->" << nodeTranslator[p.first] << ":" << nameTable[p.second] << endl;
		}
	}
}

string spectrum::decodingIdealSpectrum()
{
	string result;
	vector<vector<int>> paths = allPaths(0, adjList.size() - 1);
	for (auto& p : paths)
	{
		vector<int> spectro = peptideToSpectrum(pathToPeptide(p));
		if (spectro == linearSpectrum)
		{
			return pathToPeptide(p);
		}
	}
}

string spectrum::peptideSequencingProblem()
{
	string peptide;
	vector<vector<int>> paths = allPaths(0, adjList.size() - 1);
	int max = -9999;
	for (auto& p : paths)
	{
		int temp = pathWeight(p);
		if (temp > max)
		{
			max = temp;
			peptide = pathToPeptide(p);
		}
	}


	return peptide;
}

vector<vector<int>> spectrum::allPaths(int u, int v)
{
	vector<vector<int>> paths;
	vector<int> currentPath;
	vector<bool> visited(adjList.size(), false);
	DFS(u, v, currentPath, paths, visited);

	return paths;
}

void spectrum::DFS(int& u, int& v, vector<int> currentPath, vector<vector<int>>& simplePaths, vector<bool> visited)
{
	if (visited[u] == true)
	{
		return;
	}

	visited[u] = true;
	currentPath.push_back(u);
	if (u == v)
	{
		simplePaths.push_back(currentPath);
		visited[u] = false;
		currentPath.pop_back();
		return;
	}

	vector<int> next = getOutNeighbors(u);
	for (auto& n : next)
	{
		DFS(n, v, currentPath, simplePaths, visited);
	}

	currentPath.pop_back();
	visited[u] = false;
}

vector<int> spectrum::getOutNeighbors(int u)
{
	vector<int> out;
	for (auto& e : adjList[u])
	{
		out.push_back(e.first);
	}

	return out;
}

string spectrum::pathToPeptide(vector<int>& path)
{
	string result = "";
	for (int i = 0; i < path.size() - 1; i++)
	{
		int from = path[i];
		int to = path[i + 1];
		for (auto& p : adjList[from])
		{
			if (p.first == to)
			{
				string amino = nameTable[p.second];
				result += amino;
				break;
			}
		}
	}

	return result;
}

vector<int> spectrum::peptideToSpectrum(string peptide)
{
	vector<int> spectram;
	//Prefixes
	for (int i = 1; i <= peptide.size(); i++)
	{
		string pre = peptide.substr(0, i);
		spectram.push_back(massFragment(pre));
	}
	//Suffixes
	for (int i = 1; i < peptide.size(); i++)
	{
		string pos = peptide.substr(i);
		spectram.push_back(massFragment(pos));
	}
	//Order
	sort(spectram.begin(), spectram.end());

	return spectram;
}

int spectrum::massFragment(string fragment)
{
	int mass = 0;
	for (int i = 0; i < fragment.size(); i++)
	{
		string temp = fragment.substr(i, 1);
		mass += massTable[temp];
	}

	return mass;
}

vector<int> spectrum::peptideIntoVector(string peptide)
{
	vector<int> result;
	vector<int> masses;
	for (int i = 1; i <= peptide.size(); i++)
	{
		string pre = peptide.substr(0, i);
		masses.push_back(massFragment(pre));
	}

	int max = masses.back();
	for (int i = 0; i < max; i++)
	{
		vector<int>::iterator it = find(masses.begin(), masses.end(), i + 1);
		if (it != masses.end())
		{
			result.push_back(1);
		}
		else
		{
			result.push_back(0);
		}
	}
	return result;
}

string spectrum::vectorIntoPeptide(vector<int>& peptiVector)
{
	vector<int> prefMass;
	for (int i = 0; i < peptiVector.size(); i++)
	{
		if (peptiVector[i] != 0)
		{
			prefMass.push_back(i + 1);
		}
	}

	string result = nameTable[prefMass[0]];
	for (int i = 0; i < prefMass.size() - 1; i++)
	{
		int maso = prefMass[i + 1] - prefMass[i];
		result += nameTable[maso];
	}

	return result;
}

int spectrum::pathWeight(vector<int>& path)
{
	int weight = 0;
	for (auto& x : path)
	{
		weight += nodeTranslator[x];
	}
	return weight;
}
