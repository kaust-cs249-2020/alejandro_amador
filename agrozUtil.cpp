#include "agrozUtil.h"

vector<vector<string>> readAdjacency(string name)
{
	vector<vector<string>> result;
	ifstream file(name);
	string temp;
	while (getline(file, temp))
	{
		vector<string> row;
		int lastPos = temp.find("-");
		string nodeA = temp.substr(0, lastPos - 1);
		row.push_back(nodeA);
		int nextComma = temp.find(",", lastPos);
		lastPos = lastPos + 3;
		string newNode;
		while (nextComma != string::npos)
		{
			int length = nextComma - lastPos;
			newNode = temp.substr(lastPos, length);
			row.push_back(newNode);
			lastPos = nextComma + 1;
			nextComma = temp.find(",", lastPos);
		}
		newNode = temp.substr(lastPos);
		row.push_back(newNode);
		result.push_back(row);
	}
	file.close();
	return result;
}

void printEuler(vector<string> cycle)
{
	int size = cycle.size();
	for (int i = 0; i < size-1; i++)
	{
		cout << cycle[i] << "->";
	}
	cout << cycle.back() << endl;
}

void printKDmer(vector<vector<string>> matrix)
{
	int col = matrix[0].size();
	int row = matrix.size();

	for (int i = 0; i < row; i++)
	{
		cout << "(" << matrix[i][0] << "|" << matrix[i][1] << ") ";
	}
	cout << endl;
}

vector<vector<string>> readKDmer(string filename)
{
	vector<vector<string>> result;
	ifstream file(filename);
	string temp;
	while (getline(file, temp))
	{
		vector<string> row;
		int division = temp.find("|");
		string kmerA = temp.substr(0, division);
		string kmerB = temp.substr(division + 1);
		row.push_back(kmerA);
		row.push_back(kmerB);
		result.push_back(row);
	}
	file.close();
	return result;
}

void printPeptideSequence(vector<vector<int>> peptides)
{
	for (int i = 0; i < peptides.size(); i++)
	{
		for (int j = 0; j < peptides[i].size() - 1; j++)
		{
			cout << peptides[i][j] << "-";
		}
		cout << peptides[i].back() << " ";
	}
	cout << endl;
}
