#include "mathStuff.h"

vector<vector<double>> readVectorsDouble(string name)
{
	vector<vector<double>> res;
	ifstream file;
	file.open(name + ".txt");
	strstream s;
	string line;
	while (getline(file, line))
	{
		vector<double> newVec;
		s.clear();
		s << line;
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
				newVec.push_back(stod(temp));
			}
		}

		res.push_back(newVec);
	}

	return res;
}

vector<int> readVectorInt(string name)
{
	vector<int> res;
	ifstream file;
	file.open(name + ".txt");
	strstream s;
	string line;
	while (getline(file, line))
	{
		vector<int> newVec;
		s.clear();
		s << line;
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
				newVec.push_back(stoi(temp));
			}
		}

		res = newVec;
	}

	return res;
}

vector<vector<double>> readMatrixDouble(string name)
{
	vector<vector<double>> res;
	ifstream file;
	file.open(name + ".txt");
	strstream s;
	string line;
	while (getline(file, line))
	{
		vector<double> newVec;
		s.clear();
		s << line;
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
				newVec.push_back(stod(temp));
			}
		}

		res.push_back(newVec);
	}

	return res;
}
