#pragma once

#include <iostream>
#include <vector>
#include <map>
#include <fstream>
#include <string>
#include <strstream>
#include <iomanip>
#include <algorithm>  

using namespace std;

//Prints matrix
template<typename T>
void printMatrix(vector<vector<T>> x);

//Prints map
template<typename T, typename S>
void printMap(map<T, S> x);

//Prints vector
template<typename T>
void printVector(vector<T> x);

//Given a matrix it returns column i
template<typename T>
vector<T> getColumn(vector<vector<T>>& x, int i);

//Given a matrix it deletes column i
template<typename T>
void eraseColumn(vector<vector<T>>& x, int pos);

//Reads vectors of doubles
vector<vector<double>> readVectorsDouble(string name);

//Read vector of int
vector<int> readVectorInt(string name);

//Reads matrix of double
vector<vector<double>> readMatrixDouble(string name);

template<typename T>
inline void printMatrix(vector<vector<T>> x)
{
	int size = x.size();
	for (int i = 0; i < size; i++)
	{
		printVector(x[i]);
	}
	cout << endl;
}

template<typename T, typename S>
inline void printMap(map<T, S> x)
{
	for (auto& d : x)
	{
		cout << d.first << "," << d.second << endl;
	}

	cout << endl;
}

template<typename T>
inline void printVector(vector<T> x)
{
	for (int i = 0; i < x.size() - 1; i++)
	{
		cout << x[i] << " ";
	}

	cout << x.back() << endl;
}

template<typename T>
inline vector<T> getColumn(vector<vector<T>>& x, int i)
{
	vector<T> col;
	for (int k = 0; k < x.size(); k++)
	{
		col.push_back(x[k][i]);
	}

	return col;
}

template<typename T>
inline void eraseColumn(vector<vector<T>>& x, int pos)
{
	for (int i = 0; i < x.size(); i++)
	{
		x[i].erase(x[i].begin() + pos);
	}
}


