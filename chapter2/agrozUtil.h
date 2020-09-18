#pragma once
#include<iostream>
#include<vector>
#include<fstream>

using namespace std;

//Prints a vector
template <typename T>
void printVector(vector<T> vec)
{
	for (auto s : vec)
	{
		cout << s << " ";
	}
	cout << endl;
}

//Prints a vector with \n
template <typename T>
void printVectorRows(vector<T> vec)
{
	for (auto s : vec)
	{
		cout << s << endl;
	}
	cout << endl;
}


//Saves a vector into file
template <typename T>
void saveVector(vector<T> vec, string name)
{
	ofstream outdata;
	outdata.open(name+".txt"); 
	if (!outdata) 
	{ 
		cerr << "Error: file could not be opened" << endl;
		exit(1);
	}

	for (auto s : vec)
	{
		outdata << s << " ";
	}
	cout << endl;

	outdata << endl;
	outdata.close();
}

//Reads matrix from text without spaces (each row is a matrix row)
template <typename T>
vector<vector<T>> readMatrix(string name, int rows, int columns, T example)
{
	ifstream file(name);
	vector<vector<T>> result(columns, vector<T>(4));
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < columns; j++)
		{
			file >> result[j][i];
		}
	}
	file.close();

	return result;
}

//Prints matrix
template<typename T>
void printMatrix(vector<vector<T>> matrix)
{
	int n = matrix[0].size();
	int m = matrix.size();

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			cout << matrix[j][i] << " ";
		}
		cout << endl;
	}
	cout << endl;
}