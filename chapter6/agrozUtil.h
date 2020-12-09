#pragma once

#include<iostream>
#include<fstream>
#include<sstream>
#include<string>
#include<vector>
#include<algorithm>
#include<map>
#include<numeric>
#include<ctime>
#include <cstdlib>

using namespace std;

//Prints a vector
template <typename T>
void printVector(vector<T> vec);

//Prints a vector with \n
template <typename T>
void printVectorRows(vector<T> vec);

//Saves a vector into file
template <typename T>
void saveVector(vector<T> vec, string name);

//Reads matrix from text without spaces (each row is a matrix row)
template <typename T>
vector<vector<T>> readMatrix(string name, int rows, int columns, T example);

//Prints a matrix
template<typename T>
void printMatrix(vector<vector<T>> matrix);

//Prints matrix of adjacent vertices
template<typename T>
void printAdjacent(vector<vector<T>> matrix);

//Reads adjacency matrix
vector<vector<string>> readAdjacency(string name);

//Prints an euler cycle
void printEuler(vector<string> cycle);

//Prints (k,d)-mer in the book's format
void printKDmer(vector<vector<string>> matrix);

//Reads gapped patterns in the book's format
vector<vector<string>> readKDmer(string filename);

//Prints a sequence of pepties in book's format
void printPeptideSequence(vector<vector<int>> peptides);

//Prints a genome chapter 7
void printGenome(vector<vector<int>> gen);

/////////////////////////////////////////////////

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
	outdata.open(name + ".txt");
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
	vector<vector<T>> result(rows, vector<T>(columns));
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < columns; j++)
		{
			file >> result[i][j];
		}
	}
	file.close();

	return result;
}

//Prints matrix
template<typename T>
void printMatrix(vector<vector<T>> matrix)
{
	int col = matrix[0].size();
	int row = matrix.size();

	for (int i = 0; i < row; i++)
	{
		for (int j = 0; j < col; j++)
		{
			cout << matrix[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;
}

//Prints matrix of adjacent vertices
template<typename T>
void printAdjacent(vector<vector<T>> matrix)
{
	int row = matrix.size();

	for (int i = 0; i < row; i++)
	{
		int col = matrix[i].size();
		if (col > 1)
		{
			cout << matrix[i][0] << " -> ";
			for (int j = 1; j < col - 1; j++)
			{
				cout << matrix[i][j] << ",";
			}
			cout << matrix[i][col - 1];
			cout << endl;
		}
	}
	cout << endl;
}