#pragma once
#include<iostream>
#include<vector>
#include<fstream>

using namespace std;

//Prints a vector
template <typename T>
void printVector(vector<T> vec)
{
	vector<T>::iterator it;
	for (it = vec.begin(); it != vec.end(); it++)
	{
		cout << *it << " ";
	}
	cout << endl;
}

//Saves a vector into file
template <typename T>
void savesVector(vector<T> vec, string name)
{
	ofstream outdata;
	outdata.open(name+".txt"); 
	if (!outdata) 
	{ 
		cerr << "Error: file could not be opened" << endl;
		exit(1);
	}

	vector<T>::iterator it;
	for (it = vec.begin(); it != vec.end(); it++)
	{
		outdata << *it << " ";
	}
	outdata << endl;
	outdata.close();
}
