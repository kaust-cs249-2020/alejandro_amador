#pragma once
#include "mathStuff.h"


//Extra functions
//Given set of points and centers, returns max distance
double maxDataCenterDist(vector<vector<double>>& data, vector<vector<double>>& centers);
//Given data point and centers it returns distance
double dataPointCenterDist(vector<double>& point, vector<vector<double>>& centers);
//Euclidean distance of two points
double euclideanDistance(vector<double>& x, vector<double>& y);
//FarthestFirstTraversal
vector<vector<double>> farthestFirstTraversal(vector<vector<double>> data, double k);
//Computes distorsion of data and centers
double dataCenterDistorsion(vector<vector<double>>& point, vector<vector<double>>& centers);
//Computes center of gravity of dataset
vector<double> centerOfGravity(vector<vector<double>>& data);
//Lloyd algorithm
vector<vector<double>> lloydAlgorithm(vector<vector<double>>& data, int k);