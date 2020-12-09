#include "clusters.h"


double maxDataCenterDist(vector<vector<double>>& data, vector<vector<double>>& centers)
{
	double result = 0;
	for (auto& x : data)
	{
		double temp = dataPointCenterDist(x, centers);
		if (temp > result)
		{
			result = temp;
		}
	}

	return result;
}

double dataPointCenterDist(vector<double>& point, vector<vector<double>>& centers)
{
	double result = 99999;
	for (auto& x : centers)
	{
		double temp = euclideanDistance(point, x);
		if (temp < result)
		{
			result = temp;
		}
	}

	return result;
}

double euclideanDistance(vector<double>& x, vector<double>& y)
{
	double result = 0;
	for (int i = 0; i < x.size(); i++)
	{
		result += (x[i] - y[i])*(x[i] - y[i]);
	}
	result = sqrt(result);

	return result;
}

vector<vector<double>> farthestFirstTraversal(vector<vector<double>> data, double k)
{
	vector<vector<double>> dato(data);
	vector<vector<double>> centers = { dato[0] };
	dato.erase(dato.begin());
	while (centers.size() < k)
	{
		int maximizer = 0;
		double max = dataPointCenterDist(dato[0], centers);
		for (int i = 1; i < dato.size(); i++)
		{
			double temp = dataPointCenterDist(dato[i], centers);
			if (temp > max)
			{
				max = temp;
				maximizer = i;
			}
		}

		centers.push_back(dato[maximizer]);
		dato.erase(dato.begin() + maximizer);
	}

	return centers;
}

double dataCenterDistorsion(vector<vector<double>>& data, vector<vector<double>>& centers)
{
	int size = data.size();
	double result = 0;
	for (int i = 0; i < size; i++)
	{
		result += pow(dataPointCenterDist(data[i], centers), 2);
	}

	result /= size;

	return result;
}

vector<double> centerOfGravity(vector<vector<double>>& data)
{
	vector<double> gravity = data[0];
	for (int i = 1; i < data.size(); i++)
	{
		for (int j = 0; j < gravity.size(); j++)
		{
			gravity[j] += data[i][j];
		}
	}

	int size = data.size();
	for (auto& x : gravity)
	{
		x /= size;
	}

	return gravity;
}

vector<vector<double>> lloydAlgorithm(vector<vector<double>>& data, int k)
{
	vector<vector<double>> centers;
	for (int i = 0; i < k; i++)
	{
		centers.push_back(data[i]);
	}
	vector<vector<double>> lastCenters(centers);
	while (true)
	{
		vector<vector<vector<double>>> clusters(k);
		for (auto& x : data)
		{
			double max = 99999;
			int closer = -1;
			for (int i = 0; i < k; i++)
			{
				double temp = euclideanDistance(x, centers[i]);
				if (temp < max)
				{
					max = temp;
					closer = i;
				}
			}
			clusters[closer].push_back(x);
		}

		for (int i = 0; i < k; i++)
		{
			vector<double> temp = centerOfGravity(clusters[i]);
			centers[i] = temp;
		}

		if (lastCenters == centers)
		{
			break;
		}
		else
		{
			lastCenters = centers;
		}
	}

	return centers;
}
