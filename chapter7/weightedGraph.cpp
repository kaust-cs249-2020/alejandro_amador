#include "weightedGraph.h"

//Helper function to read file of weighted matrix
void readNodeWeight(string row, int& nodeA, int& nodeB, int& weight);
//Finds out if node is already in a path
bool nodeInPath(vector<int>& path, int node);
//Finds two leaves such that condition for phylo algorithm holds
void findLeavesAdditive(int& i, int& k, int n, vector<vector<float>>& D);
//Given a matrix it deletes row and column n
void removeRowColumn(int n, vector<vector<float>>& D);
//Constructs an n-dim array with elements 0,1,....,n-1
vector<vector<int>> initCluster(int n);
//Given two clusters it computes its distance
float clusterDistance(vector<int>& c1, vector<int>& c2, vector<vector<float>>& D);
//Given two clusters it computes its min distance
float clusterDistanceMin(vector<int>& c1, vector<int>& c2, vector<vector<float>>& D);
//Given a collection of clusters, it returns the two closest ones as indices.
vector<int> findClosestClusters(vector<vector<int>>& clusters, vector<vector<float>>& D, float& lowDist);
//Merge two clusters within a collection of clusters given its indices
void clusterMerger(vector<vector<int>>& clusters, vector<int>& orders, vector<int>& indices, int& t1, int& t2);
//Helper function, given distance matrix D removes row and column corresponding to node u
void removeNodeFromMatrix(vector<vector<float>>& D, int u);
//Helper function, given distance matrix D removes rows and columns corresponding to clusters Ci and Cj
void removeClustersFromMatrix(vector<vector<float>>& D, vector<int>& A, vector<int>& B);
//Helper function, adds column/row to D for the new cluster (it assumes the new one is the last in clusters)
void addClusterToMatrix(vector<vector<float>>& D, vector<vector<float>>& Dtemp, vector<vector<int>>& clusters, int& t1, int& t2);
//Helper function, given two clusters it merges them in increasing order
vector<int> mergePairClusters(vector<int>& A, vector<int>& B);
//Fixes indices in cluster merge/remove for UPMAG
void fixIndexUPGMA(vector<vector<int>>& clusters, vector<int>&A, vector<int>& B, vector<vector<float>>& D, map<int, int>& transform, int pop);
//Inits map from Dspace to Tspace (all nodes are one entry of D)
map<int, int> initMapSpace(int n);
//Erases key with a given value from map
void eraseKeyByValue(map<int, int>& mapi, int val);
//Erase and translate the map
void eraseAndTranslate(map<int, int>& mapi, int val);
//Finds new cluster distance
float modifiedClusterDistance(vector<int>& c1, vector<int>& c2, vector<vector<float>>& D, int ord1, int ord2);
//Constructs neighbor-joining matrix given a matrix
vector<vector<float>> neighborJoiningMatrix(vector<vector<float>>& D);
//Computes total distance of a node given a matrix
float totalDistance(vector<vector<float>>& D, int pos);
//Finds minimum element not in the diag
void minNonDiag(vector<vector<float>>& D, int& a, int& b);
//Computes delta function given matrix and entries
float deltaFromMatrix(vector<vector<float>>& D, int a, int b);
//Matrix transformatio, helper for neighbor joining algorithm
void transformNeighborMat(vector<vector<float>>& D, int i, int j);
//Map transformation, helper for neighbor joining algorithm
void transformNeighborMap(map<int,int>& mapi, int i, int j, int pop);

weightedGraph::weightedGraph(string name)
{
	vector<vector<int>> nodesHolder;
	int maxNode = 0;

	ifstream file(name + ".txt");
	string temp;
	while (getline(file, temp))
	{
		int tempA, tempB, tempWeight;
		readNodeWeight(temp, tempA, tempB, tempWeight);
		vector<int> hold = { tempA, tempB, tempWeight };
		nodesHolder.push_back(hold);

		maxNode = max(maxNode, max(tempA, tempB));
	}

	adjList = vector<vector<float>>(maxNode + 1, vector<float>(maxNode + 1, -1));
	for (auto& x : nodesHolder)
	{
		int tempA = x[0];
		int tempB = x[1];
		int tempWeight = x[2];

		adjList[tempA][tempB] = tempWeight;
	}

	file.close();
}


weightedGraph::weightedGraph(int n)
{
	adjList = vector<vector<float>>(n, vector<float>(n, -1));
}

void weightedGraph::print()
{
	printMatrix(adjList);
}

vector<vector<float>> weightedGraph::distMat()
{
	vector<int> leaves;
	for (int i = 0; i < adjList.size(); i++)
	{
		if (isNodeLeaf(i))
		{
			leaves.push_back(i);
		}
	}

	int size = leaves.size();
	vector<vector<float>> result = vector<vector<float>>(size, vector<float>(size, 0));
	for (int i = 0; i < size-1; i++)
	{
		for (int j = i + 1; j < size; j++)
		{
			int distance = pathDistanceNodes(i, j);
			result[i][j] = distance;
		}
	}

	//Matrix transpose
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < i; j++)
		{
			result[i][j] = result[j][i];
		}
	}

	return result;
}

void weightedGraph::printBookFormat()
{
	int size = adjList.size();
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			if (adjList[i][j] >= 0)
			{
				cout << i << "->" << j << ":" << adjList[i][j] << endl;
			}
		}
	}
}

void weightedGraph::limbAppender(int v, int n, int A, int B, int DA, int DB, int limb)
{
	int size = adjList.size();
	if (v < size)
	{
		adjList[v][n] = limb;
		adjList[n][v] = limb;
	}
	else
	{
		for (int i = 0; i < size; i++)
		{
			adjList[i].push_back(-1);
		}

		vector<float> newRow(size + 1, -1);
		adjList.push_back(newRow);

		adjList[A][B] = -1;
		adjList[B][A] = -1;
		adjList[A][v] = DA;
		adjList[v][A] = DA;
		adjList[B][v] = DB;
		adjList[v][B] = DB;
		adjList[v][n] = limb;
		adjList[n][v] = limb;
	}
}

weightedGraph additivePhylogeny(vector<vector<float>>& D, int& pop)
{
	int n = D.size();
	if (n == 2)
	{
		vector<vector<float>> g(pop, vector<float>(pop, -1));
		g[0][1] = D[0][1];
		g[1][0] = D[1][0];
		weightedGraph T(g);
		return T;
	}

	int limb = limbLength(D, n - 1);
	for (int j = 0; j < n - 1; j++)
	{
		D[j][n - 1] -= limb;
		D[n - 1][j] = D[j][n - 1];
	}

	int i, k;
	findLeavesAdditive(i, k, n - 1, D);
	int x = D[i][n - 1];
	removeRowColumn(n-1, D);
	weightedGraph T = additivePhylogeny(D, pop);
	int limitA, limitB, limDistA, limDistB;
	int v = T.findNodeAtDistance(i, k, x, pop, limitA, limitB, limDistA, limDistB);
	T.limbAppender(v, n-1, limitA, limitB, limDistA, limDistB, limb);
	return T;
}

weightedGraph UPGMA(vector<vector<float>>& D)
{
	int n = D.size();
	vector<vector<int>> clusters = initCluster(n);
	vector<int> clustDist(n, 1);
	vector<float> ages(n, 0);
	weightedGraph T(n);
	int pop = n;
	float lowDist;
	map<int, int> matrixToTree = initMapSpace(n);
	while (clusters.size() > 1)
	{
		vector<int> lowest = findClosestClusters(clusters, D, lowDist);
		vector<int> C1 = clusters[lowest[0]];
		vector<int> C2 = clusters[lowest[1]];
		int tempOrd1, tempOrd2;
		clusterMerger(clusters, clustDist, lowest, tempOrd1, tempOrd2);
		T.pushNode();
		pop++;
		T.connectNodeToCluster(pop - 1, C1, matrixToTree);
		T.connectNodeToCluster(pop - 1, C2, matrixToTree);
		ages.push_back(float(lowDist) / 2);
		vector<vector<float>> Dtemp = D;
		removeClustersFromMatrix(D, C1, C2);
		addClusterToMatrix(D, Dtemp, clusters, tempOrd1, tempOrd2);
		fixIndexUPGMA(clusters, C1, C2, D, matrixToTree, pop);
	}

	T.setEdgesFromAge(ages);
	return T;
}

vector<vector<int>> UPGMAcluster(vector<vector<float>>& D)
{
	int n = D.size();
	vector<vector<int>> clusters = initCluster(n);
	vector<int> clustDist(n, 1);
	vector<float> ages(n, 0);
	weightedGraph T(n);
	int pop = n;
	float lowDist;
	map<int, int> matrixToTree = initMapSpace(n);
	vector<vector<int>> realClust;
	while (clusters.size() > 1)
	{
		vector<int> lowest = findClosestClusters(clusters, D, lowDist);
		vector<int> C1 = clusters[lowest[0]];
		vector<int> C2 = clusters[lowest[1]];
		int tempOrd1, tempOrd2;
		clusterMerger(clusters, clustDist, lowest, tempOrd1, tempOrd2);
		T.pushNode();
		pop++;
		T.connectNodeToCluster(pop - 1, C1, matrixToTree);
		T.connectNodeToCluster(pop - 1, C2, matrixToTree);
		ages.push_back(float(lowDist) / 2);
		vector<vector<float>> Dtemp = D;
		removeClustersFromMatrix(D, C1, C2);
		addClusterToMatrix(D, Dtemp, clusters, tempOrd1, tempOrd2);
		//Tracks real clusters
		realClust.push_back(clusters.back());
		//
		fixIndexUPGMA(clusters, C1, C2, D, matrixToTree, pop);
	}

	T.setEdgesFromAge(ages);
	return realClust;
}

weightedGraph neighborJoining(vector<vector<float>>& D, int pop, map<int, int> dictio)
{
	map<int, int> mapi = dictio;
	int n = D.size();
	if (pop == -1)
	{
		pop = n;
		mapi = initMapSpace(n);
		dictio = initMapSpace(n);
	}
	if (n == 2)
	{
		weightedGraph T(pop);
		T.connectNodes(pop-2, pop-1, D[0][1]);
		return T;
	}

	vector<vector<float>> dAster = neighborJoiningMatrix(D);
	int i, j;
	minNonDiag(dAster, i, j);
	float delta = deltaFromMatrix(D, i, j);
	float limbI = (D[i][j] + delta) / 2;
	float limbJ = (D[i][j] - delta) / 2;
	transformNeighborMat(D, i, j);
	transformNeighborMap(mapi, i, j, pop);
	pop++;
	weightedGraph T = neighborJoining(D, pop, mapi);
	int m = mapi[mapi.size()-1];
	T.connectNodes(m, dictio[i], limbI);
	T.connectNodes(m, dictio[j], limbJ);
	return T;

}

vector<vector<int>> transformUPGMA(vector<vector<int>>& backs, int n)
{
	vector<vector<int>> clusters = initCluster(n);
	vector<vector<int>> realBacks;
	for (int i = 0; i < backs.size(); i++)
	{
		int a = backs[i][0];
		int b = backs[i][1];
		if (a > b) swap(a, b);
		vector<int> tempA = clusters[a];
		vector<int> tempB = clusters[b];
		clusters.erase(clusters.begin() + b);
		clusters.erase(clusters.begin() + a);
		tempA.insert(tempA.end(), tempB.begin(), tempB.end());
		clusters.push_back(tempA);
		realBacks.push_back(tempA);
	}
	return realBacks;
}

int weightedGraph::findNodeAtDistance(int i, int k, int d, int& pop, int& limA, int& limB, int& limDA, int& limDB)
{
	vector<int> path = pathNodesTree(i, k);
	int temp = 0;
	for (int j = 0; j < path.size() - 1; j++)
	{
		int A = path[j];
		int B = path[j + 1];
		temp += adjList[A][B];
		if (temp == d)
		{
			return B;
		}
		else if (temp > d)
		{
			limA = A;
			limB = B;
			limDA = d - temp + adjList[A][B];
			limDB = temp - d;
			break;
		}
	}

	pop++;
	return pop-1;
}

vector<int> weightedGraph::pathNodesTree(int a, int b)
{
	if (a == b)
	{
		return vector<int>{a};
	}

	vector<int> path;
	path.push_back(a);
	queue<vector<int>> q;
	q.push(path);

	while (!q.empty())
	{
		path = q.front();
		q.pop();

		int lastNode = path.back();
		if (lastNode == b)
		{
			return path;
		}

		vector<int> neighbors = adjacentsToNode(lastNode);
		for (auto& x : neighbors)
		{
			if (!nodeInPath(path, x))
			{
				vector<int> newPath(path.begin(), path.end());
				newPath.push_back(x);
				q.push(newPath);
			}
		}
	}

	vector<int> empty;
	return empty;
}

vector<int> weightedGraph::adjacentsToNode(int node)
{
	vector<int> adj;
	vector<float> nodeRow = adjList[node];
	int size = nodeRow.size();
	for (int i = 0; i < size; i++)
	{
		if (nodeRow[i] >= 0)
		{
			adj.push_back(i);
		}
	}

	return adj;
}

int weightedGraph::pathWeight(vector<int> path)
{
	if (path.empty())
	{
		return 0;
	}

	int weight = 0;
	for (int i = 0; i < path.size() - 1; i++)
	{
		int u = path[i];
		int v = path[i + 1];
		weight += adjList[u][v];
	}

	return weight;
}

vector<int> weightedGraph::inNeighborsNode(int node)
{
	vector<int> inN;
	for (int i = 0; i < adjList.size(); i++)
	{
		if (adjList[i][node] >= 0)
		{
			inN.push_back(i);
		}
	}

	return inN;
}

vector<int> weightedGraph::outNeighborsNode(int node)
{
	vector<int> outN;
	for (int i = 0; i < adjList.size(); i++)
	{
		if (adjList[node][i] >= 0)
		{
			outN.push_back(i);
		}
	}

	return outN;
}

bool weightedGraph::isNodeLeaf(int node)
{
	vector<int> outN = outNeighborsNode(node);
	if (outN.size() == 1)
	{
		return true;
	}
	else
	{
		return false;
	}
}

int weightedGraph::pathDistanceNodes(int a, int b)
{
	vector<int> path = pathNodesTree(a, b);
	int result = pathWeight(path);
	return result;
}

vector<int> weightedGraph::leavesFinder()
{
	vector<int> leaves;
	for (int i = 0; i < adjList.size(); i++)
	{
		if (isNodeLeaf(i))
		{
			leaves.push_back(i);
		}
	}
	return leaves;
}

int weightedGraph::pushNode()
{
	int size = adjList.size();
	for (int i = 0; i < size; i++)
	{
		adjList[i].push_back(-1);
	}

	vector<float> newNode(size + 1, -1);
	adjList.push_back(newNode);

	return size;
}

void weightedGraph::connectNodes(int u, int v, float dist)
{
	adjList[u][v] = dist;
	adjList[v][u] = dist;
}

void weightedGraph::connectNodeToCluster(int u, vector<int>& cluster, map<int, int>& mapi)
{
	for (auto& v : cluster)
	{
		int realV = mapi[v];
		connectNodes(u, realV, -2);
	}
}

void weightedGraph::setEdgesFromAge(vector<float>& ages)
{
	int size = adjList.size();
	for (int i = 0; i < size - 1; i++)
	{
		for (int j = i + 1; j < size; j++)
		{
			if (adjList[i][j] == -2)
			{
				adjList[i][j] = abs(ages[i] - ages[j]);
				adjList[j][i] = abs(ages[i] - ages[j]);
			}
		}
	}
}

void readNodeWeight(string row, int & nodeA, int & nodeB, int & weight)
{
	int pos1 = row.find('-');
	int pos2 = row.find('>');
	int pos3 = row.find(':');

	nodeA = stoi(row.substr(0,pos1));
	nodeB = stoi(row.substr(pos2 + 1, pos3-pos2));
	weight = stoi(row.substr(pos3 + 1));
}


bool nodeInPath(vector<int>& path, int node)
{
	for (auto& x : path)
	{
		if (node == x)
		{
			return true;
		}
	}

	return false;
}

int limbLength(vector<vector<float>>& distMatrix, int j)
{
	int length = 9999999;
	int size = distMatrix.size();
	for (int i = 0; i < size; i++)
	{
		for (int k = 0; k < size; k++)
		{
			if (i != j && k != j)
			{
				int temp = (distMatrix[i][j] + distMatrix[j][k] - distMatrix[i][k]) / 2;
				if (temp < length)
				{
					length = temp;
				}
			}
		}
	}
	return length;
}

vector<vector<float>> loadDistMatrix(string name, int n)
{
	vector<vector<float>> distMatrix = vector<vector<float>>(n, vector<float>(n, 0));
	ifstream file(name + ".txt");
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			file >> distMatrix[i][j];
		}
	}

	return distMatrix;
}


void findLeavesAdditive(int& i, int& k, int n, vector<vector<float>>& D)
{
	int size = D.size();
	for (int a = 0; a < size; a++)
	{
		for (int b = 0; b < size; b++)
		{
			int temp = D[a][b] - D[a][n] - D[n][b];
			if (temp == 0)
			{
				i = a;
				k = b;
				return;
			}
		}
	}
}

void removeRowColumn(int n, vector<vector<float>>& D)
{
	for (int i = 0; i < D.size(); i++)
	{
		D[i].erase(D[i].begin() + n);
	}

	D.erase(D.begin() + n);
}

vector<vector<int>> initCluster(int n)
{
	vector<vector<int>> result;
	for (int i = 0; i < n; i++)
	{
		vector<int> temp = { i };
		result.push_back(temp);
	}
	return result;
}

float clusterDistance(vector<int>& c1, vector<int>& c2, vector<vector<float>>& D)
{
	float S1 = c1.size();
	float S2 = c2.size();
	float temp = 0;
	for (int i = 0; i < S1; i++)
	{
		for (int j = 0; j < S2; j++)
		{
			int a = c1[i];
			int b = c2[j];
			temp += D[a][b];
		}
	}

	float distance = temp / (S1*S2);
	return distance;
}

float clusterDistanceMin(vector<int>& c1, vector<int>& c2, vector<vector<float>>& D)
{
	float S1 = c1.size();
	float S2 = c2.size();
	float distance = 99999;
	for (int i = 0; i < S1; i++)
	{
		for (int j = 0; j < S2; j++)
		{
			int a = c1[i];
			int b = c2[j];
			float temp = D[a][b];
			if (temp < distance)
			{
				distance = temp;
			}
		}
	}

	return distance;
}

vector<int> findClosestClusters(vector<vector<int>>& clusters, vector<vector<float>>& D, float& lowDist)
{
	int a, b;
    float minDist = 999999;
	int size = clusters.size();

	for (int i = 0; i < size - 1; i++)
	{
		vector<int> C1 = clusters[i];
		for (int j = i + 1; j < size; j++)
		{
			vector<int> C2 = clusters[j];
			float temp = clusterDistance(C1, C2, D);
			if (temp < minDist)
			{
				minDist = temp;
				a = i;
				b = j;
			}
		}
	}

	lowDist = minDist;
	vector<int> lowest = { a,b };
	return lowest;
}

void clusterMerger(vector<vector<int>>& clusters, vector<int>& orders, vector<int>& indices, int& t1, int& t2)
{
	int a;
	int b;
	if (indices[0] < indices[1])
	{
		a = indices[0];
		b = indices[1];
	}
	else
	{
		b = indices[0];
		a = indices[1];
	}
	vector<int> newCluster = mergePairClusters(clusters[a], clusters[b]);
	int newOrder = orders[a] + orders[b];
	t1 = orders[a];
	t2 = orders[b];

	clusters.erase(clusters.begin() + b);
	clusters.erase(clusters.begin() + a);
	clusters.push_back(newCluster);

	orders.erase(orders.begin() + b);
	orders.erase(orders.begin() + a);
	orders.push_back(newOrder);
}

void removeNodeFromMatrix(vector<vector<float>>& D, int u)
{
	D.erase(D.begin() + u);
	for (int i = 0; i < D.size(); i++)
	{
		D[i].erase(D[i].begin() + u);
	}
}

void removeClustersFromMatrix(vector<vector<float>>& D, vector<int>& A, vector<int>& B)
{
	vector<int> temp = mergePairClusters(A, B);
	int size = temp.size();
	for (int i = size - 1; i >= 0; i--)
	{
		int x = temp[i];
		removeNodeFromMatrix(D, x);
	}
}

void addClusterToMatrix(vector<vector<float>>& D, vector<vector<float>>& Dtemp, vector<vector<int>>& clusters, int& t1, int& t2)
{
	vector<int> Cnew = clusters.back();
	int size = D.size();
	vector<float> newRow(size + 1, 0);
	D.push_back(newRow);

	for (int i = 0; i < size; i++)
	{
		float dist = modifiedClusterDistance(clusters[i], Cnew, Dtemp, t1, t2);
		D.back()[i] = dist;
		D[i].push_back(dist);
	}
}

vector<int> mergePairClusters(vector<int>& A, vector<int>& B)
{
	vector<int> C(A);
	C.insert(C.end(), B.begin(), B.end());
	sort(C.begin(), C.end());
	return C;
}

void fixIndexUPGMA(vector<vector<int>>& clusters, vector<int>& A, vector<int>& B, vector<vector<float>>& D, map<int,int>& transform, int pop)
{
	vector<int> temp = mergePairClusters(A, B);
	int sizeClus = clusters.size();
	int sizeTree = D.size();
	for (int i = 0; i < sizeClus - 1; i++)
	{
		for (auto& m : clusters[i])
		{
			int tempM = m;
			for (auto& c : temp)
			{
				if (tempM > c)
				{
					m--;
				}
				else
				{
					break;
				}
			}
		}
	}

	clusters.back() = { sizeTree - 1 };

	//Fixes map
	for (int i=temp.size()-1; i>=0; i--)
	{
		eraseAndTranslate(transform, temp[i]);
	}

	transform[transform.size()] = pop - 1;
}

map<int, int> initMapSpace(int n)
{
	map<int, int> result;
	for (int i = 0; i < n; i++)
	{
		result.insert({ i,i });
	}
	return result;
}

void eraseKeyByValue(map<int, int>& mapi, int val)
{
	int pos;
	for (auto& x : mapi)
	{
		if (x.second == val)
		{
			pos = x.first;
		}
	}

	mapi.erase(pos);
}

void eraseAndTranslate(map<int, int>& mapi, int val)
{
	int last = mapi.size() - 1;
	for (int i = val + 1; i <= last; i++)
	{
		mapi[i - 1] = mapi[i];
	}

	mapi.erase(mapi.size()-1);
}

float modifiedClusterDistance(vector<int>& c1, vector<int>& c2, vector<vector<float>>& D, int ord1, int ord2)
{
	int m = c1[0];
	int i = c2[0];
	int j = c2[1];
	float distance = (ord1 * D[m][i] + ord2 * D[m][j]) / (ord1 + ord2);
	return distance;
}

vector<vector<float>> neighborJoiningMatrix(vector<vector<float>>& D)
{
	vector<vector<float>> result(D);
	int size = D.size();
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			if (i != j)
			{
				result[i][j] = (size - 2)*D[i][j] - totalDistance(D, i) - totalDistance(D, j);
			}
		}
	}

	return result;
}

float totalDistance(vector<vector<float>>& D, int pos)
{
	float dist = 0;
	for (int i = 0; i < D.size(); i++)
	{
		dist += D[i][pos];
	}

	return dist;
}

void minNonDiag(vector<vector<float>>& D, int& a, int&b)
{
	float min = 99999;
	int size = D.size();
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			if (i != j)
			{
				float temp = D[i][j];
				if (temp < min)
				{
					min = temp;
					a = i;
					b = j;
				}
			}
		}
	}
}

float deltaFromMatrix(vector<vector<float>>& D, int a, int b)
{
	float result = (totalDistance(D, a) - totalDistance(D, b)) / (D.size() - 2);
	return result;
}

void transformNeighborMat(vector<vector<float>>& D, int i, int j)
{
	//Adds column/row
	int size = D.size();
	vector<float> newRow(size + 1,0);
	for (int k = 0; k < size; k++)
	{
		float temp = (D[k][i] + D[k][j] - D[i][j]) / 2;
		newRow[k] = temp;
		D[k].push_back(temp);
	}
	D.push_back(newRow);

	//Now deletes rows columns i,j
	if (i > j)
	{
		swap(i, j);
	}

	D.erase(D.begin() + j);
	D.erase(D.begin() + i);
	for (int k = 0; k < D.size(); k++)
	{
		D[k].erase(D[k].begin() + j);
		D[k].erase(D[k].begin() + i);
	}
}

void transformNeighborMap(map<int, int>& mapi, int i, int j, int pop)
{

	//Makes sure i<j
	if (i > j)
	{
		swap(i, j);
	}
	//Fixes map
	eraseAndTranslate(mapi, j);
	eraseAndTranslate(mapi, i);
	mapi[mapi.size()] = pop;
}
