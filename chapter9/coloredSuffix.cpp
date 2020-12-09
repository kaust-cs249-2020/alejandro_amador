#include "coloredSuffix.h"

//Given a node's list of outneighbors it shows wether a given node is a neighbor, returns its position on the outArray
int isNodeNeighbors(vector<pair<pair<int, int>, int>>& vecinos, int a);
//Helper function, checks whether a trie node is already on the tree
int isNodeOnTrees(map<int, int>& mapi, int u);
//Gets ripe color based on children
nodeColor colorOfMyChildren(vector<int>& children, map<int, nodeColor>& canvas);
//Converts color to string for printing
string colorToStr(nodeColor col);

coloredSuffix::coloredSuffix(string text1, string text2)
{
	string text = text1 + "#" + text2 + "$";

	//Saves sentences for future reference
	sentence = text;
	sentence1 = text1;
	sentence2 = text2;

	//Temporal map for handling leaves
	// tree - trie
	map<int, int> mapi;
	mapi[0] = 0;

	//Book's algorithm
	modiTrie trie(text);
	vector<vector<int>> nonBran = trie.maximalNonBranchingPaths();
	addNode();

	for (auto& p : nonBran)
	{
		int position = trie.getEdgePos(p[0], p[1]);
		int length = p.size() - 1;

		int a = isNodeOnTrees(mapi, p[0]);
		int b = isNodeOnTrees(mapi, p.back());
		pair<int, int> edge = { position, length };

		if (a != -1)
		{
			if (b != -1)
			{
				connectNodes(a, b, edge);
			}
			else
			{
				int newNode = addNode();
				connectNodes(a, newNode, edge);
				mapi[newNode] = p.back();
			}
		}
		else
		{
			if (b != -1)
			{
				int newA = addNode();
				connectNodes(newA, b, edge);
				mapi[newA] = p[0];

			}
			else
			{
				int newA = addNode();
				int newB = addNode();
				connectNodes(newA, newB, edge);
				mapi[newA] = p[0];
				mapi[newB] = p.back();
			}
		}
	}

	//Adds leaves stuff
	map<int, int> leafWas = trie.getLeafLabel();

	int size = adjList.size();
	for (int i = 0; i < size; i++)
	{
		if (isNodeLeaf(i))
		{
			int oldPos = mapi[i];
			leafLabel[i] = leafWas[oldPos];
		}
	}


	giveMeColor();

}

coloredSuffix::coloredSuffix(string name)
{
	ifstream file(name + ".txt");
	string temp;
	vector<vector<int>> adjTemp;
	while (getline(file, temp))
	{
		if (temp.size() != 1)
		{
			//Finds if information is color or connection
			int type = temp.find(':');
			if (type == string::npos)
			{
				int newNode = addNode();
				vector<int> tempNeighbors;
				//Finds if empty neighbors
				if (temp.find('{') != string::npos)
				{
					adjTemp.push_back({});
				}
				else
				{
					int limA = temp.find('>');
					int limB = temp.find(',');
					int neighbor = stoi(temp.substr(limA + 2, limB - limA - 2));
					tempNeighbors.push_back(neighbor);
					limA = limB + 1;
					limB = temp.find(',', limA);
					while (limB != string::npos)
					{
						neighbor = stoi(temp.substr(limA, limB - limA));
						tempNeighbors.push_back(neighbor);
						limA = limB + 1;
						limB = temp.find(',', limA);
					}

					neighbor = stoi(temp.substr(limA));
					tempNeighbors.push_back(neighbor);
					adjTemp.push_back(tempNeighbors);
				}
			}
			else
			{
				int node = stoi(temp.substr(0, type));
				nodeColor nodeCol;
				string colStr = temp.substr(type + 2);
				if (colStr == "blue")
				{
					nodeCol = blue;
				}
				else
				{
					nodeCol = red;
				}
				nodeCanvas[node] = nodeCol;
			}
		}

	}
	file.close();

	pair<int, int> dummy = { 0,0 };
	for (int i = 0; i < adjTemp.size(); i++)
	{
		for (auto& v : adjTemp[i])
		{
			connectNodes(i, v, dummy);
		}
	}

	//Colors
	vector<int> ripes = findRipeNodes();
	while (!ripes.empty())
	{
		for (auto& r : ripes)
		{
			vector<int> children = getOutNeighbors(r);
			nodeCanvas[r] = colorOfMyChildren(children, nodeCanvas);
		}

		ripes = findRipeNodes();
	}
}

void coloredSuffix::printColors()
{
	for (auto& x : nodeCanvas)
	{
		cout << x.first << ": " << colorToStr(x.second) << endl;
	}
}

void coloredSuffix::giveMeColor()
{
	//First color the leaves of color and everything else gray
	int size = adjList.size();
	for (int i = 0; i < size; i++)
	{
		if (isNodeLeaf(i))
		{
			int leaf = i;
			int pos = adjListBack[leaf][0].first.first;
			int length = adjListBack[leaf][0].first.second;
			string leafEdge = sentence.substr(pos, length);
			nodeColor temp = findMyColor(leafEdge);
			nodeCanvas[leaf] = temp;
		}
		else
		{
			nodeCanvas[i] = gray;
		}
	}

	//Now does the algorithm for coloring everything else
	vector<int> ripes = findRipeNodes();
	while (!ripes.empty())
	{
		for (auto& r : ripes)
		{
			vector<int> children = getOutNeighbors(r);
			nodeCanvas[r] = colorOfMyChildren(children, nodeCanvas);
		}

		ripes = findRipeNodes();
	}
}

string coloredSuffix::longestSharedSubstringProblem()
{
	//Initialize visit
	map<int, bool> visited;
	vector<string> list;
	vector<nodeColor> listCol;
	string soFar = "";

	for (int i = 0; i < adjList.size(); i++)
	{
		visited[i] = false;
	}

	//Calls DFS
	stringsColorDFS(0, visited, list, listCol, soFar);

	//Looks for longest purple
	string king;
	int kingSize = 0;
	for (int i=0; i<list.size(); i++)
	{
		string x = list[i];
		nodeColor xCol = listCol[i];
		int xSize = x.size();
		if (xSize > kingSize && xCol==purple)
		{
			king = x;
			kingSize = xSize;
		}
	}

	return king;
}

string coloredSuffix::shortestNonSharedSubstringProblem()
{
	//Initialize visit
	map<int, bool> visited;
	vector<string> list;
	string soFar = "";

	for (int i = 0; i < adjList.size(); i++)
	{
		visited[i] = false;
	}

	//Calls DFS
	stringsBlueDFS(0, visited, list, soFar);

	//Looks for the shortest
	string king;
	int kingSize = 99999;
	for (auto& x : list)
	{
		int xSize = x.size();
		if (xSize < kingSize)
		{
			king = x;
			kingSize = xSize;
		}
	}

	return king;
}

void coloredSuffix::stringsColorDFS(int x, map<int, bool>& visited, vector<string>& list, vector<nodeColor>& colList, string & soFar)
{
	visited[x] = true;
	vector<int> outNeighbors = getOutNeighbors(x);

	//If all neighbors are leaves then appends the string and stops
	if (areAllLeaves(x))
	{
		list.push_back(soFar);
		colList.push_back(nodeCanvas[x]);
	}
	//Otherwise continue adventure
	else
	{
		for (auto& w : outNeighbors)
		{
			if (visited[w] == false)
			{
				//Only advance if node is not leaf
				if (!isNodeLeaf(w))
				{
					int pos = isNodeNeighbors(adjList[x], w);
					int begin = adjList[x][pos].first.first;
					int length = adjList[x][pos].first.second;
					string add = soFar + sentence.substr(begin, length);
					stringsColorDFS(w, visited, list, colList, add);
				}
			}
		}
	}
}

void coloredSuffix::stringsBlueDFS(int x, map<int, bool>& visited, vector<string>& list, string & soFar)
{
	visited[x] = true;
	vector<int> outNeighbors = getOutNeighbors(x);

	//If x is a blue leaf then just appends substring of text1
	if (isNodeLeaf(x))
	{
		if (nodeCanvas[x] == blue)
		{
			int pos = soFar.find('#');
			list.push_back(soFar.substr(0, pos));
		}
	}
	//If not a leaf then continues adventure
	else
	{
		for (auto& w : outNeighbors)
		{
			if (visited[w] == false)
			{
				int pos = isNodeNeighbors(adjList[x], w);
				int begin = adjList[x][pos].first.first;
				int length = adjList[x][pos].first.second;
				string add = soFar + sentence.substr(begin, length);
				stringsBlueDFS(w, visited, list, add);
			}
		}
	}
}

nodeColor coloredSuffix::findMyColor(string x)
{
	nodeColor result = red;
	int size1 = sentence1.size();
	int localSize = x.size();

	for (int i = 0; i < size1; i++)
	{
		string temp = sentence.substr(i, localSize);
		if (temp == x)
		{
			result = blue;
			break;
		}
	}

	return result;
}

vector<int> coloredSuffix::findRipeNodes()
{
	vector<int> ripes;
	for (int i = 0; i < adjList.size(); i++)
	{
		if (nodeCanvas[i] == gray)
		{
			vector<int> children = getOutNeighbors(i);
			bool isRipe = true;
			for (auto& c : children)
			{
				if (nodeCanvas[c] == gray)
				{
					isRipe = false;
					break;
				}
			}

			if (isRipe)
			{
				ripes.push_back(i);
			}
		}
	}

	return ripes;
}

int isNodeNeighbors(vector<pair<pair<int, int>, int>>& vecinos, int a)
{
	int pos = -1;
	int index = 0;
	for (auto& x : vecinos)
	{
		if (x.second == a)
		{
			pos = index;
			break;
		}

		index++;
	}

	return pos;
}

int isNodeOnTrees(map<int, int>& mapi, int u)
{
	int pos = -1;
	int index = 0;
	for (auto& x : mapi)
	{
		if (x.second == u)
		{
			pos = x.first;
			break;
		}

		index++;
	}

	return pos;
}

nodeColor colorOfMyChildren(vector<int>& children, map<int, nodeColor>& canvas)
{
	nodeColor col = canvas[children[0]];
	for (int i = 1; i < children.size(); i++)
	{
		nodeColor temp = canvas[children[i]];
		if (temp != col)
		{
			col = purple;
			break;
		}
	}

	return col;
}

string colorToStr(nodeColor col)
{
	switch (col)
	{
	case gray:
	{
		return "gray";
	}
	case blue:
	{
		return "blue";
	}
	case red:
	{
		return "red";
	}
	case purple:
	{
		return "purple";
	}
	}
}