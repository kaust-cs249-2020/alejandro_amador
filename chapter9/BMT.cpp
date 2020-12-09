#include "BMT.h"

//Helper function, finds if char in a map
bool isCharHere(map<char, int>& mapi, char a);
//Checks if two strings have a d-distance below bound
bool dNorm(string& a, string& b, int d);

string transformBW(string& text)
{
	int size = text.size();

	//Gets all cyclic rotations
	vector<string> rotations;
	rotations.push_back(text);
	string last = text;
	string next;
	for (int i = 1; i < size; i++)
	{
		next = last;
		next.pop_back();
		next.insert(next.begin(), text.at(size - i));
		rotations.push_back(next);
		last = next;
	}

	//Orders them to get M
	sort(rotations.begin(), rotations.end(), lexiCompare);

	//Gets last column
	string result;
	for (auto& x : rotations)
	{
		result.push_back(x.back());
	}

	return result;
}

string inverseBW(string & text)
{
	int size = text.size();
	//Creates first column
	vector<pair<char, int>> first;
	string temp = text;
	sort(temp.begin(), temp.end());
	char lastChar = temp.at(0);
	first.push_back({ lastChar, 1 });
	for (int i = 1; i < size; i++)
	{
		char newChar = temp.at(i);
		if (newChar == lastChar)
		{
			int lastNum = first[i - 1].second;
			first.push_back({ newChar, lastNum + 1 });
		}
		else
		{
			first.push_back({ newChar, 1 });
		}

		lastChar = newChar;
	}

	//Creates last column
	vector<pair<char, int>> last;
	map<char, int> currentNum;
	for (int i = 0; i < size; i++)
	{
		char tempChar = text.at(i);
		if (isCharHere(currentNum, tempChar))
		{
			currentNum[tempChar] += 1;
			last.push_back({ tempChar, currentNum[tempChar] });
		}
		else
		{
			currentNum[tempChar] = 1;
			last.push_back({ tempChar, currentNum[tempChar] });
		}
	}

	//Actual algorithm
	pair<char, int> stopPair = last[0];
	pair<char, int> nextPair = first[0];
	string result = "";
	while (nextPair != stopPair)
	{
		vector<pair<char, int>>::iterator it = find(last.begin(), last.end(), nextPair);
		int pos = it - last.begin();
		nextPair = first[pos];
		result.push_back(nextPair.first);
	}
	result += "$";

	return result;
}

bool isCharHere(map<char, int>& mapi, char a)
{
	map<char, int>::iterator it = mapi.find(a);
	if (it == mapi.end())
	{
		return false;
	}
	else
	{
		return true;
	}
}

BWT::BWT(string text)
{
	int size = text.size();
	//Creates first column
	string temp = text;
	sort(temp.begin(), temp.end());
	char lastChar = temp.at(0);
	first.push_back({ lastChar, 1 });
	for (int i = 1; i < size; i++)
	{
		char newChar = temp.at(i);
		if (newChar == lastChar)
		{
			int lastNum = first[i - 1].second;
			first.push_back({ newChar, lastNum + 1 });
		}
		else
		{
			first.push_back({ newChar, 1 });
		}

		lastChar = newChar;
	}

	//Creates last column
	map<char, int> currentNum;
	for (int i = 0; i < size; i++)
	{
		char tempChar = text.at(i);
		if (isCharHere(currentNum, tempChar))
		{
			currentNum[tempChar] += 1;
			last.push_back({ tempChar, currentNum[tempChar] });
		}
		else
		{
			currentNum[tempChar] = 1;
			last.push_back({ tempChar, currentNum[tempChar] });
		}
	}

	//Creates firsToLast map
	for (int i = 0; i < size; i++)
	{
		pair<char, int> tempo = last[i];
		vector<pair<char, int>>::iterator it = find(first.begin(), first.end(), tempo);
		int pos = it - first.begin();
		lastToFirst[i] = pos;
	}

	//Creates count matrix
	map<char, int> firsto;
	char tempChar;
	for (int i = 0; i < size; i++)
	{
		tempChar = text.at(i);
		firsto[tempChar] = 0;
	}
	count.push_back(firsto);

	for (int m = 1; m <= last.size(); m++)
	{
		tempChar = last[m - 1].first;
		firsto[tempChar] += 1;
		count.push_back(firsto);
	}


	//Creates first appereance
	char lastoChar = '&';
	for (int k = 0; k < first.size(); k++)
	{
		char tempoChar = first[k].first;
		if (tempoChar != lastoChar)
		{
			firstOcurrence[tempoChar] = k;
			lastoChar = tempoChar;
		}
	}

	sentence = inverseBW(text);
	sufijo = suffixArray(sentence);


}

int BWT::BWMatching(string pattern)
{
	int top = 0;
	int bottom = last.size() - 1;
	int topIndx, botIndx;
	while (top <= bottom)
	{
		if (!pattern.empty())
		{
			char symbol = pattern.back();
			pattern.pop_back();
			if (checkSymbolHere(top, bottom, topIndx, botIndx, symbol))
			{
				top = lastToFirst[topIndx];
				bottom = lastToFirst[botIndx];
			}
			else
			{
				return 0;
			}
		}
		else
		{
			return bottom - top + 1;
		}
	}
}

vector<int> BWT::multipleBWMatching(vector<string>& patterns)
{
	vector<int> list;
	for (auto& s : patterns)
	{
		int tempo = BWMatching(s);
		list.push_back(tempo);
	}

	return list;
}

int BWT::betterBWMatching(string pattern)
{
	int top = 0;
	int bottom = last.size() - 1;
	while (top <= bottom)
	{
		if (!pattern.empty())
		{
			char symbol = pattern.back();
			pattern.pop_back();
			if (checkSymbolHereNoPos(top, bottom, symbol))
			{
				top = firstOcurrence[symbol] + count[top][symbol];
				bottom = firstOcurrence[symbol] + count[bottom + 1][symbol] - 1;
			}
			else
			{
				return 0;
			}
		}
		else
		{
			return bottom - top + 1;
		}
	}
}

vector<int> BWT::multipleBetterBWM(vector<string>& patterns)
{
	vector<int> list;
	for (auto& s : patterns)
	{
		int tempo = betterBWMatching(s);
		list.push_back(tempo);
	}

	return list;
}

vector<int> BWT::betterBWMPositions(string& pattern)
{
	vector<int> lexOrder = sufijo.getLexOrder();

	int top = 0;
	int bottom = last.size() - 1;
	while (top <= bottom)
	{
		if (!pattern.empty())
		{
			char symbol = pattern.back();
			pattern.pop_back();
			if (checkSymbolHereNoPos(top, bottom, symbol))
			{
				top = firstOcurrence[symbol] + count[top][symbol];
				bottom = firstOcurrence[symbol] + count[bottom + 1][symbol] - 1;
			}
			else
			{
				return {};
			}
		}
		else
		{
			vector<int> positions;
			for (int i = top; i <= bottom; i++)
			{
				positions.push_back(lexOrder[i]);
			}

			return positions;
		}
	}
}

vector<int> BWT::multipleBWMPositions(vector<string>& patterns)
{
	vector<int> result;
	for (auto& x : patterns)
	{
		vector<int> xPos = betterBWMPositions(x);
		result.insert(result.end(), xPos.begin(), xPos.end());
	}

	return result;
}

vector<int> BWT::naiveMismatch(string & pattern, int d)
{
	vector<int> result;

	//First makes seeds
	vector<pair<string,int>> seeds;
	int n = pattern.size();
	int k = floor(float(n) / (d + 1));
	for (int i = 0; i < d + 1; i++)
	{
		int index = i * k;
		string temp = pattern.substr(index, k);
		seeds.push_back({ temp,index });
	}

	//Seed comparison
	for (auto& s : seeds)
	{
		vector<int> matches = betterBWMPositions(s.first);
		for (auto& m : matches)
		{
			int indexSub = m - s.second;
			string subTemp = sentence.substr(indexSub, n);
			if (dNorm(pattern, subTemp, d))
			{
				result.push_back(indexSub);
			}
		}
	}

	//Sort and eliminate duplicates
	sort(result.begin(), result.end());
	result.erase(unique(result.begin(), result.end()), result.end());

	return result;
}

vector<int> BWT::multepleNaiveMismatch(vector<string>& patterns, int d)
{
	vector<int> result;

	int count = patterns.size();
	for (auto& s : patterns)
	{
		cout << count << endl;
		count--;
		vector<int> temp = naiveMismatch(s, d);
		result.insert(result.end(), temp.begin(), temp.end());
	}

	//Sort 
	sort(result.begin(), result.end());

	return result;
}

bool BWT::checkSymbolHere(int top, int bot, int & topIndex, int & botIndex, char & symbol)
{
	bool result = false;
	topIndex = -1;
	for (int i = top; i <= bot; i++)
	{
		if (last[i].first == symbol)
		{
			result = true;
			if (topIndex == -1)
			{
				topIndex = i;
				botIndex = i;
			}
			else
			{
				botIndex = i;
			}
		}
	}

	return result;
}

bool BWT::checkSymbolHereNoPos(int top, int bot, char & symbol)
{
	bool result = false;
	for (int i = top; i <= bot; i++)
	{
		if (last[i].first == symbol)
		{
			return true;
		}
	}

	return result;
}

bool dNorm(string& a, string& b, int d)
{
	int localDist = 0;
	for (int i = 0; i < a.size(); i++)
	{
		if (a.at(i) != b.at(i))
		{
			localDist++;
			if (localDist > d)
			{
				return false;
			}
		}
	}

	return true;
}
