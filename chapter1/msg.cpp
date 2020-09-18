#include "msg.h"
#include "agrozUtil.h"

int maxMap(map<string, int> userMap);
char nucleotideConversor(char nucleotide);
map<string, int> frequencyTable(int k, string subText);
int skewLocal(char nucleotide, int previous);
int hammingDistance(string st1, string st2);
vector<string> neighbors(string pattern, int d);
string suffix(string pattern);
string reverseComplement(string text);

msg::msg(string filename)
{
	ifstream file(filename);
	string temp;
	while (getline(file, temp)) 
	{
		m_inputString = m_inputString + temp;
	}
}

void msg::print()
{
	cout << m_inputString << endl;
}

int msg::findFrequency(string pattern)
{
	int appears = 0;
	int sizeFull = m_inputString.size();
	int sizeMsg = pattern.size();
	for (int i = 0; i <= (sizeFull - sizeMsg); i++)
	{
		string subString = m_inputString.substr(i,sizeMsg);
		if (subString == pattern)
		{
			appears++;
		}
	}

	return appears;
}

int msg::findMismatchFrequency(string pattern, int distance)
{
	int appears = 0;
	int sizeFull = m_inputString.size();
	int sizeMsg = pattern.size();
	for (int i = 0; i <= (sizeFull - sizeMsg); i++)
	{
		string subString = m_inputString.substr(i, sizeMsg);
		int hamming = hammingDistance(subString, pattern);
		if (hamming <= distance)
		{
			appears++;
		}
	}

	return appears;
}

vector<string> msg::frequentWords(int k)
{
	vector<string> frequentKmers;
	vector<int> count;
	for (int i = 0; i < m_inputString.size() - k ; i++)
	{
		string pattern = m_inputString.substr(i, k);
		count.push_back(findFrequency(pattern));
	}
	int maxCount = *max_element(count.begin(), count.end());
	for (int i = 0; i < m_inputString.size() - k; i++)
	{
		if (count[i] == maxCount)
		{
			string maybe = m_inputString.substr(i, k);
			vector<string>::iterator it;
			it = find(frequentKmers.begin(), frequentKmers.end(), maybe);
			if (it == frequentKmers.end())
				frequentKmers.push_back(maybe);
		}
	}

	return frequentKmers;
}

vector<string> msg::frequentWordsWithMismatch(int k, int d)
{
	vector<string> patterns;
	map<string, int> freqMap;
	int n = m_inputString.size();
	for (int i = 0; i <= (n - k); i++)
	{
		string maybe = m_inputString.substr(i, k);
		vector<string> vecinos = neighbors(maybe, d);
		vector<string>::iterator it;
		for (it = vecinos.begin(); it != vecinos.end(); it++)
		{
			if (freqMap.find(*it) == freqMap.end())
			{
				freqMap.insert(make_pair(*it, 1));
			}
			else
			{
				freqMap[*it]++;
			}
		}
	}
	int maxValue = maxMap(freqMap);
	map<string, int>::iterator it;
	for (it = freqMap.begin(); it != freqMap.end(); it++)
	{
		if (it->second == maxValue)
		{
			patterns.push_back(it->first);
		}
	}

	return patterns;
}

vector<string> msg::betterFrequentWords(int k)
{
	vector<string> frequentKmers;
	map<string, int> freqMap = frequencyTable(k, m_inputString);
	int maxValue = maxMap(freqMap);
	map<string, int>::iterator it;
	for (it = freqMap.begin(); it != freqMap.end(); it++)
	{
		if (it->second == maxValue)
		{
			frequentKmers.push_back(it->first);
		}
	}
	return frequentKmers;
}

vector<string> msg::frequentWordsMRC(int k, int d)
{
	vector<string> patterns;
	map<string, int> freqMap;
	int n = m_inputString.size();
	for (int i = 0; i <= (n - k); i++)
	{
		cout << i - (n-k) << endl;
		string maybe = m_inputString.substr(i, k);
		string maybeReverse = reverseComplement(maybe);
		vector<string> vecinos = neighbors(maybe, d);
		vector<string> vecinosReverse = neighbors(maybeReverse, d);
		vector<string>::iterator it;
		for (it = vecinos.begin(); it != vecinos.end(); it++)
		{
			if (freqMap.find(*it) == freqMap.end())
			{
				freqMap.insert(make_pair(*it, 1));
			}
			else
			{
				freqMap[*it]++;
			}
		}
		for (it = vecinosReverse.begin(); it != vecinosReverse.end(); it++)
		{
			if (freqMap.find(*it) == freqMap.end())
			{
				freqMap.insert(make_pair(*it, 1));
			}
			else
			{
				freqMap[*it]++;
			}
		}
	}
	int maxValue = maxMap(freqMap);
	map<string, int>::iterator it;
	for (it = freqMap.begin(); it != freqMap.end(); it++)
	{
		if (it->second == maxValue)
		{
			patterns.push_back(it->first);
		}
	}

	return patterns;
}

map<string, int> frequencyTable(int k, string subText)
{
	map<string, int> freqMap;
	int n = subText.size();
	for (int i = 0; i <= n - k ; i++)
	{
		string maybe = subText.substr(i, k);
		if (freqMap.find(maybe) == freqMap.end())
		{
			freqMap.insert(make_pair(maybe, 1));
		}
		else
		{
			freqMap[maybe]++;
		}
	}

	return freqMap;
}

string reverseComplement(string text)
{
	string reverseComp = "";
	string::iterator it;
	for (it = text.begin(); it != text.end(); it++)
	{
		reverseComp = nucleotideConversor(*it) + reverseComp;
	}
	return reverseComp;
}

vector<int> msg::patternMatching(string pattern)
{
	vector<int> positions;
	int sizeFull = m_inputString.size();
	int sizeMsg = pattern.size();
	for (int i = 0; i <= (sizeFull - sizeMsg); i++)
	{
		string subString = m_inputString.substr(i, sizeMsg);
		if (subString == pattern)
		{
			positions.push_back(i);
		}
	}

	return positions;
}

vector<int> msg::approxPatternMatching(string pattern, int distance)
{
	vector<int> positions;
	int sizeFull = m_inputString.size();
	int sizeMsg = pattern.size();
	for (int i = 0; i <= (sizeFull - sizeMsg); i++)
	{
		string maybe = m_inputString.substr(i, sizeMsg);
		int hamming = hammingDistance(pattern, maybe);
		if (hamming <= distance)
		{
			positions.push_back(i);
		}
	}

	return positions;
}

vector<string> msg::findClumps(int k, int L, int t)
{
	vector<string> patterns;
	int n = m_inputString.size();
	for (int i = 0; i <= n - L; i++)
	{
		cout << i << endl;
		string window = m_inputString.substr(i,L);
		map<string, int> freqMap = frequencyTable(k, window);
		map<string, int>::iterator it;
		for (it = freqMap.begin(); it != freqMap.end(); it++)
		{
			if (it->second >= t)
			{
				vector<string>::iterator it2;
				it2 = find(patterns.begin(), patterns.end(), it->first);
				if (it2 == patterns.end())
				{
					patterns.push_back(it->first);
				}
			}
		}
	}
	return patterns;
}

vector<int> msg::frequencyPlot(string pattern, int d)
{
	int sizeFull = m_inputString.size();
	int sizeMsg = pattern.size();
	vector<int> counter(sizeFull,0);
	string patternReverse = reverseComplement(pattern);
	vector<string> vecinos = neighbors(pattern, d);
	vector<string> vecinosReverse = neighbors(patternReverse, d);
	vector<string> vecinosAll(vecinos);
	vecinosAll.insert(vecinosAll.end(),vecinosReverse.begin(), vecinosReverse.end());
	vector<string>::iterator it;

	for (int i = 0; i <= (sizeFull - sizeMsg); i++)
	{
		string maybe = m_inputString.substr(i, sizeMsg);
		for (it = vecinosAll.begin(); it != vecinosAll.end(); it++)
		{
			if (maybe == *it)
			{
				counter[i] += 2;
				break;
			}
		}
	}

	return counter;
}

vector<int> msg::skewDiagram()
{
	vector<int> values = {0};
	int last = 0;
	string::iterator it;
	for (it=m_inputString.begin(); it != m_inputString.end(); it++)
	{
		int newValue = skewLocal(*it, last);
		values.push_back(newValue);
		last = newValue;
	}
	return values;
}

vector<int> msg::skewMinPos()
{
	vector<int> positions;
	vector<int> skew = skewDiagram();
	int minValue = *min_element(skew.begin(), skew.end());
	int n = skew.size();
	for (int i = 0; i < n; i++)
	{
		if (skew[i] == minValue)
			positions.push_back(i);
	}
	return positions;
}

msg::~msg()
{
}

int maxMap(map<string, int> userMap)
{
	map<string, int>::iterator it = userMap.begin();
	int maxValue = it->second;
	it++;
	for(it; it!=userMap.end(); it++)
	{
		int maybeValue = it->second;
		if (maybeValue > maxValue)
		{
			maxValue = maybeValue;
		}
	}
	return maxValue;
}

char nucleotideConversor(char nucleotide)
{
	char converted;
	if (nucleotide == 'A')
		converted = 'T';
	else if (nucleotide == 'G')
		converted = 'C';
	else if (nucleotide == 'T')
		converted = 'A';
	else if (nucleotide == 'C')
		converted = 'G';
	else
	{
		cout << nucleotide << endl;
		cout << "You have some weird DNA" << endl;
		exit(-1);
	}
	return converted;
}

int skewLocal(char nucleotide, int previous)
{
	int value;
	switch (nucleotide)
	{
	case 'G':
		value = previous + 1;
		break;
	case 'C':
		value = previous - 1;
		break;
	default:
		value = previous;
		break;
	}
	return value;
}

int hammingDistance(string st1, string st2)
{
	int distance = 0;
	int n = st1.size();
	for (int i = 0; i < n; i++)
	{
		if (st1.at(i) != st2.at(i))
			distance++;
	}
	return distance;
}

vector<string> neighbors(string pattern, int d)
{
	vector<string> neighborhood;
	if (d == 0)
	{
		neighborhood = { pattern };
		return neighborhood;
	}
	if (pattern.size() == 1)
	{
		neighborhood = {"A", "C", "G", "T"};
		return neighborhood;
	}
	vector<string> suffixNeighbors;
	vector<string>::iterator text;
	suffixNeighbors = neighbors(suffix(pattern), d);
	for (text = suffixNeighbors.begin(); text != suffixNeighbors.end(); text++)
	{
		if (hammingDistance(suffix(pattern), *text) < d)
		{
			neighborhood.push_back("A" + *text);
			neighborhood.push_back("C" + *text);
			neighborhood.push_back("G" + *text);
			neighborhood.push_back("T" + *text);
		}

		else
		{
			neighborhood.push_back(pattern[0] + *text);
		}
	}

	return neighborhood;
}

string suffix(string pattern)
{
	string subPattern = pattern.substr(1);
	return subPattern;
}
