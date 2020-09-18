#include "motif.h"

vector<string> neighbors(string pattern, int d);
string suffix(string pattern);
int hammingDistance(string st1, string st2);
void removeDuplicates(vector<string>& vec);
vector<string> allStrings(int k);
float probGivenProfile(string kmer, vector<vector<float>> profile);
int nucleoToInt(char nucleo);
vector<vector<float>> formProfile(vector<string> motifs);
vector<vector<float>> formProfilePartial(vector<string> motifs, int except);
void profileRow(string local, vector<vector<float>>& profile, float& weight);
int motifScore(vector<string> motif);
vector<vector<int>> formCount(vector<string> motifs);
void countRow(string local, vector<vector<int>>& counts);
int randomWeight(vector<float> weight);
float randomFloat(float b);
string consent(vector<vector<int>> counts);

motif::motif(string filename)
{
	ifstream file(filename);
	string temp;
	while (getline(file, temp))
	{
		m_dna.push_back(temp);
	}
	file.close();
}

string motif::consensus(int k, algorithm option, int runs)
{
	string result;
	vector<string> motifs;
	if (option==RANDOMIZED_MOTIF_SEARCH)
	{
		motifs = randomizedMotifSearch(k, runs);
	}
	else if (option == GIBBS_SAMPLING)
	{
		motifs = gibbsSampler(k, 20, runs);
	}
	else
	{
		result = medianString(k);
		return result;
	}

	result = consent(formCount(motifs));
	return result;
}

vector<string> motif::greedyMotifSearch(int k)
{
	vector<string> bestMotifs = initialBestMotifs(k);
	int t = m_dna.size();
	int numberKmers = m_dna[0].size() - k;
	for (int i = 0; i <= numberKmers; i++)
	{
		vector<string> motifs;
		string temp = m_dna[0].substr(i, k);
		motifs.push_back(temp);
		for (int j = 1; j < t; j++)
		{
			vector<vector<float>> profile = formProfile(motifs);
			string newMotif = profileMostProbable(profile, j);
			motifs.push_back(newMotif);
		}

		if (motifScore(motifs) < motifScore(bestMotifs))
		{
			bestMotifs = motifs;
		}
	}
	return bestMotifs;
}

vector<string> motif::randomizedMotifSearch(int k, int runs)
{
	vector<string> metaBestMotifs = initialRandomMotifs(k);
	for (int i = 0; i < runs; i++)
	{
		vector<string> motifs = initialRandomMotifs(k);
		vector<string> bestMotifs = motifs;
		while (true)
		{
			vector<vector<float>> profile = formProfile(motifs);
			motifs = motifsMostProbable(profile);
			if (motifScore(motifs) < motifScore(bestMotifs))
			{
				bestMotifs = motifs;
			}
			else
			{
				break;
			}
		}

		if (motifScore(bestMotifs) < motifScore(metaBestMotifs))
		{
			metaBestMotifs = bestMotifs;
		}
	}

	return metaBestMotifs;
}

vector<string> motif::gibbsSampler(int k, int pop, int runs)
{
	srand((unsigned int)time(NULL));

	int t = m_dna.size();
	vector<string> metaBestMotifs = initialRandomMotifs(k);
	for (int p = 0; p < pop; p++)
	{
		vector<string> motifs = initialRandomMotifs(k);
		vector<string> bestMotifs = motifs;
		for (int i = 0; i < runs; i++)
		{
			int except = rand()%t;
			vector<vector<float>> profile = formProfilePartial(motifs, except);
			string motifNew = profileRandomly(profile, except);
			if (motifScore(motifs) < motifScore(bestMotifs))
			{
				bestMotifs = motifs;
			}
		}

		if (motifScore(bestMotifs) < motifScore(metaBestMotifs))
		{
			metaBestMotifs = bestMotifs;
		}
	}

	return metaBestMotifs;
}

vector<string> motif::initialRandomMotifs(int k)
{
	vector<string> result;
	int max = m_dna[0].size() - k + 1;
	for (auto s : m_dna)
	{
		int index = rand() % max;
		string temp = s.substr(index, k);
		result.push_back(temp);
	}
	return result;
	return vector<string>();
}

vector<string> motif::motifsMostProbable(vector<vector<float>> profile)
{
	vector<string> result;
	int size = m_dna.size();
	for (int i = 0; i < size; i++)
	{
		string temp = profileMostProbable(profile, i);
		result.push_back(temp);
	}
	return result;
}

vector<string> motif::initialBestMotifs(int k)
{
	vector<string> result;
	for (auto s : m_dna)
	{
		string temp = s.substr(0, k);
		result.push_back(temp);
	}
	return result;
}

string motif::profileMostProbable(vector<vector<float>> profile, int position)
{
	int k = profile.size();
	string result;
	float maxProb = -1;
	int size = m_dna[position].size();
	for (int i = 0; i <= size - k; i++)
	{
		string temp = m_dna[position].substr(i, k);
		float tempProb = probGivenProfile(temp, profile);
		if (tempProb > maxProb)
		{
			maxProb = tempProb;
			result = temp;
		}
	}
	return result;
}

string motif::profileRandomly(vector<vector<float>> profile, int position)
{
	string result;
	int k = profile.size();
	int size = m_dna[position].size();
	vector<float> weight;
	for (int i = 0; i <= size - k; i++)
	{
		string temp = m_dna[position].substr(i, k);
		float tempProb = probGivenProfile(temp, profile);
		weight.push_back(tempProb);
	}

	int neo = randomWeight(weight);
	result = m_dna[position].substr(neo, k);
	return result;
}

string motif::orderLine()
{
	string order;
	for (auto s : m_dna)
	{
		order += s;
	}
	return order;
}

string motif::medianString(int k)
{
	int distance = m_dna.size() * k;
	string median;
	vector<string> patterns = allStrings(k);
	for (auto s : patterns)
	{
		int localDistance = DistanceBetweenPatternAndStrings(s);
		if (localDistance < distance)
		{
			distance = localDistance;
			median = s;
		}
	}
	return median;
}

vector<string> motif::motifEnumeration(int k, int d)
{
	vector<string> motifs;
	string ordered = orderLine();
	int fullSize = ordered.size();
	for (int i = 0; i < fullSize - k; i++)
	{
		string localKmer = ordered.substr(i, k);
		vector<string> localNeighbors = neighbors(localKmer, d);
		for (auto s : localNeighbors)
		{
			bool appears = appearsMismatch(s, d);
			if (appears)
			{
				motifs.push_back(s);
			}
		}
	}

	removeDuplicates(motifs);
	return motifs;
}

void motif::print()
{
	cout << "The strings are:\n";
	for (auto s : m_dna)
	{
		cout << s << endl;
	}
}

int motif::DistanceBetweenPatternAndStrings(string pattern)
{
	int size = pattern.size();
	int distance = 0;
	for (auto s : m_dna)
	{
		int hamming = size + 1;
		for (int i = 0; i < s.size() - size; i++)
		{
			string localKmer = s.substr(i, size);
			int maybe = hammingDistance(pattern, localKmer);
			if (maybe < hamming)
			{
				hamming = maybe;
			}
		}

		distance += hamming;
	}
	return distance;
}

bool motif::appearsMismatch(string pattern, int distance)
{
	for (auto dna : m_dna)
	{
		bool localAppear = false;
		int localSize = dna.size();
		int sizeMsg = pattern.size();
		for (int i = 0; i <= (localSize - sizeMsg); i++)
		{
			string subString = dna.substr(i, sizeMsg);
			int hamming = hammingDistance(subString, pattern);
			if (hamming <= distance)
			{
				localAppear = true;
				break;
			}
		}

		if (not localAppear)
		{
			return false;
		}
	}

	return true;
}

int motif::numberAppearsMismatch(string pattern, int distance)
{
	int count = 0;
	int sizeMsg = pattern.size();
	int size = m_dna[0].size() - sizeMsg;
	for (auto dna : m_dna)
	{
		for (int i = 0; i <= size; i++)
		{
			string subString = dna.substr(i, sizeMsg);
			int hamming = hammingDistance(subString, pattern);
			if (hamming <= distance)
			{
				count++;
			}
		}
	}

	return count;
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
		neighborhood = { "A", "C", "G", "T" };
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

void removeDuplicates(vector<string>& vec)
{
	sort(vec.begin(), vec.end());
	vec.erase(unique(vec.begin(), vec.end()), vec.end());
}

vector<string> allStrings(int k)
{
	string temp;
	for (int i = 0; i < k; i++)
	{
		temp += "A";
	}
	vector<string> result = neighbors(temp, k);
	return result;
}

float probGivenProfile(string kmer, vector<vector<float>> profile)
{
	float result = 1;
	int size = kmer.size();
	for (int i = 0; i < size; i++)
	{
		int nucleo = nucleoToInt(kmer[i]);
		result *= profile[i][nucleo];
	}
	return result;
}


int nucleoToInt(char nucleo)
{
	int result;
	if (nucleo == 'A')
	{
		result = 0;
	}
	else if (nucleo == 'C')
	{
		result = 1;
	}
	else if (nucleo == 'G')
	{
		result = 2;
	}
	else
	{
		result = 3;
	}

	return result;
}

vector<vector<float>> formProfile(vector<string> motifs)
{
	vector<vector<float>> result(motifs[0].size(), vector<float>(4, 1));
	float numberMotifs = motifs.size();
	for (int i = 0; i < numberMotifs; i++)
	{
		profileRow(motifs[i], result, numberMotifs);
	}
	return result;
}


void profileRow(string local, vector<vector<float>>& profile, float& weight)
{
	int size = local.size();
	for (int i = 0; i < size; i++)
	{
		char temp = local.at(i);
		profile[i][nucleoToInt(temp)] += 1.0/weight;
	}
}

int motifScore(vector<string> motif)
{
	vector<vector<int>> profile = formCount(motif);
	int score = 0;
	int size = motif.size();

	for (auto c : profile)
	{
		int localScore = size - *max_element(c.begin(), c.end());
		score += localScore;
	}

	return score;
}


vector<vector<int>> formCount(vector<string> motifs)
{
	vector<vector<int>> result(motifs[0].size(), vector<int>(4));
	float numberMotifs = motifs.size();
	for (int i = 0; i < numberMotifs; i++)
	{
		countRow(motifs[i], result);
	}
	return result;
}

void countRow(string local, vector<vector<int>>& counts)
{
	int size = local.size();
	for (int i = 0; i < size; i++)
	{
		char temp = local.at(i);
		counts[i][nucleoToInt(temp)]++;
	}
}

int randomWeight(vector<float> weight)
{
	float sum = accumulate(weight.begin(), weight.end(), 0);
	float rnd = randomFloat(sum);
	int size = weight.size();
	for (int i = 0; i < size; i++)
	{
		if (rnd < weight[i])
		{
			return i;
		}
		else
		{
			rnd -= weight[i];
		}
	}
}

float randomFloat(float b)
{
	float r = float(rand()) / float((RAND_MAX)) * b;
	return r;
}

vector<vector<float>> formProfilePartial(vector<string> motifs, int except)
{
	vector<vector<float>> result(motifs[0].size(), vector<float>(4, 1));
	float numberMotifs = motifs.size() - 1;
	for (int i = 0; i < numberMotifs; i++)
	{
		if (i != except)
		{
			profileRow(motifs[i], result, numberMotifs);
		}
	}
	return result;
}

string consent(vector<vector<int>> counts)
{
	string result;
	for (auto s : counts)
	{
		vector<int>::iterator it = max_element(s.begin(), s.end());
		int index = distance(s.begin(), it);
		if (index == 0)
		{
			result += "A";
		}
		else if (index == 1)
		{
			result += "C";
		}
		else if (index == 2)
		{
			result += "G";
		}
		else
		{
			result += "T";
		}
	}
	return result;
}