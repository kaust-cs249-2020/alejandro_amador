#include "amino.h"

struct golfist
{
	vector<string> name;
	int score;
};

char nucleotideConversor(char nucleotide);
string reverseComplement(string text);
string dnaToRna(string dna);
int findElementVector(vector<vector<string>>& set, vector<string> element);
void sortLeaders(vector<vector<string>>& leaders, vector<int>& points);
void removeElements(vector<vector<string>>& set, int index);
bool compareLinearCircular(vector<string> p1, vector<string> p2);
bool compareAlreadyCircular(vector<string> p, vector<vector<string>> kings);
bool compareAlreadyIncluded(vector<string> p, vector<vector<string>> kings);
void pushScore(vector<string> p, int score, vector<vector<string>>& kings);


//Constructor
amino::amino(string filename, dataType type)
{
	switch (type)
	{
	case dna:
	{
		loadTranslation();
		loadMasses();
		ifstream file(filename + ".txt");
		string temp;
		while (getline(file, temp))
		{
			m_rna.push_back(temp);
		}
		file.close();
		break;
	}
	case spectrum:
	{
		loadTranslation();
		loadMasses();
		ifstream file(filename + ".txt");
		while (!file.fail())
		{
			int temp;
			file >> temp;
			m_spectrum.push_back(temp);
		}
		file.close();
		m_spectrum.pop_back();
		break;
	}
	case totalAlphabet:
	{
		loadTranslation();
		loadAllMass();
		ifstream file(filename + ".txt");
		while (!file.fail())
		{
			int temp;
			file >> temp;
			m_spectrum.push_back(temp);
		}
		file.close();
		m_spectrum.pop_back();
		break;
	}
	default:
	{
		break;
	}
	}
}

//Loads translation file 3-mer to aminos
void amino::loadTranslation()
{
	ifstream file("translation.txt");
	string temp;
	while (getline(file, temp))
	{
		if (temp.size() == 3)
		{
			m_trans.push_back({ temp, "*" });
		}
		else
		{
			string rna = temp.substr(0, 3);
			string aminoAcid = temp.substr(4);
			m_trans.push_back({ rna, aminoAcid });
		}
	}
	file.close();
}

void amino::loadMasses()
{
	ifstream file("integerMass.txt");
	string temp;
	while (getline(file, temp))
	{
		string amino = temp.substr(0, 1);
		int mass = stoi(temp.substr(2));
		m_masses.insert({ amino, mass });
	}
	file.close();
}

void amino::loadAllMass()
{
	for (int i = 57; i <= 200; i++)
	{
		string fakeName = to_string(i);
		m_masses.insert({ fakeName, i });
	}
}

//Translate a 3-mer into aminoacid
string amino::codonToAmino(string codon)
{
	string amino;
	for (int i = 0; i < 64; i++)
	{
		string check = m_trans[i][0];
		if (check == codon)
		{
			return m_trans[i][1];
		}
	}
}

//Converts string of RNA into translation aminoacid
string amino::rnaToAmino(string rna)
{
	string translation ="";
	int size = rna.size() / 3;
	for (int i = 0; i < size; i++)
	{
		string codon = rna.substr(3 * i, 3);
		string amino = codonToAmino(codon);
		if (amino != "*")
		{
			translation += amino;
		}
	}
	return translation;
}


//Counts multiplicity of each aminoacid
void amino::multiplicityCount()
{
	vector<vector<string>> result;
	for (int i = 0; i < m_trans.size(); i++)
	{
		string temp = m_trans[i][1];
		bool appears = false;
		for (int j = 0; j < result.size(); j++)
		{
			if (result[j][0] == temp)
			{
				appears = true;
				break;
			}
		}
		if (!appears)
		{
			int counter = 1;
			for (int k = i+1; k < m_trans.size(); k++)
			{
				if (m_trans[k][1] == temp)
				{
					counter++;
				}
			}
			vector<string> lol = { temp, to_string(counter) };
			result.push_back(lol);
		}
	}

	printMatrix(result);
}

vector<string> amino::substringsEncoding(string aminos)
{
	vector<string> result;
	int localSize = 3 * aminos.size();
	int size = m_rna[0].size() - localSize;
	for (int i = 0; i < size; i++)
	{
		string temp = m_rna[0].substr(i, localSize);
		string tempReverse = reverseComplement(temp);
		string check1 = rnaToAmino(dnaToRna(temp));
		string check2 = rnaToAmino(dnaToRna(tempReverse));
		if (check1 == aminos || check2 == aminos)
		{
			result.push_back(temp);
		}
	}
	return result;
}


//Returns linear spectrum of a peptide
vector<int> amino::linearSpectrum(vector<string> peptide)
{
	vector<int> prefixMass(peptide.size()+1,0);
	for (int i = 0; i < peptide.size(); i++)
	{
		string amino = peptide[i];
		int aminoMass = m_masses[amino];
		prefixMass[i + 1] = prefixMass[i] + aminoMass;
	}
	vector<int> linearSpectrum = { 0 };
	for (int i = 0; i < peptide.size(); i++)
	{
		for (int j = i + 1; j <= peptide.size(); j++)
		{
			int temp = prefixMass[j] - prefixMass[i];
			linearSpectrum.push_back(temp);
		}
	}

	sort(linearSpectrum.begin(), linearSpectrum.end());
	return linearSpectrum;
}

//Returns cyclic spectrum of a peptide
vector<int> amino::cyclicSpectrum(vector<string> peptide)
{
	vector<int> prefixMass(peptide.size() + 1, 0);
	for (int i = 0; i < peptide.size(); i++)
	{
		string amino = peptide[i];
		int aminoMass = m_masses[amino];
		prefixMass[i + 1] = prefixMass[i] + aminoMass;
	}

	int peptideMass = prefixMass.back();
	vector<int> cyclicSpectrum = { 0 };

	for (int i = 0; i < peptide.size(); i++)
	{
		for (int j = i + 1; j <= peptide.size(); j++)
		{
			int temp = prefixMass[j] - prefixMass[i];
			cyclicSpectrum.push_back(temp);
			if (i > 0 and j < peptide.size())
			{
				int tempTwo = peptideMass - temp;
				cyclicSpectrum.push_back(tempTwo);
			}
		}
	}

	sort(cyclicSpectrum.begin(), cyclicSpectrum.end());
	return cyclicSpectrum;
}

bool amino::isConsistent(vector<int> subset)
{
	sort(subset.begin(), subset.end());
	return includes(m_spectrum.begin(), m_spectrum.end(), subset.begin(), subset.end());
}

long long int amino::numberOfPeptides(int m)
{
	vector<long long int> counter(m + 1, 0);
	for (const auto& x : m_masses)
	{
		int tempMass = x.second;
		counter[tempMass]++;
	}

	for (int i = 0; i <= m; i++)
	{
		for (const auto& x : m_masses)
		{
			long long int temp = i - x.second;
			if (temp > 0)
			{
				counter[i] += counter[temp];
			}
		}
	}

	long long int result = counter.back();
	return result;
}

int amino::peptideMass(vector<string> peptide)
{
	long long int mass = 0;
	for (int i = 0; i < peptide.size(); i++)
	{
		string amino = peptide[i];
		mass += m_masses[amino];
	}
	return mass;
}

void amino::expand(vector<vector<string>>& peptides)
{
	int originalSize = peptides.size();

	if (peptides[0][0] == "")
	{
		for (int i = 0; i < originalSize; i++)
		{
			vector<string> temp = peptides[0];
			for (auto& s : m_masses)
			{
				if (s.first != "Q" && s.first != "L")
				{
					vector<string> branch = {s.first};
					peptides.push_back(branch);
				}
			}
			peptides.erase(peptides.begin());
		}
	}
	else
	{
		for (int i = 0; i < originalSize; i++)
		{
			vector<string> temp = peptides[0];
			for (auto& s : m_masses)
			{
				if (s.first != "Q" && s.first != "L")
				{
					vector<string> branch = temp;
					branch.push_back(s.first);
					peptides.push_back(branch);
				}
			}
			peptides.erase(peptides.begin());
		}
	}
}

vector<vector<string>> amino::cycloPeptideSequencing()
{
	int parentMass = m_spectrum.back();
	vector<vector<string>> candidatePeptides = { {""} };
	vector<vector<string>> finalPeptides;
	while (candidatePeptides.size() != 0)
	{
		expand(candidatePeptides);
		int originalSize = candidatePeptides.size();
		int index = 0;
		while (index < originalSize)
		{
			vector<string> localPep = candidatePeptides[index];
			vector<int> localSpectrum = linearSpectrum(localPep);
			if (peptideMass(localPep) == parentMass)
			{
				int localPos = findElementVector(finalPeptides, localPep);
				vector<int> localCyclic = cyclicSpectrum(localPep);
				if (localCyclic == m_spectrum && localPos == -1)
				{
					finalPeptides.push_back(localPep);
				}
				candidatePeptides.erase(candidatePeptides.begin() + index);
				originalSize--;
			}

			else if (!isConsistent(localSpectrum))
			{
				candidatePeptides.erase(candidatePeptides.begin() + index);
				originalSize--;
			}
			else
			{
				index++;
			}
		}
	}

	return finalPeptides;
}

vector<string> amino::convolutionCycloPeptideSequencing(int M, int N)
{
	vector<string> king;
	extendedAlphabet(M);
	king = leaderBoardCycloPeptideSequencing(N);
	return king;
}

vector<string> amino::leaderBoardCycloPeptideSequencing(int N)
{
	vector<vector<string>> cumbiaKings;
	int parentMass = m_spectrum.back();
	vector<vector<string>> leaderBoard = { {""} };
	vector<string> king = {""};
	while (leaderBoard.size() != 0)
	{
		expand(leaderBoard);
		int originalSize = leaderBoard.size();
		int index = 0;
		while (index < originalSize)
		{
			vector<string> localPep = leaderBoard[index];
			if (peptideMass(localPep) == parentMass)
			{
				pushScore(localPep, cumbiaKings);
				if (score(localPep) > score(king))
				{
					cout << score(localPep) << endl;
					king = localPep;
				}


				index++;
			}
			else if (peptideMass(localPep) > parentMass)
			{
				leaderBoard.erase(leaderBoard.begin() + index);
				originalSize--;
			}
			else
			{
				index++;
			}
		}

		trim(leaderBoard, N);
	}

	vector<vector<string>> realKings;
	for (int i = 0; i < 86; i++)
	{
		realKings.insert(realKings.begin(), cumbiaKings[cumbiaKings.size() - 1 - i]);
	}
	vector<vector<int>> lo = peptidesToMass(realKings);
	printPeptideSequence(lo);
	cout << realKings.size() << endl;
	cout << cumbiaKings.size() << endl;
	for (auto s : realKings)
	{
		cout << score(s) << endl;
	}
	return king;
}

vector<vector<int>> amino::peptidesToMass(vector<vector<string>> peptides)
{
	vector<vector<int>> result;
	for (int i = 0; i < peptides.size(); i++)
	{
		vector<int> temp = peptideToMass(peptides[i]);
		result.push_back(temp);
	}
	return result;
}

vector<int> amino::peptideToMass(vector<string> peptide)
{
	vector<int> result;
	for (int i = 0; i < peptide.size(); i++)
	{
		string temp = peptide[i];
		for (auto&s : m_masses)
		{
			if (s.first == temp)
			{
				result.push_back(s.second);
				break;
			}
		}
	}
	return result;
}

int amino::score(vector<string> peptide)
{
	vector<string> fatality = { "" };
	if (peptide == fatality)
	{
		return 0;
	}
	else
	{
		int score = 0;
		vector<int> localCyclic = cyclicSpectrum(peptide);
		vector<int> localCopy = m_spectrum;
		for (auto& s : localCyclic)
		{
			vector<int>::iterator it = find(localCopy.begin(), localCopy.end(), s);
			if (it != localCopy.end())
			{
				score++;
				localCopy.erase(it);
			}
		}
		return score;
	}
}

int amino::linearScore(vector<string> peptide)
{
	vector<string> fatality = { "" };
	if (peptide == fatality)
	{
		return 0;
	}
	else
	{
		int score = 0;
		vector<int> localCyclic = linearSpectrum(peptide);
		vector<int> localCopy = m_spectrum;
		for (auto& s : localCyclic)
		{
			vector<int>::iterator it = find(localCopy.begin(), localCopy.end(), s);
			if (it != localCopy.end())
			{
				score++;
				localCopy.erase(it);
			}
		}
		return score;
	}
}

void amino::trim(vector<vector<string>>& leaderboard, int N)
{
	int sizeBoard = leaderboard.size();
	vector<int> linearScores;
	for (int i = 0; i < sizeBoard; i++)
	{
		vector<string> pept = leaderboard[i];
		int peptScore = linearScore(pept);
		linearScores.push_back(peptScore);
	}

	sortLeaders(leaderboard, linearScores);
	for (int j = N; j < sizeBoard; j++)
	{
		if (linearScores[j] < linearScores[N-1])
		{
			removeElements(leaderboard, j);
			return;
		}
	}

	return;

}

vector<int> amino::convolution()
{
	vector<int> result;
	for (int i = 0; i < m_spectrum.size(); i++)
	{
		int localMass = m_spectrum[i];
		for (int j = 0; j < m_spectrum.size(); j++)
		{
			int temp = localMass - m_spectrum[j];
			if (temp > 0)
			{
				result.push_back(temp);
			}
		}
	}

	sort(result.begin(), result.end());
	return result;
}

void amino::extendedAlphabet(int M)
{
	//First makes convolution, checks elements with repetitions and sort them
	struct counting
	{
		int mass;
		int times;
	};

	vector<int> convo = convolution();
	int index = 0;
	int lastMass = convo[0];
	vector<counting> counter;
	counting first = { lastMass, 1 };
	counter.push_back(first);
	for (int i = 1; i < convo.size(); i++)
	{
		int tempMass = convo[i];
		if (tempMass == lastMass)
		{
			counter[index].times++;
		}
		else
		{
			counting newValue = { tempMass, 1 };
			counter.push_back(newValue);
			lastMass = tempMass;
			index++;
		}
	}

	sort(counter.begin(), counter.end(), [](auto const &a, auto const &b) { return a.times > b.times; });

	//Creates new alphabet with most frequent masses.
	m_masses.clear();
	int idx = 0;
	int top = 0;
	while(top < M)
	{
		int mass = counter[idx].mass;
		if (mass >= 57 && mass <= 200)
		{
			string name = "P" + to_string(mass);//to_string(top+1);
			m_masses.insert({ name, mass });
			idx++;
			top++;
		}
		else
		{
			idx++;
		}
	}

}

void amino::pushScore(vector<string> p, vector<vector<string>>& kings)
{
	int newScore = score(p);
	for (int i = 0; i < kings.size(); i++)
	{
		int temp = score(kings[i]);
		if (newScore <= temp)
		{
			kings.insert(kings.begin() + i, p);
			return;
		}
	}

	kings.push_back(p);
}

//Reverse complement of a DNA string
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


//For reverse complement
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

//Transforms dna sequence into rna sequence
string dnaToRna(string dna)
{
	string rna;
	for (auto c : dna)
	{
		if (c == 'T')
		{
			rna += 'U';
		}
		else
		{
			rna += c;
		}
	}
	return rna;
}

int findElementVector(vector<vector<string>>& set, vector<string> element)
{
	vector<vector<string>>::iterator it = find(set.begin(), set.end(), element);
	if (it == set.end())
	{
		return -1;
	}
	else
	{
		int temp = it - set.begin();
		return temp;
	}
}

void sortLeaders(vector<vector<string>>& leaders, vector<int>& points)
{
	vector<golfist> leaderList;
	for (int i = 0; i < leaders.size(); i++)
	{
		golfist temp = { leaders[i], points[i] };
		leaderList.push_back(temp);
	}

	sort(leaderList.begin(), leaderList.end(),[](auto const &a, auto const &b) { return a.score > b.score; });
	for (int i = 0; i < leaders.size(); i++)
	{
		leaders[i] = leaderList[i].name;
		points[i] = leaderList[i].score;
	}
}

void removeElements(vector<vector<string>>& set, int index)
{
	set.erase(set.begin() + index, set.end());
}

bool compareLinearCircular(vector<string> p1, vector<string> p2)
{
	int size = p1.size();
	int sizeTemp = p2.size();
	if (size != sizeTemp)
	{
		return false;
	}

	vector<string> t1 = p1;
	vector<string> t2 = p2;
	for (int i = 0; i < size; i++)
	{
		string temp = t1[i];
		vector<string>::iterator it = find(t2.begin(), t2.end(), temp);
		if (it == t2.end())
		{
			return false;
		}
		else
		{
			t2.erase(it);
		}
	}

	return true;
}

bool compareAlreadyCircular(vector<string> p, vector<vector<string>> kings)
{
	for (int i = 0; i < kings.size(); i++)
	{
		bool temp = compareLinearCircular(p, kings[i]);
		if (temp)
		{
			return true;
		}
	}

	return false;
}

bool compareAlreadyIncluded(vector<string> p, vector<vector<string>> kings)
{
	for (int i = 0; i < kings.size(); i++)
	{
		bool temp = kings[i] == p;
		if (temp)
		{
			return true;
		}
	}

	return false;
}
