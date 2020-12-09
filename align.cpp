#include "align.h"

vector<vector<int>> loadScoreMatrix();
int maxInGraph(vector<vector<int>> graph, string& x);
int fittingMax(vector<vector<int>> graph, string& x);
int overlapMax(vector<vector<int>> graph, string& x);
void checkJump(string check, int&x, int& y);
void checkJumpFitting(string check, int&x, int& y);
int max3(int a, int b, int c);
int max4(int a, int b, int c, int d);
int case7(int a, int b, int c, int d, int e, int f, int g, multiSteps& type);
string reverseStr(string str);
void initMulti(vector<vector<vector<multiSteps>>>& s, int n, int m, int l);

alignment::alignment(string name, problemType type)
{
	if (type == coin)
	{
		ifstream file(name+".txt");
		string temp;
		getline(file, temp);
		int firstPos = temp.find(",");
		int secondPos = temp.find(",", firstPos + 1);
		int num = stoi(temp.substr(0, firstPos));
		m_data.push_back(num);
		while (true)
		{
			if (secondPos == string::npos)
			{
				num = stoi(temp.substr(firstPos + 1));
				m_data.push_back(num);
				break;
			}
			else
			{
				num = stoi(temp.substr(firstPos + 1, secondPos - firstPos));
				m_data.push_back(num);
				firstPos = secondPos;
				secondPos = temp.find(",", firstPos + 1);
			}
		}

		m_scoreMatrix = loadScoreMatrix();
	}
}

int alignment::matchScore(char x)
{
	string table = "ACDEFGHIKLMNPQRSTVWY";
	int index = table.find(x);
	int result = m_scoreMatrix[index][index];
	return 1;//result;
}

int alignment::mismatchScore(char x, char y)
{
	string table = "ACDEFGHIKLMNPQRSTVWY";
	int a = table.find(x);
	int b = table.find(y);
	int result = m_scoreMatrix[a][b];
	return -1;//result;
}

int alignment::DPChange(int money)
{
	int infty = 99999;
	vector<int> minNumCoins = { 0 };
	for (int m = 1; m <= money; m++)
	{
		minNumCoins.push_back(infty);
		for (int i = 0; i < m_data.size(); i++)
		{
			if (m >= m_data[i])
			{
				if (minNumCoins[m - m_data[i]] + 1 < minNumCoins[m])
				{
					minNumCoins[m] = minNumCoins[m - m_data[i]] + 1;
				}
			}
		}
	}

	return minNumCoins.back();
}

int alignment::ManhattanTourist(int n, int m)
{
	vector<vector<int>> down = readMatrix("down.txt", n, m + 1, 0);
	vector<vector<int>> right = readMatrix("right.txt", n + 1, m, 0);
	vector<vector<int>> s(n + 1, vector<int>(m + 1, 0));
	for (int i = 1; i <= n; i++)
	{
		s[i][0] = s[i - 1][0] + down[i - 1][0];
	}
	for (int j = 1; j <= m; j++)
	{
		s[0][j] = s[0][j - 1] + right[0][j - 1];
	}
	for (int i = 1; i <= n; i++)
	{
		for (int j = 1; j <= m; j++)
		{
			int temp1 = s[i - 1][j] + down[i - 1][j];
			int temp2 = s[i][j - 1] + right[i][j - 1];
			s[i][j] = max(temp1, temp2);
		}
	}

	return s[n][m];
}

vector<vector<steps>> alignment::LCSBackTrack(string v, string w)
{
	int n = v.size();
	int m = w.size();
	vector<vector<steps>> backtrack(n+1, vector<steps>(m+1, down));
	vector<vector<int>> s(n + 1, vector<int>(m + 1, 0));
	for (int i = 1; i <= n; i++)
	{
		for (int j = 1; j <= m; j++)
		{
			int match = 0;
			if (v.at(i - 1) == w.at(j - 1))
			{
				match = 1;
			}
			s[i][j] = max(s[i - 1][j], max(s[i][j - 1], s[i - 1][j - 1] + match));
			if (s[i][j] == s[i - 1][j])
			{
				backtrack[i][j] = down;
			}
			else if (s[i][j] == s[i][j - 1])
			{
				backtrack[i][j] = rights;
			}
			else if (s[i][j] == s[i - 1][j - 1] + match)
			{
				backtrack[i][j] = diagonal;
			}
		}
	}

	return backtrack;
}

string alignment::outputLCS(vector<vector<steps>> backtrack, string v, int i, int j)
{
	if (i == 0 || j == 0)
	{
		return "";
	}
	else
	{

		if (backtrack[i][j] == down)
		{
			return outputLCS(backtrack, v, i - 1, j);
		}
		else if (backtrack[i][j] == rights)
		{
			return outputLCS(backtrack, v, i, j - 1);
		}
		else if(backtrack[i][j] == diagonal)
		{
			return outputLCS(backtrack, v, i - 1, j - 1) + v.at(i-1);
		}
	}
}

vector<vector<steps>> alignment::generalBackTrack(string v, string w, int sigma, int& score)
{
	int n = v.size();
	int m = w.size();
	vector<vector<steps>> backtrack(n + 1, vector<steps>(m + 1, down));
	vector<vector<int>> s(n + 1, vector<int>(m + 1, 0));
	for (int i = 1; i <= n; i++)
	{
		s[i][0] = s[i - 1][0] - sigma;
		backtrack[i][0] = down;
	}
	for (int j = 1; j <= m; j++)
	{
		s[0][j] = s[0][j - 1] - sigma;
		backtrack[0][j] = rights;
	}

	for (int i = 1; i <= n; i++)
	{
		for (int j = 1; j <= m; j++)
		{
			char vChar = v.at(i - 1);
			char wChar = w.at(j - 1);
			int match;
			if (vChar==wChar)
			{
				match = matchScore(vChar);
			}
			else
			{
				match = mismatchScore(vChar, wChar);
			}
			s[i][j] = max(s[i - 1][j] - sigma, max(s[i][j - 1] - sigma, s[i - 1][j - 1] + match));
			if (s[i][j] == s[i - 1][j] - sigma)
			{
				backtrack[i][j] = down;
				cout << "ONE" << endl;
			}
			else if (s[i][j] == s[i][j - 1] - sigma)
			{
				backtrack[i][j] = rights;
				cout << "Two" << endl;
			}
			else if (s[i][j] == s[i - 1][j - 1] + match)
			{
				backtrack[i][j] = diagonal;
				cout << "three" << endl;
			}
		}
	}

	printMatrix(s);
	score =  s[n][m];
	return backtrack;
}

vector<string> alignment::generalLCS(vector<vector<steps>> backtrack, string v, string w, int i, int j)
{
	vector<string> result;
	if (i == 0 && j == 0)
	{
		result =  { "" ,"" };
	}
	else
	{

		if (backtrack[i][j] == down)
		{
			result =  generalLCS(backtrack, v, w, i - 1, j);
			result[0] += v.substr(i - 1,1);
			result[1] += "-";
		}
		else if (backtrack[i][j] == rights)
		{
			result = generalLCS(backtrack, v, w, i, j - 1);
			result[0] += "-";
			result[1] += w.substr(j - 1, 1);
		}
		else if (backtrack[i][j] == diagonal)
		{
			result = generalLCS(backtrack, v, w, i - 1, j - 1);
			result[0] += v.substr(i - 1, 1);
			result[1] += w.substr(j - 1, 1);
		}
	}

	return result;
}

vector<vector<vector<multiSteps>>> alignment::multiDimBackTrack(string u, string v, string w, int & score)
{
	int n = u.size();
	int m = v.size();
	int l = w.size();
	vector<vector<vector<multiSteps>>> back(n + 1, vector<vector<multiSteps>>(m + 1, vector<multiSteps>(l + 1, moveX)));
	vector<vector<vector<int>>> s(n + 1, vector<vector<int>>(m + 1, vector<int>(l + 1, 0)));
	initMulti(back, n, m, l);

	for (int i = 1; i <= n; i++)
	{
		for (int j = 1; j <= m; j++)
		{
			for (int k = 1; k <= l; k++)
			{
				char uChar = u.at(i - 1);
				char vChar = v.at(j - 1);
				char wChar = w.at(k - 1);
				int match = 0;
				if (vChar == uChar && uChar == wChar)
				{
					match = 1;
				}

				multiSteps tempCase;
				s[i][j][k] = case7(s[i - 1][j][k], s[i][j - 1][k], s[i][j][k - 1], s[i - 1][j - 1][k],
					s[i - 1][j][k - 1], s[i][j - 1][k - 1], s[i - 1][j - 1][k - 1] + match, tempCase);

				back[i][j][k] = tempCase;
			}
		}
	}

	score = s[n][m][l];
	return back;
}

vector<string> alignment::multiDimLCS(vector<vector<vector<multiSteps>>> back, string u, string v, string w, int i, int j, int k)
{
	vector<string> result;
	if (i == 0 && j == 0 && k==0)
	{
		result = { "" ,"", "" };
	}
	else
	{
		switch (back[i][j][k])
		{
		case(moveX):
		{
				result = multiDimLCS(back, u, v, w, i - 1, j, k);
				result[0] += u.substr(i - 1, 1);
				result[1] += "-";
				result[2] += "-";
				break;
		}
		case(moveY):
		{
			result = multiDimLCS(back, u, v, w, i, j - 1, k);
			result[0] += "-";
			result[1] += v.substr(j - 1, 1);
			result[2] += "-";
			break;
		}
		case(moveZ):
		{
			result = multiDimLCS(back, u, v, w, i, j, k - 1);
			result[0] += "-";
			result[1] += "-";
			result[2] += w.substr(k - 1, 1);
			break;
		}
		case(diagX):
		{
			result = multiDimLCS(back, u, v, w, i, j - 1, k - 1);
			result[0] += "-";
			result[1] += v.substr(j - 1, 1);
			result[2] += w.substr(k - 1, 1);
			break;
		}
		case(diagY):
		{
			result = multiDimLCS(back, u, v, w, i - 1, j, k - 1);
			result[0] += u.substr(i - 1, 1);
			result[1] += "-";
			result[2] += w.substr(k - 1, 1);
			break;
		}
		case(diagZ):
		{
			result = multiDimLCS(back, u, v, w, i - 1, j - 1, k);
			result[0] += u.substr(i - 1, 1);
			result[1] += v.substr(j - 1, 1);
			result[2] += "-";
			break;
		}
		case(superDiag):
		{
			result = multiDimLCS(back, u, v, w, i - 1, j - 1, k - 1);
			result[0] += u.substr(i - 1, 1);
			result[1] += v.substr(j - 1, 1);
			result[2] += w.substr(k - 1, 1);
			break;
		}
		}
	}

	return result;
}

vector<vector<string>> alignment::taxiBackTrack(string v, string w, int sigma)
{
	int n = v.size();
	int m = w.size();
	vector<vector<string>> backtrack(n + 1, vector<string>(m + 1, "jumpStart"));
	vector<vector<int>> s(n + 1, vector<int>(m + 1, 0));

	for (int i = 1; i <= n; i++)
	{
		for (int j = 1; j <= m; j++)
		{
			char vChar = v.at(i - 1);
			char wChar = w.at(j - 1);
			int match;
			if (vChar == wChar)
			{
				match = matchScore(vChar);
			}
			else
			{
				match = mismatchScore(vChar, wChar);
			}
			s[i][j] = max(s[i - 1][j] - sigma, max(s[i][j - 1] - sigma, max(s[i - 1][j - 1] + match, 0)));
			if (s[i][j] == s[i - 1][j] - sigma)
			{
				backtrack[i][j] = "down";
			}
			else if (s[i][j] == s[i][j - 1] - sigma)
			{
				backtrack[i][j] = "rights";
			}
			else if (s[i][j] == s[i - 1][j - 1] + match)
			{
				backtrack[i][j] = "diagonal";
			}
			else if (s[i][j] == 0)
			{
				backtrack[i][j] = "jumpStart";
			}
		}
	}


	int tempMax = s[n][m];
	string comesFrom;
	int newMax = maxInGraph(s, comesFrom);
	if (newMax > tempMax)
	{
		s[n][m] = newMax;
		backtrack[n][m] = comesFrom;
	}
	cout << s[n][m] << endl;
	return backtrack;
}

vector<string> alignment::taxiLCS(vector<vector<string>> backtrack, string v, string w, int i, int j)
{
	vector<string> result;
	if (i == 0 && j == 0)
	{
		result = { "" ,"" };
	}
	else
	{

		if (backtrack[i][j] == "down")
		{
			result = taxiLCS(backtrack, v, w, i - 1, j);
			result[0] += v.substr(i - 1, 1);
			result[1] += "-";
		}
		else if (backtrack[i][j] == "rights")
		{
			result = taxiLCS(backtrack, v, w, i, j - 1);
			result[0] += "-";
			result[1] += w.substr(j - 1, 1);
		}
		else if (backtrack[i][j] == "diagonal")
		{
			result = taxiLCS(backtrack, v, w, i - 1, j - 1);
			result[0] += v.substr(i - 1, 1);
			result[1] += w.substr(j - 1, 1);
		}
		else if (backtrack[i][j] == "jumpStart")
		{
			result = { "" ,"" };
		}
		else 
		{
			int newI, newJ;
			checkJump(backtrack[i][j], newI, newJ);
			result = taxiLCS(backtrack, v, w, newI, newJ);
		}

	}

	return result;
}

vector<vector<string>> alignment::fittingBackTrack(string v, string w, int sigma, int& score)
{
	int n = v.size();
	int m = w.size();
	vector<vector<string>> backtrack(n + 1, vector<string>(m + 1, "down"));
	vector<vector<int>> s(n + 1, vector<int>(m + 1, 0));
	for (int i = 1; i <= n; i++)
	{
		s[i][0] = s[i - 1][0];
		backtrack[i][0] = "down";
	}
	for (int j = 1; j <= m; j++)
	{
		s[0][j] = s[0][j - 1] - sigma;
		backtrack[0][j] = "rights";
	}

	for (int i = 1; i <= n; i++)
	{
		for (int j = 1; j <= m; j++)
		{
			char vChar = v.at(i - 1);
			char wChar = w.at(j - 1);
			int match;
			if (vChar == wChar)
			{
				match = matchScore(vChar);
			}
			else
			{
				match = mismatchScore(vChar, wChar);
			}
			s[i][j] = max(s[i - 1][j] - sigma, max(s[i][j - 1] - sigma, s[i - 1][j - 1] + match));
			if (s[i][j] == s[i - 1][j] - sigma)
			{
				backtrack[i][j] = "down";
			}
			else if (s[i][j] == s[i][j - 1] - sigma)
			{
				backtrack[i][j] = "rights";
			}
			else if (s[i][j] == s[i - 1][j - 1] + match)
			{
				backtrack[i][j] = "diagonal";
			}
		}
	}

	string comesFrom;
	int tempMax = s[n][m];
	int newMax = fittingMax(s, comesFrom);
	if (newMax > tempMax)
	{
		s[n][m] = newMax;
		backtrack[n][m] = comesFrom;
	}

	score = s[n][m];
	return backtrack;
}

vector<string> alignment::fittingString(vector<vector<string>> backtrack, string v, string w, int i, int j)
{
	vector<string> result;
	if (i == 0 && j == 0)
	{
		result = { "" ,"" };
	}
	else
	{

		if (backtrack[i][j] == "down")
		{
			result = fittingString(backtrack, v, w, i - 1, j);
			if (j != 0)
			{
				result[0] += v.substr(i - 1, 1);
				result[1] += "-";
			}
		}
		else if (backtrack[i][j] == "rights")
		{
			result = fittingString(backtrack, v, w, i, j - 1);
			result[0] += "-";
			result[1] += w.substr(j - 1, 1);
		}
		else if (backtrack[i][j] == "diagonal")
		{
			result = fittingString(backtrack, v, w, i - 1, j - 1);
			result[0] += v.substr(i - 1, 1);
			result[1] += w.substr(j - 1, 1);
		}
		else
		{
			int newI, newJ;
			checkJumpFitting(backtrack[i][j], newI, newJ);
			result = fittingString(backtrack, v, w, newI, newJ);
		}
	}

	return result;
}

vector<vector<string>> alignment::overlapBackTrack(string v, string w, int sigma, int & score)
{
	int n = v.size();
	int m = w.size();
	vector<vector<string>> backtrack(n + 1, vector<string>(m + 1, "down"));
	vector<vector<int>> s(n + 1, vector<int>(m + 1, 0));
	for (int i = 1; i <= n; i++)
	{
		backtrack[i][0] = "down";
	}
	for (int j = 1; j <= m; j++)
	{
		backtrack[0][j] = "rights";
	}

	for (int i = 1; i <= n; i++)
	{
		for (int j = 1; j <= m; j++)
		{
			char vChar = v.at(i - 1);
			char wChar = w.at(j - 1);
			int match;
			if (vChar == wChar)
			{
				match = 1;
			}
			else
			{
				match = -2;
			}
			s[i][j] = max3(s[i - 1][j] - sigma, s[i][j - 1] - sigma, s[i - 1][j - 1] + match);
			if (s[i][j] == s[i - 1][j - 1] + match)
			{
			backtrack[i][j] = "diagonal";
			}
			else if (s[i][j] == s[i - 1][j] - sigma)
			{
				backtrack[i][j] = "down";
			}
			else if (s[i][j] == s[i][j - 1] - sigma)
			{
				backtrack[i][j] = "rights";
			}

		}
	}

	string comesFrom;
	int tempMax = s[n][m];
	int newMax = overlapMax(s, comesFrom);
	if (newMax > tempMax)
	{
		s[n][m] = newMax;
		backtrack[n][m] = comesFrom;
	}

	score = s[n][m];
	return backtrack;
}

vector<string> alignment::overlapAllignment(vector<vector<string>> backtrack, string v, string w, int i, int j)
{
	vector<string> result;
	if (i == 0 && j == 0)
	{
		result = { "" ,"" };
	}
	else
	{

		if (backtrack[i][j] == "down")
		{
			result = overlapAllignment(backtrack, v, w, i - 1, j);
			if (j >= 0)
			{
				result[0] += v.substr(i - 1, 1);
				result[1] += "-";
			}
		}
		else if (backtrack[i][j] == "rights")
		{
			result = overlapAllignment(backtrack, v, w, i, j - 1);
			if (i != 0)
			{
				result[0] += "-";
				result[1] += w.substr(j - 1, 1);
			}
		}
		else if (backtrack[i][j] == "diagonal")
		{
			result = overlapAllignment(backtrack, v, w, i - 1, j - 1);
			result[0] += v.substr(i - 1, 1);
			result[1] += w.substr(j - 1, 1);
		}
		else if (backtrack[i][j] == "jumpStart")
		{
			result = overlapAllignment(backtrack, v, w, i - 1, 0);
		}
		else
		{
			int newI, newJ;
			checkJump(backtrack[i][j], newI, newJ);
			result = overlapAllignment(backtrack, v, w, newI, newJ);
		}

	}

	return result;
}

vector<vector<vector<string>>> alignment::gapBackTrack(string v, string w, int sigma, int epsi, int & score, int & level)
{
	int n = v.size();
	int m = w.size();
	vector<vector<string>> backTop(n + 1, vector<string>(m + 1, "right"));
	vector<vector<string>> backMid(n + 1, vector<string>(m + 1, "diagonal"));
	vector<vector<string>> backBot(n + 1, vector<string>(m + 1, "down"));
	vector<vector<int>> sTop(n + 1, vector<int>(m + 1, 0));
	vector<vector<int>> sMid(n + 1, vector<int>(m + 1, 0));
	vector<vector<int>> sBot(n + 1, vector<int>(m + 1, 0));


	for (int i = 1; i <= n; i++)
	{
		for (int j = 1; j <= m; j++)
		{
			char vChar = v.at(i - 1);
			char wChar = w.at(j - 1);
			int match;
			if (vChar == wChar)
			{
				match = matchScore(vChar);
			}
			else
			{
				match = mismatchScore(vChar, wChar);
			}

			sBot[i][j] = max(sBot[i - 1][j] - epsi, sMid[i - 1][j] - sigma);
			sTop[i][j] = max(sTop[i][j - 1] - epsi, sMid[i][j - 1] - sigma);
			sMid[i][j] = max3(sBot[i][j], sMid[i - 1][j - 1] + match, sTop[i][j]);

			//Bot score
			if (sBot[i][j] == sBot[i - 1][j] - epsi)
			{
				backBot[i][j] = "bot";
			}
			else if (sBot[i][j] == sMid[i - 1][j] - sigma)
			{
				backBot[i][j] = "mid";
			}

			//Top score
			if (sTop[i][j] == sTop[i][j - 1] - epsi)
			{
				backTop[i][j] = "top";
			}
			else if (sTop[i][j] == sMid[i][j - 1] - sigma)
			{
				backTop[i][j] = "mid";
			}

			//Mid score
			if (sMid[i][j] == sBot[i][j])
			{
				backMid[i][j] = "bot";
			}
			else if (sMid[i][j] == sMid[i - 1][j - 1] + match)
			{
				backMid[i][j] = "mid";
			}
			else if (sMid[i][j] == sTop[i][j])
			{
				backMid[i][j] = "top";
			}

		}
	}

	score = max3(sMid[n][m], sBot[n][m], sTop[n][m]);
	if (score == sMid[n][m])
	{
		level = 1;
	}
	else if (score == sBot[n][m])
	{
		level = 0;
	}
	else if (score == sTop[n][m])
	{
		level = 2;
	}

	vector<vector<vector<string>>> multiScore = { backBot, backMid, backTop };
	return multiScore;
}

vector<string> alignment::gapAllignment(vector<vector<vector<string>>> backtrack, string v, string w, int i, int j, int k)
{
	vector<string> result;
	if (i == 0 && j == 0)
	{
		result = { "" ,"" };
	}
	else
	{
		if (k == 0)
		{
			if (backtrack[k][i][j] == "bot")
			{
				result = gapAllignment(backtrack, v, w, i - 1, j, 0);
				result[0] += v.substr(i - 1, 1);
				result[1] += "-";
			}
			else if (backtrack[k][i][j] == "mid")
			{
				result = gapAllignment(backtrack, v, w, i - 1, j, 1);
				result[0] += v.substr(i - 1, 1);
				result[1] += "-";
			}
		}
		else if (k == 1)
		{
			if (backtrack[k][i][j] == "bot")
			{
				result = gapAllignment(backtrack, v, w, i, j, 0);
			}
			else if (backtrack[k][i][j] == "mid")
			{
				result = gapAllignment(backtrack, v, w, i - 1, j - 1, 1);
				result[0] += v.substr(i - 1, 1);
				result[1] += w.substr(j - 1, 1);
			}
			else if (backtrack[k][i][j] == "top")
			{
				result = gapAllignment(backtrack, v, w, i, j, 2);
			}
		}
		else if (k == 2)
		{

			if (backtrack[k][i][j] == "top")
			{
				result = gapAllignment(backtrack, v, w, i, j - 1, 2);
				result[0] += "-";
				result[1] += w.substr(j - 1, 1);
			}
			else if (backtrack[k][i][j] == "mid")
			{
				result = gapAllignment(backtrack, v, w, i, j - 1, 1);
				result[0] += "-";
				result[1] += w.substr(j - 1, 1);
			}
		}

	}

	return result;
}

steps alignment::middleEdge(string realV, string realW, int top, int bottom, int left, int right, int sigma, int& point)
{
	string v = realV.substr(top, bottom - top);
	string w = realW.substr(left, right - left);
	steps result;
	int n = v.size();
	int m = w.size();
	int midCol = floor(m/2);
	vector<vector<int>> sNormal(n + 1, vector<int>(midCol + 1, 0));
	vector<vector<int>> sReverse(n + 1, vector<int>(m - midCol + 1, 0));
	vector<steps> backtrack(n + 1, rights);
	string vRev = reverseStr(v);
	string wRev = reverseStr(w);


	for (int i = 1; i <= n; i++)
	{
		sNormal[i][0] = sNormal[i - 1][0] - sigma;
		sReverse[i][0] = sReverse[i - 1][0] - sigma;
	}
	for (int j = 1; j <= midCol; j++)
	{
		sNormal[0][j] = sNormal[0][j - 1] - sigma;
	}
	for (int j = 1; j <= m - midCol; j++)
	{
		sReverse[0][j] = sReverse[0][j - 1] - sigma;
	}

	for (int i = 1; i <= n; i++)
	{
		for (int j = 1; j <= midCol; j++)
		{
			char vChar = v.at(i - 1);
			char wChar = w.at(j - 1);
			int match;
			if (vChar == wChar)
			{
				match = matchScore(vChar);
			}
			else
			{
				match = mismatchScore(vChar, wChar);
			}
			sNormal[i][j] = max3(sNormal[i - 1][j] - sigma, sNormal[i][j - 1] - sigma, sNormal[i - 1][j - 1] + match);
		}

		for (int j = 1; j <= m - midCol; j++)
		{
			char vChar = vRev.at(i - 1);
			char wChar = wRev.at(j - 1);
			int match;
			if (vChar == wChar)
			{
				match = matchScore(vChar);
			}
			else
			{
				match = mismatchScore(vChar, wChar);
			}
			sReverse[i][j] = max3(sReverse[i - 1][j] - sigma, sReverse[i][j - 1] - sigma, sReverse[i - 1][j - 1] + match);
			if (j == m - midCol)
			{
				if (sReverse[i][j] == sReverse[i - 1][j - 1] + match)
				{
					backtrack[i] = diagonal;
				}
				else if (sReverse[i][j] == sReverse[i - 1][j] - sigma)
				{
					backtrack[i] = down;
				}
				else if (sReverse[i][j] == sReverse[i][j - 1] - sigma)
				{
					backtrack[i] = rights;
				}
			}
		}
	}

	int posX = 0;
	int score = sNormal[0][midCol] + sReverse[n][m - midCol];
	for (int i = 1; i <= n; i++)
	{
		int temp = sNormal[i][midCol] + sReverse[n - i][m - midCol];
		if (temp > score)
		{
			score = temp;
			posX = i;
		}
	}

	int nextX, nextY;
	result = backtrack[n - posX];

	point += posX;
	return result;
}

void alignment::linearSpaceAlignment(string v, string w, int top, int bottom, int left, int right, int sigma, vector<steps>& res)
{
	if (left == right)
	{
		for (int i = 0; i < (bottom - top); i++)
		{
			res.push_back(down);
		}
		return;
	}

	if (top == bottom)
	{
		for (int i = 0; i < (right - left); i++)
		{
			res.push_back(rights);
		}
		return;
	}

	int middle = floor((left + right) / 2);
	int midNode = top;
	steps midEdge = middleEdge(v, w, top, bottom, left, right, sigma, midNode);
	linearSpaceAlignment(v, w, top, midNode, left, middle, sigma, res);
	res.push_back(midEdge);

	if (midEdge == rights || midEdge == diagonal)
	{
		middle++;
	}
	if (midEdge == down || midEdge == diagonal)
	{
		midNode++;
	}

	linearSpaceAlignment(v, w, midNode, bottom, middle, right, sigma, res);

}

vector<string> alignment::pathTostrings(string str1, string str2, vector<steps> path, int sigma, int & score)
{
	score = 0;
	vector<string> result(2, "");
	int x = 0;
	int y = 0;
	for (int i = 0; i < path.size(); i++)
	{
		switch (path[i])
		{
		case(down):
		{
			result[0] += str1.substr(x, 1);
			result[1] += "-";
			x++;
			score -= sigma;
			break;
		}
		case(rights):
		{
			result[0] += "-";
			result[1] += str2.substr(y, 1);
			y++;
			score -= sigma;
			break;
		}
		case(diagonal):
		{
			result[0] += str1.substr(x, 1);
			result[1] += str2.substr(y, 1);

			char vChar = str1.at(x);
			char wChar = str2.at(y);
			int match;
			if (vChar == wChar)
			{
				match = matchScore(vChar);
			}
			else
			{
				match = mismatchScore(vChar, wChar);
			}

			score += match;

			x++;
			y++;
			break;
		}
		}
	}
	return result;
}

int alignment::editDistance(string str1, string str2)
{
	{
		int len1 = str1.length();
		int len2 = str2.length();

		// Create a DP array to memoize result 
		// of previous computations 
		vector<vector<int>> DP(2, vector<int>(len1 + 1, 0));

		// Base condition when second string 
		// is empty then we remove all characters 
		for (int i = 0; i <= len1; i++)
			DP[0][i] = i;

		// Start filling the DP 
		// This loop run for every 
		// character in second string 
		for (int i = 1; i <= len2; i++) {
			// This loop compares the char from 
			// second string with first string 
			// characters 
			for (int j = 0; j <= len1; j++) {
				// if first string is empty then 
				// we have to perform add character 
				// operation to get second string 
				if (j == 0)
					DP[i % 2][j] = i;

				// if character from both string 
				// is same then we do not perform any 
				// operation . here i % 2 is for bound 
				// the row number. 
				else if (str1[j - 1] == str2[i - 1]) {
					DP[i % 2][j] = DP[(i - 1) % 2][j - 1];
				}

				// if character from both string is 
				// not same then we take the minimum 
				// from three specified operation 
				else {
					DP[i % 2][j] = 1 + min(DP[(i - 1) % 2][j],
						min(DP[i % 2][j - 1],
							DP[(i - 1) % 2][j - 1]));
				}
			}
		}

		// after complete fill the DP array 
		// if the len2 is even then we end 
		// up in the 0th row else we end up 
		// in the 1th row so we take len2 % 2 
		// to get row 
		return DP[len2 % 2][len1];
	}
}

void alignment::cutHairString(vector<string>& match, int& score)
{
	int size = match[0].size() - 1;
	for (int i = size; i > -1; i--)
	{
		char temp1 = match[0].at(i);
		char temp2 = match[1].at(i);
		if (temp1 == '-' || temp2 == '-')
		{
			score++;
			match[0].pop_back();
			match[1].pop_back();
		}
		else
		{
			break;
		}
	}

	for (int i = 0; i < size; i++)
	{
		char temp1 = match[0].at(i);
		char temp2 = match[1].at(i);
		if (temp1 == '-' || temp2 == '-')
		{
			score++;
			match[0].erase(match[0].begin());
			match[1].erase(match[1].begin());
		}
		else
		{
			break;
		}
	}
}


vector<vector<int>> loadScoreMatrix()
{
	int size = 20;
	vector<vector<int>> result(size, vector<int>(size, 0));
	ifstream file("scoreMatrix.txt");
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			file >> result[i][j];
		}
	}
	file.close();

	return result;
}

int maxInGraph(vector<vector<int>> graph, string& x)
{
	int n = graph.size();
	int m = graph[0].size();
	int result = 0;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			int temp = graph[i][j];
			if (temp > result)
			{
				x = to_string(i) + "," + to_string(j);
				result = temp;
			}
		}
	}

	return result;
}


void checkJump(string check, int&x, int& y)
{
	int pos = check.find(",");
	x = stoi(check.substr(0, pos));
	y = stoi(check.substr(pos + 1));
}

int max3(int a, int b, int c)
{
	return max(a, max(b, c));
}

int max4(int a, int b, int c, int d)
{
	return max(max3(a, b, c), d);
}

string reverseStr(string str)
{
	string rev = string(str.rbegin(), str.rend());
	return rev;
}

int fittingMax(vector<vector<int>> graph, string& x)
{
	int n = graph.size();
	int m = graph[0].size();
	int result = 0;
	for (int i = 0; i < n; i++)
	{
		int temp = graph[i].back();
		if (temp > result)
		{
			x = to_string(i) + "," + to_string(m-1);
			result = temp;
		}
	}

	return result;
}

void checkJumpFitting(string check, int&x, int& y)
{
	int pos = check.find(",");
	x = stoi(check.substr(0, pos));
	y = stoi(check.substr(pos + 1));
}

int overlapMax(vector<vector<int>> graph, string& x)
{
	int n = graph.size();
	int m = graph[0].size();
	int result = 0;
	for (int i = 0; i < n; i++)
	{
		int temp = graph[i].back();
		if (temp >= result)
		{
			x = to_string(i) + "," + to_string(m - 1);
			result = temp;
		}
	}

	return result;
}

int case7(int a, int b, int c, int d, int e, int f, int g, multiSteps& type)
{
	int temp1 = max3(a, b, c);
	int temp2 = max3(d, e, f);
	int result = max3(temp1, temp2, g);

	if (result == a)
	{
		type = moveX;
	}
	else if (result == b)
	{
		type = moveY;
	}
	else if (result == c)
	{
		type = moveZ;
	}
	else if (result == d)
	{
		type = diagZ;
	}
	else if (result == e)
	{
		type = diagY;
	}
	else if (result == f)
	{
		type = diagX;
	}
	else if (result == g)
	{
		type = superDiag;
	}

	return result;
}

void initMulti(vector<vector<vector<multiSteps>>>& s, int n, int m, int l)
{
	for (int i = 0; i <= n; i++)
	{
		s[i][0][0] = moveX;
	}

	for (int j = 0; j <= m; j++)
	{
		s[0][j][0] = moveY;
	}

	for (int k = 0; k <= l; k++)
	{
		s[0][0][k] = moveZ;
	}


	for (int i = 1; i <= n; i++)
	{
		for (int j = 1; j <= m; j++)
		{
			s[i][j][0] = diagZ;
		}
	}

	for (int j = 1; j <= m; j++)
	{
		for (int k = 1; k <= l; k++)
		{
			s[0][j][k] = diagX;
		}
	}

	for (int i = 1; i <= n; i++)
	{
		for (int k = 1; k <= l; k++)
		{
			s[i][0][k] = diagY;
		}
	}
}