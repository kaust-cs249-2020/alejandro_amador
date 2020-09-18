#include "motif.h"
#include "agrozUtil.h"

int main()
{
	motif test("test.txt");
	string consen = test.consensus(12, MEDIAN_STRING);
	cout << consen << endl;
	return 0;
}