#include<iostream>
#include "msg.h"
#include "agrozUtil.h"

int main()
{
	msg mensaje("prueba.txt");
	vector<string> lala = mensaje.frequentWordsMRC(10, 3);
	printVector(lala);
	return 0;
}