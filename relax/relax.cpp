#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <windows.h>
#include "Functional.h"
#include <iostream>

int main()
{
	system("Color 06");
	setlocale(LC_ALL, "Ukrainian");
	int N = 7, p = 1;
	vector<vector<double>> A; vector<double> b;
	FillMatrix(b,A,N,p);
	OutputMatrix(A,N);
	for (int i(0); i < N; i++)
		cout << b[i] << endl;
	cout << endl;
	double w = 0.8,toler=0.0001;
	int count(0);
	vector<double>x=SOR(A,b,N,w,toler,count);
	cout << endl;
	cout<<"Кiлькiсть iтерацiй:"<<endl;
	cout << count;
	vector<double> temp = MultiplyMatrixVector(A, x, N);
	vector<double> r = Substruct(temp, b, N);
	cout << endl << endl;
	cout << "Нев'язка: " << std::endl;
	for (int i(0); i < N; i++) {
		cout << scientific << r[i] << endl;
	}
}


