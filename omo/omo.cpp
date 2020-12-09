#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <windows.h>
#include "Functional.h"
using namespace std;


int main()
{
	system("Color 06");
	setlocale(LC_ALL, "Ukrainian");
	int N = 7, p = 1, L = 0; 
	vector<vector<double>> Matrix; 
	vector<vector<double>> A; vector<double> b;
	Solve(Matrix,A,b,N,p,L);
	vector<vector<double>> Inv = Inverse(Matrix,N);
	vector<vector<double>> Identity=MultiplyMatrix(A,Inv,N);
	cout << "Добуток заданої матриці та оберненої до неї:" << endl;
	for (int i(0); i < N; i++) {
		for (int j(0); j < N; j++)
			cout << fixed << setprecision(10) << setw(20) << Identity[i][j] << setw(20);
		cout << endl;
	}
	cout << endl;
	cout << "Число обумовленості: " << endl;
	cout << Cond(A,Inv,N);
}
