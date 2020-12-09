#pragma once
#include <iostream>
#include <vector>
using namespace std;


void FillMatrix(vector<double>& b, vector<vector<double>>& A, const int& N, const int& p) {
	A.resize(N); b.resize(N);
	for (int i(0); i < N; i++) {
		A[i].resize(N);
	}
	for (int i(0); i < N; i++) {
		for (int j(0); j < N; j++) {
			A[i][j] = 35 / pow((5 * (i + 1) * (j + 1) + pow((i + 1) * (j + 1), p + 1)), p - 0.5);
			
		}
		b[i] = N - i;
	}
	for (int i(0); i < N; i++) {
		A[i][i] += 10;
	}
}
void OutputMatrix(vector<vector<double>>& Matrix, const int& N) {
	for (int i(0); i < N; i++) {
		for (int j(0); j < N; j++)
			cout << fixed << setprecision(10) << Matrix[i][j] << "\t";
		cout << endl;
	}
	cout << endl;
}

double Norm(vector<double> x, vector<double> y, const int& N) {
	double res = 0,temp;
	for (int i(0); i < N; i++) {
		temp = abs(y[i] - x[i]);
		if (temp > res)res = temp;
	}
	return res;
}

vector<double> SOR(vector<vector<double>> A, vector<double> b,const int& N, double& w, double& toler,int& count) {
	vector<double> x_prev, x_next; x_prev.resize(N); x_next.resize(N);
	double sum1 = 0, sum2 = 0; int counter(0);
	while (true) {
		counter++;
		for (int i(0); i < N; i++) {
			sum1 = sum2 = 0;
			for (int j(0); j < i; j++) {
				sum1 += A[i][j] * x_next[j];
			}
			for (int j(i+1); j < N; j++) {
				sum2 += A[i][j] * x_prev[j];
			}
			x_next[i] = (1-w)*x_prev[i] + (w / A[i][i]) * (b[i] - sum1 - sum2);
		}
		if (Norm(x_next, x_prev, N) <= toler) { count = counter; break; }
		x_prev = x_next;
	}
	cout << endl;
	cout << "Розв'язок: " << endl;
	for (int i(0); i < N; i++) {
		cout << x_next[i] << endl;
	}
	return x_next;
	
}


vector<double> MultiplyMatrixVector(vector<vector<double>> Matrix, vector<double> x, const int& N) {
	vector<double> res; res.resize(N);
	for (int i(0); i < N; i++) {
		for (int j(0); j < N; j++) {
			res[i] += x[j] * Matrix[i][j];
		}
	}
	return res;
}

vector<double> Substruct(vector<double> x, vector<double> y, const int& N) {
	vector<double> res; res.resize(N);
	for (int i(0); i < N; i++) {
		res[i] = x[i] - y[i];
	}
	return res;
}