#pragma once
using namespace std;
#include <vector>

void FillMatrix1(vector<double>&b, vector<vector<double>>& A,vector<vector<double>> &Matrix, const int& N, const int& p) {
	Matrix.resize(N);
	A.resize(N); b.resize(N);
	for (int i(0); i < N; i++) {
		Matrix[i].resize(N + 1);
		A[i].resize(N);
	}
	for (int i(0); i < N; i++) {
		for (int j(0); j < N; j++) {
			Matrix[i][j] = 35 / pow((5 * (i + 1) * (j + 1) + pow((i + 1) * (j + 1), p + 1)), p - 0.5);
			A[i][j] = Matrix[i][j];
		}
		Matrix[i][N] = N - i; b[i] = Matrix[i][N];
	}
}

void OutputMatrix(vector<vector<double>>&Matrix, const int& N) {
	for (int i(0); i < N; i++) {
		for (int j(0); j < N+1; j++)
			cout << fixed << setprecision(10) << Matrix[i][j] << "\t";
		cout << endl;
	}
	cout << endl;
}

vector<vector<double>> ForwardElimination(vector<vector<double>> Matrix, const int& N, int&L) {
	for (int j = 0; j < N - 1; j++)
	{
		for (int i = j + 1; i < N; i++)
		{
			double temp = Matrix[i][j] / Matrix[j][j];

			for (int k = 0; k < N + 1; k++)
				Matrix[i][k] -= Matrix[j][k] * temp;
		}
		double max = Matrix[j + 1][j + 1]; int index = j + 1;
		for (int k = j + 1; k < N; k++) {
			if (Matrix[k][j + 1] > max) { max = Matrix[k][j + 1]; index = k; }
		}
		//OutputMatrix(Matrix, N);
		vector<double> temp; temp.resize(N + 1);
		if (index != j + 1) {
			L++;
			for (int i(0); i < N + 1; i++) {
				temp[i] = Matrix[index][i];
			}
			for (int i(0); i < N + 1; i++) {
				Matrix[index][i] = Matrix[j+1][i];
			}
			for (int i(0); i < N + 1; i++) {
				Matrix[j + 1][i] = temp[i];
			}
		}
		
		//OutputMatrix(Matrix, N);
	}
	//OutputMatrix(Matrix,N);
	return Matrix;
}

vector<double> BackSubstitution(vector<vector<double>> Matrix, const int& N) {
	vector<double> x; x.resize(N);
	for (int i = N - 1; i >= 0; i--)
	{
		double s = 0;
		for (int j = i + 1; j < N; j++)
			s += Matrix[i][j] * x[j];
		x[i] = (Matrix[i][N] - s) / Matrix[i][i];
	}
	return x;
}

vector<double> MultiplyMatrixVector(vector<vector<double>> Matrix, vector<double> x,const int& N) {
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

double det(vector<vector<double>>&Matrix, const int& N, const int& L) {
	double res = 1;
	for (int i(0); i < N; i++) {
		res *= Matrix[i][i];
	}
	res *= pow(-1, L);
	return res;
}

void Solve(vector<vector<double>>& Matrix, vector<vector<double>>& A, vector<double>& b, const int&N,const int&p, int&L) {
	vector<double> x; x.resize(N);
	FillMatrix1(b, A, Matrix, N, p);
	OutputMatrix(Matrix, N);
	vector<vector<double>> UpperTriangle = ForwardElimination(Matrix, N, L);
	x = BackSubstitution(UpperTriangle, N);
	cout << endl << endl;
	cout << "Розв'язок: " << std::endl;
	for (int i(0); i < N; i++) {
		cout << setw(20) <<x[i] << endl;
	}
	vector<double> temp = MultiplyMatrixVector(A, x, N);
	vector<double> r = Substruct(temp, b, N);
	cout << endl << endl;
	cout << "Нев'язка: " << std::endl;
	for (int i(0); i < N; i++) {
		cout << scientific << r[i] << endl;
	}
	cout << endl << endl;
	double D = det(UpperTriangle, N, L);
	cout << "Визначник: " << std::endl;
	cout << scientific << D << endl << endl;
}

vector<vector<double>> Inverse(vector<vector<double>> A, const int& N) {
	vector<vector<double>> Inv,Temp; vector<double> b;
	Inv.resize(N); Temp.resize(N);
	b.resize(N);
	int omg(0);
	for (int i(0); i < N; i++) {
		Inv[i].resize(N);
		Temp[i].resize(N + 1);
		Temp[i] = A[i];
	}
	for (int g(0); g < N; g++) {
		for (int i(0); i < N; i++) {
			Temp[i] = A[i];
		}
		if(g-1>=0)
		b[g - 1] = 0;
		b[g] = 1;
		for (int h(0); h < N; h++) {
			Temp[h][N] = b[h];
		}
		Temp = ForwardElimination(Temp, N, omg);
		vector<double> x = BackSubstitution(Temp,N);
		for (int i(0); i < N; i++) {
			Inv[i][g] = x[i];
		}
	}
	cout << endl;
	cout << "Оберена матриця :" << endl;
	for (int i(0); i < N; i++) {
		for (int j(0); j < N; j++)
			cout << fixed << setprecision(6)<< setw(20) << Inv[i][j] << setw(20);
		cout << endl;
	}
	cout << endl;
	return Inv;
}



vector<vector<double>> MultiplyMatrix(vector<vector<double>>&A, vector<vector<double>>&B, const int& N) {
	vector<vector<double>> H;
	H.resize(N); 
	for (int i(0); i < N; i++) {
		H[i].resize(N);
	}
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			H[i][j] = 0;
			for (int t = 0; t < N; t++) {
				H[i][j] += A[i][t] * B[t][j];
			}
		}
	}
	return H;
}

double MaxSumByRow(vector<vector<double>>& A, const int& N) {
	double res = 0,temp;
	for (int i(0); i < N; i++) {
		temp = 0;
		for (int j(0); j < N; j++) {
			temp += abs(A[i][j]);
		}
		if (temp > res)
			res = temp;

	}
	return res;
}

double Cond(vector<vector<double>>& A, vector<vector<double>>& B, const int& N) {

	return MaxSumByRow(A, N) * MaxSumByRow(B, N);
}