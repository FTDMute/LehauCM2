#include <iostream>
#include <conio.h>
#include <stdlib.h>
#include <cmath>
#include <math.h>
#include <vector>
#include <cmath>
#include <fstream>
#include <random>
#include "Source.h"
using namespace std;
double Function(double x) {
	return pow(x, 3) * sqrt(1 - pow(x, 2));
}
double FunctionNI(double x,double a = 0, double b = 1) {
	double t = (b + a) / 2 + x * (b - a) / 2;
	return Function(t);
}
double d2F(double x) {
	double res = (12 * pow(x,5) - 19 * pow(x,3) + 6 * x) / (sqrt(1 - pow(x,2)) * (1 - pow(x,2)));
	return res;
}
void FileImput(vector <vector<double>> x)
{
	string s;
	cin >> s;
	ofstream FileS(s, ios_base::trunc);
	for (int i = 0; i < x[1].size(); i++)
	{
		FileS << x[0][i] << "  " << x[1][i] << endl;
	}
}
void MakeScript(int x)
{
	string s;
	cin >> s;
	ofstream FileS(s, ios_base::trunc);
	for (int i = 0; i < x; i++)
	{
		FileS << " set xrange [0:1] set ticslevel 0plot " << s << " u 1:2 w p " << endl;
	}
}
vector <vector<double>> MakeGrid(int number, double start, double end) {
	vector <vector<double>> grid(2);
	double step = (end - start) / number;
	for (int i = 0; i < number; i++) {
		grid[0].push_back(start + step);
		grid[1].push_back(Function(grid[i][0]));
	}
	return grid;
}
vector <vector<double>> Rectangular(int number) {
	vector <double> x;
	vector <vector<double>> errorForce(2);
	double answer = 0;
	for (int j = 10; j <= number; j += 10) {
		double step = (1. / j);
		x.push_back(0);
		answer = 0;
		for (int i = 1; i < j; i++) {
			x.push_back(0 + i*step);
			answer += (Function((x[i - 1] + x[i]) / 2.) * step);
		}
		errorForce[0].push_back(j);
		errorForce[1].push_back((2. / 15.) - answer);
		x.clear();
	}
	return errorForce;
}
vector <vector<double>> Trapeze(int number) {
	vector <double> x;
	vector <vector<double>> errorForce(2);
	vector <vector<double>> trapezeResults(2);
	double answer;
	for (int j = 10; j <= number; j += 10) {
		double step = (1. / j);
		x.push_back(0);
		answer = 0;
		for (int i = 1; i < j; i++) {
			x.push_back(0 + i * step);
			answer += step*((Function(x[i-1]) + Function(x[i]))/2);
		}
		errorForce[0].push_back(j);
		errorForce[1].push_back((2. / 15.) - answer);
		x.clear();
	}
	cout << "imhere" << endl;
	return errorForce;
	
}
double TrapezeEval(int number) {
	double step = 1.0 / number;
	vector <double> x;
	double answer=0;
	x.push_back(0);
	for (double j = step; j <= 1; j += step) {
		x.push_back(0 + j);
	}
	cout << x.size() << " problem?" << endl;
	for (int i = 1; i < x.size(); i++) {
		answer += step * ((Function(x[i - 1]) + Function(x[i])) / 2);
	}
	cout << "--------------------------" << endl;
	return answer;
}
double TrapezeEvalForOptT(vector <double> x, int rn) {
	double step = 1./rn;
	double answer = 0.;
	for (int i = 1; i < x.size(); i++) {
		answer += (x[i]-x[i-1]) * ((Function(x[i - 1]) + Function(x[i])) / 2);
	}
	return answer;
}
vector <vector<double>> Simpson(int number) {
	vector <double> x;
	vector <vector<double>> errorForce(2);
	double answer = 0;
	for (int j = 10; j <= number; j += 10) {
		double step = (1. / j);
		x.push_back(0);
		answer = (0);
		for (int i = 1; i < j; i++) {
			x.push_back(0 + i * step);
			answer += step*(Function(x[i-1])+4*Function((x[i-1]+x[i])/2)+Function(x[i]))/6;
		}
		errorForce[0].push_back(j);
		errorForce[1].push_back((2. / 15.) - answer);
		x.clear();
	}
	return errorForce;
}
double EitkenT(int number,double q) {
	vector <double> x; 
	x.push_back(0);
	vector <double> answers_for_T_3;
	double answer;
	double step = 1./number;
	int k = 1./step;
	for (double i = 0; i < k; i++) {
			x.push_back(0 + i * step);
	}
	for (int j = 0; j < 3; j++) {
		answer = 0;
		for (double i = 1; i < x.size(); i++) {
			answer += step * ((Function(x[i-1]) + Function(x[i])) / 2.);
		}
		answers_for_T_3.push_back(answer);
		step = step * q;
		k = 1. / step;
		x.clear();
		x.push_back(0);
		for (double i = 0; i < k; i++) {
			x.push_back(0 + i * step);
		}
	}
	double p = (log((answers_for_T_3[2] - answers_for_T_3[1]) / (answers_for_T_3[1] - answers_for_T_3[0]))) / (log(q));
	return p;
}
double Runge(int number,double p) {
	vector <double> x;
	x.push_back(0);
	vector <double> answers_for_T_3;
	double answer;
	double result;
	double step = 1. / number;
	int ki = 1. / step;
	for (double i = 0; i < ki; i++) {
		x.push_back(0 + i * step);
	}
	double k;
	for (int j = 0; j < 2; j++) {
		answer = 0;
		for (int i = 0; i < number; i ++) {
			cout << "im here" << endl;
			x.push_back(0 + (i+1) * step);
			answer += step * ((Function(x[i]) + Function(x[i + 1])) / 2.);
		}
		step = step * 0.5;
		answers_for_T_3.push_back(answer);
		step = 1. / number;
		ki = 1. / step;
		x.clear();
		x.push_back(0);
		for (double i = 0; i < ki; i++) {
			x.push_back(0 + i * step);
		}
	}
	k = step / (step*2);
	double delta = pow(k, p) / (pow(k, p) - 1);
	result = delta * answers_for_T_3[0] + (1 - delta) * answers_for_T_3[1];
	return result;
}
double Romberg(int number) {
	int n = 5;
	int k = 1;
	vector<vector<double>> results(n);
	vector<vector<double>> res(2);
	for (int i = 0; i < n; i++) {
		results[i].resize(i + 1);
	}
	for (int i = 0; i < n; i++) {
		cout << k << endl;
		k += 1;
		results[i][0] = TrapezeEval(pow(2,i)*number);
	}
	for (int i = 1; i < n; i++) {
		for (int j = 1; j <= i; j++) {
			if (j == n - 1)
				results[i][j] = results[i][j - 1] + (results[i][j - 1] - results[i - 1][j - 1]) / (pow(2, j) - 1);
			else 
				results[i][j] = results[i][j - 1] + (results[i][j - 1] - results[i - 1][j - 1]) / (pow(2, j+1) - 1);
		}
	}
	return results[n - 1][n - 1];
}
double Lejandr(int count) {
	switch (count)
	{
	case 0:
		return 2.0 * FunctionNI(0);
	case 1:
		return FunctionNI(-0.5773502692) + FunctionNI(0.5773502692);
	case 2:
		return 0.5555555556 * FunctionNI(-0.7745966692) + 0.8888888888 * FunctionNI(0) + 0.5555555556 * FunctionNI(0.7745966692);
	case 3:
		return 0.3478547451 * FunctionNI(-0.8611363115) + 0.6521451549 * FunctionNI(-0.3399810436) + 0.6521451549 * FunctionNI(0.3399810436) + 0.3478547451 * FunctionNI(0.8611363115);
	case 4:
		return 0.2369268851 * FunctionNI(-0.9061798459) + 0.4786286705 * FunctionNI(-0.5384693101) + 0.5688888888 * FunctionNI(0) + 0.2369268851 * FunctionNI(0.9061798459) + 0.4786286705 * FunctionNI(0.5384693101);
	case 5:
		return 0.1713244924 * FunctionNI(-0.9324695142) + 0.3607615730 * FunctionNI(-0.6612093684) + 0.4679139346 * FunctionNI(-0.2386191861) + 0.1713244924 * FunctionNI(0.9324695142) + 0.3607615730 * FunctionNI(0.6612093684) + 0.4679139346 * FunctionNI(0.2386191861);
	case 6:
		return 0.12948496 * (FunctionNI(0.94910791) + FunctionNI(-0.94910791)) + 0.27970540 * (FunctionNI(0.74153119) + FunctionNI(-0.74153119)) + 0.38183006 * (FunctionNI(0.40584515) + FunctionNI(-0.40584515)) + 0.41795918 * FunctionNI(0);
	case 7:
		return 0.10122854 * (FunctionNI(0.96028986) + FunctionNI(-0.96028986)) + 0.22238104 * (FunctionNI(0.79666648) + FunctionNI(-0.79666648)) + 0.31370664 * (FunctionNI(0.52553242) + FunctionNI(-0.52553242)) + 0.36268378 * (FunctionNI(0.18343464) + FunctionNI(-0.18343464));
	}
}
vector <vector<double>> Gauss(int number) {
	vector <double> x;
	vector <vector<double>> errorForce(2);
	double answer = 0;
	double step = 1. / number;
	for (double i = 0; i < 1; i += step){
		x.push_back(i);
	}
	for (double i = 0; i < number; i ++) {
		answer = Lejandr(i)/2;
		cout << abs((2. / 15.) - answer) << endl;
		errorForce[0].push_back(i);
		errorForce[1].push_back(abs((2. / 15.) - answer));
	}
	cout << "imhere" << endl;
	return errorForce;
}
double A(double a, double b) {
	double h = (b - a) / 100; //10000
	double res = 0.;
	for (int i = 0; i < 100; i++) {;
		res = max(res, abs(d2F(a + i * h)));
	}
	return res;
}
double opt_trapezoid(int n, int rn) {
	double step = 1. / rn;
	double Asum = 0.0;
	int trn = 0;
	vector<double> vecA;
	double res = 0.0;
	for (int i = 0; i < rn; i++) {
		double c = cbrt(A(i * step, (i + 1) * step));
		vecA.push_back(c);
		Asum += c;
	}
	for (int i = 0; i < rn; i++) {
		vector<double> X;
		int nl = int(vecA[i] * n / Asum);
		trn += nl;
		if (nl == 0)  return 0;
		double h = step / nl;
		for (int j = 0; j <= nl; j++) {
			X.push_back(step * i + h * j);
		}
		ofstream FileS("OptTrap.dat", ios_base::trunc);
		for (int i = 0; i < X.size(); i++)
		{
			FileS << X[i] << "  " << 1 << endl;
		}
		res += TrapezeEvalForOptT(X,rn);
	}
	return res;
}
double MonteCarlo(int n)
{
	double a, b, result;
	a = 0, b = 1;
	auto gen = default_random_engine();
	auto dis = uniform_real_distribution<>(a,b);
	double S = 0.;
	for (int i = 0; i < n; i++)
	{
		S += Function(dis(gen));
	}
	result = ((b-a) / (double)n) * S;
	return result;
}
void Task1()
{
	int number;
	cin >> number;
	double start;
	cin >> start;
	double end;
	cin >> end;
	vector <vector<double>> error;
	for (int g = 0; g < 2; g++) {
		if (g == 0) {
			error = Rectangular(number);
		}
		if (g==1) {
			error = Trapeze(number);
		}

		for (int i = 0; i < error[1].size(); i++) {
			if (error[1][i] < 1e-5) {
				cout << "needed error claimed on " << error[0][i] << "  Pogreshnost  " << error[1][i] << endl;
				break;
			}
		}
		FileImput(error);
		error.clear();
	}
}
void Task2()
{
	int number;
	cin >> number;
	double start;
	cin >> start;
	double end;
	cin >> end;
	vector <vector<double>> error;
	error = Simpson(number);
	for (int i = 0; i < error[1].size(); i++) {
		if (error[1][i] < 1e-5) {
			cout << "needed error claimed on " << error[0][i] << "  Pogreshnost  " << error[1][i] << endl;
			break;
		}
	}
	FileImput(error);
}
void Task3()
{
	int number_forRomberg;
	cin >> number_forRomberg;
	int number;
	cin >> number;
	double start;
	cin >> start;
	double end;
	cin >> end;
	double result;
	double p = EitkenT(number, 0.5);
	cout << Runge(number, p) << endl;
	result = Romberg(number_forRomberg);
	cout << result << endl;
}
void Task4()
{
	int number;
	cin >> number;
	vector <vector<double>> error;
	error = Gauss(number);
	for (int i = 0; i < error[1].size(); i++) {
		if (abs(error[1][i]) < 2e-3) {
			cout << "needed error claimed on " << error[0][i] << "  Pogreshnost  " << error[1][i] << endl;
			break;
		}
	}
	FileImput(error);
}
void Task5()
{
	int n = 78;
	int q = 2;
	int p = 0;
	int s;
	cout << "show while? (1-Y, 0-N) " << endl;
	cin >> s;
	switch (s) {
	case 1: {
		n = 500;
		while (q < n) {
			int i = n;
			double error = 0;
			while (error <= 1e-2 && i > q) {
				i--;
				error = abs(2. / 15. - opt_trapezoid(i, q));
				p = q;
				cout << i << "    " << q << "    " << error << "    " << endl;
			}
			n = i + 1;
			q++;
		}
	}
	case 0: {
		double error = abs(2. / 15. - opt_trapezoid(77, 2));
		cout << n << "    " << q << " " << error << endl;
	}
	}
}
void Task6()
{
	int n;
	cin >> n;
	double m = MonteCarlo(n);
	cout << m << endl;
}
int main() {
	int Task;
	cout << "choose task" << endl;
	cin >> Task;
	switch (Task) {
	case 1: {
		Task1();
		break;
	}
	case 2: {
		Task2();
		break;
	}
	case 3: {
		Task3();
		break;
	}
	case 4: {
		Task4();
		break;
	}
	case 5: {
		Task5();
	}
	case 6: {
		Task6();
	}
	}
}