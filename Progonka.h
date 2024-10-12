#pragma once
#include <vector>
#include <string>
#include <iostream>

using namespace std;

class Progonka
{
private:
	int n;
	vector<double> a;
	vector<double> b;
	vector<double> c;
	vector<double> d;
	vector<double> e;
	vector<double> f;
	vector<double> alpha;
	vector<double> beta;
	vector<double> gama;
	vector<double> x;

public:

	Progonka(int n_);// дл проверки
	Progonka(vector<double>& a_, vector<double>& b_, vector<double>& c_, vector<double>& f_);
	Progonka(vector<double>& a_, vector<double>& b_, vector<double>& c_, vector<double>& d_, vector<double>& e_, vector<double>& f_);
	vector<double> SolveThreeDiagonal();
	vector<double> SolveFiveDiagonal();
};
