// Master's work.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//
#include <iostream>
#include <thread>
#include <fstream>
#include <vector>
#include <windows.h>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include <random>
#include <chrono>
#include <algorithm>
#include <map>
#include <exception>
#include <string>
#include <math.h>
#include "GnuPlot.h"
#include "Splayn.h"

using namespace std;

double pi = 3.141592654;

enum Metal { Au, Al, Cu, Ni };
enum TypeBeam { Gauss, Vortex };

template <typename T>
T convertString(const string& s) {
	return T{ s };
}

template<> 
double convertString<double>(const string& s) {
	return stod(s);
}

template <typename T> 
map<string, T> parseFile(string filename) {
	fstream file(filename);
	string tmp_line;
	map<string, T> result;

	while (file >> tmp_line) {
		if (tmp_line.front() == ';')
			continue;

		std::string::size_type equal_sign_at = tmp_line.find("=");
		if (equal_sign_at == std::string::npos)
			throw runtime_error("No '=' sign in a string");

		string key = tmp_line.substr(0, equal_sign_at);
		string value = tmp_line.substr(equal_sign_at + 1, tmp_line.length());
		result[key] = convertString<T>(value);
	}

	file.close();

	return result;
}

struct Parametrs {
	double t0;
	double r0;
	double kabs;
	double P0;
	double T00;
	double beta;
};

struct Melting
{
	Metal mt;
	double Q_fusion;
	double T_melting;
	double Density;
	double DensityLiquid;
	double gi;
	double ge;
	double u0;
	double u0_Liquid;
};

struct Interval
{
	Interval(int begin, int end) : begin{ begin }, end{ end } {}
	int begin;
	int end;
};

string ConvertNumToStringdouble(double i)
{
	string tmp_num_line;

	int length = snprintf(NULL, 0, "%lf", i);
	char* str = new char[length + 1];
	snprintf(str, length + 1, "%lf", i);
	for (int j = 0; j < length; j++)
	{
		tmp_num_line.push_back(str[j]); // номер столбца
	}

	return tmp_num_line;
}

struct Point3D // структура необходимая для определение точек, где произошел разрыв
{
	Point3D(int index_x, int index_y, int index_z) : index_x{ index_x }, index_y{ index_y }, index_z{ index_z } {}
	int index_x;
	int index_y;
	int index_z;
};

bool operator ==(const Point3D& p1, const Point3D& p2)
{
	return p1.index_x == p2.index_x && p1.index_y == p2.index_y && p1.index_z == p2.index_z;
}

bool comp_x(const Point3D& pt1, const Point3D& pt2)
{
	return (pt1.index_x < pt2.index_x);
}

bool comp_y(const Point3D& pt1, const Point3D& pt2)
{
	return (pt1.index_y < pt2.index_y);
}

bool comp_z(const Point3D& pt1, const Point3D& pt2)
{
	return (pt1.index_z < pt2.index_z);
}

void MySort_Point3D_y(vector<Point3D>& v)
{
	vector<Point3D> copy_v, v_for_partial_sort;
	vector<Point3D>::iterator it;
	copy_v = v;
	v.erase(v.begin(), v.end());
	v.clear();
	vector<Point3D>::iterator begin_pos, end_pos, temp;
	begin_pos = copy_v.begin();
	int count = 0;
	for (it = copy_v.begin(); it != copy_v.end(); it++)// цикл по точкам с которым сравниваем
	{
		if (((*it).index_z == (*(it + 1)).index_z) && (it != copy_v.end()))
		{
			count++;
		}
		else
		{
			end_pos = it;
			Point3D tmp(0, 0, 0);
			for (temp = begin_pos; temp <= end_pos; temp++)
			{
				v_for_partial_sort.push_back(tmp);
			}
			copy(begin_pos, end_pos + 1, v_for_partial_sort.begin());
			sort(v_for_partial_sort.begin(), v_for_partial_sort.end(), comp_y);
			copy(v_for_partial_sort.begin(), v_for_partial_sort.end(), v.end());

			for (temp = v_for_partial_sort.begin(); temp != v_for_partial_sort.end(); temp++)
			{
				v.push_back((*temp));
			}

			begin_pos = (it + 1);
			count = 0;
			v_for_partial_sort.clear();
			v_for_partial_sort.erase(v_for_partial_sort.begin(), v_for_partial_sort.end());
		}
	}
}

void MySort_Point3D_z(vector<Point3D> & v)
{
	vector<Point3D> copy_v, v_for_partial_sort;
	vector<Point3D>::iterator it;
	copy_v = v;
	v.clear();
	v.erase(v.begin(), v.end());
	vector<Point3D>::iterator begin_pos, end_pos, temp;
	begin_pos = copy_v.begin();
	int count = 0;
	for (it = copy_v.begin(); it != copy_v.end(); it++)// цикл по точкам с которым сравниваем
	{
		if (((*it).index_x == (*(it + 1)).index_x) && ((*it).index_y == (*(it + 1)).index_y) && (it != copy_v.end()))
		{
			count++;
		}
		else
		{
			end_pos = it;
			Point3D tmp(0, 0, 0);
			for (temp = begin_pos; temp <= end_pos; temp++)
			{
				v_for_partial_sort.push_back(tmp);
			}
			copy(begin_pos, end_pos + 1, v_for_partial_sort.begin());
			sort(v_for_partial_sort.begin(), v_for_partial_sort.end(), comp_z);
			copy(v_for_partial_sort.begin(), v_for_partial_sort.end(), v.end());

			for (temp = v_for_partial_sort.begin(); temp != v_for_partial_sort.end(); temp++)
			{
				v.push_back((*temp));
			}

			begin_pos = (it + 1);
			count = 0;
			v_for_partial_sort.clear();
			v_for_partial_sort.erase(v_for_partial_sort.begin(), v_for_partial_sort.end());
		}
	}
}

void MySort_Point3D_x(vector<Point3D> & v)
{
	vector<Point3D> copy_v, v_for_partial_sort;
	vector<Point3D>::iterator it;
	copy_v = v;
	v.clear();
	v.erase(v.begin(), v.end());
	vector<Point3D>::iterator begin_pos, end_pos, temp;
	begin_pos = copy_v.begin();
	int count = 0;
	for (it = copy_v.begin(); it != copy_v.end(); it++)// цикл по точкам с которым сравниваем
	{
		if (((*it).index_z == (*(it + 1)).index_z) && ((*it).index_y == (*(it + 1)).index_y) && (it != copy_v.end()))
		{
			count++;
		}
		else
		{
			end_pos = it;
			Point3D tmp(0, 0, 0);
			for (temp = begin_pos; temp <= end_pos; temp++)
			{
				v_for_partial_sort.push_back(tmp);
			}
			copy(begin_pos, end_pos + 1, v_for_partial_sort.begin());
			sort(v_for_partial_sort.begin(), v_for_partial_sort.end(), comp_x);
			copy(v_for_partial_sort.begin(), v_for_partial_sort.end(), v.end());

			for (temp = v_for_partial_sort.begin(); temp != v_for_partial_sort.end(); temp++)
			{
				v.push_back((*temp));
			}

			begin_pos = (it + 1);
			count = 0;
			v_for_partial_sort.clear();
			v_for_partial_sort.erase(v_for_partial_sort.begin(), v_for_partial_sort.end());
		}
	}
}

void My_unique(vector<Point3D> & v)
{
	vector<Point3D>::iterator it1_comp = v.begin();
	vector<Point3D>::iterator it2;

	for (it1_comp = v.begin(); it1_comp != v.end(); it1_comp++)// цикл по точкам с которым сравниваем
	{
		for (it2 = it1_comp + 1; it2 != v.end(); it2++)// перебор по всем точкам остальным
		{
			if ((*it1_comp) == (*it2))
			{
				v.erase(it2);
				it1_comp--;
				it2--;
			}
		}
	}
}

void ReadFile(vector<string> nametxtfile, double*** tmpe0, double*** tmpe1, double*** tmpe2, double*** tmpi0, double*** tmpi1, double*** tmpi2, double*** e1, double*** e2, double*** a1x, double*** a2x, double*** a1y, double*** a2y, double*** a1z, double*** a2z, double*** b1x, double*** b2x, double*** b1y, double*** b2y, double*** b1z, double*** b2z, double*** V2, int& Nx, int& Ny, int& Nz)
{
	vector<ifstream> file;
	file.reserve(nametxtfile.size());
	for (int i = 0; i < nametxtfile.size(); i++)
	{
		file.emplace_back(ifstream{ nametxtfile[i] }); // создаие потока вывода в файл и самого файла дл записи туда данных и ео открытие file.open(filename);
	}

	for (int numfile = 0; numfile < nametxtfile.size(); numfile++)
	{
		for (int i = 0; i < Nx; i++)
		{
			for (int j = 0; j < Ny; j++)
			{
				for (int k = 0; k < Nz; k++)
				{
					switch (numfile)
					{
					case 0:
						file[numfile] >> tmpe0[i][j][k];
						break;
					case 1:
						file[numfile] >> tmpe1[i][j][k];
						break;
					case 2:
						file[numfile] >> tmpe2[i][j][k];
						break;
					case 3:
						file[numfile] >> tmpi0[i][j][k];
						break;
					case 4:
						file[numfile] >> tmpi1[i][j][k];
						break;
					case 5:
						file[numfile] >> tmpi2[i][j][k];
						break;
					case 6:
						file[numfile] >> e1[i][j][k];
						break;
					case 7:
						file[numfile] >> e2[i][j][k];
						break;
					case 8:
						file[numfile] >> a1x[i][j][k];
						break;
					case 9:
						file[numfile] >> a2x[i][j][k];
						break;
					case 10:
						file[numfile] >> a1y[i][j][k];
						break;
					case 11:
						file[numfile] >> a2y[i][j][k];
						break;
					case 12:
						file[numfile] >> a1z[i][j][k];
						break;
					case 13:
						file[numfile] >> a2z[i][j][k];
						break;
					case 14:
						file[numfile] >> b1x[i][j][k];
						break;
					case 15:
						file[numfile] >> b2x[i][j][k];
						break;
					case 16:
						file[numfile] >> b1y[i][j][k];
						break;
					case 17:
						file[numfile] >> b2y[i][j][k];
						break;
					case 18:
						file[numfile] >> b1z[i][j][k];
						break;
					case 19:
						file[numfile] >> b2z[i][j][k];
						break;
					case 20:
						file[numfile] >> V2[i][j][k];
						break;
					}
				}
			}
		}
	}
}

void WriteToFile(vector<string> nametxtfile, double*** tmpe0, double*** tmpe1, double*** tmpe2, double*** tmpi0, double*** tmpi1, double*** tmpi2, double*** e1, double*** e2, double*** a1x, double*** a2x, double*** a1y, double*** a2y, double*** a1z, double*** a2z, double*** b1x, double*** b2x, double*** b1y, double*** b2y, double*** b1z, double*** b2z, double*** V2, int& Nx, int& Ny, int& Nz)
{
	vector<FILE*> gnuplotPipe;
	vector<ofstream> file;
	//vector<string> filename;
	//string filename;
	FILE* gnuplotPipe_tmp;

	file.reserve(nametxtfile.size());

	for (int i = 0; i < nametxtfile.size(); i++)
	{
		//string file_name = nametxtfile[i];
	//	filename.push_back(file_name);
		file.emplace_back(ofstream{ nametxtfile[i] }); // создаие потока вывода в файл и самого файла дл записи туда данных и ео открытие file.open(filename);
		gnuplotPipe_tmp = new FILE; // создание файла дл гнуплота чтобы записывать туда комады дл гнуплота
		gnuplotPipe.push_back(gnuplotPipe_tmp); // вектор файлов
	}

	for (int numfile = 0; numfile < nametxtfile.size(); numfile++)
	{
		for (int i = 0; i < Nx; i++)
		{
			for (int j = 0; j < Ny; j++)
			{
				for (int k = 0; k < Nz; k++)
				{
					switch (numfile)
					{
					case 0:
						file[numfile] << tmpe0[i][j][k] << "  ";
						break;
					case 1:
						file[numfile] << tmpe1[i][j][k] << "  ";
						break;
					case 2:
						file[numfile] << tmpe2[i][j][k] << "  ";
						break;
					case 3:
						file[numfile] << tmpi0[i][j][k] << "  ";
						break;
					case 4:
						file[numfile] << tmpi1[i][j][k] << "  ";
						break;
					case 5:
						file[numfile] << tmpi2[i][j][k] << "  ";
						break;
					case 6:
						file[numfile] << e1[i][j][k] << "  ";
						break;
					case 7:
						file[numfile] << e2[i][j][k] << "  ";
						break;
					case 8:
						file[numfile] << a1x[i][j][k] << "  ";
						break;
					case 9:
						file[numfile] << a2x[i][j][k] << "  ";
						break;
					case 10:
						file[numfile] << a1y[i][j][k] << "  ";
						break;
					case 11:
						file[numfile] << a2y[i][j][k] << "  ";
						break;
					case 12:
						file[numfile] << a1z[i][j][k] << "  ";
						break;
					case 13:
						file[numfile] << a2z[i][j][k] << "  ";
						break;
					case 14:
						file[numfile] << b1x[i][j][k] << "  ";
						break;
					case 15:
						file[numfile] << b2x[i][j][k] << "  ";
						break;
					case 16:
						file[numfile] << b1y[i][j][k] << "  ";
						break;
					case 17:
						file[numfile] << b2y[i][j][k] << "  ";
						break;
					case 18:
						file[numfile] << b1z[i][j][k] << "  ";
						break;
					case 19:
						file[numfile] << b2z[i][j][k] << "  ";
						break;
					case 20:
						file[numfile] << V2[i][j][k] << "  ";
						break;
					}
				}
				//file[numfile] << endl;
			}
			//file[numfile] << endl;
		}
	}
}

double My_function(double x) // от 0 до п/2
{
	return (pow(x, 4) * exp(x)) / (pow(exp(x) - 1, 2));
}

double delta_function(double Tl, double Tm, double delta)
{
	return (1 / (sqrt(2 * pi) * delta)) * exp(-pow(Tl - Tm, 2) / (2 * delta * delta));
}

//Splayn Calculation_Interpolation(string fiename)
//{
//	vector<Point> new_points;
//	ifstream in(fiename); // окрываем файл для чтения
//
//	VecD X;
//	VecD Y;
//
//	if (in.is_open())
//	{
//		double x, y;
//		while (in >> x >> y)
//		{
//			new_points.push_back(Point{ x, y });
//		}
//	}
//	in.close();
//
//	for (int i = 0; i < new_points.size(); i++)
//	{
//		X.push_back(new_points[i].x);
//		Y.push_back(new_points[i].y);
//	}
//
//	Splayn spl;
//
//	spl.Interpolate(X, Y, 1);
//	return spl;
//}

//Splayn Calculation_Interpolation(double*** Tei, int& Nz_heat, double& dz_heat, int fix_x, int fix_y)
//{
//	VecD X;
//	VecD Y;
//
//	//ofstream fout("zTei.txt");
//
//	for (int i = 0; i < Nz_heat; i++)
//	{
//		X.push_back(i * dz_heat);
//		Y.push_back(Tei[fix_x][fix_y][i]);
//		//fout << X[i] << "    " << Y[i] << endl;
//	}
//
//	//system("pause");
//	//for (int i = 0; i <= Nz_heat; i++)
//	//{
//	//	cout << Tei[fix_x][fix_y][i] << "   ";
//	//}
//	//fout.close();
//	//system("pause");
//
//	Splayn spl;
//
//	spl.Interpolate(X, Y, 1);
//	return spl;
//}

double Dependence_k_e_on_T(Metal mt, double T_e) // температура в К
{
	double k_b = 1.38e-23; // постоянная Больцмана
	int i;
	//   Au, Al, Cu, Ni
	//double n_a[4] = { 5.90e+28, 18.10e+28, 8.45e+28, 9.13e+28 };// 18.10e+28;//18.10e+28;//8.45e+28; // концентрация атомов
	double T_D[4] = { 165, 428, 340, 375 };//428;//340; // температура Дебая (К)
	double T_l = 300;
	double v_F[4] = { 1.40e+6, 1.98e+6, 1.57e+6, 2.04e+6 }; // скорость Ферми (м/с)
	double E_F[4] = { 5.53 * 1.6e-19, 11.70 * 1.6e-19, 7.03 * 1.6e-19, 11.67 * 1.6e-19 }; //5.53*1.6e-19; //11.70 * 1.6e-19; // энергия Ферми (Дж)
	double Y[4] = { 68, 91.2, 97, 46 };// 1065;//91.2; // (Дж/м3 К2)
	double A[4] = { 1.18e+7, 0.376e+7, 1.28e+7, 0.59e+7 };//0.376e+7; // (1/(с К2))
	double B[4] = { 1.25e+11, 3.9e+11, 1.23e+11, 1.4e+11 }; // 3.9e+11; // (1/(с К))

	switch (mt)
	{
	case Au:
		i = 0;
		break;
	case Al:
		i = 1;
		break;
	case Cu:
		i = 2;
		break;
	case Ni:
		i = 3;
		break;
	}

	double bb = (B[i] * k_b) / (A[i] * E_F[i]);
	double K = (k_b * pow(v_F[i], 2) * Y[i] / 3) / (0.147 * A[i] * E_F[i]);
	double v_i = k_b * T_l / E_F[i];
	double v_e = k_b * T_e / E_F[i];
	double k_e = K * v_e * (pow(v_e * v_e + 0.16, 5 / 4) * (v_e * v_e + 0.44)) / (pow(v_e * v_e + 0.092, 1 / 2) * (v_e * v_e + bb * v_i));

	return k_e;
}

double Dependence_C_l_on_T(Metal mt, double T_l)
{
	int i;
	//   Au, Al, Cu, Ni
	//double n_a[4] = { 5.90e+28, 18.10e+28, 8.45e+28, 9.13e+28 };// 18.10e+28;//18.10e+28;//8.45e+28; // концентрация атомов
	double n_a[4] = { 5.90e+28, 6.02e+28, 8.45e+28, 9.13e+28 };// 18.10e+28;//18.10e+28;//8.45e+28; // концентрация атомов
	double T_D[4] = { 165, 428, 340, 375 };//428;//340; // температура Дебая (К)
	double k_b = 1.38e-23; // постоянная Больцмана

	switch (mt)
	{
	case Au:
		i = 0;
		break;
	case Al:
		i = 1;
		break;
	case Cu:
		i = 2;
		break;
	case Ni:
		i = 3;
		break;
	}

	double a = 0;// 1e-5;
	double b = T_D[i] / T_l;
	double eps = 1.e-7;
	double result = 0;

	//result = opt_h(a, b, eps);
	result = ((b - a) / 2) * (My_function(((b + a) / 2) - ((b - a) / 2) * (1 / (sqrt(3)))) + My_function(((b + a) / 2) + ((b - a) / 2) * (1 / (sqrt(3)))));
	double C_l = 9 * n_a[i] * k_b * pow(T_l / T_D[i], 3) * result;

	return C_l;
}

vector<double> proced(double*** T, int fixNx, int fixNy, int totalNzacoust, int totalNzheat)
{
	vector<double> newT;
	for (int i = 0; i < totalNzheat - 1; i++)
	{
		//double temph = (T[fixNx][fixNy][i] + T[fixNx][fixNy][i + 1]) / (totalNzacoust / totalNzheat); 
		double temph = (T[fixNx][fixNy][i + 1] - T[fixNx][fixNy][i]) / (totalNzacoust / totalNzheat);
		for (int j = 0; j < (totalNzacoust / totalNzheat); j++)
		{
			newT.push_back(T[fixNx][fixNy][i] + j * temph);
		}
	}
	newT.push_back(T[fixNx][fixNy][totalNzheat - 1]);
	return newT;
}

void GRU1(double*** e2, int& Nx, int& Ny)
{ // Boundary condition of
  // GRU 1

	for (int k = 0; k < Nx; k++)
	{
		for (int kk = 0; kk < Ny; kk++)
		{
			e2[k][kk][0] = 1e-16;
		}
	}
}

void GRU2(double*** b2x, double*** b2y, double*** b2z, int& Nx, int& Ny, int& Nz)
{// Boundary condition of
 // GRU 2

	for (int k = 0; k < Nx; k++)
	{
		for (int kk = 0; kk < Ny; kk++)
		{
			b2x[k][kk][Nz - 1] = 1e-16;
			b2y[k][kk][Nz - 1] = 1e-16;
			b2z[k][kk][Nz - 1] = 1e-16;
		}
	}

	for (int kk = 0; kk < Ny; kk++)
	{
		for (int kkk = 0; kkk < Nz; kkk++)
		{
			b2x[0][kk][kkk] = 1e-16;
			b2y[0][kk][kkk] = 1e-16;
			b2z[0][kk][kkk] = 1e-16;
			b2x[Nx - 1][kk][kkk] = 1e-16;
			b2y[Nx - 1][kk][kkk] = 1e-16;
			b2z[Nx - 1][kk][kkk] = 1e-16;
		}
	}

	for (int k = 0; k < Nx; k++)
	{
		for (int kkk = 0; kkk < Nz; kkk++)
		{
			b2x[k][0][kkk] = 1e-16;
			b2y[k][0][kkk] = 1e-16;
			b2z[k][0][kkk] = 1e-16;
			b2x[k][Ny - 1][kkk] = 1e-16;
			b2y[k][Ny - 1][kkk] = 1e-16;
			b2z[k][Ny - 1][kkk] = 1e-16;
		}
	}
}

void GaussBeam(double*** F, int& Nx, int& Ny, int& Nz, double& dx, double& dy, double& dz, double& dt, int& n, double& beta)
{
	for (int i = 0; i < Nx; i++)
	{
		for (int j = 0; j < Ny; j++)
		{
			for (int k = 0; k < Nz; k++)
			{
				F[i][j][k] = exp(-dx * dx * (i - Nx / 2) * (i - Nx / 2)) * exp(-dy * dy * (j - Ny / 2) * (j - Ny / 2)) * exp(-dz * k) * dt * n * exp(-beta * dt * n);
			}
		}
	}
}

void VortexBeam(double*** F, int& Nx, int& Ny, int& Nz, double& dx, double& dy, double& dz, double& dt, int& n, double& beta, int m)
{
	double h11 = 0, h12 = 0, hf1 = 0;
	for (int i = 0; i < Nx; i++)
	{
		for (int j = 0; j < Ny; j++)
		{
			for (int k = 0; k < Nz; k++)
			{
				h11 = dx * (i - Nx / 2);
				h12 = dy * (j - Ny / 2);
				hf1 = atan(h12 / h11);
				if (h11 == 0 && h12 == 0)
				{
					hf1 = 0;
				}

				if (h11 == 0)
				{
					hf1 = 0;
				}

				F[i][j][k] = (pow(sqrt(pow(h11 * h11 + h12 * h12, abs(m))) * cos(m * hf1) * (exp(-dx * dx * (i - Nx / 2) * (i - Nx / 2)) * exp(-dy * dy * (j - Ny / 2) * (j - Ny / 2))), 2) +
					pow(sqrt(pow(h11 * h11 + h12 * h12, abs(m))) * sin(m * hf1) * (exp(-dx * dx * (i - Nx / 2) * (i - Nx / 2)) * exp(-dy * dy * (j - Ny / 2) * (j - Ny / 2))), 2)) * exp(-dz * k) * dt * n * exp(-beta * dt * n);
			}
		}
	}
}

void Calculation00(Metal mt, Parametrs param, Splayn spl_C_e_on_T, Splayn spl_G_e_on_T, /*double*** F,*/ TypeBeam tbeam, double*** tmpe0, double*** tmpe1, double*** tmpi0, double*** tmpi1, int& Nx, int& Ny, int& Nz, double& dx, double& dy, double& dz, double& dt, /*double& A1,*/ double& A2, /*double& B1, double& CC1, double& CC2, */int& n, double& beta)
{// Calculation of heat conduction by explicit schem

	//if (tbeam == Gauss)
	//{
	//	//GaussBeam(F, Nx, Ny, Nz, dx, dy, dz, dt, n, beta);
	//}

	//if (tbeam == Vortex)
	//{
	//	//VortexBeam(F, Nx, Ny, Nz, dx, dy, dz, dt, n, beta, 1);
	//}

	double A111 = 0, B111 = 0, A1i111 = 0, CC1111 = 0, CC2111 = 0, FFF = 0;

	//for (int i = 2; i < Nx - 1; i++)
	for (int i = 1; i < Nx - 1; i++)
	{
		//for (int j = 2; j < Ny - 1; j++)
		for (int j = 1; j < Ny - 1; j++)
		{
			//for (int k = 2; k < Nz - 1; k++)
			for (int k = 1; k < Nz - 1; k++)
			{
				A111 = (Dependence_k_e_on_T(mt, param.T00 * tmpe0[i][j][k]) / (100.)) * param.t0 / ((spl_C_e_on_T.GetY(param.T00 * tmpe0[i][j][k]) / (100. * 100. * 100.)) * param.r0 * param.r0);
				B111 = param.kabs * param.P0 * param.t0 * beta / ((spl_C_e_on_T.GetY(param.T00 * tmpe0[i][j][k]) / (100. * 100. * 100.)) * param.T00);
				A1i111 = (Dependence_k_e_on_T(mt, 300.) / (100. * 100.)) * param.t0 / ((Dependence_C_l_on_T(mt, param.T00 * tmpi0[i][j][k]) / (100. * 100. * 100.)) * param.r0 * param.r0);
				FFF = exp(-dx * dx * (i - Nx / 2) * (i - Nx / 2)) * exp(-dy * dy * (j - Ny / 2) * (j - Ny / 2)) * exp(-dz * k) * dt * n * exp(-beta * dt * n);

				tmpe1[i][j][k] = tmpe0[i][j][k] +
					dt * A111 * ((tmpe0[i + 1][j][k] + tmpe0[i - 1][j][k] - 2 * tmpe0[i][j][k]) / (dx * dx) + (tmpe0[i][j + 1][k] + tmpe0[i][j - 1][k] - 2 * tmpe0[i][j][k]) / (dy * dy) + A2 * (tmpe0[i][j][k - 1] + tmpe0[i][j][k + 1] - 2 * tmpe0[i][j][k]) / (dz * dz))
					+ dt * B111 * /*F[i][j][k]*/ FFF;// -dt * CC1 * (tmpe0[i][j][k]);

													//tmpi1[i][j][k] = dt * CC2 * tmpe0[i][j][k];
													//tmpi1[i][j][k] = tmpi1[i][j][k];// +dt * CC2 * tmpe0[i][j][k];
				tmpi1[i][j][k] = tmpi0[i][j][k] +
					dt * A1i111 * ((tmpi0[i + 1][j][k] + tmpi0[i - 1][j][k] - 2 * tmpi0[i][j][k]) / (dx * dx) + (tmpi0[i][j + 1][k] + tmpi0[i][j - 1][k] - 2 * tmpi0[i][j][k]) / (dy * dy) + A2 * (tmpi0[i][j][k - 1] + tmpi0[i][j][k + 1] - 2 * tmpi0[i][j][k]) / (dz * dz));

			}
		}
	}

	for (int i = 0; i < Nx; i++)
	{
		for (int j = 0; j < Ny; j++)
		{
			tmpe1[i][j][Nz - 1] = 1.0;// 1e-16;
			tmpi1[i][j][Nz - 1] = 1.0;// 1e-16;
			tmpe1[i][j][0] = tmpe1[i][j][1];
			tmpi1[i][j][0] = tmpi1[i][j][1];
		}
	}

	for (int j = 0; j < Ny; j++)
	{
		for (int k = 0; k < Nz; k++)
		{
			tmpe1[0][j][k] = tmpe1[1][j][k];
			tmpi1[0][j][k] = tmpi1[1][j][k];
			tmpe1[Nx - 1][j][k] = tmpe1[Nx - 2][j][k];
			tmpi1[Nx - 1][j][k] = tmpi1[Nx - 2][j][k];
		}
	}

	for (int i = 0; i < Nx; i++)
	{
		for (int k = 0; k < Nz; k++)
		{
			tmpe1[i][0][k] = tmpe1[i][1][k];
			tmpi1[i][0][k] = tmpi1[i][1][k];
			tmpe1[i][Ny - 1][k] = tmpe1[i][Ny - 2][k];
			tmpi1[i][Ny - 1][k] = tmpi1[i][Ny - 2][k];
		}
	}
}

void fun2(int& Nx, int& Ny, int& Nz, double& dx, double& dy, double& dz, double& dt, double*** tmpe0, double*** tmpe1, double*** tmpe2, double*** tmpi0, double*** tmpi1, double*** tmpi2, double& A2,/*double*** A1, double*** A1i, double& A2, double*** B1, double*** CC1, double*** CC2,*/ /*double*** F,*/ /*double*** copy_tmpe1, double*** copy_tmpi1,*/ int& begin, int& end, Parametrs& param, vector<Melting>& Melt_metal, Metal mt, /*double*** Ci, double*** gamma,*/ bool& melting, Splayn& spl_G_e_on_T, Splayn spl_C_e_on_T, double beta, int n, Splayn spl_C_l_on_T)
{
	double A111 = 0, B111 = 0, A1i111 = 0, CC1111 = 0, CC2111 = 0, FFF = 0;
	for (int k = begin; k < end; k++)
		//for (int i = 1; i < Nx - 1; i++)
		//for (int k = 1; k < Nz - 1; k++)
	{
		for (int j = 1; j < Ny - 1; j++)
		{
			for (int i = 1; i < Nx - 1; i++)
				//for (int k = 1; k < Nz - 1; k++)
			{
				A1i111 = ((Dependence_k_e_on_T(mt, 300.) / (100. * 100.)) * param.t0 / ((Dependence_C_l_on_T(mt, param.T00 * tmpi1[i][j][k]) / (100. * 100. * 100.)) * param.r0 * param.r0));
				A111 = ((Dependence_k_e_on_T(mt, param.T00 * tmpe1[i][j][k]) / (100.)) * param.t0 / ((spl_C_e_on_T.GetY(param.T00 * tmpe1[i][j][k]) / (100. * 100. * 100.)) * param.r0 * param.r0));
				B111 = (param.kabs * param.P0 * param.t0 * beta / ((spl_C_e_on_T.GetY(param.T00 * tmpe1[i][j][k]) / (100. * 100. * 100.)) * param.T00));
				CC1111 = ((spl_G_e_on_T.GetY(param.T00 * tmpe1[i][j][k]) / (100. * 100. * 100.)) * param.t0 / ((spl_C_e_on_T.GetY(param.T00 * tmpe1[i][j][k]) / (100. * 100. * 100.))));
				CC2111 = ((spl_G_e_on_T.GetY(param.T00 * tmpe1[i][j][k]) / (100. * 100. * 100.)) * param.t0 / ((Dependence_C_l_on_T(mt, param.T00 * tmpi1[i][j][k]) / (100. * 100. * 100.))));
				FFF = exp(-dx * dx * (i - Nx / 2) * (i - Nx / 2)) * exp(-dy * dy * (j - Ny / 2) * (j - Ny / 2)) * exp(-dz * k) * dt * n * exp(-beta * dt * n);

				tmpe2[i][j][k] = tmpe0[i][j][k] * (1 - 2 * dt * /*A1[i][j][k]*/ A111 / (dx * dx) - 2 * dt * /*A1[i][j][k]*/ A111 / (dy * dy) - 2 * dt * /*A1[i][j][k]*/ A111 * A2 / (dz * dz))
					+ 2 * dt */* A1[i][j][k]*/ A111 * ((tmpe1[i + 1][j][k] + tmpe1[i - 1][j][k]) / (dx * dx) + (tmpe1[i][j + 1][k] + tmpe1[i][j - 1][k]) / (dy * dy) + A2 * (tmpe1[i][j][k - 1] + tmpe1[i][j][k + 1]) / (dz * dz))
					+ 2 * dt * B111/*B1[i][j][k]*/ * /*F[i][j][k]*/ FFF - 2 * dt * /*CC1[i][j][k]*/ CC1111 * (tmpe1[i][j][k] - tmpi1[i][j][k]);

				tmpe2[i][j][k] = tmpe2[i][j][k] / (1 + 2 * dt * /*A1[i][j][k]*/ A111 / (dx * dx) + 2 * dt */* A1[i][j][k]*/ A111 / (dy * dy) + 2 * dt * /*A1[i][j][k] */ A111 * A2 / (dz * dz));

				double delta = 1.;
				if ((tmpi1[i][j][k] * param.T00) < (Melt_metal[mt].T_melting - delta))
				{

					tmpi2[i][j][k] = tmpi0[i][j][k] * (1 - 2 * dt * /*A1i[i][j][k]*/ A1i111 / (dx * dx) - 2 * dt * /*A1i[i][j][k]*/ A1i111 / (dy * dy) - 2 * dt * /*A1i[i][j][k]*/ A1i111 * A2 / (dz * dz))
						+ 2 * dt * /*A1i[i][j][k]*/ A1i111 * ((tmpi1[i + 1][j][k] + tmpi1[i - 1][j][k]) / (dx * dx) + (tmpi1[i][j + 1][k] + tmpi1[i][j - 1][k]) / (dy * dy) + A2 * (tmpi1[i][j][k - 1] + tmpi1[i][j][k + 1]) / (dz * dz))
						+ 2 * dt * /*CC2[i][j][k]*/ CC2111 * (tmpe1[i][j][k] - tmpi1[i][j][k]);

					tmpi2[i][j][k] = tmpi2[i][j][k] / (1 + 2 * dt * /*A1i[i][j][k]*/ A1i111 / (dx * dx) + 2 * dt * /*A1i[i][j][k]*/ A1i111 / (dy * dy) + 2 * dt * /*A1i[i][j][k]*/ A1i111 * A2 / (dz * dz));

				}

				if ((tmpi1[i][j][k] * param.T00) >= (Melt_metal[mt].T_melting - delta) && (tmpi1[i][j][k] * param.T00) <= (Melt_metal[mt].T_melting + delta))
				{
					//melting = true;
					// старый способ расплава
					/*A1i111 = (Dependence_k_e_on_T(mt, 300.) / (100. * 100.)) * param.t0 / ((Dependence_C_l_on_T(mt, param.T00 * tmpi1[i][j][k]) / (100. * 100. * 100.) + (Melt_metal[mt].Q_fusion * Melt_metal[mt].Density) / (2 * delta)) * param.r0 * param.r0);
					CC2111 = (spl_G_e_on_T.GetY(param.T00 * tmpe1[i][j][k]) / (100. * 100. * 100.)) * param.t0 / ((Dependence_C_l_on_T(mt, param.T00 * tmpi1[i][j][k]) / (100. * 100. * 100.) + (Melt_metal[mt].Q_fusion * Melt_metal[mt].Density) / (2 * delta)));*/

					/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

					// новый способ расплава
					A1i111 = (Dependence_k_e_on_T(mt, 300.) / (100. * 100.)) * param.t0 / ((/*Dependence_C_l_on_T(mt, param.T00 * tmpi1[i][j][k]) / (100. * 100. * 100.) +*/ spl_C_l_on_T.GetY(param.T00 * tmpi1[i][j][k])) * param.r0 * param.r0);
					CC2111 = (spl_G_e_on_T.GetY(param.T00 * tmpe1[i][j][k]) / (100. * 100. * 100.)) * param.t0 / ((/*Dependence_C_l_on_T(mt, param.T00 * tmpi1[i][j][k]) / (100. * 100. * 100.) + */ spl_C_l_on_T.GetY(param.T00 * tmpi1[i][j][k])));

					////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

					tmpi2[i][j][k] = tmpi0[i][j][k] * (1 - 2 * dt * /*A1i[i][j][k]*/ A1i111 / (dx * dx) - 2 * dt * /*A1i[i][j][k]*/ A1i111 / (dy * dy) - 2 * dt * /*A1i[i][j][k]*/ A1i111 * A2 / (dz * dz))
						+ 2 * dt */* A1i[i][j][k]*/ A1i111 * ((tmpi1[i + 1][j][k] + tmpi1[i - 1][j][k]) / (dx * dx) + (tmpi1[i][j + 1][k] + tmpi1[i][j - 1][k]) / (dy * dy) + A2 * (tmpi1[i][j][k - 1] + tmpi1[i][j][k + 1]) / (dz * dz))
						+ 2 * dt * /*CC2[i][j][k]*/ CC2111 * (tmpe1[i][j][k] - tmpi1[i][j][k]);

					tmpi2[i][j][k] = tmpi2[i][j][k] / (1 + 2 * dt * /*A1i[i][j][k]*/ A1i111 / (dx * dx) + 2 * dt * /*A1i[i][j][k]*/ A1i111 / (dy * dy) + 2 * dt * /*A1i[i][j][k]*/ A1i111 * A2 / (dz * dz));

				}

				if ((tmpi1[i][j][k] * param.T00) > (Melt_metal[mt].T_melting + delta))
				{
					A1i111 = (Dependence_k_e_on_T(mt, 300.) / (100. * 100.)) * param.t0 / ((Dependence_C_l_on_T(mt, param.T00 * tmpi1[i][j][k]) / (100. * 100. * 100.)) * param.r0 * param.r0);
					CC2111 = (spl_G_e_on_T.GetY(param.T00 * tmpe1[i][j][k]) / (100. * 100. * 100.)) * param.t0 / ((Dependence_C_l_on_T(mt, param.T00 * tmpi1[i][j][k]) / (100. * 100. * 100.)));

					/*A1i111 = (Dependence_k_e_on_T(mt, 300.) / (100. * 100.)) * param.t0 / ((spl_C_l_on_T.GetY(param.T00 * tmpi1[i][j][k])) * param.r0 * param.r0);
					CC2111 = (spl_G_e_on_T.GetY(param.T00 * tmpe1[i][j][k]) / (100. * 100. * 100.)) * param.t0 / ((spl_C_l_on_T.GetY(param.T00 * tmpi1[i][j][k])));*/

					tmpi2[i][j][k] = tmpi0[i][j][k] * (1 - 2 * dt * /*A1i[i][j][k]*/ A1i111 / (dx * dx) - 2 * dt * /*A1i[i][j][k]*/ A1i111 / (dy * dy) - 2 * dt * /*A1i[i][j][k]*/ A1i111 * A2 / (dz * dz))
						+ 2 * dt * /*A1i[i][j][k]*/ A1i111 * ((tmpi1[i + 1][j][k] + tmpi1[i - 1][j][k]) / (dx * dx) + (tmpi1[i][j + 1][k] + tmpi1[i][j - 1][k]) / (dy * dy) + A2 * (tmpi1[i][j][k - 1] + tmpi1[i][j][k + 1]) / (dz * dz))
						+ 2 * dt * /*CC2[i][j][k]*/ CC2111 * (tmpe1[i][j][k] - tmpi1[i][j][k]);

					tmpi2[i][j][k] = tmpi2[i][j][k] / (1 + 2 * dt * /*A1i[i][j][k]*/ A1i111 / (dx * dx) + 2 * dt * /*A1i[i][j][k]*/ A1i111 / (dy * dy) + 2 * dt * /*A1i[i][j][k]*/ A1i111 * A2 / (dz * dz));
				}
			}
		}
	}
}

void fun21(int& Nx, int& Ny, int& Nz, double*** tmpe2, double*** tmpi1, double*** tmpi2, /*double** fict_masive,*/ double& dz, Parametrs& param, Metal mt)
{
	double a4, a5, T_air = 300, eps = 0.1, sigma = 5.67e-12; // Вт/см2 К4 //???????????????
	double h_c = 0.0015; // Вт/cm2 K
	double f1;
	for (int i = 0; i < Nx; i++)
	{
		for (int j = 0; j < Ny; j++)
		{
			tmpe2[i][j][Nz - 1] = 1.0;// 1e-16;
			tmpi2[i][j][Nz - 1] = 1.0;// 1e-16;
			tmpe2[i][j][0] = tmpe2[i][j][1];
			tmpi2[i][j][0] = tmpi2[i][j][1];
			//эта строчка заменяются на следующие
		   //a4 = h_c / (param.kabs * (Dependence_k_e_on_T(mt, 300.) / (100. * 100.))); // kl = 1% от ke (при комнат темп)
		   //a5 = eps * sigma * pow(param.T00, 3) / (param.kabs * (Dependence_k_e_on_T(mt, 300.) / (100. * 100.))); //  kl = 1% от ke (при комнат темп)
		   //f1 = a4 * (tmpi1[i][j][0] - 1.) + a5 * (pow(tmpi1[i][j][0], 4) - 1.);// ГУ-2-го рода
		   //fict_masive[i][j] = tmpi2[i][j][1] - 2 * dz * f1; // С использованием центральной разностной аппроксимации градиента на поверхности z=0 находим значение температуры во внешнем узле Q(i,j,0,2).
		   //tmpi2[i][j][0] = (fict_masive[i][j] + tmpi2[i][j][1]) / 2; // Значение температуры на поверхности раздела сред (узел k=1) находим как среднее между внешним узлом Q(i,j,0,2) и внутренним узлом Q(i,j,2,2

		}
	}
}

void fun22(int& Nx, int& Ny, int& Nz, double*** tmpe2, double*** tmpi2)
{
	for (int j = 0; j < Ny; j++)
	{
		for (int k = 0; k < Nz; k++)
		{
			tmpe2[0][j][k] = tmpe2[1][j][k];
			tmpi2[0][j][k] = tmpi2[1][j][k];
			tmpe2[Nx - 1][j][k] = tmpe2[Nx - 2][j][k];
			tmpi2[Nx - 1][j][k] = tmpi2[Nx - 2][j][k];
		}
	}
}

void fun23(int& Nx, int& Ny, int& Nz, double*** tmpe2, double*** tmpi2)
{
	for (int i = 0; i < Nx; i++)
	{
		for (int k = 0; k < Nz; k++)
		{
			tmpe2[i][0][k] = tmpe2[i][1][k];
			tmpi2[i][0][k] = tmpi2[i][1][k];
			tmpe2[i][Ny - 1][k] = tmpe2[i][Ny - 2][k];
			tmpi2[i][Ny - 1][k] = tmpi2[i][Ny - 2][k];
		}
	}
}

void Calculation0(Metal mt, Parametrs& param, Splayn& spl_C_e_on_T, Splayn& spl_G_e_on_T,/* double*** F,*/ TypeBeam tbeam, double*** tmpe0, double*** tmpe1, double*** tmpe2, double*** tmpi0, double*** tmpi1, double*** tmpi2,/* double& A1,*/ double& A2, /*double& B1, double& CC1, double& CC2,*/ double& dx, double& dy, double& dz, double& dt, int& Nx, int& Ny, int& Nz, int& n, double& beta, vector<Interval>& new_interv,/* double** fict_masive,*/ vector<Melting>& Melt_metal, bool& melting, Splayn spl_C_l_on_T)
{
	//if (tbeam == Gauss)
	//{
	//	//GaussBeam(F, Nx, Ny, Nz, dx, dy, dz, dt, n, beta);
	//}

	//if (tbeam == Vortex)
	//{
	////	VortexBeam(F, Nx, Ny, Nz, dx, dy, dz, dt, n, beta, 1);
	//}

	//fun2(ref(Nx), ref(Ny), ref(Nz), ref(dx), ref(dy), ref(dz), ref(dt), tmpe0, tmpe1, tmpe2, tmpi0, tmpi1, tmpi2, ref(A2), /* A1, A1i, ref(A2), B1, CC1, CC2, F,*/ /*copy_tmpe1, copy_tmpi1,*/ ref(new_interv[0].begin), ref(new_interv[0].end), ref(param), ref(Melt_metal), mt, /*Ci, gamma,*/ ref(melting), ref(spl_G_e_on_T), spl_C_e_on_T, beta, n);

	std::thread t26(fun2, ref(Nx), ref(Ny), ref(Nz), ref(dx), ref(dy), ref(dz), ref(dt), tmpe0, tmpe1, tmpe2, tmpi0, tmpi1, tmpi2, ref(A2), /* A1, A1i, ref(A2), B1, CC1, CC2, F,*/ /*copy_tmpe1, copy_tmpi1,*/ ref(new_interv[0].begin), ref(new_interv[0].end), ref(param), ref(Melt_metal), mt, /*Ci, gamma,*/ ref(melting), ref(spl_G_e_on_T), spl_C_e_on_T, beta, n, spl_C_l_on_T);
	std::thread t27(fun2, ref(Nx), ref(Ny), ref(Nz), ref(dx), ref(dy), ref(dz), ref(dt), tmpe0, tmpe1, tmpe2, tmpi0, tmpi1, tmpi2, ref(A2), /* A1, A1i, ref(A2), B1, CC1, CC2,*/ /*F,*//*copy_tmpe1, copy_tmpi1,*/ ref(new_interv[1].begin), ref(new_interv[1].end), ref(param), ref(Melt_metal), mt, /*Ci, gamma,*/ ref(melting), ref(spl_G_e_on_T), spl_C_e_on_T, beta, n, spl_C_l_on_T);
	std::thread t28(fun2, ref(Nx), ref(Ny), ref(Nz), ref(dx), ref(dy), ref(dz), ref(dt), tmpe0, tmpe1, tmpe2, tmpi0, tmpi1, tmpi2, ref(A2), /* A1, A1i, ref(A2), B1, CC1, CC2,*//* F,*/ /*copy_tmpe1, copy_tmpi1,*/ ref(new_interv[2].begin), ref(new_interv[2].end), ref(param), ref(Melt_metal), mt, /*Ci, gamma,*/ ref(melting), ref(spl_G_e_on_T), spl_C_e_on_T, beta, n, spl_C_l_on_T);
	std::thread t29(fun2, ref(Nx), ref(Ny), ref(Nz), ref(dx), ref(dy), ref(dz), ref(dt), tmpe0, tmpe1, tmpe2, tmpi0, tmpi1, tmpi2, ref(A2), /* A1, A1i, ref(A2), B1, CC1, CC2,*//* F,*/ /*copy_tmpe1, copy_tmpi1,*/ ref(new_interv[3].begin), ref(new_interv[3].end), ref(param), ref(Melt_metal), mt, /*Ci, gamma,*/ ref(melting), ref(spl_G_e_on_T), spl_C_e_on_T, beta, n, spl_C_l_on_T);
	std::thread t30(fun2, ref(Nx), ref(Ny), ref(Nz), ref(dx), ref(dy), ref(dz), ref(dt), tmpe0, tmpe1, tmpe2, tmpi0, tmpi1, tmpi2, ref(A2), /* A1, A1i, ref(A2), B1, CC1, CC2,*//* F,*/ /*copy_tmpe1, copy_tmpi1,*/ ref(new_interv[4].begin), ref(new_interv[4].end), ref(param), ref(Melt_metal), mt, /*Ci, gamma,*/ ref(melting), ref(spl_G_e_on_T), spl_C_e_on_T, beta, n, spl_C_l_on_T);
	std::thread t31(fun2, ref(Nx), ref(Ny), ref(Nz), ref(dx), ref(dy), ref(dz), ref(dt), tmpe0, tmpe1, tmpe2, tmpi0, tmpi1, tmpi2, ref(A2), /* A1, A1i, ref(A2), B1, CC1, CC2,*/ /*F,*/ /*copy_tmpe1, copy_tmpi1,*/ ref(new_interv[5].begin), ref(new_interv[5].end), ref(param), ref(Melt_metal), mt, /*Ci, gamma,*/ ref(melting), ref(spl_G_e_on_T), spl_C_e_on_T, beta, n, spl_C_l_on_T);
	std::thread t32(fun2, ref(Nx), ref(Ny), ref(Nz), ref(dx), ref(dy), ref(dz), ref(dt), tmpe0, tmpe1, tmpe2, tmpi0, tmpi1, tmpi2, ref(A2), /* A1, A1i, ref(A2), B1, CC1, CC2,*/ /*F, *//*copy_tmpe1, copy_tmpi1,*/ ref(new_interv[6].begin), ref(new_interv[6].end), ref(param), ref(Melt_metal), mt, /*Ci, gamma,*/ ref(melting), ref(spl_G_e_on_T), spl_C_e_on_T, beta, n, spl_C_l_on_T);
	std::thread t33(fun2, ref(Nx), ref(Ny), ref(Nz), ref(dx), ref(dy), ref(dz), ref(dt), tmpe0, tmpe1, tmpe2, tmpi0, tmpi1, tmpi2, ref(A2), /* A1, A1i, ref(A2), B1, CC1, CC2,*//* F,*/ /*copy_tmpe1, copy_tmpi1,*/ ref(new_interv[7].begin), ref(new_interv[7].end), ref(param), ref(Melt_metal), mt, /*Ci, gamma,*/ ref(melting), ref(spl_G_e_on_T), spl_C_e_on_T, beta, n, spl_C_l_on_T);
	std::thread t34(fun2, ref(Nx), ref(Ny), ref(Nz), ref(dx), ref(dy), ref(dz), ref(dt), tmpe0, tmpe1, tmpe2, tmpi0, tmpi1, tmpi2, ref(A2), /* A1, A1i, ref(A2), B1, CC1, CC2,*/ /*F,*/ /*copy_tmpe1, copy_tmpi1,*/ ref(new_interv[8].begin), ref(new_interv[8].end), ref(param), ref(Melt_metal), mt, /*Ci, gamma,*/ ref(melting), ref(spl_G_e_on_T), spl_C_e_on_T, beta, n, spl_C_l_on_T);
	std::thread t35(fun2, ref(Nx), ref(Ny), ref(Nz), ref(dx), ref(dy), ref(dz), ref(dt), tmpe0, tmpe1, tmpe2, tmpi0, tmpi1, tmpi2, ref(A2), /* A1, A1i, ref(A2), B1, CC1, CC2,*/  /*F, *//*copy_tmpe1, copy_tmpi1,*/ ref(new_interv[9].begin), ref(new_interv[9].end), ref(param), ref(Melt_metal), mt,/* Ci, gamma,*/ ref(melting), ref(spl_G_e_on_T), spl_C_e_on_T, beta, n, spl_C_l_on_T);

	t26.join();
	t27.join();
	t28.join(); ////////////////////////////////////
	t29.join();
	t30.join();
	t31.join();
	t32.join(); ////////////////////////////////////
	t33.join();
	t34.join();
	t35.join(); ////////////////////////////////////

	thread t51(fun21, ref(Nx), ref(Ny), ref(Nz), tmpe2, tmpi1, tmpi2,/* fict_masive,*/ ref(dz), ref(param), mt);
	thread t52(fun22, ref(Nx), ref(Ny), ref(Nz), tmpe2, tmpi2);
	thread t53(fun23, ref(Nx), ref(Ny), ref(Nz), tmpe2, tmpi2);

	t51.join();
	t52.join();
	t53.join(); ////////////////////////////////////
}

void Calculation1b2y(int& Nx, int& Ny, int& Nz, double& dx, double& dy, double& dz, double& dt, double& CC0, double*** e1, double*** a1x, double*** a1z, double*** b1y, double*** b2y, int& begin, int& end)
{
	/*for (int i = 0; i < Nx; i++)*/
	for (int l = begin; l < end; l++)
		//for (int i = 0; i < Nx; i++)
	//for (int l = 0; l < Nz; l++)
	{
		for (int j = 0; j < Ny; j++)
		{
			for (int i = 0; i < Nx; i++)
				//for (int l = 0; l < Nz; l++)
			{
				if ((i - 1) < 0 || (i + 1) > (Nx - 1) || (j - 1) < 0 || (j + 1) > (Ny - 1) || (l - 1) < 0 || (l + 1) > (Nz - 1))
				{
					//b2y[i][j][l] = b1y[i][j][l] + dt*CC0*(e1[i][j][l] - tmp)*((tmp - a1x[i][j][l])*(tmp - a1z[i][j][l]) - (tmp - a1x[i][j][l])*(tmp - a1z[i][j][l])) / (dx*dy*dz);
					//b2y[i][j][l] = b2y[i][j][l] + dt*CC0*(e1[i][j][l] - tmp)*((tmp - a1x[i][j][l])*(tmp - a1z[i][j][l]) - (tmp - a1x[i][j][l])*(tmp - a1z[i][j][l])) / (dx*dy*dz);
					//b2y[i][j][l] = b2y[i][j][l] + dt*CC0*(e1[i][j][l] - tmp)*((tmp - a1x[i][j][l])*(tmp - a1z[i][j][l]) - (tmp - a1x[i][j][l])*(tmp - a1z[i][j][l])) / (dx*dy*dz);
				}
				else
				{
					b2y[i][j][l] = b1y[i][j][l] + dt * CC0 * (e1[i][j][l] - e1[i - 1][j][l]) * ((a1x[i][j + 1][l] - a1x[i][j][l]) * (a1z[i][j][l + 1] - a1z[i][j][l]) - (a1x[i][j][l + 1] - a1x[i][j][l]) * (a1z[i][j + 1][l] - a1z[i][j][l])) / (dx * dy * dz);
					b2y[i][j][l] = b2y[i][j][l] + dt * CC0 * (e1[i][j][l] - e1[i][j - 1][l]) * ((a1x[i][j][l + 1] - a1x[i][j][l]) * (a1z[i + 1][j][l] - a1z[i][j][l]) - (a1x[i + 1][j][l] - a1x[i][j][l]) * (a1z[i][j][l + 1] - a1z[i][j][l])) / (dx * dy * dz);
					b2y[i][j][l] = b2y[i][j][l] + dt * CC0 * (e1[i][j][l] - e1[i][j][l - 1]) * ((a1x[i + 1][j][l] - a1x[i][j][l]) * (a1z[i][j + 1][l] - a1z[i][j][l]) - (a1x[i][j + 1][l] - a1x[i][j][l]) * (a1z[i + 1][j][l] - a1z[i][j][l])) / (dx * dy * dz);
				}
			}
		}
	}
}

void Calculation1b2z(int& Nx, int& Ny, int& Nz, double& dx, double& dy, double& dz, double& dt, double*** e1, double*** a1x, double*** a1y, double*** b1z, double*** b2z, int& begin, int& end)
{
	/*for (int i = 0; i < Nx; i++)*/
	for (int l = begin; l < end; l++)
		//for (int i = 0; i < Nx; i++)
	{
		for (int j = 0; j < Ny; j++)
		{
			for (int i = 0; i < Nx; i++)
				//for (int l = 0; l < Nz; l++)
			{
				if ((i - 1) < 0 || (i + 1) > (Nx - 1) || (j - 1) < 0 || (j + 1) > (Ny - 1) || (l - 1) < 0 || (l + 1) > (Nz - 1))
				{
					//b2z[i][j][l] = b1z[i][j][l] + dt*(e1[i][j][l] - tmp)*((tmp - a1x[i][j][l])*(tmp - a1y[i][j][l]) - (tmp - a1x[i][j][l])*(tmp - a1y[i][j][l])) / (dx*dy*dz);
					//b2z[i][j][l] = b2z[i][j][l] + dt*(e1[i][j][l] - tmp)*((tmp - a1x[i][j][l])*(tmp - a1y[i][j][l]) - (tmp - a1x[i][j][l])*(tmp - a1y[i][j][l])) / (dx*dy*dz);
					//b2z[i][j][l] = b2z[i][j][l] + dt*(e1[i][j][l] - tmp)*((tmp - a1x[i][j][l])*(tmp - a1y[i][j][l]) - (tmp - a1x[i][j][l])*(tmp - a1y[i][j][l])) / (dx*dy*dz);
				}
				else
				{
					b2z[i][j][l] = b1z[i][j][l] + dt * (e1[i][j][l] - e1[i - 1][j][l]) * ((a1x[i][j][l + 1] - a1x[i][j][l]) * (a1y[i][j + 1][l] - a1y[i][j][l]) - (a1x[i][j + 1][l] - a1x[i][j][l]) * (a1y[i][j][l + 1] - a1y[i][j][l])) / (dx * dy * dz);
					b2z[i][j][l] = b2z[i][j][l] + dt * (e1[i][j][l] - e1[i][j - 1][l]) * ((a1x[i + 1][j][l] - a1x[i][j][l]) * (a1y[i][j][l + 1] - a1y[i][j][l]) - (a1x[i][j][l + 1] - a1x[i][j][l]) * (a1y[i + 1][j][l] - a1y[i][j][l])) / (dx * dy * dz);
					b2z[i][j][l] = b2z[i][j][l] + dt * (e1[i][j][l] - e1[i][j][l - 1]) * ((a1x[i][j + 1][l] - a1x[i][j][l]) * (a1y[i + 1][j][l] - a1y[i][j][l]) - (a1x[i + 1][j][l] - a1x[i][j][l]) * (a1y[i][j + 1][l] - a1y[i][j][l])) / (dx * dy * dz);
				}
			}
		}
	}
}

void Calculation1b2xPart1(int& Ny, int& Nz, double& dx, double& dy, double& dz, double& dt, double*** e1, double*** a1y, double*** a1z, double*** b1x, double*** b2x, double& CC0, int& begin, int& end)
{
	for (int l = begin; l < end; l++)
		//for (int j = 0; j < Ny; j++)//????????????????????????????????????????????????????????????????????????
	{// Или тут рассчитываются только внутренние узлы???????????????????????????????????И почему тут часть кода рабоатет???????????
		//for (int l = 0; l < Nz; l++)
		for (int j = 0; j < Ny; j++)
		{
			if ((j - 1) < 0 || (j + 1) > (Ny - 1) || (l - 1) < 0 || (l + 1) > (Nz - 1))
			{
				//cout << j << "  " << l << endl;
				//b2x[0][j][l] = b1x[0][j][l] + dt*CC0*(e1[1][j][l] - e1[0][j][l])*(-(tmp - a1y[0][j][l])*(tmp - a1z[0][j][l]) + (tmp - a1y[0][j][l])*(tmp - a1z[0][j][l])) / (dx*dy*dz);
				//b2x[0][j][l] = b2x[0][j][l] + dt*CC0*(e1[0][j][l] - tmp)*((a1y[1][j][l] - a1y[0][j][l])*(tmp - a1z[0][j][l]) - (tmp - a1y[0][j][l])*(tmp - a1z[0][j][l])) / (dx*dy*dz);
				//b2x[0][j][l] = b2x[0][j][l] + dt*CC0*(e1[0][j][l] - tmp)*((tmp - a1y[0][j][l])*(a1z[1][j][l] - a1z[0][j][l]) - (tmp - a1y[0][j][l])*(tmp - a1z[0][j][l])) / (dx*dy*dz);
			}
			else {
				//cout << j << "  " << l << endl;
				b2x[0][j][l] = b1x[0][j][l] + dt * CC0 * (e1[1][j][l] - e1[0][j][l]) * (-(a1y[0][j + 1][l] - a1y[0][j][l]) * (a1z[0][j][l + 1] - a1z[0][j][l]) + (a1y[0][j][l + 1] - a1y[0][j][l]) * (a1z[0][j + 1][l] - a1z[0][j][l])) / (dx * dy * dz);
				b2x[0][j][l] = b2x[0][j][l] + dt * CC0 * (e1[0][j][l] - e1[0][j - 1][l]) * ((a1y[1][j][l] - a1y[0][j][l]) * (a1z[0][j][l + 1] - a1z[0][j][l]) - (a1y[0][j][l + 1] - a1y[0][j][l]) * (a1z[0][j][l + 1] - a1z[0][j][l])) / (dx * dy * dz);
				b2x[0][j][l] = b2x[0][j][l] + dt * CC0 * (e1[0][j][l] - e1[0][j][l - 1]) * ((a1y[0][j + 1][l] - a1y[0][j][l]) * (a1z[1][j][l] - a1z[0][j][l]) - (a1y[0][j + 1][l] - a1y[0][j][l]) * (a1z[0][j + 1][l] - a1z[0][j][l])) / (dx * dy * dz);
			}
		}
	}
}

void Calculation1b2xPart2(int& Nx, int& Ny, int& Nz, double& dx, double& dy, double& dz, double& dt, double*** e1, double*** a1y, double*** a1z, double*** b1x, double*** b2x, double& CC0, int& begin, int& end)
{
	for (int l = begin; l < end; l++)
		//for (int i = 1; i < Nx - 1; i++)
	{
		for (int j = 0; j < Ny; j++)
		{
			for (int i = 1; i < Nx - 1; i++)
				//for (int l = 0; l < Nz; l++)
			{
				if ((j - 1) < 0 || (j + 1) > (Ny - 1) || (l - 1) < 0 || (l + 1) > (Nz - 1))
				{
					//b2x[i][j][l] = b1x[i][j][l] + dt*CC0*(e1[i][j][l] - tmp)*(-(tmp - a1y[i][j][l])*(tmp - a1z[i][j][l]) + (tmp - a1y[i][j][l])*(tmp - a1z[i][j][l])) / (dx*dy*dz);
					//b2x[i][j][l] = b2x[i][j][l] + dt*CC0*(e1[i][j][l] - tmp)*((tmp - a1y[i][j][l])*(tmp - a1z[i][j][l]) - (tmp - a1y[i][j][l])*(tmp - a1z[i][j][l])) / (dx*dy*dz);
					//b2x[i][j][l] = b2x[i][j][l] + dt*CC0*(e1[i][j][l] - tmp)*((tmp - a1y[i][j][l])*(tmp - a1z[i][j][l]) - (tmp - a1y[i][j][l])*(tmp - a1z[i][j][l])) / (dx*dy*dz);
				}
				else
				{
					b2x[i][j][l] = b1x[i][j][l] + dt * CC0 * (e1[i][j][l] - e1[i - 1][j][l]) * (-(a1y[i][j + 1][l] - a1y[i][j][l]) * (a1z[i][j][l + 1] - a1z[i][j][l]) + (a1y[i][j][l + 1] - a1y[i][j][l]) * (a1z[i][j + 1][l] - a1z[i][j][l])) / (dx * dy * dz);
					b2x[i][j][l] = b2x[i][j][l] + dt * CC0 * (e1[i][j][l] - e1[i][j - 1][l]) * ((a1y[i + 1][j][l] - a1y[i][j][l]) * (a1z[i][j][l + 1] - a1z[i][j][l]) - (a1y[i][j][l + 1] - a1y[i][j][l]) * (a1z[i][j][l + 1] - a1z[i][j][l])) / (dx * dy * dz);
					b2x[i][j][l] = b2x[i][j][l] + dt * CC0 * (e1[i][j][l] - e1[i][j][l - 1]) * ((a1y[i][j + 1][l] - a1y[i][j][l]) * (a1z[i + 1][j][l] - a1z[i][j][l]) - (a1y[i][j + 1][l] - a1y[i][j][l]) * (a1z[i][j + 1][l] - a1z[i][j][l])) / (dx * dy * dz);
				}
			}
		}
	}
}

void EqMotio(double*** a1x, double*** a2x, double*** b1x, double*** b2x, double*** a1z, double*** b1y, double*** b2y, double*** b1z, double*** b2z, double*** a1y, double*** e1, double& dx, double& dy, double& dz, double& dt, double& CC0, int& Nx, int& Ny, int& Nz)
{

	/*std::thread t1(Calculation1b2y, ref(Nx), ref(Ny), ref(Nz), ref(dx), ref(dy), ref(dz), ref(dt), ref(CC0), e1, a1x, a1z, b1y, b2y, ref(new_interv_z[0].begin), ref(new_interv_z[0].end));
	std::thread t2(Calculation1b2y, ref(Nx), ref(Ny), ref(Nz), ref(dx), ref(dy), ref(dz), ref(dt), ref(CC0), e1, a1x, a1z, b1y, b2y, ref(new_interv_z[1].begin), ref(new_interv_z[1].end));
	std::thread t3(Calculation1b2y, ref(Nx), ref(Ny), ref(Nz), ref(dx), ref(dy), ref(dz), ref(dt), ref(CC0), e1, a1x, a1z, b1y, b2y, ref(new_interv_z[2].begin), ref(new_interv_z[2].end));
	std::thread t4(Calculation1b2y, ref(Nx), ref(Ny), ref(Nz), ref(dx), ref(dy), ref(dz), ref(dt), ref(CC0), e1, a1x, a1z, b1y, b2y, ref(new_interv_z[3].begin), ref(new_interv_z[3].end));
	std::thread t5(Calculation1b2y, ref(Nx), ref(Ny), ref(Nz), ref(dx), ref(dy), ref(dz), ref(dt), ref(CC0), e1, a1x, a1z, b1y, b2y, ref(new_interv_z[4].begin), ref(new_interv_z[4].end));
	std::thread t6(Calculation1b2y, ref(Nx), ref(Ny), ref(Nz), ref(dx), ref(dy), ref(dz), ref(dt), ref(CC0), e1, a1x, a1z, b1y, b2y, ref(new_interv_z[5].begin), ref(new_interv_z[5].end));
	std::thread t7(Calculation1b2y, ref(Nx), ref(Ny), ref(Nz), ref(dx), ref(dy), ref(dz), ref(dt), ref(CC0), e1, a1x, a1z, b1y, b2y, ref(new_interv_z[6].begin), ref(new_interv_z[6].end));
	std::thread t8(Calculation1b2y, ref(Nx), ref(Ny), ref(Nz), ref(dx), ref(dy), ref(dz), ref(dt), ref(CC0), e1, a1x, a1z, b1y, b2y, ref(new_interv_z[7].begin), ref(new_interv_z[7].end));
	std::thread t9(Calculation1b2y, ref(Nx), ref(Ny), ref(Nz), ref(dx), ref(dy), ref(dz), ref(dt), ref(CC0), e1, a1x, a1z, b1y, b2y, ref(new_interv_z[8].begin), ref(new_interv_z[8].end));
	std::thread t10(Calculation1b2y, ref(Nx), ref(Ny), ref(Nz), ref(dx), ref(dy), ref(dz), ref(dt), ref(CC0), e1, a1x, a1z, b1y, b2y, ref(new_interv_z[9].begin), ref(new_interv_z[9].end));

	t1.join();
	t2.join();
	t3.join();
	t4.join();
	t5.join();
	t6.join();
	t7.join();
	t8.join();
	t9.join();
	t10.join();*/


	for (int j = 0; j < Ny; j++)//????????????????????????????????????????????????????????????????????????
	{// Или тут рассчитываются только внутренние узлы???????????????????????????????????И почему тут часть кода рабоатет???????????
		for (int l = 0; l < Nz; l++)
		{
			if ((j - 1) < 0 || (j + 1) > (Ny - 1) || (l - 1) < 0 || (l + 1) > (Nz - 1))
			{
				//cout << j << "  " << l << endl;
				//b2x[0][j][l] = b1x[0][j][l] + dt*CC0*(e1[1][j][l] - e1[0][j][l])*(-(tmp - a1y[0][j][l])*(tmp - a1z[0][j][l]) + (tmp - a1y[0][j][l])*(tmp - a1z[0][j][l])) / (dx*dy*dz);
				//b2x[0][j][l] = b2x[0][j][l] + dt*CC0*(e1[0][j][l] - tmp)*((a1y[1][j][l] - a1y[0][j][l])*(tmp - a1z[0][j][l]) - (tmp - a1y[0][j][l])*(tmp - a1z[0][j][l])) / (dx*dy*dz);
				//b2x[0][j][l] = b2x[0][j][l] + dt*CC0*(e1[0][j][l] - tmp)*((tmp - a1y[0][j][l])*(a1z[1][j][l] - a1z[0][j][l]) - (tmp - a1y[0][j][l])*(tmp - a1z[0][j][l])) / (dx*dy*dz);
			}
			else {
				//cout << j << "  " << l << endl;
				b2x[0][j][l] = b1x[0][j][l] + dt * CC0 * (e1[1][j][l] - e1[0][j][l]) * (-(a1y[0][j + 1][l] - a1y[0][j][l]) * (a1z[0][j][l + 1] - a1z[0][j][l]) + (a1y[0][j][l + 1] - a1y[0][j][l]) * (a1z[0][j + 1][l] - a1z[0][j][l])) / (dx * dy * dz);
				b2x[0][j][l] = b2x[0][j][l] + dt * CC0 * (e1[0][j][l] - e1[0][j - 1][l]) * ((a1y[1][j][l] - a1y[0][j][l]) * (a1z[0][j][l + 1] - a1z[0][j][l]) - (a1y[0][j][l + 1] - a1y[0][j][l]) * (a1z[0][j][l + 1] - a1z[0][j][l])) / (dx * dy * dz);
				b2x[0][j][l] = b2x[0][j][l] + dt * CC0 * (e1[0][j][l] - e1[0][j][l - 1]) * ((a1y[0][j + 1][l] - a1y[0][j][l]) * (a1z[1][j][l] - a1z[0][j][l]) - (a1y[0][j + 1][l] - a1y[0][j][l]) * (a1z[0][j + 1][l] - a1z[0][j][l])) / (dx * dy * dz);

			}
		}
	}

	for (int i = 1; i < Nx - 1; i++)
	{
		for (int j = 0; j < Ny; j++)
		{
			for (int l = 0; l < Nz; l++)
			{
				if ((j - 1) < 0 || (j + 1) > (Ny - 1) || (l - 1) < 0 || (l + 1) > (Nz - 1))
				{
					//b2x[i][j][l] = b1x[i][j][l] + dt*CC0*(e1[i][j][l] - tmp)*(-(tmp - a1y[i][j][l])*(tmp - a1z[i][j][l]) + (tmp - a1y[i][j][l])*(tmp - a1z[i][j][l])) / (dx*dy*dz);
					//b2x[i][j][l] = b2x[i][j][l] + dt*CC0*(e1[i][j][l] - tmp)*((tmp - a1y[i][j][l])*(tmp - a1z[i][j][l]) - (tmp - a1y[i][j][l])*(tmp - a1z[i][j][l])) / (dx*dy*dz);
					//b2x[i][j][l] = b2x[i][j][l] + dt*CC0*(e1[i][j][l] - tmp)*((tmp - a1y[i][j][l])*(tmp - a1z[i][j][l]) - (tmp - a1y[i][j][l])*(tmp - a1z[i][j][l])) / (dx*dy*dz);
				}
				else
				{
					b2x[i][j][l] = b1x[i][j][l] + dt * CC0 * (e1[i][j][l] - e1[i - 1][j][l]) * (-(a1y[i][j + 1][l] - a1y[i][j][l]) * (a1z[i][j][l + 1] - a1z[i][j][l]) + (a1y[i][j][l + 1] - a1y[i][j][l]) * (a1z[i][j + 1][l] - a1z[i][j][l])) / (dx * dy * dz);
					b2x[i][j][l] = b2x[i][j][l] + dt * CC0 * (e1[i][j][l] - e1[i][j - 1][l]) * ((a1y[i + 1][j][l] - a1y[i][j][l]) * (a1z[i][j][l + 1] - a1z[i][j][l]) - (a1y[i][j][l + 1] - a1y[i][j][l]) * (a1z[i][j][l + 1] - a1z[i][j][l])) / (dx * dy * dz);
					b2x[i][j][l] = b2x[i][j][l] + dt * CC0 * (e1[i][j][l] - e1[i][j][l - 1]) * ((a1y[i][j + 1][l] - a1y[i][j][l]) * (a1z[i + 1][j][l] - a1z[i][j][l]) - (a1y[i][j + 1][l] - a1y[i][j][l]) * (a1z[i][j + 1][l] - a1z[i][j][l])) / (dx * dy * dz);
				}
			}
		}
	}

	//std::thread t111(Calculation1b2xPart1, ref(Ny), ref(Nz), ref(dx), ref(dy), ref(dz), ref(dt), e1, a1y, a1z, b1x, b2x, ref(CC0), ref(new_interv_z[0].begin), ref(new_interv_z[0].end));
	//std::thread t212(Calculation1b2xPart1, ref(Ny), ref(Nz), ref(dx), ref(dy), ref(dz), ref(dt), e1, a1y, a1z, b1x, b2x, ref(CC0), ref(new_interv_z[1].begin), ref(new_interv_z[1].end));
	//std::thread t313(Calculation1b2xPart1, ref(Ny), ref(Nz), ref(dx), ref(dy), ref(dz), ref(dt), e1, a1y, a1z, b1x, b2x, ref(CC0), ref(new_interv_z[2].begin), ref(new_interv_z[2].end));
	//std::thread t414(Calculation1b2xPart1, ref(Ny), ref(Nz), ref(dx), ref(dy), ref(dz), ref(dt), e1, a1y, a1z, b1x, b2x, ref(CC0), ref(new_interv_z[3].begin), ref(new_interv_z[3].end));
	//std::thread t515(Calculation1b2xPart1, ref(Ny), ref(Nz), ref(dx), ref(dy), ref(dz), ref(dt), e1, a1y, a1z, b1x, b2x, ref(CC0), ref(new_interv_z[4].begin), ref(new_interv_z[4].end));
	//std::thread t616(Calculation1b2xPart1, ref(Ny), ref(Nz), ref(dx), ref(dy), ref(dz), ref(dt), e1, a1y, a1z, b1x, b2x, ref(CC0), ref(new_interv_z[5].begin), ref(new_interv_z[5].end));
	//std::thread t717(Calculation1b2xPart1, ref(Ny), ref(Nz), ref(dx), ref(dy), ref(dz), ref(dt), e1, a1y, a1z, b1x, b2x, ref(CC0), ref(new_interv_z[6].begin), ref(new_interv_z[6].end));
	//std::thread t818(Calculation1b2xPart1, ref(Ny), ref(Nz), ref(dx), ref(dy), ref(dz), ref(dt), e1, a1y, a1z, b1x, b2x, ref(CC0), ref(new_interv_z[7].begin), ref(new_interv_z[7].end));
	//std::thread t919(Calculation1b2xPart1, ref(Ny), ref(Nz), ref(dx), ref(dy), ref(dz), ref(dt), e1, a1y, a1z, b1x, b2x, ref(CC0), ref(new_interv_z[8].begin), ref(new_interv_z[8].end));
	//std::thread t1010(Calculation1b2xPart1, ref(Ny), ref(Nz), ref(dx), ref(dy), ref(dz), ref(dt), e1, a1y, a1z, b1x, b2x, ref(CC0), ref(new_interv_z[9].begin), ref(new_interv_z[9].end));

	//t111.join();
	//t212.join();
	//t313.join();
	//t414.join();
	//t515.join();
	//t616.join();
	//t717.join();
	//t818.join();
	//t919.join();
	//t1010.join();


	//std::thread t121(Calculation1b2xPart2, ref(Nx), ref(Ny), ref(Nz), ref(dx), ref(dy), ref(dz), ref(dt), e1, a1y, a1z, b1x, b2x, ref(CC0), ref(new_interv_z[0].begin), ref(new_interv_z[0].end));
	//std::thread t221(Calculation1b2xPart2, ref(Nx), ref(Ny), ref(Nz), ref(dx), ref(dy), ref(dz), ref(dt), e1, a1y, a1z, b1x, b2x, ref(CC0), ref(new_interv_z[1].begin), ref(new_interv_z[1].end));
	//std::thread t321(Calculation1b2xPart2, ref(Nx), ref(Ny), ref(Nz), ref(dx), ref(dy), ref(dz), ref(dt), e1, a1y, a1z, b1x, b2x, ref(CC0), ref(new_interv_z[2].begin), ref(new_interv_z[2].end));
	//std::thread t421(Calculation1b2xPart2, ref(Nx), ref(Ny), ref(Nz), ref(dx), ref(dy), ref(dz), ref(dt), e1, a1y, a1z, b1x, b2x, ref(CC0), ref(new_interv_z[3].begin), ref(new_interv_z[3].end));
	//std::thread t521(Calculation1b2xPart2, ref(Nx), ref(Ny), ref(Nz), ref(dx), ref(dy), ref(dz), ref(dt), e1, a1y, a1z, b1x, b2x, ref(CC0), ref(new_interv_z[4].begin), ref(new_interv_z[4].end));
	//std::thread t621(Calculation1b2xPart2, ref(Nx), ref(Ny), ref(Nz), ref(dx), ref(dy), ref(dz), ref(dt), e1, a1y, a1z, b1x, b2x, ref(CC0), ref(new_interv_z[5].begin), ref(new_interv_z[5].end));
	//std::thread t721(Calculation1b2xPart2, ref(Nx), ref(Ny), ref(Nz), ref(dx), ref(dy), ref(dz), ref(dt), e1, a1y, a1z, b1x, b2x, ref(CC0), ref(new_interv_z[6].begin), ref(new_interv_z[6].end));
	//std::thread t821(Calculation1b2xPart2, ref(Nx), ref(Ny), ref(Nz), ref(dx), ref(dy), ref(dz), ref(dt), e1, a1y, a1z, b1x, b2x, ref(CC0), ref(new_interv_z[7].begin), ref(new_interv_z[7].end));
	//std::thread t921(Calculation1b2xPart2, ref(Nx), ref(Ny), ref(Nz), ref(dx), ref(dy), ref(dz), ref(dt), e1, a1y, a1z, b1x, b2x, ref(CC0), ref(new_interv_z[8].begin), ref(new_interv_z[8].end));
	//std::thread t1021(Calculation1b2xPart2, ref(Nx), ref(Ny), ref(Nz), ref(dx), ref(dy), ref(dz), ref(dt), e1, a1y, a1z, b1x, b2x, ref(CC0), ref(new_interv_z[9].begin), ref(new_interv_z[9].end));

	//t121.join();
	//t221.join();
	//t321.join();
	//t421.join();
	//t521.join();
	//t621.join();
	//t721.join();
	//t821.join();
	//t921.join();
	//t1021.join();

	//Calculation1b2y(ref(Nx), ref(Ny), ref(Nz), ref(dx), ref(dy), ref(dz), ref(dt), ref(CC0), e1, a1x, a1z, b1y, /*b2y,*/ ref(new_interv_z[0].begin), ref(new_interv_z[0].end));


	//Calculation1b2z(ref(Nx), ref(Ny), ref(Nz), ref(dx), ref(dy), ref(dz), ref(dt), e1, a1x, a1y, b1z, /*b2z,*/ ref(new_interv_z[0].begin), ref(new_interv_z[0].end));

	/*std::thread t11(Calculation1b2z, ref(Nx), ref(Ny), ref(Nz), ref(dx), ref(dy), ref(dz), ref(dt), e1, a1x, a1y, b1z, b2z, ref(new_interv_z[0].begin), ref(new_interv_z[0].end));
	std::thread t12(Calculation1b2z, ref(Nx), ref(Ny), ref(Nz), ref(dx), ref(dy), ref(dz), ref(dt), e1, a1x, a1y, b1z, b2z, ref(new_interv_z[1].begin), ref(new_interv_z[1].end));
	std::thread t13(Calculation1b2z, ref(Nx), ref(Ny), ref(Nz), ref(dx), ref(dy), ref(dz), ref(dt), e1, a1x, a1y, b1z, b2z, ref(new_interv_z[2].begin), ref(new_interv_z[2].end));
	std::thread t14(Calculation1b2z, ref(Nx), ref(Ny), ref(Nz), ref(dx), ref(dy), ref(dz), ref(dt), e1, a1x, a1y, b1z, b2z, ref(new_interv_z[3].begin), ref(new_interv_z[3].end));
	std::thread t15(Calculation1b2z, ref(Nx), ref(Ny), ref(Nz), ref(dx), ref(dy), ref(dz), ref(dt), e1, a1x, a1y, b1z, b2z, ref(new_interv_z[4].begin), ref(new_interv_z[4].end));
	std::thread t16(Calculation1b2z, ref(Nx), ref(Ny), ref(Nz), ref(dx), ref(dy), ref(dz), ref(dt), e1, a1x, a1y, b1z, b2z, ref(new_interv_z[5].begin), ref(new_interv_z[5].end));
	std::thread t17(Calculation1b2z, ref(Nx), ref(Ny), ref(Nz), ref(dx), ref(dy), ref(dz), ref(dt), e1, a1x, a1y, b1z, b2z, ref(new_interv_z[6].begin), ref(new_interv_z[6].end));
	std::thread t18(Calculation1b2z, ref(Nx), ref(Ny), ref(Nz), ref(dx), ref(dy), ref(dz), ref(dt), e1, a1x, a1y, b1z, b2z, ref(new_interv_z[7].begin), ref(new_interv_z[7].end));
	std::thread t19(Calculation1b2z, ref(Nx), ref(Ny), ref(Nz), ref(dx), ref(dy), ref(dz), ref(dt), e1, a1x, a1y, b1z, b2z, ref(new_interv_z[8].begin), ref(new_interv_z[8].end));
	std::thread t20(Calculation1b2z, ref(Nx), ref(Ny), ref(Nz), ref(dx), ref(dy), ref(dz), ref(dt), e1, a1x, a1y, b1z, b2z, ref(new_interv_z[9].begin), ref(new_interv_z[9].end));*/

	//t11.join();
	//t12.join();
	//t13.join();
	//t14.join();
	//t15.join();
	//t16.join();
	//t17.join();
	//t18.join();
	//t19.join();
	//t20.join();



























	
}

void Calculation1(double*** a1x, double*** a2x, double*** b1x, double*** b2x, double*** a1z, double*** b1y, double*** b2y, double*** b1z, double*** b2z, double*** a1y, double*** e1, double& dx, double& dy, double& dz, double& dt, double& CC0, int& Nx, int& Ny, int& Nz, vector<Interval>& new_interv_z)
{
	// Calculation 1
	// solving the equation of motion

	for (int j = 0; j < Ny; j++)//????????????????????????????????????????????????????????????????????????
	{// Или тут рассчитываются только внутренние узлы???????????????????????????????????И почему тут часть кода рабоатет???????????
		for (int l = 0; l < Nz; l++)
		{
			if ((j - 1) < 0 || (j + 1) > (Ny - 1) || (l - 1) < 0 || (l + 1) > (Nz - 1))
			{
				//cout << j << "  " << l << endl;
				//b2x[0][j][l] = b1x[0][j][l] + dt*CC0*(e1[1][j][l] - e1[0][j][l])*(-(tmp - a1y[0][j][l])*(tmp - a1z[0][j][l]) + (tmp - a1y[0][j][l])*(tmp - a1z[0][j][l])) / (dx*dy*dz);
				//b2x[0][j][l] = b2x[0][j][l] + dt*CC0*(e1[0][j][l] - tmp)*((a1y[1][j][l] - a1y[0][j][l])*(tmp - a1z[0][j][l]) - (tmp - a1y[0][j][l])*(tmp - a1z[0][j][l])) / (dx*dy*dz);
				//b2x[0][j][l] = b2x[0][j][l] + dt*CC0*(e1[0][j][l] - tmp)*((tmp - a1y[0][j][l])*(a1z[1][j][l] - a1z[0][j][l]) - (tmp - a1y[0][j][l])*(tmp - a1z[0][j][l])) / (dx*dy*dz);
			}
			else {
				//cout << j << "  " << l << endl;
				b2x[0][j][l] = b1x[0][j][l] + dt * CC0 * (e1[1][j][l] - e1[0][j][l]) * (-(a1y[0][j + 1][l] - a1y[0][j][l]) * (a1z[0][j][l + 1] - a1z[0][j][l]) + (a1y[0][j][l + 1] - a1y[0][j][l]) * (a1z[0][j + 1][l] - a1z[0][j][l])) / (dx * dy * dz);
				b2x[0][j][l] = b2x[0][j][l] + dt * CC0 * (e1[0][j][l] - e1[0][j - 1][l]) * ((a1y[1][j][l] - a1y[0][j][l]) * (a1z[0][j][l + 1] - a1z[0][j][l]) - (a1y[0][j][l + 1] - a1y[0][j][l]) * (a1z[0][j][l + 1] - a1z[0][j][l])) / (dx * dy * dz);
				b2x[0][j][l] = b2x[0][j][l] + dt * CC0 * (e1[0][j][l] - e1[0][j][l - 1]) * ((a1y[0][j + 1][l] - a1y[0][j][l]) * (a1z[1][j][l] - a1z[0][j][l]) - (a1y[0][j + 1][l] - a1y[0][j][l]) * (a1z[0][j + 1][l] - a1z[0][j][l])) / (dx * dy * dz);

			}
		}
	}

	for (int i = 1; i < Nx - 1; i++)
	{
		for (int j = 0; j < Ny; j++)
		{
			for (int l = 0; l < Nz; l++)
			{
				if ((j - 1) < 0 || (j + 1) > (Ny - 1) || (l - 1) < 0 || (l + 1) > (Nz - 1))
				{
					//b2x[i][j][l] = b1x[i][j][l] + dt*CC0*(e1[i][j][l] - tmp)*(-(tmp - a1y[i][j][l])*(tmp - a1z[i][j][l]) + (tmp - a1y[i][j][l])*(tmp - a1z[i][j][l])) / (dx*dy*dz);
					//b2x[i][j][l] = b2x[i][j][l] + dt*CC0*(e1[i][j][l] - tmp)*((tmp - a1y[i][j][l])*(tmp - a1z[i][j][l]) - (tmp - a1y[i][j][l])*(tmp - a1z[i][j][l])) / (dx*dy*dz);
					//b2x[i][j][l] = b2x[i][j][l] + dt*CC0*(e1[i][j][l] - tmp)*((tmp - a1y[i][j][l])*(tmp - a1z[i][j][l]) - (tmp - a1y[i][j][l])*(tmp - a1z[i][j][l])) / (dx*dy*dz);
				}
				else
				{
					b2x[i][j][l] = b1x[i][j][l] + dt * CC0 * (e1[i][j][l] - e1[i - 1][j][l]) * (-(a1y[i][j + 1][l] - a1y[i][j][l]) * (a1z[i][j][l + 1] - a1z[i][j][l]) + (a1y[i][j][l + 1] - a1y[i][j][l]) * (a1z[i][j + 1][l] - a1z[i][j][l])) / (dx * dy * dz);
					b2x[i][j][l] = b2x[i][j][l] + dt * CC0 * (e1[i][j][l] - e1[i][j - 1][l]) * ((a1y[i + 1][j][l] - a1y[i][j][l]) * (a1z[i][j][l + 1] - a1z[i][j][l]) - (a1y[i][j][l + 1] - a1y[i][j][l]) * (a1z[i][j][l + 1] - a1z[i][j][l])) / (dx * dy * dz);
					b2x[i][j][l] = b2x[i][j][l] + dt * CC0 * (e1[i][j][l] - e1[i][j][l - 1]) * ((a1y[i][j + 1][l] - a1y[i][j][l]) * (a1z[i + 1][j][l] - a1z[i][j][l]) - (a1y[i][j + 1][l] - a1y[i][j][l]) * (a1z[i][j + 1][l] - a1z[i][j][l])) / (dx * dy * dz);
				}
			}
		}
	}

	for (int i = 0; i < Nx; i++)
	{
		for (int j = 0; j < Ny; j++)
		{
			for (int l = 0; l < Nz; l++)
			{
				if ((i - 1) < 0 || (i + 1) > (Nx - 1) || (j - 1) < 0 || (j + 1) > (Ny - 1) || (l - 1) < 0 || (l + 1) > (Nz - 1))
				{
					//b2y[i][j][l] = b1y[i][j][l] + dt*CC0*(e1[i][j][l] - tmp)*((tmp - a1x[i][j][l])*(tmp - a1z[i][j][l]) - (tmp - a1x[i][j][l])*(tmp - a1z[i][j][l])) / (dx*dy*dz);
					//b2y[i][j][l] = b2y[i][j][l] + dt*CC0*(e1[i][j][l] - tmp)*((tmp - a1x[i][j][l])*(tmp - a1z[i][j][l]) - (tmp - a1x[i][j][l])*(tmp - a1z[i][j][l])) / (dx*dy*dz);
					//b2y[i][j][l] = b2y[i][j][l] + dt*CC0*(e1[i][j][l] - tmp)*((tmp - a1x[i][j][l])*(tmp - a1z[i][j][l]) - (tmp - a1x[i][j][l])*(tmp - a1z[i][j][l])) / (dx*dy*dz);
				}
				else
				{
					b2y[i][j][l] = b1y[i][j][l] + dt * CC0 * (e1[i][j][l] - e1[i - 1][j][l]) * ((a1x[i][j + 1][l] - a1x[i][j][l]) * (a1z[i][j][l + 1] - a1z[i][j][l]) - (a1x[i][j][l + 1] - a1x[i][j][l]) * (a1z[i][j + 1][l] - a1z[i][j][l])) / (dx * dy * dz);
					b2y[i][j][l] = b2y[i][j][l] + dt * CC0 * (e1[i][j][l] - e1[i][j - 1][l]) * ((a1x[i][j][l + 1] - a1x[i][j][l]) * (a1z[i + 1][j][l] - a1z[i][j][l]) - (a1x[i + 1][j][l] - a1x[i][j][l]) * (a1z[i][j][l + 1] - a1z[i][j][l])) / (dx * dy * dz);
					b2y[i][j][l] = b2y[i][j][l] + dt * CC0 * (e1[i][j][l] - e1[i][j][l - 1]) * ((a1x[i + 1][j][l] - a1x[i][j][l]) * (a1z[i][j + 1][l] - a1z[i][j][l]) - (a1x[i][j + 1][l] - a1x[i][j][l]) * (a1z[i + 1][j][l] - a1z[i][j][l])) / (dx * dy * dz);
				}
			}
		}
	}

	for (int i = 0; i < Nx; i++)
	{
		for (int j = 0; j < Ny; j++)
		{
			for (int l = 0; l < Nz; l++)
			{
				if ((i - 1) < 0 || (i + 1) > (Nx - 1) || (j - 1) < 0 || (j + 1) > (Ny - 1) || (l - 1) < 0 || (l + 1) > (Nz - 1))
				{
					//b2z[i][j][l] = b1z[i][j][l] + dt*(e1[i][j][l] - tmp)*((tmp - a1x[i][j][l])*(tmp - a1y[i][j][l]) - (tmp - a1x[i][j][l])*(tmp - a1y[i][j][l])) / (dx*dy*dz);
					//b2z[i][j][l] = b2z[i][j][l] + dt*(e1[i][j][l] - tmp)*((tmp - a1x[i][j][l])*(tmp - a1y[i][j][l]) - (tmp - a1x[i][j][l])*(tmp - a1y[i][j][l])) / (dx*dy*dz);
					//b2z[i][j][l] = b2z[i][j][l] + dt*(e1[i][j][l] - tmp)*((tmp - a1x[i][j][l])*(tmp - a1y[i][j][l]) - (tmp - a1x[i][j][l])*(tmp - a1y[i][j][l])) / (dx*dy*dz);
				}
				else
				{
					b2z[i][j][l] = b1z[i][j][l] + dt * (e1[i][j][l] - e1[i - 1][j][l]) * ((a1x[i][j][l + 1] - a1x[i][j][l]) * (a1y[i][j + 1][l] - a1y[i][j][l]) - (a1x[i][j + 1][l] - a1x[i][j][l]) * (a1y[i][j][l + 1] - a1y[i][j][l])) / (dx * dy * dz);
					b2z[i][j][l] = b2z[i][j][l] + dt * (e1[i][j][l] - e1[i][j - 1][l]) * ((a1x[i + 1][j][l] - a1x[i][j][l]) * (a1y[i][j][l + 1] - a1y[i][j][l]) - (a1x[i][j][l + 1] - a1x[i][j][l]) * (a1y[i + 1][j][l] - a1y[i][j][l])) / (dx * dy * dz);
					b2z[i][j][l] = b2z[i][j][l] + dt * (e1[i][j][l] - e1[i][j][l - 1]) * ((a1x[i][j + 1][l] - a1x[i][j][l]) * (a1y[i + 1][j][l] - a1y[i][j][l]) - (a1x[i + 1][j][l] - a1x[i][j][l]) * (a1y[i][j + 1][l] - a1y[i][j][l])) / (dx * dy * dz);
				}
			}
		}
	}
}


void a2xyz(int& Nx, int& Ny, int& Nz, double& dt, double& CC0, double*** a1x, double*** a1y, double*** a1z, double*** a2x, double*** a2y, double*** a2z, double*** b2x, double*** b2y, double*** b2z, int& begin, int& end)
{
	/*for (int i = 0; i < Nx; i++)*/
	for (int l = begin; l < end; l++)
		//for (int l = 0; l < Nz; l++)
	{
		for (int j = 0; j < Ny; j++)
		{
			for (int i = 0; i < Nx; i++)
				//for (int l = 0; l < Nz; l++)
			{
				//calculation of the change in the Euler coordinates of Lagrangian particles
				a2x[i][j][l] = a1x[i][j][l] + dt * b2x[i][j][l] / CC0;
				a2y[i][j][l] = a1y[i][j][l] + dt * b2y[i][j][l] / CC0;
				a2z[i][j][l] = a1z[i][j][l] + dt * b2z[i][j][l];
			}
		}
	}
}

void funV2(int& Nx, int& Ny, int& Nz, double& dx, double& dy, double& dz, double*** a1x, double*** a1y, double*** a1z, double*** a2x, double*** a2y, double*** a2z, double*** V2, int& begin, int& end)
{
	for (int l = begin; l < end; l++)
		//for (int i = 1; i < Nx - 1; i++)
	{
		for (int j = 1; j < Ny - 1; j++)
		{
			//for (int l = 1; l < Nz - 1; l++)
			for (int i = 1; i < Nx - 1; i++)
			{
				// Calculation of the continuity equation
				V2[i][j][l] = (a2x[i + 1][j][l] - a2x[i][j][l]) * ((a2y[i][j + 1][l] - a2y[i][j][l]) * (a2z[i][j][l + 1] - a2z[i][j][l]) - (a2y[i][j][l + 1] - a2y[i][j][l]) * (a2z[i][j + 1][l] - a2z[i][j][l])) / (dx * dy * dz);
				V2[i][j][l] = V2[i][j][l] - (a2y[i + 1][j][l] - a2y[i][j][l]) * ((a2x[i][j + 1][l] - a2x[i][j][l]) * (a2z[i][j][l + 1] - a2z[i][j][l]) - (a2x[i][j][l + 1] - a2x[i][j][l]) * (a2z[i][j + 1][l] - a2z[i][j][l])) / (dx * dy * dz);
				V2[i][j][l] = V2[i][j][l] + (a2z[i + 1][j][l] - a2z[i][j][l]) * ((a2x[i][j + 1][l] - a2x[i][j][l]) * (a2y[i][j][l + 1] - a2y[i][j][l]) - (a2x[i][j][l + 1] - a2x[i][j][l]) * (a2y[i][j + 1][l] - a2y[i][j][l])) / (dx * dy * dz);
			}
		}
	}
}

void Calculation2(double*** a1x, double*** a1y, double*** a1z, double*** a2x, double*** a2y, double*** a2z, double*** b2x, double*** b2y, double*** b2z, double*** V2, double*** e2, double*** tmpe2, double*** tmpi2, double& dx, double& dy, double& dz, double& dt, double& CC0, int& Nx, int& Ny, int& Nz, double& V0, vector<Melting>& Melt_metal, Metal mt, Parametrs& param, Splayn& spl_C_e_on_T, Splayn& spl_Te, Splayn& spl_Ti, double*** Te_acoustic, double*** Ti_acoustic, vector<Interval>& new_interv_z, Splayn spl_C_l_on_T)
{// Calculation of 
 // Calculatuion 2

	for (int i = 0; i < Nx; i++)
	{
		for (int j = 0; j < Ny; j++)
		{
			for (int l = 0; l < Nz; l++)
			{
				//calculation of the change in the Euler coordinates of Lagrangian particles
				a2x[i][j][l] = a1x[i][j][l] + dt * b2x[i][j][l] / CC0;
				a2y[i][j][l] = a1y[i][j][l] + dt * b2y[i][j][l] / CC0;
				a2z[i][j][l] = a1z[i][j][l] + dt * b2z[i][j][l];
			}
		}
	}

	for (int i = 1; i < Nx - 1; i++)
	{
		for (int j = 1; j < Ny - 1; j++)
		{
			for (int l = 1; l < Nz - 1; l++)
			{
				// Calculation of the continuity equation
				V2[i][j][l] = (a2x[i + 1][j][l] - a2x[i][j][l]) * ((a2y[i][j + 1][l] - a2y[i][j][l]) * (a2z[i][j][l + 1] - a2z[i][j][l]) - (a2y[i][j][l + 1] - a2y[i][j][l]) * (a2z[i][j + 1][l] - a2z[i][j][l])) / (dx * dy * dz);
				V2[i][j][l] = V2[i][j][l] - (a2y[i + 1][j][l] - a2y[i][j][l]) * ((a2x[i][j + 1][l] - a2x[i][j][l]) * (a2z[i][j][l + 1] - a2z[i][j][l]) - (a2x[i][j][l + 1] - a2x[i][j][l]) * (a2z[i][j + 1][l] - a2z[i][j][l])) / (dx * dy * dz);
				V2[i][j][l] = V2[i][j][l] + (a2z[i + 1][j][l] - a2z[i][j][l]) * ((a2x[i][j + 1][l] - a2x[i][j][l]) * (a2y[i][j][l + 1] - a2y[i][j][l]) - (a2x[i][j][l + 1] - a2x[i][j][l]) * (a2y[i][j + 1][l] - a2y[i][j][l])) / (dx * dy * dz);
			}
		}
	}

	for (int j = 0; j < Nx; j++)
	{
		for (int l = 0; l < Ny; l++)
		{
			V2[j][l][Nz - 1] = V0;
		}
	}

	for (int j = 0; j < Ny; j++)
	{
		for (int l = 0; l < Nz; l++)
		{
			e2[0][j][l] = 1e-16;
		}
	}


	// этот блок циклов необходим для того, чтобы распаралелить по оси z
	for (int i = 0; i < Nx; i++)
	{
		for (int j = 0; j < Ny; j++)
		{
			spl_Te.InterpolateFast(1, tmpe2, i, j);
			spl_Ti.InterpolateFast(1, tmpi2, i, j);
			for (int k = 0; k < Nz; k++)
			{
				Te_acoustic[i][j][k] = spl_Te.GetY(k * dz);
				Ti_acoustic[i][j][k] = spl_Ti.GetY(k * dz);
			}
		}
	}

	//for (int i = 0; i < Nx; i++) // Nz_heat, dz_heat - не нужны
	//{
	//	for (int j = 0; j < Ny; j++)// циклы по узлам акустики
	//	{
	//		for (int k = 1; k < Nz; k++)//l=1
	//		{
	//			// почему тако вид ур-я состояния???
	//			//   e2[i,j,l]:=(1-V2[i,j,l])+at1i*(tmpi2[i,j,l])/(V2[i,j,l])+at1e*(tmpe2[i,j,l]);

	//			///e2[i][j][k] = (1 - V2[i][j][k]) + at1i * (tmpi2[i][j][k]) + at1e * (tmpe2[i][j][k]);

	//			e2[i][j][k] = (1. - V2[i][j][k]) +
	//				(Melt_metal[mt].gi * param.T00 * (Dependence_C_l_on_T(mt, param.T00 * Ti_acoustic[i][j][k]) / (100. * 100. * 100.)) / (Melt_metal[mt].Density * Melt_metal[mt].Density * pow(Melt_metal[mt].u0, 2) * V0 * 1e-6)) *
	//				(Ti_acoustic[i][j][k] - 1.) / V2[i][j][k] + (Melt_metal[mt].ge * param.T00 * (spl_C_e_on_T.GetY(param.T00 * Te_acoustic[i][j][k]) / (100. * 100. * 100.)) / (Melt_metal[mt].Density * pow(Melt_metal[mt].u0, 2) * 1e-6)) *
	//				(Te_acoustic[i][j][k] - 1.);
	//		}
	//	}
	//}

	double delta = 1.;

	for (int i = 0; i < Nx; i++) // Nz_heat, dz_heat - не нужны
	{
		for (int j = 0; j < Ny; j++)// циклы по узлам акустики
		{
			for (int k = 1; k < Nz; k++)//l=1
			{
				// добавить плотность и скорость звука в зависимоси от фазы (тверд, жидк) 
				if ((Ti_acoustic[i][j][k] * param.T00) < (Melt_metal[mt].T_melting - delta))
				{
					
						e2[i][j][k] = (1. - V2[i][j][k]) +
							(Melt_metal[mt].gi * param.T00 * (Dependence_C_l_on_T(mt, param.T00 * Ti_acoustic[i][j][k]) / (100. * 100. * 100.)) / (Melt_metal[mt].Density * Melt_metal[mt].Density * pow(Melt_metal[mt].u0, 2) * V0 * 1e-6)) *
							(Ti_acoustic[i][j][k] - 1.) / V2[i][j][k] + (Melt_metal[mt].ge * param.T00 * (spl_C_e_on_T.GetY(param.T00 * Te_acoustic[i][j][k]) / (100. * 100. * 100.)) / (Melt_metal[mt].Density * pow(Melt_metal[mt].u0, 2) * 1e-6)) *
							(Te_acoustic[i][j][k] - 1.);
				}

				if ((Ti_acoustic[i][j][k] * param.T00) >= (Melt_metal[mt].T_melting - delta) && (Ti_acoustic[i][j][k] * param.T00) <= (Melt_metal[mt].T_melting + delta))
				{
					//Ci[J/cm3 K] = ((Dependence_C_l_on_T(mt, param.T00 * tmpi2[i][j][k]) / (100. * 100. * 100.) + (Melt_metal[mt].Q_fusion * Melt_metal[mt].Density) / (2 * delta)))
					// spl_C_l_on_T.GetY(param.T00 * Ti_acoustic[i][j][k]) 
					// spl_Ti.GetY(k * dz) = tmpi2

					/*e2[i][j][k] = (1. - V2[i][j][k]) +
						(Melt_metal[mt].gi * param.T00 * (((((Dependence_C_l_on_T(mt, param.T00* tmpi2[i][j][k]) / (100. * 100. * 100.) + (Melt_metal[mt].Q_fusion * Melt_metal[mt].Density) / (2 * delta)))))) / (pow((Melt_metal[mt].Density + Melt_metal[mt].DensityLiquid) / 2, 2) * pow((Melt_metal[mt].u0 + Melt_metal[mt].u0_Liquid) / 2, 2) * V0 * 1e-6)) *
						(Ti_acoustic[i][j][k] - 1.) / V2[i][j][k] + (Melt_metal[mt].ge * param.T00 * (spl_C_e_on_T.GetY(param.T00 * Te_acoustic[i][j][k]) / (100. * 100. * 100.)) / (((Melt_metal[mt].Density + Melt_metal[mt].DensityLiquid) / 2) * pow((Melt_metal[mt].u0 + Melt_metal[mt].u0_Liquid) / 2, 2) * 1e-6)) *
						(Te_acoustic[i][j][k] - 1.);*/

					e2[i][j][k] = (1. - V2[i][j][k]) +
						(Melt_metal[mt].gi * param.T00 * (((( (spl_C_l_on_T.GetY(param.T00* Ti_acoustic[i][j][k])))))) / (pow((Melt_metal[mt].Density + Melt_metal[mt].DensityLiquid) / 2, 2) * pow((Melt_metal[mt].u0 + Melt_metal[mt].u0_Liquid) / 2, 2) * V0 * 1e-6)) *
						(Ti_acoustic[i][j][k] - 1.) / V2[i][j][k] + (Melt_metal[mt].ge * param.T00 * (spl_C_e_on_T.GetY(param.T00 * Te_acoustic[i][j][k]) / (100. * 100. * 100.)) / (((Melt_metal[mt].Density + Melt_metal[mt].DensityLiquid) / 2) * pow((Melt_metal[mt].u0 + Melt_metal[mt].u0_Liquid) / 2, 2) * 1e-6)) *
						(Te_acoustic[i][j][k] - 1.);
				}

				if ((Ti_acoustic[i][j][k] * param.T00) > (Melt_metal[mt].T_melting + delta))
				{
					e2[i][j][k] = (1. - V2[i][j][k]) +
						(Melt_metal[mt].gi * param.T00 * ((Dependence_C_l_on_T(mt, param.T00 * Ti_acoustic[i][j][k]) / (100. * 100. * 100.))) / (Melt_metal[mt].DensityLiquid * Melt_metal[mt].DensityLiquid * pow(Melt_metal[mt].u0_Liquid, 2) * V0 * 1e-6)) *
						(Ti_acoustic[i][j][k] - 1.) / V2[i][j][k] + (Melt_metal[mt].ge * param.T00 * (spl_C_e_on_T.GetY(param.T00 * Te_acoustic[i][j][k]) / (100. * 100. * 100.)) / (Melt_metal[mt].DensityLiquid * pow(Melt_metal[mt].u0_Liquid, 2) * 1e-6)) *
						(Te_acoustic[i][j][k] - 1.);
				}

			}
		}
	}


}

void MainProcedure(Metal mt, TypeBeam tbeam, double kte_kte, double ro0, double CeCe, double CiCi, double gammagamma, double g_e, double g_i, double u00, double tptp, double P00, double*** V2, double*** a1y, double*** a2y, double*** a1z, double*** a2z, double*** a2x, double*** a1x, double*** b1x, double*** b2x, double*** b1y, double*** b2y, double*** b1z, double*** b2z, double*** e2, double*** e1, double*** tmpe0, double*** tmpe1, double*** tmpe2, double*** tmpi0, double*** tmpi1, double*** tmpi2, int Nx_heat, int  Ny_heat, int  Nz_heat, int  Nx_acoustic, int  Ny_acoustic, int  Nz_acoustic, Splayn spl_C_l_on_T, double*** Te_acoustic, double ***Ti_acoustic, string current_namefile)    // davlenie, gorizontal'naja skorost', vertikal'naja skorost', temperatura
{
	// parametrs of laser and area

	double r0 = 1e-2; // radius of light beam, cm
	double kabs = 5e+5; // absorption coeff, cm-1
	double P0 = P00;//1e+8; // input intensity, W/cm2
	double tp = tptp;// 1e-13; // pulse duration, s
	cout << " tp, s = " << tp << endl;
	cout << " tp, fs = " << tp * 1e+15 << endl;
	double xy0 = 1e-1; // transverse size, cm
	double z0 = 1e-4;  // longsize, cm
	cout << " Fluence = " << tp * P0 << endl;

	// parametrs of metals

	///double kte = kte_kte;//3.115;  // heat cond electron, W/cmK
	//double roe = roe_roe;// 4.56e-5; // electron density, g/cm3
	//double Ce = CeCe;// 1.2e+3;  // electron heat capacity, J/gK
	//double roi = roiroi;// 19.32;  // ion density, g/cm3
	//double Ci = CiCi;// 0.132;   // ion heat capacity, J/gK
	double T00 = 300;    // char temperature, K
	double u0 = u00;// 3.24e+5; // velocity of sound, cm/s
	double t0 = 1 / (kabs * u0);  // char time, s
	double beta = t0 / tp;  // parameter of pulse
	///double gamma = gammagamma;// 0.25e+11;  //  gamma el-phonon, W/cm3K
	double ge = g_e; 
	double gi = g_i;

	cout << endl;

	///double A1 = kte * t0 / (Ce * r0 * r0);
	double A2 = kabs * kabs * r0 * r0;
	///double B1 = kabs * P0 * t0 * beta / (Ce * T00);
	///double CC1 = gamma * t0 / (Ce);
	///double CC2 = gamma * t0 / (Ci);
	Parametrs param;
	param.beta = beta;
	param.kabs = kabs;
	param.P0 = P0;
	param.r0 = r0;
	param.t0 = t0;
	param.T00 = T00;

	string metall, type_beam;
	vector<Melting> Melt_metal(4);

	switch (mt)
	{
	case Au:
		metall = "Au";
		Melt_metal[mt].mt = mt;
		Melt_metal[mt].Q_fusion = 63.7; // Дж/г
		Melt_metal[mt].T_melting = 1338.; // K
		Melt_metal[mt].Density = 19.3; // g/cm3
		Melt_metal[mt].DensityLiquid = 17.; // g/cm3
		Melt_metal[mt].ge = g_e;
		Melt_metal[mt].gi = g_i;
		Melt_metal[mt].u0 = u00;
		Melt_metal[mt].u0_Liquid = 2.567e+5;
		break;
	case Al:
		metall = "Al";
		Melt_metal[mt].mt = mt;
		Melt_metal[mt].Q_fusion = 397.; // Дж/г
		Melt_metal[mt].T_melting = 934.; // K
		Melt_metal[mt].Density = 2.7; // g/cm3
		Melt_metal[mt].ge = g_e;
		Melt_metal[mt].gi = g_i;
		Melt_metal[mt].u0 = u00;
		break;
	case Cu:
		metall = "Cu";
		Melt_metal[mt].mt = mt;
		Melt_metal[mt].Q_fusion = 209.; // Дж/г
		Melt_metal[mt].T_melting = 1358.; // K
		Melt_metal[mt].Density = 8.96; // g/cm3
		Melt_metal[mt].ge = g_e;
		Melt_metal[mt].gi = g_i;
		Melt_metal[mt].u0 = u00;
		break;
	case Ni:
		metall = "Ni";
		Melt_metal[mt].mt = mt;
		Melt_metal[mt].Q_fusion = 298.; // Дж/г
		Melt_metal[mt].T_melting = 1726.; // K
		Melt_metal[mt].Density = 8.9; // g/cm3
		Melt_metal[mt].ge = g_e;
		Melt_metal[mt].gi = g_i;
		Melt_metal[mt].u0 = u00;
		break;
	}

	switch (tbeam)
	{
	case Gauss:
		type_beam = "Gauss";
		break;
	case Vortex:
		type_beam = "Vortex";
		break;
	}

	/////////////////////////////////////////

	double V0 = 1;
	double CC0 = 1 / (kabs * r0);

	//double p0v = ro0 * 1000 * pow(u0 * 0.01, 2) * 1e-5;//202.8e+3;   // c00*c00*ro0 - normirovka davlenija bar (r0 *u0* u0) (perevod v SI i v bar)
	double p0v = ro0 * u0 * u0 * (1e-6); // 
	cout << " p0v = " << p0v << endl;
	double temp = ((Melt_metal[mt].Density + Melt_metal[mt].DensityLiquid) / 2) * pow((Melt_metal[mt].u0 + Melt_metal[mt].u0_Liquid) / 2, 2) * (1e-6);
	cout << " p0v 1 = " << temp << endl;
	temp = Melt_metal[mt].DensityLiquid * Melt_metal[mt].u0_Liquid * Melt_metal[mt].u0_Liquid * (1e-6);
	cout << " p0v = " << temp << endl;
	double sigma = 1400; // bar Au - предел прочности
	bool melting_split = false; // индикатор того, что (произошел/не произошло) плавление и откол
	bool melting = false;
	bool split = false;
	//system("pause");

	// откуда вычисляются???????
	//double at1i = gi * T00 * Ci / p0v;// 1.11* 0.132 * 300 * 19.32 / p0v;  // parameter vo vtorom (teplom) slagaemom uravnenija sostojanija  dlja Au
	//double at1e = ge * T00 * Ce / p0v;//1.5 * 300 * 4.56 * 1.2 * 0.01 / p0v;   // parameter vo vtorom (teplom) slagaemom uravnenija sostojanija  dlja Au

	//int Nx = 150;//100; // chislo uzlov po x // 100
	//int Ny = 150;// 100; // chislo uzlov po y //100
	//int Nz = 200;// 100; // chislo uzlov po z // 200
	//cout << Nx << "   " << Ny << "   " << Nz << endl;

	////parametry v uravnenijah

	//double dxy = xy0 / (r0 * Nx);
	//double dx = dxy;
	//double dy = dxy;

	//double dz;

	//if (metall == "Ni")
	//{
	//	dz = z0 * kabs / Nz / 10; // для никеля
	//}

	//if (metall == "Au")
	//{
	//	dz = z0 * kabs / Nz / 2; // для никеля
	//}


	cout << Nx_heat << "   " << Ny_heat << "   " << Nz_heat << endl;
	cout << Nx_acoustic << "   " << Ny_acoustic << "   " << Nz_acoustic << endl;

	double dxy_heat = xy0 / (r0 * Nx_heat);
	double dx_heat = dxy_heat;
	double dy_heat = dxy_heat;
	double dz_heat = z0 * kabs / Nz_heat;
	double dxy_acoustic = xy0 / (r0 * Nx_acoustic);
	double dx_acoustic = dxy_acoustic;
	double dy_acoustic = dxy_acoustic;
	double dz_acoustic = z0 * kabs / Nz_acoustic;
	cout << dx_heat << "  " << dy_heat << "  " << dz_heat << endl;
	cout << dx_acoustic << "  " << dy_acoustic << "  " << dz_acoustic << endl;
	cout << " Razmern shag x y (mkm) heat = " << 1e+4* r0* dx_heat << endl;
	cout << " Razmern shag z (mkm) heat = " << 1e+4* dz_heat / kabs << endl;
	cout << " Razmern shag x y (mkm) acoustic = " << 1e+4* r0* dx_acoustic << endl;
	cout << " Razmern shag z (mkm) acoustic = " << 1e+4* dz_acoustic / kabs;
	cout << endl << endl;
	cout << " Razmern shag x y (nm) heat = " << 1e+4* r0* dx_heat * 1000 << endl;
	cout << " Razmern shag z (nm) heat = " << 1e+4* dz_heat / kabs * 1000 << endl;
	cout << " Razmern shag x y (nm) acoustic = " << 1e+4* r0* dx_acoustic * 1000 << endl;
	cout << " Razmern shag z (nm) acoustic = " << 1e+4* dz_acoustic / kabs * 1000;
	cout << endl;


	//cout << dx << "  " << dy << "  " << dz << endl;
	//  dz:=0.25;
	double dt = 2e-4;
	int n = 1;
	double tt = 0;
	// 32 переменные double

	Splayn spl_G_e_on_T;// = Calculation_Interpolation("Ge_" + metall + "_new.txt");
	Splayn spl_C_e_on_T;// = Calculation_Interpolation("Ce_" + metall + "_new.txt");
	spl_C_e_on_T.Calculation_Interpolation("Ce_" + metall + "_new.txt");
	spl_G_e_on_T.Calculation_Interpolation("Ge_" + metall + "_new.txt");
	Splayn spl_Te;
	Splayn spl_Ti;

	// это для того, чтобы сетку координат один раз задать для интерполяции
	VecD X;
	for (int i = 0; i < Nz_heat; i++)//dz_heat
	{
		X.push_back(i * dz_heat); // если будет неравномерная сетка, то нужно одпраить в методе способ заполнения диагонали
	}

	spl_Te.SetInitialData(X, 1);
	spl_Ti.SetInitialData(X, 1);
	X.clear();

	// initial conditions

	for (int i = 0; i < Nx_heat; i++)
	{
		for (int j = 0; j < Ny_heat; j++)
		{
			for (int k = 0; k < Nz_heat; k++)
			{
				tmpe0[i][j][k] = 1.0;// 1e-16;
				tmpe1[i][j][k] = 1e-16;
				tmpe2[i][j][k] = 1e-16;
				tmpi0[i][j][k] = 1.0;// 1e-16;
				tmpi1[i][j][k] = 1.0;// 1e-16;
				tmpi2[i][j][k] = 1e-16;
			}
		}
	}

	for (int i = 0; i < Nx_acoustic; i++)
	{
		for (int j = 0; j < Ny_acoustic; j++)
		{
			for (int k = 0; k < Nz_acoustic; k++)
			{
				e1[i][j][k] = 1e-16;
				e2[i][j][k] = 1e-16;
				a1x[i][j][k] = (i - 1) * dx_acoustic;
				a2x[i][j][k] = 1e-16;
				a1y[i][j][k] = (j - 1) * dy_acoustic;
				a2y[i][j][k] = 1e-16;
				a1z[i][j][k] = (k - 1) * dz_acoustic;
				a2z[i][j][k] = 1e-16;
				b1x[i][j][k] = 1e-16;
				b2x[i][j][k] = 1e-16;
				b1y[i][j][k] = 1e-16;
				b2y[i][j][k] = 1e-16;
				b1z[i][j][k] = 1e-16;
				b2z[i][j][k] = 1e-16;
				V2[i][j][k] = V0;
			}
		}
	}
	
	//для распаралелл тепловой задачи
	vector<Interval> new_interv = {
		Interval{ 1, Nz_heat / 10 },
		Interval{ Nz_heat / 10, 2 * Nz_heat / 10 },
		Interval{ 2 * Nz_heat / 10, 3 * Nz_heat / 10 },
		Interval{ 3 * Nz_heat / 10, 4 * Nz_heat / 10 },
		Interval{ 4 * Nz_heat / 10, 5 * Nz_heat / 10 },
		Interval{ 5 * Nz_heat / 10, 6 * Nz_heat / 10 },
		Interval{ 6 * Nz_heat / 10, 7 * Nz_heat / 10 },
		Interval{ 7 * Nz_heat / 10, 8 * Nz_heat / 10 },
		Interval{ 8 * Nz_heat / 10, 9 * Nz_heat / 10 },
		Interval{ 9 * Nz_heat / 10,  Nz_heat - 1 }
	};

	vector<Interval> new_interv_acoustic_z = {
		Interval{ 1, Nz_acoustic / 10 },
		Interval{ Nz_acoustic / 10, 2 * Nz_acoustic / 10 },
		Interval{ 2 * Nz_acoustic / 10, 3 * Nz_acoustic / 10 },
		Interval{ 3 * Nz_acoustic / 10, 4 * Nz_acoustic / 10 },
		Interval{ 4 * Nz_acoustic / 10, 5 * Nz_acoustic / 10 },
		Interval{ 5 * Nz_acoustic / 10, 6 * Nz_acoustic / 10 },
		Interval{ 6 * Nz_acoustic / 10, 7 * Nz_acoustic / 10 },
		Interval{ 7 * Nz_acoustic / 10, 8 * Nz_acoustic / 10 },
		Interval{ 8 * Nz_acoustic / 10, 9 * Nz_acoustic / 10 },
		Interval{ 9 * Nz_acoustic / 10,  Nz_acoustic - 1 }
	};
	
	//Calculation00(F, tbeam, tmpe0, tmpe1, tmpi1, Nx, Ny, Nz, dx, dy, dz, dt, A1, A2, B1, CC1, CC2, n, beta);
	Calculation00(mt, param, spl_C_e_on_T, spl_G_e_on_T, /*double*** F,*/  tbeam, tmpe0, tmpe1, tmpi0, tmpi1, Nx_heat, Ny_heat, Nz_heat, dx_heat, dy_heat, dz_heat, dt, A2, n, beta);
	Calculation1(a1x, a2x, b1x, b2x, a1z, b1y, b2y, b1z, b2z, a1y, e1, dx_acoustic, dy_acoustic, dz_acoustic, dt, CC0, Nx_acoustic, Ny_acoustic, Nz_acoustic, new_interv_acoustic_z);
	GRU1(e2, Nx_acoustic, Ny_acoustic);
	GRU2(b2x, b2y, b2z, Nx_acoustic, Ny_acoustic, Nz_acoustic);
	Calculation2(a1x, a1y, a1z, a2x, a2y, a2z, b2x, b2y, b2z, V2, e2, tmpe2, tmpi2, dx_acoustic, dy_acoustic, dz_acoustic, dt, CC0, Nx_acoustic, Ny_acoustic, Nz_acoustic,  V0, Melt_metal, mt, param, spl_C_e_on_T, spl_Te, spl_Ti, Te_acoustic, Ti_acoustic, new_interv_acoustic_z, spl_C_l_on_T);

	tt = n * dt * t0 * 1e+15; // vremja v fs

	int number_plots = 24; //24
	int count_of_lines = 9;
	int npxzt = 1;
	int current_number_line_melting = 1;
	GnuPlot plt(number_plots); // объект хранит 5 плотиков
	system("pause");
	plt.SetParametrs2D(0, 3, 3, "Te,Ti", "x,mkm", "Te,Ti,P");
	plt.SetParametrs2D(1, 2, 3, "Te,Ti", "time,fs", "Te,Ti (K)");
	plt.SetParametrsOnPlotColor(2, "Te", "z,mkm", "x,mkm", (1e+4* dz_heat* Nz_heat / kabs), (1e+4* r0* dx_heat* Nx_heat));
	plt.SetParametrsOnPlotColor(3, "Ti", "z,mkm", "x,mkm", (1e+4* dz_heat* Nz_heat / kabs), (1e+4* r0* dx_heat* Nx_heat));
	plt.SetParametrsOnPlotColor(4, "P(x,N/2,z)", "z,mkm", "x,mkm", (1e+4* dz_acoustic* Nz_acoustic / kabs), (1e+4* r0* dx_acoustic* Nx_acoustic));
	plt.SetParametrsOnPlotColor(5, "P(x,y,10)", "y,mkm", "x,mkm", (1e+4* r0* dy_acoustic* Ny_acoustic), (1e+4* r0* dx_acoustic* Nx_acoustic));

	//plt.SetParametrsOnPlotColor(29, "Pac(x,N/2,z)", "z,mkm", "x,mkm", (1e+4 * dz_acoustic * Nz_acoustic / kabs), (1e+4 * r0 * dx_acoustic * Nx_acoustic)); // подобие на 4 и 5
	//plt.SetParametrsOnPlotColor(30, "Pac(x,y,10)", "y,mkm", "x,mkm", (1e+4 * r0 * dy_acoustic * Ny_acoustic), (1e+4 * r0 * dx_acoustic * Nx_acoustic));
	//plt.SetParametrsOnPlotColor(31, "Pti(x,N/2,z)", "z,mkm", "x,mkm", (1e+4 * dz_acoustic * Nz_acoustic / kabs), (1e+4 * r0 * dx_acoustic * Nx_acoustic));
	//plt.SetParametrsOnPlotColor(32, "Pti(x,y,10)", "y,mkm", "x,mkm", (1e+4 * r0 * dy_acoustic * Ny_acoustic), (1e+4 * r0 * dx_acoustic * Nx_acoustic));
	//plt.SetParametrsOnPlotColor(33, "Pte(x,N/2,z)", "z,mkm", "x,mkm", (1e+4 * dz_acoustic * Nz_acoustic / kabs), (1e+4 * r0 * dx_acoustic * Nx_acoustic));
	//plt.SetParametrsOnPlotColor(34, "Pte(x,y,10)", "y,mkm", "x,mkm", (1e+4 * r0 * dy_acoustic * Ny_acoustic), (1e+4 * r0 * dx_acoustic * Nx_acoustic)); // подобие на 4 и 5

	plt.SetParametrsOnPlotColor(6, "Te", "y,mkm", "x,mkm", (1e+4* r0* dy_heat* Ny_heat), (1e+4* r0* dx_heat* Nx_heat));
	plt.SetParametrsOnPlotColor(7, "Ti", "y,mkm", "x,mkm", (1e+4* r0* dy_heat* Ny_heat), (1e+4* r0* dx_heat* Nx_heat));
	plt.SetParametrs2D(8, count_of_lines, 3, "P,bar", "x,mkm", "P,bar");
	plt.SetParametrs2D(9, count_of_lines, 3, "P, bar", "z,mkm", "P, bar");
	plt.SetParametrs2D(10, 1, 3, "Temperature in central zone", "time,fs", "Ti, K"); // копии 1 и 10 графиков
	plt.SetParametrs2D(11, 11, 3, "Temperature T(z)", "z,mkm", "Ti, K");
	plt.SetParametrs2D(12, 11, 3, "Temperature T(x)", "x,mkm", "Ti, K");

	plt.SetParametrsOnPlotColor(13, "Melting zone (xz)", "z,mkm", "x,mkm", (1e+4* dz_heat* Nz_heat / kabs), (1e+4* r0* dx_heat* Nx_heat));
	plt.SetParametrsOnPlotColor(14, "Melting zone (xy)", "y,mkm", "x,mkm", (1e+4* r0* dy_heat* Ny_heat), (1e+4* r0* dx_heat* Nx_heat));
	plt.SetParametrs2D(15, 2, 3, "Te,Ti", "x,mkm", "Te,Ti (K)");
	plt.SetParametrs2D(16, 2, 3, "Te,Ti", "z,mkm", "Te,Ti (K)");

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	plt.SetParametrs2D(17, 1, 3, "P", "time,fs", "P (bar)");
	plt.SetParametrs2D(18, 1, 3, "P", "x,mkm", "P (bar)");
	plt.SetParametrs2D(19, 1, 3, "P", "z,mkm", "P (bar)");

	/*plt.SetParametrs2D(20, 1, 3, "Pacoustic", "time,fs", "P (bar)");
	plt.SetParametrs2D(21, 1, 3, "Pacoustic", "x,mkm", "P (bar)");
	plt.SetParametrs2D(22, 1, 3, "Pacoustic", "z,mkm", "P (bar)");

	plt.SetParametrs2D(23, 1, 3, "PTi", "time,fs", "P (bar)");
	plt.SetParametrs2D(24, 1, 3, "PTi", "x,mkm", "P (bar)");
	plt.SetParametrs2D(25, 1, 3, "PTi", "z,mkm", "P (bar)");

	plt.SetParametrs2D(26, 1, 3, "PTe", "time,fs", "P (bar)");
	plt.SetParametrs2D(27, 1, 3, "PTe", "x,mkm", "P (bar)");
	plt.SetParametrs2D(28, 1, 3, "PTe", "z,mkm", "P (bar)");*/

	//plt.SetParametrs2D(20, 1, 3, "Temperature in central zone", "time,fs", "Ti, K"); // копии 1 и 10 графиков
	//plt.SetParametrs2D(21, 2, 3, "Te,Ti", "time,fs", "Te,Ti (K)"); // копии 1 и 10 графиков
	//plt.SetParametrs2D(22, 1, 3, "Temperature in central zone", "time,fs", "Ti, K"); // копии 1 и 10 графиков
	//plt.SetParametrs2D(23, 2, 3, "Te,Ti", "time,fs", "Te,Ti (K)"); // копии 1 и 10 графиков
	//plt.SetParametrs2D(24, 1, 3, "Temperature in central zone", "time,fs", "Ti, K"); // копии 1 и 10 графиков
	//plt.SetParametrs2D(25, 2, 3, "Te,Ti", "time,fs", "Te,Ti (K)"); // копии 1 и 10 графиков
	//plt.SetParametrs2D(26, 1, 3, "Temperature in central zone", "time,fs", "Ti, K"); // копии 1 и 10 графиков


	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// возможо здесь нужно отобразить в начальный момент времени, но это не обзательно

	plt.SetGridOnPlot3D(8, Nx_acoustic, Ny_acoustic, Nz_acoustic, 1e+4* dx_acoustic* r0, 1e+4* dy_acoustic* r0, 1e+4* dz_acoustic / kabs, count_of_lines, fun_on_x);
	plt.SetGridOnPlot3D(9, Nx_acoustic, Ny_acoustic, Nz_acoustic, 1e+4* dx_acoustic* r0, 1e+4* dy_acoustic* r0, 1e+4* dz_acoustic / kabs, count_of_lines, fun_on_z);
	plt.SetGridOnPlot3D(11, Nx_heat, Ny_heat, Nz_heat, 1e+4* dx_heat* r0, 1e+4* dy_heat* r0, 1e+4* dz_heat / kabs, 11, fun_on_z);
	plt.SetGridOnPlot3D(12, Nx_heat, Ny_heat, Nz_heat, 1e+4* dx_heat* r0, 1e+4* dy_heat* r0, 1e+4* dz_heat / kabs, 11, fun_on_x);

	clock_t t;
	double sum_t = 0;
	t = clock();

	vector<double***> TeTi = { tmpe2 , tmpi2 };//????
	vector<double***> Ti = { tmpi2 };//????
	vector<double***> P = { e1 };/////
	//vector<double***> Pac = { Pacoustic };//?????
	//vector<double***> PTTi = { PTi };//????
	//vector<double***> PTTe = { PTe }; //????
	vector<double***> empty;
	vector<string> TeTiTfs = { "Te","Ti" };
	vector<string> Pfs = { "P" };
	vector<string> TiTfs = { "Ti" };
	vector<string> Melting_zone = { "Melting_zone" };
	vector<string> Legenda_melting_tmp;
	vector<string> Legenda;
	vector<string> FilesForWritting;  // имена файлов куда запишутся поля темепратур и т.д
	vector<string> FilesForReading;
	vector<int> for_close_and_open = { 2, 3, 4, 5, 6, 7, 13, 14, 15, 16, 18, 19, 21/*, 22, 24, 25, 27, 28, 29, 30, 31, 32, 33, 34 */ };
	vector<int> for_close_and_open_ac = { 0, 1, 2 };
	string namefile;
	string fileoflistnamefiles = "List.txt";
	string depth = "Depth.txt";
	vector<Point3D> points_rupture;
	vector<Point3D> points_rupture_for_plot_z;
	vector<Point3D> points_rupture_for_plot_x;
	vector<Point3D>::iterator it;
	vector<Point3D>::iterator it_z, it_x;
	points_rupture.clear();
	points_rupture.erase(points_rupture.begin(), points_rupture.end());

	string depth2 = current_namefile;// current_namefile;//"Depth2.txt"; // I0 = 0.25e+10
	//string depth23 = "Depth23.txt"; // I0 = 0.275e+10
	//string depth234 = "Depth234.txt"; // I0 = 0.3e+10
	//string depth2345 = "Depth2345.txt"; // I0 = 0.225e+10
	//string depth23456 = "Depth23456.txt"; // I0 = 0.325e+10
	//string depth234567 = "Depth234567.txt"; // I0 = 0.35e+10
	//ofstream file;
	//ifstream fin(fileoflistnamefiles);
	//ofstream fout_depth(depth);
	ofstream fout_depth2(depth2);

	vector<string> Pacteti = { "Pac.txt", "PTi.txt", "PTe.txt", "Pacx.txt", "PTix.txt", "PTex.txt" };
	GnuPlot pac(6, Pacteti);
	pac.SetParametrs2D(0, 1, 3, "Pacoustic", "z,mkm", "P (bar)");
	pac.SetParametrs2D(1, 1, 3, "PTi", "z,mkm", "P (bar)");
	pac.SetParametrs2D(2, 1, 3, "PTe", "z,mkm", "P (bar)");
	pac.SetParametrs2D(3, 1, 3, "Pacoustic", "x,mkm", "P (bar)");
	pac.SetParametrs2D(4, 1, 3, "PTi", "x,mkm", "P (bar)");
	pac.SetParametrs2D(5, 1, 3, "PTe", "x,mkm", "P (bar)");
	ofstream fout_Pac(Pacteti[0]);
	ofstream fout_PTi(Pacteti[1]);
	ofstream fout_PTe(Pacteti[2]);
	ofstream fout_Pacx(Pacteti[3]);
	ofstream fout_PTix(Pacteti[4]);
	ofstream fout_PTex(Pacteti[5]);
	ofstream fout_point("Points.txt");

	int count_depth = 1;
	bool change_time_step = true;
	cout << endl << endl;

	do // цикл по времени
	{
		//Calculation0;

		t = clock();

		Calculation0(mt, param, spl_C_e_on_T, spl_G_e_on_T, /*double*** F,*/  tbeam, tmpe0, tmpe1, tmpe2, tmpi0, tmpi1, tmpi2, A2, dx_heat, dy_heat, dz_heat, dt, Nx_heat, Ny_heat, Nz_heat, n, beta, new_interv, Melt_metal, melting, spl_C_l_on_T);
		Calculation1(a1x, a2x, b1x, b2x, a1z, b1y, b2y, b1z, b2z, a1y, e1, dx_acoustic, dy_acoustic, dz_acoustic, dt, CC0, Nx_acoustic, Ny_acoustic, Nz_acoustic, new_interv_acoustic_z);
		GRU1(e2, Nx_acoustic, Ny_acoustic);
		GRU2(b2x, b2y, b2z, Nx_acoustic, Ny_acoustic, Nz_acoustic);
		Calculation2(a1x, a1y, a1z, a2x, a2y, a2z, b2x, b2y, b2z, V2, e2, tmpe2, tmpi2, dx_acoustic, dy_acoustic, dz_acoustic, dt, CC0, Nx_acoustic, Ny_acoustic, Nz_acoustic, V0, Melt_metal, mt, param, spl_C_e_on_T, spl_Te, spl_Ti, Te_acoustic, Ti_acoustic, new_interv_acoustic_z, spl_C_l_on_T);

		cout << "It took me %d clicks (%f seconds)." << endl <<
			(std::clock() - t) / (double)CLOCKS_PER_SEC << endl;
		sum_t += (std::clock() - t) / (double)CLOCKS_PER_SEC;
		cout << " General duration calculation = " << sum_t << endl;
		cout << endl;

		double level_depth = 0.;

		tt = n * dt * t0 * 1e+15; // vremja v fs
		cout << " time, fs = " << tt << endl;


		if ((tmpi1[Nx_heat / 2][Ny_heat / 2][1] * param.T00) >= (Melt_metal[mt].T_melting + 1.)) // 2-й вариант определения макс глубины расплава
		{
			for (int i = 0; i < Nz_heat; i++)
			{
				//Depth[i][count_depth] = (tmpi1[Nx / 2][Ny / 2][i] * param.T00);
				if ((tmpi1[Nx_heat / 2][Ny_heat / 2][i] * param.T00) >= (Melt_metal[mt].T_melting + 1.))
				{
					level_depth = (1e+4 * dz_heat * i / kabs);
				}
			}
			count_depth++;

			fout_depth2 << level_depth << endl;
		}

		plt.SetDataOnPlot3D(1, Nx_heat, Ny_heat, Nz_heat, 1e+4 * dx_heat * r0, 1e+4 * dy_heat * r0, 1e+4 * dz_heat / kabs, NULL, T00, Nx_heat / 2, Ny_heat / 2, 1, 2, 0, tt, TeTi, fun_on_t);
		plt.SetDataOnPlot3D(10, Nx_heat, Ny_heat, Nz_heat, 1e+4 * dx_heat * r0, 1e+4 * dy_heat * r0, 1e+4 * dz_heat / kabs, NULL, T00, Nx_heat / 2, Ny_heat / 2, 1, 1, 0, tt, Ti, fun_on_t);
		/*plt.SetDataOnPlot3D(17, Nx, Ny, Nz, 1e+4* dx* r0, 1e+4* dy* r0, 1e+4* dz / kabs, NULL, p0v, Nx / 2, Ny / 2, 1, 1, 0, tt, P, fun_on_t);
		plt.SetDataOnPlot3D(20, Nx, Ny, Nz, 1e+4* dx* r0, 1e+4* dy* r0, 1e+4* dz / kabs, NULL, 1., Nx / 2, Ny / 2, 1, 1, 0, tt, Pac, fun_on_t);
		plt.SetDataOnPlot3D(23, Nx, Ny, Nz, 1e+4* dx* r0, 1e+4* dy* r0, 1e+4* dz / kabs, NULL, p0v, Nx / 2, Ny / 2, 1, 1, 0, tt, PTTi, fun_on_t);
		plt.SetDataOnPlot3D(26, Nx, Ny, Nz, 1e+4* dx* r0, 1e+4* dy* r0, 1e+4* dz / kabs, NULL, p0v, Nx / 2, Ny / 2, 1, 1, 0, tt, PTTe, fun_on_t);*/
		plt.SetDataOnPlot3D(17, Nx_acoustic, Ny_acoustic, Nz_acoustic, 1e+4 * dx_acoustic * r0, 1e+4 * dy_acoustic * r0, 1e+4 * dz_acoustic / kabs, NULL, p0v, Nx_acoustic / 2, Ny_acoustic / 2, 0, 1, 0, tt, P, fun_on_t);
		/*plt.SetDataOnPlot3D(20, Nx_acoustic, Ny_acoustic, Nz_acoustic, 1e+4 * dx_acoustic * r0, 1e+4 * dy_acoustic * r0, 1e+4 * dz_acoustic / kabs, NULL, p0v, Nx_acoustic / 2, Ny_acoustic / 2, 0, 1, 0, tt, Pac, fun_on_t);
		plt.SetDataOnPlot3D(23, Nx_acoustic, Ny_acoustic, Nz_acoustic, 1e+4 * dx_acoustic * r0, 1e+4 * dy_acoustic * r0, 1e+4 * dz_acoustic / kabs, NULL, p0v, Nx_acoustic / 2, Ny_acoustic / 2, 0, 1, 0, tt, PTTi, fun_on_t);
		plt.SetDataOnPlot3D(26, Nx_acoustic, Ny_acoustic, Nz_acoustic, 1e+4 * dx_acoustic * r0, 1e+4 * dy_acoustic * r0, 1e+4 * dz_acoustic / kabs, NULL, p0v, Nx_acoustic / 2, Ny_acoustic / 2, 0, 1, 0, tt, PTTe, fun_on_t);*/

		/*if (melting)
		if ((tmpi1[Nx / 2][Ny / 2][1] * param.T00) >= (Melt_metal[mt].T_melting + 1.)) // Зона плавления
		{
			if ((int(ceil(tt)) % 100 == 0) || (int(trunc(tt)) % 100 == 0)) //
			{
				for (int i = 0; i < Nx; i++)
				{
					for (int k = 0; k < Nz; k++)
					{
						//if ((tmpi1[i][Ny / 2][k] * param.T00) >= (Melt_metal[mt].T_melting - 1.) && (tmpi1[i][Ny / 2][k] * param.T00) <= (Melt_metal[mt].T_melting + 1.))
						if ((tmpi1[i][Ny / 2][k] * param.T00) >= (Melt_metal[mt].T_melting + 1.))
						{
							massive_melting_tmp1[i][Ny / 2][k] = tmpi1[i][Ny / 2][k];
						}
					}
				}

				for (int i = 0; i < Nx; i++)
				{
					for (int j = 0; j < Ny; j++)
					{
						if ((tmpi1[i][j][1] * param.T00) >= (Melt_metal[mt].T_melting + 1.))
						{
							massive_melting_tmp2[i][j][1] = tmpi1[i][j][1];
						}
					}
				}

				plt.SetDataOnPlotColor3D(13, Nx, Ny, Nz, 1e+4 * dx * r0, 1e+4 * dy * r0, 1e+4 * dz / kabs, massive_melting_tmp1, T00, Ny / 2, xz);
				plt.SetDataOnPlotColor3D(14, Nx, Ny, Nz, 1e+4 * dx * r0, 1e+4 * dy * r0, 1e+4 * dz / kabs, massive_melting_tmp2, T00, 1, xy);

				namefile.clear();
				namefile = "Melting zone (zx) " + metall + "," + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs";
				plt.ShowDataOnPlotColor(13, namefile, true);
				Sleep(1000);

				namefile.clear();
				namefile = "Melting zone (yx) " + metall + "," + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs";
				plt.ShowDataOnPlotColor(14, namefile, true);
				Sleep(1000);

				massive_melting_tmp1 = massive_melting_null;
				massive_melting_tmp2 = massive_melting_null;

				plt.Close_and_open_files_for_replot(for_close_and_open);
			}
		}*/

		if ((tmpi1[Nx_heat / 2][Ny_heat / 2][1] * param.T00) >= (Melt_metal[mt].T_melting + 1.)) // Профили температур ионной решетки
		{
			if ((int(ceil(tt)) % 1000 == 0 || int(trunc(tt)) % 1000 == 0) && current_number_line_melting <= 10)
			{
				plt.SetDataOnPlot3D(12, Nx_heat, Ny_heat, Nz_heat, 1e+4 * dx_heat * param.r0, 1e+4 * dy_heat * param.r0, 1e+4 * dz_heat / param.kabs, tmpi1, param.T00, NULL, Ny_heat / 2, 1, 12, current_number_line_melting, NULL, empty, fun_on_x);
				plt.SetDataOnPlot3D(11, Nx_heat, Ny_heat, Nz_heat, 1e+4 * dx_heat * param.r0, 1e+4 * dy_heat * param.r0, 1e+4 * dz_heat / param.kabs, tmpi1, param.T00, Nx_heat / 2, Ny_heat / 2, NULL, 12, current_number_line_melting, NULL, empty, fun_on_z);
				current_number_line_melting++;
				Legenda_melting_tmp.push_back(ConvertNumToStringdouble(tt) + " fs");
			}
		}


		//!!!!!!!!!!
		if ((int)tt == 500 || (int)tt == 1000 || (int)tt == 2000 || (int)tt == 3000 || (int)tt == 4000 || (int)tt == 5000 || (int)tt == 6000 || (int)tt == 7000 || (int)tt == 8000 || (int)tt == 9000
			|| (int)tt == 10000 || (int)tt == 11000 || (int)tt == 12000 || (int)tt == 13000 || (int)tt == 14000 || (int)tt == 15000 || (int)tt == 16000 || (int)tt == 17000 || (int)tt == 18000 || (int)tt == 19000 || (int)tt == 20000)
		//if (int(ceil(tt)) == 500 || int(trunc(tt)) == 500 || int(ceil(tt)) % 1000 == 0 || int(trunc(tt)) % 1000 == 0)
			//|| (int(ceil(tt)) % 50 == 0) || (int(trunc(tt)) % 50 == 0)|| (int(ceil(tt)) == 70) || (int(trunc(tt)) == 70))
		{//!!!!!!!!!!!!!!!!!!!!!!!!!!	tp * 1e+15
			//plt.Close_and_open_files_for_replot(for_close_and_open);
			namefile.clear();
			namefile = "Te,Ti (t,fms) " + metall + "," + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs, 1";
			plt.ShowDataOnPlot2D(1, 2, TeTiTfs, namefile, true);
			namefile.clear();
			namefile = "Ti (t,fms) " + metall + "," + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs, 1";
			plt.ShowDataOnPlot2D(10, 1, TiTfs, namefile, true);
			namefile.clear();
			namefile = "P (t,fms) " + metall + "," + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs, 1";
			plt.ShowDataOnPlot2D(17, 1, Pfs, namefile, true);
			/*namefile.clear();
			namefile = "Pac (t,fms) " + metall + "," + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs, 1";
			plt.ShowDataOnPlot2D(20, 1, Pfs, namefile, true);
			namefile.clear();
			namefile = "PTi (t,fms) " + metall + "," + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs, 1";
			plt.ShowDataOnPlot2D(23, 1, Pfs, namefile, true);
			namefile.clear();
			namefile = "PTe (t,fms) " + metall + "," + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs, 1";
			plt.ShowDataOnPlot2D(26, 1, Pfs, namefile, true);*/

			plt.SetDataOnPlotColor3D(2, Nx_heat, Ny_heat, Nz_heat, 1e+4 * dx_heat * r0, 1e+4 * dy_heat * r0, 1e+4 * dz_heat / kabs, tmpe2, T00, Ny_heat / 2, xz);
			plt.SetDataOnPlotColor3D(3, Nx_heat, Ny_heat, Nz_heat, 1e+4 * dx_heat * r0, 1e+4 * dy_heat * r0, 1e+4 * dz_heat / kabs, tmpi2, T00, Ny_heat / 2, xz);
			plt.SetDataOnPlotColor3D(4, Nx_acoustic, Ny_acoustic, Nz_acoustic, 1e+4 * dx_acoustic * r0, 1e+4 * dy_acoustic * r0, 1e+4 * dz_acoustic / kabs, e1, p0v, Ny_acoustic / 2, xz);
			plt.SetDataOnPlotColor3D(5, Nx_acoustic, Ny_acoustic, Nz_acoustic, 1e+4 * dx_acoustic * r0, 1e+4 * dy_acoustic * r0, 1e+4 * dz_acoustic / kabs, e1, p0v, 10, xy);

			/*plt.SetDataOnPlotColor3D(29, Nx_acoustic, Ny_acoustic, Nz_acoustic, 1e+4 * dx_acoustic * r0, 1e+4 * dy_acoustic * r0, 1e+4 * dz_acoustic / kabs, Pacoustic, p0v, Ny_acoustic / 2, xz);
			plt.SetDataOnPlotColor3D(30, Nx_acoustic, Ny_acoustic, Nz_acoustic, 1e+4 * dx_acoustic * r0, 1e+4 * dy_acoustic * r0, 1e+4 * dz_acoustic / kabs, Pacoustic, p0v, 10, xy);
			plt.SetDataOnPlotColor3D(31, Nx_acoustic, Ny_acoustic, Nz_acoustic, 1e+4 * dx_acoustic * r0, 1e+4 * dy_acoustic * r0, 1e+4 * dz_acoustic / kabs, PTi, p0v, Ny_acoustic / 2, xz);
			plt.SetDataOnPlotColor3D(32, Nx_acoustic, Ny_acoustic, Nz_acoustic, 1e+4 * dx_acoustic * r0, 1e+4 * dy_acoustic * r0, 1e+4 * dz_acoustic / kabs, PTi, p0v, 10, xy);
			plt.SetDataOnPlotColor3D(33, Nx_acoustic, Ny_acoustic, Nz_acoustic, 1e+4 * dx_acoustic * r0, 1e+4 * dy_acoustic * r0, 1e+4 * dz_acoustic / kabs, PTe, p0v, Ny_acoustic / 2, xz);
			plt.SetDataOnPlotColor3D(34, Nx_acoustic, Ny_acoustic, Nz_acoustic, 1e+4 * dx_acoustic * r0, 1e+4 * dy_acoustic * r0, 1e+4 * dz_acoustic / kabs, PTe, p0v, 10, xy);
*/
			plt.SetDataOnPlotColor3D(6, Nx_heat, Ny_heat, Nz_heat, 1e+4 * dx_heat * r0, 1e+4 * dy_heat * r0, 1e+4 * dz_heat / kabs, tmpe2, T00, 1, xy);
			plt.SetDataOnPlotColor3D(7, Nx_heat, Ny_heat, Nz_heat, 1e+4 * dx_heat * r0, 1e+4 * dy_heat * r0, 1e+4 * dz_heat / kabs, tmpi2, T00, 1, xy);
			if (npxzt <= count_of_lines)
			{
				plt.SetDataOnPlot3D(8, Nx_acoustic, Ny_acoustic, Nz_acoustic, 1e+4 * dx_acoustic * r0, 1e+4 * dy_acoustic * r0, 1e+4 * dz_acoustic / kabs, e1, p0v, NULL, Ny_acoustic / 2, 10, count_of_lines, npxzt, NULL, empty, fun_on_x);
				plt.SetDataOnPlot3D(9, Nx_acoustic, Ny_acoustic, Nz_acoustic, 1e+4 * dx_acoustic * r0, 1e+4 * dy_acoustic * r0, 1e+4 * dz_acoustic / kabs, e1, p0v, Nx_acoustic / 2, Ny_acoustic / 2, NULL, count_of_lines, npxzt, NULL, empty, fun_on_z);
				npxzt++;
				Legenda.push_back(ConvertNumToStringdouble(tt) + " fs");
			}
			namefile.clear();
			namefile = "file Te(x, N div 2, z) " + metall + "," + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs";
			plt.ShowDataOnPlotColor(2, namefile, true);
			Sleep(1000);

			namefile.clear();
			namefile = "file Ti(x, N div 2, z) " + metall + "," + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs";
			plt.ShowDataOnPlotColor(3, namefile, true);
			Sleep(1000);

			namefile.clear();
			namefile = "file P (x,N div 2,z) " + metall + "," + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs";
			plt.ShowDataOnPlotColor(4, namefile, true);
			Sleep(1000);

			namefile.clear();
			namefile = "file P (x,y,10) " + metall + "," + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs";
			plt.ShowDataOnPlotColor(5, namefile, true);
			Sleep(1000);


			/*namefile.clear();
			namefile = "file Pac (x,N div 2,z) " + metall + "," + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs";
			plt.ShowDataOnPlotColor(29, namefile, true);
			Sleep(1000);

			namefile.clear();
			namefile = "file Pac (x,y,10) " + metall + "," + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs";
			plt.ShowDataOnPlotColor(30, namefile, true);
			Sleep(1000);

			namefile.clear();
			namefile = "file Pti (x,N div 2,z) " + metall + "," + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs";
			plt.ShowDataOnPlotColor(31, namefile, true);
			Sleep(1000);

			namefile.clear();
			namefile = "file Pti (x,y,10) " + metall + "," + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs";
			plt.ShowDataOnPlotColor(32, namefile, true);
			Sleep(1000);

			namefile.clear();
			namefile = "file Pte (x,N div 2,z) " + metall + "," + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs";
			plt.ShowDataOnPlotColor(33, namefile, true);
			Sleep(1000);

			namefile.clear();
			namefile = "file Pte (x,y,10) " + metall + "," + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs";
			plt.ShowDataOnPlotColor(34, namefile, true);
			Sleep(1000);*/

			namefile.clear();
			namefile = "file Te (x,y,1) " + metall + "," + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs";
			plt.ShowDataOnPlotColor(6, namefile, true);
			Sleep(1000);

			namefile.clear();
			namefile = "file Ti (x,Ny,1) " + metall + "," + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs";
			plt.ShowDataOnPlotColor(7, namefile, true);
			Sleep(1000);

			namefile.clear();

			plt.SetGridOnPlot3D(15, Nx_heat, Ny_heat, Nz_heat, 1e+4 * dx_heat * r0, 1e+4 * dy_heat * r0, 1e+4 * dz_heat / kabs, 2, fun_on_x);
			plt.SetGridOnPlot3D(16, Nx_heat, Ny_heat, Nz_heat, 1e+4 * dx_heat * r0, 1e+4 * dy_heat * r0, 1e+4 * dz_heat / kabs, 2, fun_on_z);

			plt.SetGridOnPlot3D(18, Nx_acoustic, Ny_acoustic, Nz_acoustic, 1e+4 * dx_acoustic * r0, 1e+4 * dy_acoustic * r0, 1e+4 * dz_acoustic / kabs, 1, fun_on_x);
			plt.SetGridOnPlot3D(19, Nx_acoustic, Ny_acoustic, Nz_acoustic, 1e+4 * dx_acoustic * r0, 1e+4 * dy_acoustic * r0, 1e+4 * dz_acoustic / kabs, 1, fun_on_z);
			/*plt.SetGridOnPlot3D(21, Nx_acoustic, Ny_acoustic, Nz_acoustic, 1e+4 * dx_acoustic * r0, 1e+4 * dy_acoustic * r0, 1e+4 * dz_acoustic / kabs, 1, fun_on_x);
			plt.SetGridOnPlot3D(22, Nx_acoustic, Ny_acoustic, Nz_acoustic, 1e+4 * dx_acoustic * r0, 1e+4 * dy_acoustic * r0, 1e+4 * dz_acoustic / kabs, 1, fun_on_z);
			plt.SetGridOnPlot3D(24, Nx_acoustic, Ny_acoustic, Nz_acoustic, 1e+4 * dx_acoustic * r0, 1e+4 * dy_acoustic * r0, 1e+4 * dz_acoustic / kabs, 1, fun_on_x);
			plt.SetGridOnPlot3D(25, Nx_acoustic, Ny_acoustic, Nz_acoustic, 1e+4 * dx_acoustic * r0, 1e+4 * dy_acoustic * r0, 1e+4 * dz_acoustic / kabs, 1, fun_on_z);
			plt.SetGridOnPlot3D(27, Nx_acoustic, Ny_acoustic, Nz_acoustic, 1e+4 * dx_acoustic * r0, 1e+4 * dy_acoustic * r0, 1e+4 * dz_acoustic / kabs, 1, fun_on_x);
			plt.SetGridOnPlot3D(28, Nx_acoustic, Ny_acoustic, Nz_acoustic, 1e+4 * dx_acoustic * r0, 1e+4 * dy_acoustic * r0, 1e+4 * dz_acoustic / kabs, 1, fun_on_z);*/

			plt.SetDataOnPlot3D(15, Nx_heat, Ny_heat, Nz_heat, 1e+4 * dx_heat * r0, 1e+4 * dy_heat * r0, 1e+4 * dz_heat / kabs, tmpe2, T00, NULL, Ny_heat / 2, 1, 2, 1, NULL, empty, fun_on_x);
			plt.SetDataOnPlot3D(15, Nx_heat, Ny_heat, Nz_heat, 1e+4 * dx_heat * r0, 1e+4 * dy_heat * r0, 1e+4 * dz_heat / kabs, tmpi2, T00, NULL, Ny_heat / 2, 1, 2, 2, NULL, empty, fun_on_x);
			plt.SetDataOnPlot3D(16, Nx_heat, Ny_heat, Nz_heat, 1e+4 * dx_heat * r0, 1e+4 * dy_heat * r0, 1e+4 * dz_heat / kabs, tmpe2, T00, Nx_heat / 2, Ny_heat / 2, NULL, 2, 1, NULL, empty, fun_on_z);
			plt.SetDataOnPlot3D(16, Nx_heat, Ny_heat, Nz_heat, 1e+4 * dx_heat * r0, 1e+4 * dy_heat * r0, 1e+4 * dz_heat / kabs, tmpi2, T00, Nx_heat / 2, Ny_heat / 2, NULL, 2, 2, NULL, empty, fun_on_z);

			plt.SetDataOnPlot3D(18, Nx_acoustic, Ny_acoustic, Nz_acoustic, 1e+4 * dx_acoustic * r0, 1e+4 * dy_acoustic * r0, 1e+4 * dz_acoustic / kabs, e1, p0v, NULL, Ny_acoustic / 2, 1, 1, 1, NULL, empty, fun_on_x);
			plt.SetDataOnPlot3D(19, Nx_acoustic, Ny_acoustic, Nz_acoustic, 1e+4 * dx_acoustic * r0, 1e+4 * dy_acoustic * r0, 1e+4 * dz_acoustic / kabs, e1, p0v, Nx_acoustic / 2, Ny_acoustic / 2, NULL, 1, 1, NULL, empty, fun_on_z);

			/*fout_Pac.open("Pac.txt");
			fout_PTi.open("PTi.txt");
			fout_PTe.open("PTe.txt");*/
			/*spl_Te.Calculation_InterpolationFast(tmpe2, Nz_heat, dz_heat, Nx_acoustic / 2, Ny_acoustic / 2);
			spl_Ti.Calculation_InterpolationFast(tmpi2, Nz_heat, dz_heat, Nx_acoustic / 2, Ny_acoustic / 2);*/

			
			// ЗАПИСЬ ОТДЕЛЬНЫХ СЛАГАЕМЫХ УРАВНЕНИЯ СОСТОЯНИЯ МИ-ГРЮНАЙЗЕНА 
			spl_Te.InterpolateFast(1, tmpe2, Nx_acoustic / 2, Ny_acoustic / 2);
			spl_Ti.InterpolateFast(1, tmpi2, Nx_acoustic / 2, Ny_acoustic / 2);
			// может значения Т перезаписать в матрицу под сетку акустики?
			// если есть изменения в акустике не забывать менять там где идет расчет

			if (!points_rupture.empty())
			{
				for (it = points_rupture.begin(); it != points_rupture.end(); it++)// цикл по точкам где произошел разрыв
				{
					if ((*it).index_x == (Nx_acoustic / 2) && (*it).index_y == (Ny_acoustic / 2))
					{
						points_rupture_for_plot_z.push_back((*it));
					}

					if ((*it).index_y == (Ny_acoustic / 2) && (*it).index_z == 1)
					{
						points_rupture_for_plot_x.push_back((*it));
					}
				}
				sort(points_rupture_for_plot_z.begin(), points_rupture_for_plot_z.end(), comp_z);
				it_z = points_rupture_for_plot_z.begin();
				it_x = points_rupture_for_plot_x.begin();

				cout << " Points z for plotting " << endl;
				for (it = points_rupture_for_plot_z.begin(); it != points_rupture_for_plot_z.end(); it++)// цикл по точкам где произошел разрыв
				{
					cout << (*it).index_x << "   " << (*it).index_y << "   " << (*it).index_z << endl;
				}

				cout << endl << endl;
				cout << " Points x for plotting " << endl;

				for (it = points_rupture_for_plot_x.begin(); it != points_rupture_for_plot_x.end(); it++)// цикл по точкам где произошел разрыв
				{
					cout << (*it).index_x << "   " << (*it).index_y << "   " << (*it).index_z << endl;
				}

				cout << endl;
			}

			for (int k = 0; k < Nz_acoustic; k++)// цикл по z
			{
				//if (points_rupture.empty())
				{
					fout_Pac << 1e+4 * dz_acoustic * k / kabs << "   " << p0v * (1. - V2[Nx_acoustic / 2][Ny_acoustic / 2][k]) << endl;

					if ((spl_Ti.GetY(k * dz_acoustic) * param.T00) < (Melt_metal[mt].T_melting - 1.))
					{
						/*if ((spl_Ti.GetY(k * dz_acoustic) * param.T00) >= 1336.)
						{
							fout_PTi << 1e+4 * dz_acoustic * k / kabs << "   " << p0v * ((Melt_metal[mt].gi * param.T00 * (spl_C_l_on_T.GetY(param.T00 * spl_Ti.GetY(k * dz_acoustic))) / (Melt_metal[mt].Density * Melt_metal[mt].Density * pow(Melt_metal[mt].u0, 2) * V0 * 1e-6)) *
								(spl_Ti.GetY(k * dz_acoustic) - 1.) / V2[Nx_acoustic / 2][Ny_acoustic / 2][k]) << endl;
						}
						else*/
						//{
							fout_PTi << 1e+4 * dz_acoustic * k / kabs << "   " << p0v * ((Melt_metal[mt].gi * param.T00 * (Dependence_C_l_on_T(mt, param.T00 * spl_Ti.GetY(k * dz_acoustic)) / (100. * 100. * 100.)) / (Melt_metal[mt].Density * Melt_metal[mt].Density * pow(Melt_metal[mt].u0, 2) * V0 * 1e-6)) *
								(spl_Ti.GetY(k * dz_acoustic) - 1.) / V2[Nx_acoustic / 2][Ny_acoustic / 2][k]) << endl;
						//}

						fout_PTe << 1e+4 * dz_acoustic * k / kabs << "   " << p0v * ((Melt_metal[mt].ge * param.T00 * (spl_C_e_on_T.GetY(param.T00 * spl_Te.GetY(k * dz_acoustic)) / (100. * 100. * 100.)) / (Melt_metal[mt].Density * pow(Melt_metal[mt].u0, 2) * 1e-6)) *
							(spl_Te.GetY(k * dz_acoustic) - 1.)) << endl;
					}

					if ((spl_Ti.GetY(k * dz_acoustic) * param.T00) >= (Melt_metal[mt].T_melting - 1.) && (spl_Ti.GetY(k * dz_acoustic) * param.T00) <= (Melt_metal[mt].T_melting + 1.))
					{
						//if ((spl_Ti.GetY(k * dz_acoustic) * param.T00) >= 1336.)
						//{
						// новый способ плавления
							fout_PTi << 1e+4 * dz_acoustic * k / kabs << "   " << p0v * ((Melt_metal[mt].gi * param.T00 * ((( /*Dependence_C_l_on_T(mt, param.T00* spl_Ti.GetY(k* dz_acoustic)) / (100. * 100. * 100.) + (Melt_metal[mt].Q_fusion * Melt_metal[mt].Density) / (2 * 1.0)*/  spl_C_l_on_T.GetY(param.T00 * spl_Ti.GetY(k * dz_acoustic))))) / (pow((Melt_metal[mt].Density + Melt_metal[mt].DensityLiquid) / 2, 2) * pow((Melt_metal[mt].u0 + Melt_metal[mt].u0_Liquid) / 2, 2) * V0 * 1e-6)) *
								(spl_Ti.GetY(k * dz_acoustic) - 1.) / V2[Nx_acoustic / 2][Ny_acoustic / 2][k]) << endl;

						// старый способ плавления
							/*fout_PTi << 1e+4 * dz_acoustic * k / kabs << "   " << p0v * ((Melt_metal[mt].gi * param.T00 * (((Dependence_C_l_on_T(mt, param.T00 * spl_Ti.GetY(k * dz_acoustic)) / (100. * 100. * 100.) + (Melt_metal[mt].Q_fusion * Melt_metal[mt].Density) / (2 * 1.0)))) / (pow((Melt_metal[mt].Density + Melt_metal[mt].DensityLiquid) / 2, 2) * pow((Melt_metal[mt].u0 + Melt_metal[mt].u0_Liquid) / 2, 2) * V0 * 1e-6)) *
								(spl_Ti.GetY(k * dz_acoustic) - 1.) / V2[Nx_acoustic / 2][Ny_acoustic / 2][k]) << endl;*/
						//}

						fout_PTe << 1e+4 * dz_acoustic * k / kabs << "   " << p0v * ((Melt_metal[mt].ge * param.T00 * (spl_C_e_on_T.GetY(param.T00 * spl_Te.GetY(k * dz_acoustic)) / (100. * 100. * 100.)) / (((Melt_metal[mt].Density + Melt_metal[mt].DensityLiquid) / 2) * pow((Melt_metal[mt].u0 + Melt_metal[mt].u0_Liquid) / 2, 2) * 1e-6)) *
							(spl_Te.GetY(k * dz_acoustic) - 1.)) << endl;
					}

					// ДОБАВИТЬ ИФ СООТВ ЖИДКОЙ ФАЗЕ

					if ((spl_Ti.GetY(k * dz_acoustic) * param.T00) > (Melt_metal[mt].T_melting + 1.))
					{
						fout_PTi << 1e+4* dz_acoustic* k / kabs << "   " << p0v * ((Melt_metal[mt].gi* param.T00* (Dependence_C_l_on_T(mt, param.T00* spl_Ti.GetY(k* dz_acoustic)) / (100. * 100. * 100.)) / (Melt_metal[mt].DensityLiquid * Melt_metal[mt].DensityLiquid * pow(Melt_metal[mt].u0_Liquid, 2) * V0 * 1e-6))*
							(spl_Ti.GetY(k* dz_acoustic) - 1.) / V2[Nx_acoustic / 2][Ny_acoustic / 2][k]) << endl;
						//}

						fout_PTe << 1e+4 * dz_acoustic * k / kabs << "   " << p0v * ((Melt_metal[mt].ge * param.T00 * (spl_C_e_on_T.GetY(param.T00 * spl_Te.GetY(k * dz_acoustic)) / (100. * 100. * 100.)) / (Melt_metal[mt].DensityLiquid * pow(Melt_metal[mt].u0_Liquid, 2) * 1e-6)) *
							(spl_Te.GetY(k * dz_acoustic) - 1.)) << endl;
					}
				}
				//else
				//{
				//	if ((*it_z).index_x == (Nx_acoustic / 2) && (*it_z).index_y == (Ny_acoustic / 2) && (*it_z).index_z == k)
				//	{
				//		fout_Pac << 1e+4 * dz_acoustic * k / kabs << "   " << 0. << endl;
				//		fout_PTi << 1e+4 * dz_acoustic * k / kabs << "   " << 0. << endl;
				//		fout_PTe << 1e+4 * dz_acoustic * k / kabs << "   " << 0. << endl;
				//		it_z++;
				//	}
				//	else
				//	{
				//		fout_Pac << 1e+4 * dz_acoustic * k / kabs << "   " << p0v * (1. - V2[Nx_acoustic / 2][Ny_acoustic / 2][k]) << endl;

				//		if ((spl_Ti.GetY(k * dz_acoustic) * param.T00) < (Melt_metal[mt].T_melting - 1.))
				//		{
				//			if ((spl_Ti.GetY(k * dz_acoustic) * param.T00) >= 1336.)
				//			{
				//				fout_PTi << 1e+4 * dz_acoustic * k / kabs << "   " << p0v * ((Melt_metal[mt].gi * param.T00 * (spl_C_l_on_T.GetY(param.T00 * spl_Ti.GetY(k * dz_acoustic))) / (Melt_metal[mt].Density * Melt_metal[mt].Density * pow(Melt_metal[mt].u0, 2) * V0 * 1e-6)) *
				//					(spl_Ti.GetY(k * dz_acoustic) - 1.) / V2[Nx_acoustic / 2][Ny_acoustic / 2][k]) << endl;
				//			}
				//			else
				//			{
				//				fout_PTi << 1e+4 * dz_acoustic * k / kabs << "   " << p0v * ((Melt_metal[mt].gi * param.T00 * (Dependence_C_l_on_T(mt, param.T00 * spl_Ti.GetY(k * dz_acoustic)) / (100. * 100. * 100.)) / (Melt_metal[mt].Density * Melt_metal[mt].Density * pow(Melt_metal[mt].u0, 2) * V0 * 1e-6)) *
				//					(spl_Ti.GetY(k * dz_acoustic) - 1.) / V2[Nx_acoustic / 2][Ny_acoustic / 2][k]) << endl;
				//			}

				//			fout_PTe << 1e+4 * dz_acoustic * k / kabs << "   " << p0v * ((Melt_metal[mt].ge * param.T00 * (spl_C_e_on_T.GetY(param.T00 * spl_Te.GetY(k * dz_acoustic)) / (100. * 100. * 100.)) / (Melt_metal[mt].Density * pow(Melt_metal[mt].u0, 2) * 1e-6)) *
				//				(spl_Te.GetY(k * dz_acoustic) - 1.)) << endl;
				//		}

				//		if ((spl_Ti.GetY(k * dz_acoustic) * param.T00) >= (Melt_metal[mt].T_melting - 1.) && (spl_Ti.GetY(k * dz_acoustic) * param.T00) <= (Melt_metal[mt].T_melting + 1.))
				//		{
				//			if ((spl_Ti.GetY(k * dz_acoustic) * param.T00) >= 1336.)
				//			{
				//				fout_PTi << 1e+4 * dz_acoustic * k / kabs << "   " << p0v * ((Melt_metal[mt].gi * param.T00 * ((( /*Dependence_C_l_on_T(mt, param.T00* spl_Ti.GetY(k* dz_acoustic)) / (100. * 100. * 100.) + (Melt_metal[mt].Q_fusion * Melt_metal[mt].Density) / (2 * 1.0)*/
				//					spl_C_l_on_T.GetY(param.T00 * spl_Ti.GetY(k * dz_acoustic))))) / (pow((Melt_metal[mt].Density + Melt_metal[mt].DensityLiquid) / 2, 2) * pow((Melt_metal[mt].u0 + Melt_metal[mt].u0_Liquid) / 2, 2) * V0 * 1e-6)) *
				//					(spl_Ti.GetY(k * dz_acoustic) - 1.) / V2[Nx_acoustic / 2][Ny_acoustic / 2][k]) << endl;

				//				/*fout_PTi << 1e+4 * dz_acoustic * k / kabs << "   " << p0v * ((Melt_metal[mt].gi * param.T00 * (((Dependence_C_l_on_T(mt, param.T00 * spl_Ti.GetY(k * dz_acoustic)) / (100. * 100. * 100.) + (Melt_metal[mt].Q_fusion * Melt_metal[mt].Density) / (2 * 1.0)))) / (pow((Melt_metal[mt].Density + Melt_metal[mt].DensityLiquid) / 2, 2) * pow((Melt_metal[mt].u0 + Melt_metal[mt].u0_Liquid) / 2, 2) * V0 * 1e-6)) *
				//					(spl_Ti.GetY(k * dz_acoustic) - 1.) / V2[Nx_acoustic / 2][Ny_acoustic / 2][k]) << endl;*/
				//			}

				//			fout_PTe << 1e+4 * dz_acoustic * k / kabs << "   " << p0v * ((Melt_metal[mt].ge * param.T00 * (spl_C_e_on_T.GetY(param.T00 * spl_Te.GetY(k * dz_acoustic)) / (100. * 100. * 100.)) / (((Melt_metal[mt].Density + Melt_metal[mt].DensityLiquid) / 2) * pow((Melt_metal[mt].u0 + Melt_metal[mt].u0_Liquid) / 2, 2) * 1e-6)) *
				//				(spl_Te.GetY(k * dz_acoustic) - 1.)) << endl;
				//		}
				//	}
				//}
			}

			for (int i = 0; i < Nx_acoustic; i++)// цикл по x
			{
				spl_Te.InterpolateFast(1, tmpe2, i, Ny_acoustic / 2);
				spl_Ti.InterpolateFast(1, tmpi2, i, Ny_acoustic / 2);
				//if (points_rupture.empty())
				{
					fout_Pacx << 1e+4 * r0 * dx_acoustic * i << "   " << p0v * (1. - V2[i][Ny_acoustic / 2][1]) << endl;

					if ((spl_Ti.GetY(1 * dz_acoustic) * param.T00) < (Melt_metal[mt].T_melting - 1.))
					{
						/*if ((spl_Ti.GetY(1 * dz_acoustic) * param.T00) >= 1336.)
						{
							fout_PTix << 1e+4 * r0 * dx_acoustic * i << "   " << p0v * ((Melt_metal[mt].gi * param.T00 * (spl_C_l_on_T.GetY(spl_Ti.GetY(1 * dz_acoustic) * param.T00)) / (Melt_metal[mt].Density * Melt_metal[mt].Density * pow(Melt_metal[mt].u0, 2) * V0 * 1e-6)) *
								(spl_Ti.GetY(1 * dz_acoustic) - 1.) / V2[i][Ny_acoustic / 2][1]) << endl;
						}
						else
						{*/
							fout_PTix << 1e+4 * r0 * dx_acoustic * i << "   " << p0v * ((Melt_metal[mt].gi * param.T00 * (Dependence_C_l_on_T(mt, param.T00 * spl_Ti.GetY(1 * dz_acoustic)) / (100. * 100. * 100.)) / (Melt_metal[mt].Density * Melt_metal[mt].Density * pow(Melt_metal[mt].u0, 2) * V0 * 1e-6)) *
								(spl_Ti.GetY(1 * dz_acoustic) - 1.) / V2[i][Ny_acoustic / 2][1]) << endl;
						//}

						fout_PTex << 1e+4 * r0 * dx_acoustic * i << "   " << p0v * ((Melt_metal[mt].ge * param.T00 * (spl_C_e_on_T.GetY(param.T00 * spl_Te.GetY(1 * dz_acoustic)) / (100. * 100. * 100.)) / (Melt_metal[mt].Density * pow(Melt_metal[mt].u0, 2) * 1e-6)) *
							(spl_Te.GetY(1 * dz_acoustic) - 1.)) << endl;
					}

					if ((spl_Ti.GetY(1 * dz_acoustic) * param.T00) >= (Melt_metal[mt].T_melting - 1.) && (spl_Ti.GetY(1 * dz_acoustic) * param.T00) <= (Melt_metal[mt].T_melting + 1.))
					{
						// новый способ плавления
						fout_PTix << 1e+4 * r0 * dx_acoustic * i << "   " << p0v * ((Melt_metal[mt].gi * param.T00 * ((spl_C_l_on_T.GetY(spl_Ti.GetY(1 * dz_acoustic) * param.T00))) / (pow((Melt_metal[mt].Density + Melt_metal[mt].DensityLiquid) / 2, 2) * pow((Melt_metal[mt].u0 + Melt_metal[mt].u0_Liquid) / 2, 2) * V0 * 1e-6)) *
							(spl_Ti.GetY(1 * dz_acoustic) - 1.) / V2[i][Ny_acoustic / 2][1]) << endl;

						// старый способ плавления
						/*fout_PTix << 1e+4 * r0 * dx_acoustic * i << "   " << p0v * ((Melt_metal[mt].gi * param.T00 * (((Dependence_C_l_on_T(mt, param.T00 * spl_Ti.GetY(1 * dz_acoustic)) / (100. * 100. * 100.) + (Melt_metal[mt].Q_fusion * Melt_metal[mt].Density) / (2 * 1.0)))) / (pow((Melt_metal[mt].Density + Melt_metal[mt].DensityLiquid) / 2, 2) * pow((Melt_metal[mt].u0 + Melt_metal[mt].u0_Liquid) / 2, 2) * V0 * 1e-6)) *
							(spl_Ti.GetY(1 * dz_acoustic) - 1.) / V2[i][Ny_acoustic / 2][1]) << endl;*/

						fout_PTex << 1e+4 * r0 * dx_acoustic * i << "   " << p0v * ((Melt_metal[mt].ge * param.T00 * (spl_C_e_on_T.GetY(param.T00 * spl_Te.GetY(1 * dz_acoustic)) / (100. * 100. * 100.)) / (((Melt_metal[mt].Density + Melt_metal[mt].DensityLiquid) / 2) * pow((Melt_metal[mt].u0 + Melt_metal[mt].u0_Liquid) / 2, 2) * 1e-6)) *
							(spl_Te.GetY(1 * dz_acoustic) - 1.)) << endl;
					}

					if ((spl_Ti.GetY(1 * dz_acoustic) * param.T00) > (Melt_metal[mt].T_melting + 1.))
					{
						fout_PTix << 1e+4 * r0 * dx_acoustic * i << "   " << p0v * ((Melt_metal[mt].gi * param.T00 * (Dependence_C_l_on_T(mt, param.T00 * spl_Ti.GetY(1 * dz_acoustic)) / (100. * 100. * 100.)) / (Melt_metal[mt].DensityLiquid * Melt_metal[mt].DensityLiquid * pow(Melt_metal[mt].u0_Liquid, 2) * V0 * 1e-6)) *
							(spl_Ti.GetY(1 * dz_acoustic) - 1.) / V2[i][Ny_acoustic / 2][1]) << endl;

						fout_PTex << 1e+4 * r0 * dx_acoustic * i << "   " << p0v * ((Melt_metal[mt].ge * param.T00 * (spl_C_e_on_T.GetY(param.T00 * spl_Te.GetY(1 * dz_acoustic)) / (100. * 100. * 100.)) / (Melt_metal[mt].DensityLiquid * pow(Melt_metal[mt].u0_Liquid, 2) * 1e-6)) *
							(spl_Te.GetY(1 * dz_acoustic) - 1.)) << endl;
					}
				}
				//else
				//{
				//	if ((*it_x).index_x == i && (*it_x).index_y == (Ny_acoustic / 2) && (*it_x).index_z == 1)
				//	{
				//		fout_Pacx << 1e+4 * r0 * dx_acoustic * i << "   " << 0. << endl;
				//		fout_PTix << 1e+4 * r0 * dx_acoustic * i << "   " << 0. << endl;
				//		fout_PTex << 1e+4 * r0 * dx_acoustic * i << "   " << 0. << endl;
				//		it_x++;
				//	}
				//	else
				//	{
				//		fout_Pacx << 1e+4 * r0 * dx_acoustic * i << "   " << p0v * (1. - V2[i][Ny_acoustic / 2][1]) << endl;

				//		if ((spl_Ti.GetY(1 * dz_acoustic) * param.T00) < (Melt_metal[mt].T_melting - 1.))
				//		{
				//			if ((spl_Ti.GetY(1 * dz_acoustic) * param.T00) >= 1336.)
				//			{
				//				fout_PTix << 1e+4 * r0 * dx_acoustic * i << "   " << p0v * ((Melt_metal[mt].gi * param.T00 * (spl_C_l_on_T.GetY(spl_Ti.GetY(1 * dz_acoustic) * param.T00)) / (Melt_metal[mt].Density * Melt_metal[mt].Density * pow(Melt_metal[mt].u0, 2) * V0 * 1e-6)) *
				//					(spl_Ti.GetY(1 * dz_acoustic) - 1.) / V2[i][Ny_acoustic / 2][1]) << endl;
				//			}
				//			else
				//			{
				//				fout_PTix << 1e+4 * r0 * dx_acoustic * i << "   " << p0v * ((Melt_metal[mt].gi * param.T00 * (Dependence_C_l_on_T(mt, param.T00 * spl_Ti.GetY(1 * dz_acoustic)) / (100. * 100. * 100.)) / (Melt_metal[mt].Density * Melt_metal[mt].Density * pow(Melt_metal[mt].u0, 2) * V0 * 1e-6)) *
				//					(spl_Ti.GetY(1 * dz_acoustic) - 1.) / V2[i][Ny_acoustic / 2][1]) << endl;
				//			}

				//			fout_PTex << 1e+4 * r0 * dx_acoustic * i << "   " << p0v * ((Melt_metal[mt].ge * param.T00 * (spl_C_e_on_T.GetY(param.T00 * spl_Te.GetY(1 * dz_acoustic)) / (100. * 100. * 100.)) / (Melt_metal[mt].Density * pow(Melt_metal[mt].u0, 2) * 1e-6)) *
				//				(spl_Te.GetY(1 * dz_acoustic) - 1.)) << endl;
				//		}

				//		if ((spl_Ti.GetY(1 * dz_acoustic) * param.T00) >= (Melt_metal[mt].T_melting - 1.) && (spl_Ti.GetY(1 * dz_acoustic) * param.T00) <= (Melt_metal[mt].T_melting + 1.))
				//		{
				//			fout_PTix << 1e+4 * r0 * dx_acoustic * i << "   " << p0v * ((Melt_metal[mt].gi * param.T00 * ((spl_C_l_on_T.GetY(spl_Ti.GetY(1 * dz_acoustic) * param.T00))) / (pow((Melt_metal[mt].Density + Melt_metal[mt].DensityLiquid) / 2, 2) * pow((Melt_metal[mt].u0 + Melt_metal[mt].u0_Liquid) / 2, 2) * V0 * 1e-6)) *
				//				(spl_Ti.GetY(1 * dz_acoustic) - 1.) / V2[i][Ny_acoustic / 2][1]) << endl;

				//			/*fout_PTix << 1e+4 * r0 * dx_acoustic * i << "   " << p0v * ((Melt_metal[mt].gi * param.T00 * (((Dependence_C_l_on_T(mt, param.T00 * spl_Ti.GetY(1 * dz_acoustic)) / (100. * 100. * 100.) + (Melt_metal[mt].Q_fusion * Melt_metal[mt].Density) / (2 * 1.0)))) / (pow((Melt_metal[mt].Density + Melt_metal[mt].DensityLiquid) / 2, 2) * pow((Melt_metal[mt].u0 + Melt_metal[mt].u0_Liquid) / 2, 2) * V0 * 1e-6)) *
				//				(spl_Ti.GetY(1 * dz_acoustic) - 1.) / V2[i][Ny_acoustic / 2][1]) << endl;*/

				//			fout_PTex << 1e+4 * r0 * dx_acoustic * i << "   " << p0v * ((Melt_metal[mt].ge * param.T00 * (spl_C_e_on_T.GetY(param.T00 * spl_Te.GetY(1 * dz_acoustic)) / (100. * 100. * 100.)) / (((Melt_metal[mt].Density + Melt_metal[mt].DensityLiquid) / 2) * pow((Melt_metal[mt].u0 + Melt_metal[mt].u0_Liquid) / 2, 2) * 1e-6)) *
				//				(spl_Te.GetY(1 * dz_acoustic) - 1.)) << endl;
				//		}
				//	}

				//}
			}

			//system("pause");

			namefile.clear();
			namefile = "Pacoustic (z,mkm) " + metall + "," + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs";
			pac.ShowDataOnPlot2D(0, 1, Pfs, namefile, true);
			Sleep(1000);
			namefile.clear();
			namefile = "PTi (z,mkm) " + metall + "," + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs";
			pac.ShowDataOnPlot2D(1, 1, Pfs, namefile, true);
			Sleep(1000);
			namefile.clear();
			namefile = "PTe (z,mkm) " + metall + "," + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs";
			pac.ShowDataOnPlot2D(2, 1, Pfs, namefile, true);
			Sleep(1000);
			namefile.clear();
			namefile = "Pacoustic (x,mkm) " + metall + "," + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs";
			pac.ShowDataOnPlot2D(3, 1, Pfs, namefile, true);
			Sleep(1000);
			namefile.clear();
			namefile = "PTi (x,mkm) " + metall + "," + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs";
			pac.ShowDataOnPlot2D(4, 1, Pfs, namefile, true);
			Sleep(1000);
			namefile.clear();
			namefile = "PTe (x,mkm) " + metall + "," + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs";
			pac.ShowDataOnPlot2D(5, 1, Pfs, namefile, true);
			Sleep(1000);

			fout_Pac.close();
			fout_PTe.close();
			fout_PTi.close();
			fout_Pacx.close();
			fout_PTex.close();
			fout_PTix.close();
			points_rupture_for_plot_z.clear();
			points_rupture_for_plot_z.erase(points_rupture_for_plot_z.begin(), points_rupture_for_plot_z.end());
			points_rupture_for_plot_x.clear();
			points_rupture_for_plot_x.erase(points_rupture_for_plot_x.begin(), points_rupture_for_plot_x.end());
			fout_Pac.open(Pacteti[0]);
			fout_PTi.open(Pacteti[1]);
			fout_PTe.open(Pacteti[2]);
			fout_Pacx.open(Pacteti[3]);
			fout_PTix.open(Pacteti[4]);
			fout_PTex.open(Pacteti[5]);

			// открытие файлов произвожу выше

			//plt.SetDataOnPlot3D(21, Nx_acoustic, Ny_acoustic, Nz_acoustic, 1e+4 * dx_acoustic * r0, 1e+4 * dy_acoustic * r0, 1e+4 * dz_acoustic / kabs, Pacoustic, p0v, NULL, Ny_acoustic / 2, 1, 1, 1, NULL, empty, fun_on_x);//??
			//plt.SetDataOnPlot3D(22, Nx_acoustic, Ny_acoustic, Nz_acoustic, 1e+4 * dx_acoustic * r0, 1e+4 * dy_acoustic * r0, 1e+4 * dz_acoustic / kabs, Pacoustic, p0v, Nx_acoustic / 2, Ny_acoustic / 2, NULL, 1, 1, NULL, empty, fun_on_z);//??
			//plt.SetDataOnPlot3D(24, Nx_acoustic, Ny_acoustic, Nz_acoustic, 1e+4 * dx_acoustic * r0, 1e+4 * dy_acoustic * r0, 1e+4 * dz_acoustic / kabs, PTi, p0v, NULL, Ny_acoustic / 2, 1, 1, 1, NULL, empty, fun_on_x);//??
			//plt.SetDataOnPlot3D(25, Nx_acoustic, Ny_acoustic, Nz_acoustic, 1e+4 * dx_acoustic * r0, 1e+4 * dy_acoustic * r0, 1e+4 * dz_acoustic / kabs, PTi, p0v, Nx_acoustic / 2, Ny_acoustic / 2, NULL, 1, 1, NULL, empty, fun_on_z);//??
			//plt.SetDataOnPlot3D(27, Nx_acoustic, Ny_acoustic, Nz_acoustic, 1e+4 * dx_acoustic * r0, 1e+4 * dy_acoustic * r0, 1e+4 * dz_acoustic / kabs, PTe, p0v, NULL, Ny_acoustic / 2, 1, 1, 1, NULL, empty, fun_on_x);//??
			//plt.SetDataOnPlot3D(28, Nx_acoustic, Ny_acoustic, Nz_acoustic, 1e+4 * dx_acoustic * r0, 1e+4 * dy_acoustic * r0, 1e+4 * dz_acoustic / kabs, PTe, p0v, Nx_acoustic / 2, Ny_acoustic / 2, NULL, 1, 1, NULL, empty, fun_on_z);//??

			namefile.clear();
			namefile = "Te,Ti (x,mkm) " + metall + "," + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs";
			plt.ShowDataOnPlot2D(15, 2, TeTiTfs, namefile, true);
			Sleep(1000);
			namefile.clear();
			namefile = "Te,Ti (z,mkm) " + metall + "," + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs";
			plt.ShowDataOnPlot2D(16, 2, TeTiTfs, namefile, true);
			Sleep(1000);

			namefile.clear();
			namefile = "P (x,mkm) " + metall + "," + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs";
			plt.ShowDataOnPlot2D(18, 1, Pfs, namefile, true);
			Sleep(1000);
			namefile.clear();
			namefile = "P (z,mkm) " + metall + "," + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs";
			plt.ShowDataOnPlot2D(19, 1, Pfs, namefile, true);
			Sleep(1000);

			/*namefile.clear();
			namefile = "Pacoustic (x,mkm) " + metall + "," + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs";
			plt.ShowDataOnPlot2D(21, 1, Pfs, namefile, true);
			Sleep(1000);
			namefile.clear();
			namefile = "Pacoustic (z,mkm) " + metall + "," + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs";
			plt.ShowDataOnPlot2D(22, 1, Pfs, namefile, true);
			Sleep(1000);
			namefile.clear();
			namefile = "PTi (x,mkm) " + metall + "," + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs";
			plt.ShowDataOnPlot2D(24, 1, Pfs, namefile, true);
			Sleep(1000);
			namefile.clear();
			namefile = "PTi (z,mkm) " + metall + "," + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs";
			plt.ShowDataOnPlot2D(25, 1, Pfs, namefile, true);
			Sleep(1000);
			namefile.clear();
			namefile = "PTe (x,mkm) " + metall + "," + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs";
			plt.ShowDataOnPlot2D(27, 1, Pfs, namefile, true);
			Sleep(1000);
			namefile.clear();
			namefile = "PTe (z,mkm) " + metall + "," + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs";
			plt.ShowDataOnPlot2D(28, 1, Pfs, namefile, true);
			Sleep(1000);
*/

			plt.Close_and_open_files_for_replot(for_close_and_open);
		}
		// next step
		n++;
		cout << " step = " << n << endl;

		// Perehod na sledujuschij shag po vremeni
		tmpe0 = tmpe1;
		tmpe1 = tmpe2;

		tmpi0 = tmpi1;
		tmpi1 = tmpi2;

		a1x = a2x;
		a1y = a2y;
		a1z = a2z;
		b1x = b2x;
		b1y = b2y;
		b1z = b2z;

		e1 = e2;

		cout << endl;
	} while (tt <= 8000/*nn < 1000*/);

	t = clock() - t;

	cout << endl << endl;
	cout << "It took me %d clicks (%f seconds)." << endl <<
		(int)t << "   " << ((double)t) / CLOCKS_PER_SEC;
	cout << endl;

	// могут приодитьс (если вдруг моментов времени меньше заявленных линий на графике, и как результат мы не запишем данные из матрицы в файл )
	if (npxzt < count_of_lines)
	{
		npxzt = count_of_lines;
		plt.SetDataOnPlot3D(8, Nx_acoustic, Ny_acoustic, Nz_acoustic, 1e+4 * dx_acoustic * r0, 1e+4 * dy_acoustic * r0, 1e+4 * dz_acoustic / kabs, e1, p0v, NULL, Ny_acoustic / 2, 10, count_of_lines, npxzt, NULL, empty, fun_on_x);
		plt.SetDataOnPlot3D(9, Nx_acoustic, Ny_acoustic, Nz_acoustic, 1e+4 * dx_acoustic * r0, 1e+4 * dy_acoustic * r0, 1e+4 * dz_acoustic / kabs, e1, p0v, Nx_acoustic / 2, Ny_acoustic / 2, NULL, count_of_lines, npxzt, NULL, empty, fun_on_z);
	}

	//double*** Field_melting;
	//Field_melting = new double** [Nx];
	//for (int i = 0; i < Nx; i++)
	//{
	//	Field_melting[i] = new double* [Ny];
	//	for (int j = 0; j < Ny; j++)
	//	{
	//		Field_melting[i][j] = new double[Nz];
	//	}
	//}

	//for (int i = 0; i < Nx; i++)
	//{
	//	for (int j = 0; j < Ny; j++)
	//	{
	//		for (int k = 0; k < Nz; k++)
	//		{
	//			Field_melting[i][j][k] = Melt_metal[mt].T_melting;// 1e-16;
	//		}
	//	}
	//}

	//if (current_number_line_melting < 15)
	//{
	//	current_number_line_melting = 11;
	//	plt.SetDataOnPlot3D(12, Nx, Ny, Nz, 1e+4 * dx * param.r0, 1e+4 * dy * param.r0, 1e+4 * dz / param.kabs, Field_melting, 1., NULL, Ny / 2, 1, 11, current_number_line_melting, NULL, empty, fun_on_x);
	//	plt.SetDataOnPlot3D(11, Nx, Ny, Nz, 1e+4 * dx * param.r0, 1e+4 * dy * param.r0, 1e+4 * dz / param.kabs, Field_melting, 1., Nx / 2, Ny / 2, NULL, 11, current_number_line_melting, NULL, empty, fun_on_z);
	//}

	//vector<string> Legenda = { "500 fs", "1000 fs", "1500 fs","2000 fs" , "3000 fs" , "4000 fs" , "5000 fs" , "6000 fs", "7000 fs", "8000 fs" };
	namefile.clear();
	namefile = "file P(x), bar " + metall + ", " + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs";
	plt.ShowDataOnPlot2D(8, count_of_lines, Legenda, namefile, true);
	namefile.clear();
	namefile = "file P(z), bar " + metall + ", " + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs";
	plt.ShowDataOnPlot2D(9, count_of_lines, Legenda, namefile, true);
	namefile.clear();
	//vector<string> Legenda_melting = { "2500 fs" , "3000 fs", "3500 fs", "4000 fs", "4500 fs", "5000 fs", "5500 fs", "6000 fs", "6500 fs", "7000 fs","T melting" };
	vector<string> Legenda_melting = Legenda_melting_tmp;
	Legenda_melting.push_back("T melting");
	namefile = "file Profile Teperature T(x)";
	//plt.ShowDataOnPlot2D(12, 11, Legenda_melting, namefile, true);
	namefile.clear();
	namefile = "file Profile Teperature T(z)";
	//plt.ShowDataOnPlot2D(11, 11, Legenda_melting, namefile, true);
	namefile.clear();

	sum_t = clock() - t;

	plt.Close_all_files_and_plots(number_plots);
}


inline int random_ab(int a, int b, std::default_random_engine& g)
{
	//return a + (b - a) * ((double)rand() / (double)RAND_MAX);
	std::uniform_int_distribution<int> distribution(a, b);
	return distribution(g);
}

int main()
{
	TypeBeam tp = Gauss;
	cout << " tp = " << tp << endl;
	cout << endl;

	int dxy = 100;
	int Nx_heat = dxy; // 100 //100; // chislo uzlov po x // 100
	int Ny_heat = dxy;// 100// 100; // chislo uzlov po y //100
	int Nz_heat = 250;//250; //!!!!!!!!!!!!!!

	int Nx_acoustic = dxy; // 100 //100; // chislo uzlov po x // 100
	int Ny_acoustic = dxy;// 100// 100; // chislo uzlov po y //100
	int Nz_acoustic = 400;// 400;// 10000; // !!!!!!

	double*** V2, *** a1y, *** a2y, *** a1z, *** a2z, *** a2x, *** a1x, *** b1x, *** b2x, *** b1y, *** b2y, *** b1z, *** b2z, *** e2, *** e1, *** tmpe0, *** tmpe1, *** tmpe2, *** tmpi0,*** tmpi1, *** tmpi2;    // davlenie, gorizontal'naja skorost', vertikal'naja skorost', temperatura
	double*** Ti_acoustic, *** Te_acoustic; // массивы для переопрелеления сетки температуры под акустику
																																																						 
	// 21 массива

	V2 = new double** [Nx_acoustic];
	a1y = new double** [Nx_acoustic]; // Euler coordinates
	a2y = new double** [Nx_acoustic]; // Euler coordinates
	a1z = new double** [Nx_acoustic]; // Euler coordinates
	a2z = new double** [Nx_acoustic]; // Euler coordinates
	a1x = new double** [Nx_acoustic]; // Euler coordinates
	a2x = new double** [Nx_acoustic]; // Euler coordinates
	b1x = new double** [Nx_acoustic]; // particle velocity
	b2x = new double** [Nx_acoustic]; // particle velocity
	b1y = new double** [Nx_acoustic]; // particle velocity
	b2y = new double** [Nx_acoustic]; // particle velocity
	b1z = new double** [Nx_acoustic]; // particle velocity
	b2z = new double** [Nx_acoustic]; // particle velocity
	e1 = new double** [Nx_acoustic]; // Pressure
	e2 = new double** [Nx_acoustic]; // Pressure
	tmpe0 = new double** [Nx_heat]; // temperature of electrons 
	tmpe1 = new double** [Nx_heat]; // temperature of electrons 
	tmpe2 = new double** [Nx_heat]; // temperature of electrons 
	tmpi0 = new double** [Nx_heat]; // temperature of ions
	tmpi1 = new double** [Nx_heat]; // temperature of ions
	tmpi2 = new double** [Nx_heat]; // temperature of ions
	Ti_acoustic = new double** [Nx_acoustic];
	Te_acoustic = new double** [Nx_acoustic];
	for (int i = 0; i < Nx_acoustic; i++)
	{
		V2[i] = new double* [Ny_acoustic];
		a1y[i] = new double* [Ny_acoustic];
		a2y[i] = new double* [Ny_acoustic];
		a1z[i] = new double* [Ny_acoustic];
		a2z[i] = new double* [Ny_acoustic];
		a1x[i] = new double* [Ny_acoustic];
		a2x[i] = new double* [Ny_acoustic];
		b1x[i] = new double* [Ny_acoustic];
		b2x[i] = new double* [Ny_acoustic];
		b1y[i] = new double* [Ny_acoustic];
		b2y[i] = new double* [Ny_acoustic];
		b1z[i] = new double* [Ny_acoustic];
		b2z[i] = new double* [Ny_acoustic];
		e1[i] = new double* [Ny_acoustic];
		e2[i] = new double* [Ny_acoustic];
		Ti_acoustic[i] = new double* [Ny_acoustic];
		Te_acoustic[i] = new double* [Ny_acoustic];
		for (int j = 0; j < Ny_acoustic; j++)
		{
			V2[i][j] = new double[Nz_acoustic];
			a1y[i][j] = new double[Nz_acoustic];
			a2y[i][j] = new double[Nz_acoustic];
			a1z[i][j] = new double[Nz_acoustic];
			a2z[i][j] = new double[Nz_acoustic];
			a1x[i][j] = new double[Nz_acoustic];
			a2x[i][j] = new double[Nz_acoustic];
			b1x[i][j] = new double[Nz_acoustic];
			b2x[i][j] = new double[Nz_acoustic];
			b1y[i][j] = new double[Nz_acoustic];
			b2y[i][j] = new double[Nz_acoustic];
			b1z[i][j] = new double[Nz_acoustic];
			b2z[i][j] = new double[Nz_acoustic];
			e1[i][j] = new double[Nz_acoustic];
			e2[i][j] = new double[Nz_acoustic];
			Ti_acoustic[i][j] = new double[Nz_acoustic];
			Te_acoustic[i][j] = new double[Nz_acoustic];
		}
	}

	for (int i = 0; i < Nx_heat; i++)
	{
		tmpe0[i] = new double* [Ny_heat];
		tmpe1[i] = new double* [Ny_heat];
		tmpe2[i] = new double* [Ny_heat];
		tmpi0[i] = new double* [Ny_heat];
		tmpi1[i] = new double* [Ny_heat];
		tmpi2[i] = new double* [Ny_heat];
		for (int j = 0; j < Ny_heat; j++)
		{
			tmpe0[i][j] = new double[Nz_heat];
			tmpe1[i][j] = new double[Nz_heat];
			tmpe2[i][j] = new double[Nz_heat];
			tmpi0[i][j] = new double[Nz_heat];
			tmpi1[i][j] = new double[Nz_heat];
			tmpi2[i][j] = new double[Nz_heat];
		}
	}

	double TTT = 300.;
	/*  отрисовка графиков */

	Splayn spl_G_e_on_T, spl_C_e_on_T; //= Calculation_Interpolation("Ge_Au_new.txt");
	spl_C_e_on_T.Calculation_Interpolation("Ce_Au_new.txt");
	spl_G_e_on_T.Calculation_Interpolation("Ge_Au_new.txt");
	cout << "C_e_on_T = " << spl_C_e_on_T.GetY(TTT) / (100. * 100. * 100.) << "   " << 5.472e-2 << endl;
	cout << "G_e_on_T = " << spl_G_e_on_T.GetY(TTT) / (100. * 100. * 100.) << "   " << 2.5e+10 << endl;
	cout << " kte = " << Dependence_k_e_on_T(Au, TTT) / (100.) << "   " << 3.115 << endl;
	cout << " Cl = " << Dependence_C_l_on_T(Al, TTT) / (100. * 100. * 100.) << "   " << 2.550 << endl;

	ofstream fout1("Ce_interp.txt");
	ofstream fout2("Ge_interp.txt");
	ofstream fout3("ke_depend.txt");
	ofstream fout4("Cl_depend.txt");

	for (double i = 300.; i < 50250.; i++)
	{
		fout1 << i << "   " << spl_C_e_on_T.GetY(i) << endl;
		fout2 << i << "   " << spl_G_e_on_T.GetY(i) << endl;
		fout3 << i << "   " << Dependence_k_e_on_T(Au, i) << endl;
	}


	/*for (double i = 300.; i < 3000.; i++)
	{
		fout4 << i << "   " << Dependence_C_l_on_T(Au, i) << endl;
	}*/
	fout4.precision(10);
	double d = 0.5;//121.;

	for (double i = 1338 - d - 2; i <= 1338. + d + 2; i += 0.5)
	{
		if (i < (1338. - d))
		{
			//////fout4 << i << "   " << (1.11 * 300. * (Dependence_C_l_on_T(Au, i) / (100. * 100. * 100.)) / (19.3 * 19.3 * pow(3.24e+5, 2) * 1e-6)) << endl;
			fout4 << i << "   " << (Dependence_C_l_on_T(Au, i) / (100. * 100. * 100.)) << endl;
		}

		if (i >= (1338. - d) && i <= (1338. + d))
		{
			///	//fout4 << i << "   " << (1.11* 300.* (((Dependence_C_l_on_T(Au, i) / (100. * 100. * 100.) + (63.7 * 19.3) / (2 * 1.)))) / (pow((19.3 + 17.) / 2, 2) * pow((3.24e+5 + 2.567e+5) / 2, 2) * 1e-6)) << endl;
			fout4 << i << "   " << (((Dependence_C_l_on_T(Au, i) / (100. * 100. * 100.) + (63.7 * 19.3) * delta_function(i, 1338., d)))) << endl;
			///	//fout4 << i << "   " << 0. << endl;
		}

		if (i > (1338. + d))
		{
			fout4 << i << "   " << (Dependence_C_l_on_T(Au, i) / (100. * 100. * 100.)) << endl;
			///////fout4 << i << "   " << (1.11 * 300.* (Dependence_C_l_on_T(Au, i) / (100. * 100. * 100.)) / (17. * 17. * pow(2.567e+5, 2) * 1e-6)) << endl;
		}
		//(Melt_metal[mt].gi * param.T00 * (Dependence_C_l_on_T(mt, param.T00 * spl_Ti.GetY(k * dz)) / (100. * 100. * 100.)) / (Melt_metal[mt].Density * Melt_metal[mt].Density * pow(Melt_metal[mt].u0, 2) * V0 * 1e-6))
	}

	//system("pause");
	Splayn spl_C_l_on_T;
	spl_C_l_on_T.Calculation_Interpolation("Cl_depend.txt");
	fout4.close();
	fout4.open("Cl_depend.txt");

	for (double i = 1338 - d - 2.; i <= 1338. + d + 2.; i += 0.001)
		//for (double i = 1287.99; i <= 1288.05; i += 0.001)
	{
		fout4 << i << "   " << spl_C_l_on_T.GetY(i) << endl;
	}

	vector<string> vv = { "Ce_interp.txt", "Ge_interp.txt", "ke_depend.txt", "Cl_depend.txt", "Dependence Hmax on Fluence.txt" };
	GnuPlot plt(5, vv); // объект хранит 5 плотиков
	plt.SetParametrs2D(0, 1, 3, "C_e(T_e)", "T_e, K", "C_e, J/m^3K");
	plt.SetParametrs2D(1, 1, 3, "G_e(T_e)", "T_e, K", "G_e, W/m^3K");
	plt.SetParametrs2D(2, 1, 3, "k_e(T_e)", "T_e, K", "k_e, W/mK");
	plt.SetParametrs2D(3, 1, 3, "C_i(T_i)", "T_i, K", "C_i, J/m^3K");
	plt.SetParametrs2D(4, 1, 3, "h max(Fluence)", "F, J/m^2", "z, nm");

	string  namefile;
	vector<string> Metal = { "1" };
	namefile = "Dependence Ce on Te";
	plt.ShowDataOnPlot2D(0, 1, Metal, namefile, true);
	namefile.clear();
	namefile = "Dependence Ge on Te";
	plt.ShowDataOnPlot2D(1, 1, Metal, namefile, true);
	namefile.clear();
	namefile = "Dependence ke on Te";
	plt.ShowDataOnPlot2D(2, 1, Metal, namefile, true);
	namefile.clear();
	namefile = "Dependence Cl on Tl";
	plt.ShowDataOnPlot2D(3, 1, Metal, namefile, true);
	namefile.clear();
	namefile = "Dependence hmax on Fluence";
	plt.ShowDataOnPlot2D(4, 1, Metal, namefile, true);
	namefile.clear();

	//fout1.close();
	//fout2.close();
	//fout3.close();
	fout4.close();

	system("pause");


	//type metal, type beam, kte, Ce,Ci,	   gamma,    ge,   gi, us	tp s	
																								//I0 = 1e+7
	// double*** V2, double*** a1y, double*** a2y, double*** a1z, double*** a2z, double*** a2x, double*** a1x, double*** b1x, double*** b2x, double*** b1y, double*** b2y, double*** b1z, double*** b2z, double*** e2, double*** e1, double*** F, double*** tmpe0, double*** tmpe1, double*** tmpe2, double*** tmpi0, double*** tmpi1, double*** tmpi2, double** Pxt, double** Pzt, int Nx, int Ny, int Nz, Splayn spl_C_l_on_T)
	MainProcedure(Au, Gauss, 3.115, 19.32, 5.472e-2, 2.550, 2.5e+10, 1.5, 1.11, 3.24e+5, 1e-13, /*0.15e+12*/ 1.5e+11, V2, a1y, a2y, a1z, a2z, a2x, a1x, b1x, b2x, b1y, b2y, b1z, b2z, e2, e1, tmpe0, tmpe1, tmpe2, tmpi0, tmpi1, tmpi2, Nx_heat, Ny_heat, Nz_heat, Nx_acoustic, Ny_acoustic, Nz_acoustic, spl_C_l_on_T, Te_acoustic, Ti_acoustic, "Depth1.txt");
	system("pause");
	//MainProcedure(Ni, Gauss, 0.900, 8.96, 2.733e-2, 3.910, 4.05e+12, 1.5, 1.83, 5.63e+5, 1e-13, 1e+7, V2, a1y, a2y, a1z, a2z, a2x, a1x, b1x, b2x, b1y, b2y, b1z, b2z, e2, e1, tmpe0, tmpe1, tmpe2, tmpi0, tmpi1, tmpi2, Nx_heat, Ny_heat, Nz_heat, Nx_acoustic, Ny_acoustic, Nz_acoustic, spl_C_l_on_T, Te_acoustic, Ti_acoustic);
	system("pause");

	/*MainProcedure(Au, Gauss, 3.115, 19.32, 5.472e-2, 2.550, 2.5e+10, 1.5, 1.11, 3.24e+5, 5e-14, 0.2e+8, V2, a1y, a2y, a1z, a2z, a2x, a1x, b1x, b2x, b1y, b2y, b1z, b2z, e2, e1, F, tmpe0, tmpe1, tmpe2, tmpi1, tmpi2, Pxt, Pzt, Nx, Ny, Nz);
	MainProcedure(Au, Gauss, 3.115, 19.32, 5.472e-2, 2.550, 2.5e+10, 1.5, 1.11, 3.24e+5, 2e-13, 0.5e+7, V2, a1y, a2y, a1z, a2z, a2x, a1x, b1x, b2x, b1y, b2y, b1z, b2z, e2, e1, F, tmpe0, tmpe1, tmpe2, tmpi1, tmpi2, Pxt, Pzt, Nx, Ny, Nz);
	MainProcedure(Au, Gauss, 3.115, 19.32, 5.472e-2, 2.550, 2.5e+10, 1.5, 1.11, 3.24e+5, 5e-13, 0.2e+7, V2, a1y, a2y, a1z, a2z, a2x, a1x, b1x, b2x, b1y, b2y, b1z, b2z, e2, e1, F, tmpe0, tmpe1, tmpe2, tmpi1, tmpi2, Pxt, Pzt, Nx, Ny, Nz);
	*///MainProcedure(Au, Gauss, 3.115, 19.32, 5.472e-2, 2.550, 2.5e+10, 1.5, 1.11, 3.24e+5, 1e-13, V2, a1y, a2y, a1z, a2z, a2x, a1x, b1x, b2x, b1y, b2y, b1z, b2z, e2, e1, F, tmpe0, tmpe1, tmpe2, tmpi1, tmpi2, Pxt, Pzt, Nx, Ny, Nz);
	//MainProcedure(Au, Gauss, 3.115, 19.32, 5.472e-2, 2.550, 2.5e+10, 1.5, 1.11, 3.24e+5, 2e-13, V2, a1y, a2y, a1z, a2z, a2x, a1x, b1x, b2x, b1y, b2y, b1z, b2z, e2, e1, F, tmpe0, tmpe1, tmpe2, tmpi1, tmpi2, Pxt, Pzt, Nx, Ny, Nz);
	system("pause");

		//type metal, type beam, kte, ro0, 	Ce,			Ci,	   gamma,  ge, gi,		us	  tp s	
//	MainProcedure(Au, Vortex, 3.115, 19.32, 5.472e-2, 2.550, 2.5e+10, 1.5, 1.11, 3.24e+5, 1e-14, V2, a1y, a2y, a1z, a2z, a2x, a1x, b1x, b2x, b1y, b2y, b1z, b2z, e2, e1, F, tmpe0, tmpe1, tmpe2, tmpi1, tmpi2, Pxt, Pzt, Nx, Ny, Nz);
	//MainProcedure(Au, Vortex, 3.115, 19.32, 5.472e-2, 2.550, 2.5e+10, 1.5, 1.11, 3.24e+5, 5e-14, V2, a1y, a2y, a1z, a2z, a2x, a1x, b1x, b2x, b1y, b2y, b1z, b2z, e2, e1, F, tmpe0, tmpe1, tmpe2, tmpi1, tmpi2, Pxt, Pzt, Nx, Ny, Nz);
	//MainProcedure(Au, Vortex, 3.115, 19.32, 5.472e-2, 2.550, 2.5e+10, 1.5, 1.11, 3.24e+5, 1e-13, V2, a1y, a2y, a1z, a2z, a2x, a1x, b1x, b2x, b1y, b2y, b1z, b2z, e2, e1, F, tmpe0, tmpe1, tmpe2, tmpi1, tmpi2, Pxt, Pzt, Nx, Ny, Nz);
	//MainProcedure(Au, Vortex, 3.115, 19.32, 5.472e-2, 2.550, 2.5e+10, 1.5, 1.11, 3.24e+5, 2e-13, V2, a1y, a2y, a1z, a2z, a2x, a1x, b1x, b2x, b1y, b2y, b1z, b2z, e2, e1, F, tmpe0, tmpe1, tmpe2, tmpi1, tmpi2, Pxt, Pzt, Nx, Ny, Nz);

	//system("pause");
	//type metal, type beam, kte,   Ce,			Ci,	   gamma,  ge,  gi,     us	   tp s	
	//MainProcedure(Ni, Gauss, 0.900, 8.96, 2.733e-2, 3.910, 4.05e+12, 1.5, 1.83, 5.63e+5, 1e-14, V2, a1y, a2y, a1z, a2z, a2x, a1x, b1x, b2x, b1y, b2y, b1z, b2z, e2, e1, F, tmpe0, tmpe1, tmpe2, tmpi1, tmpi2, Pxt, Pzt, Nx, Ny, Nz);
	//MainProcedure(Ni, Gauss, 0.900, 8.96, 2.733e-2, 3.910, 4.05e+12, 1.5, 1.83, 5.63e+5, 5e-14, V2, a1y, a2y, a1z, a2z, a2x, a1x, b1x, b2x, b1y, b2y, b1z, b2z, e2, e1, F, tmpe0, tmpe1, tmpe2, tmpi1, tmpi2, Pxt, Pzt, Nx, Ny, Nz);
	//MainProcedure(Ni, Gauss, 0.900, 8.96, 2.733e-2, 3.910, 4.05e+12, 1.5, 1.83, 5.63e+5, 1e-13, V2, a1y, a2y, a1z, a2z, a2x, a1x, b1x, b2x, b1y, b2y, b1z, b2z, e2, e1, F, tmpe0, tmpe1, tmpe2, tmpi1, tmpi2, Pxt, Pzt, Nx, Ny, Nz);
	//MainProcedure(Ni, Gauss, 0.900, 8.96, 2.733e-2, 3.910, 4.05e+12, 1.5, 1.83, 5.63e+5, 2e-13, V2, a1y, a2y, a1z, a2z, a2x, a1x, b1x, b2x, b1y, b2y, b1z, b2z, e2, e1, F, tmpe0, tmpe1, tmpe2, tmpi1, tmpi2, Pxt, Pzt, Nx, Ny, Nz);

	//type metal, type beam,  kte,   Ce,			Ci,	   gamma,  ge,  gi,     us	   tp s		
	//MainProcedure(Ni, Vortex, 0.900, 8.96, 2.733e-2, 3.910, 4.05e+12, 1.5, 1.83, 5.63e+5, 1e-14, V2, a1y, a2y, a1z, a2z, a2x, a1x, b1x, b2x, b1y, b2y, b1z, b2z, e2, e1, F, tmpe0, tmpe1, tmpe2, tmpi1, tmpi2, Pxt, Pzt, Nx, Ny, Nz);
	//MainProcedure(Ni, Vortex, 0.900, 8.96, 2.733e-2, 3.910, 4.05e+12, 1.5, 1.83, 5.63e+5, 5e-14, V2, a1y, a2y, a1z, a2z, a2x, a1x, b1x, b2x, b1y, b2y, b1z, b2z, e2, e1, F, tmpe0, tmpe1, tmpe2, tmpi1, tmpi2, Pxt, Pzt, Nx, Ny, Nz);
	//MainProcedure(Ni, Vortex, 0.900, 8.96, 2.733e-2, 3.910, 4.05e+12, 1.5, 1.83, 5.63e+5, 1e-13, V2, a1y, a2y, a1z, a2z, a2x, a1x, b1x, b2x, b1y, b2y, b1z, b2z, e2, e1, F, tmpe0, tmpe1, tmpe2, tmpi1, tmpi2, Pxt, Pzt, Nx, Ny, Nz);
	//MainProcedure(Ni, Vortex, 0.900, 8.96, 2.733e-2, 3.910, 4.05e+12, 1.5, 1.83, 5.63e+5, 2e-13, V2, a1y, a2y, a1z, a2z, a2x, a1x, b1x, b2x, b1y, b2y, b1z, b2z, e2, e1, F, tmpe0, tmpe1, tmpe2, tmpi1, tmpi2, Pxt, Pzt, Nx, Ny, Nz);	// продублировать по разным длительностям импульса
	system("pause");


	//type metal, type beam, kte,   Ce,	???		Ci,	   gamma,  ge,  gi,     us	   tp s	
	/*MainProcedure(Cu, Gauss, 3.930, 8.92, 2.733e-2, 3.568, 1.0e+11, 1.5, 1.96, 4.73e+5, 1e-14, V2, a1y, a2y, a1z, a2z, a2x, a1x, b1x, b2x, b1y, b2y, b1z, b2z, e2, e1, F, tmpe0, tmpe1, tmpe2, tmpi1, tmpi2, Pxt, Pzt, Nx, Ny, Nz);
	MainProcedure(Cu, Gauss, 3.930, 8.92, 2.733e-2, 3.568, 1.0e+11, 1.5, 1.96, 4.73e+5, 5e-14, V2, a1y, a2y, a1z, a2z, a2x, a1x, b1x, b2x, b1y, b2y, b1z, b2z, e2, e1, F, tmpe0, tmpe1, tmpe2, tmpi1, tmpi2, Pxt, Pzt, Nx, Ny, Nz);
	MainProcedure(Cu, Gauss, 3.930, 8.92, 2.733e-2, 3.568, 1.0e+11, 1.5, 1.96, 4.73e+5, 1e-13, V2, a1y, a2y, a1z, a2z, a2x, a1x, b1x, b2x, b1y, b2y, b1z, b2z, e2, e1, F, tmpe0, tmpe1, tmpe2, tmpi1, tmpi2, Pxt, Pzt, Nx, Ny, Nz);
	MainProcedure(Cu, Gauss, 3.930, 8.92, 2.733e-2, 3.568, 1.0e+11, 1.5, 1.96, 4.73e+5, 2e-13, V2, a1y, a2y, a1z, a2z, a2x, a1x, b1x, b2x, b1y, b2y, b1z, b2z, e2, e1, F, tmpe0, tmpe1, tmpe2, tmpi1, tmpi2, Pxt, Pzt, Nx, Ny, Nz);

	//type metal, type beam,  kte,   Ce,??			Ci,	   gamma,  ge,  gi,     us	   tp s		
	MainProcedure(Cu, Vortex, 3.930, 8.92, 2.733e-2, 3.568, 1.0e+11, 1.5, 1.96, 4.73e+5, 1e-14, V2, a1y, a2y, a1z, a2z, a2x, a1x, b1x, b2x, b1y, b2y, b1z, b2z, e2, e1, F, tmpe0, tmpe1, tmpe2, tmpi1, tmpi2, Pxt, Pzt, Nx, Ny, Nz);
	MainProcedure(Cu, Vortex, 3.930, 8.92, 2.733e-2, 3.568, 1.0e+11, 1.5, 1.96, 4.73e+5, 5e-14, V2, a1y, a2y, a1z, a2z, a2x, a1x, b1x, b2x, b1y, b2y, b1z, b2z, e2, e1, F, tmpe0, tmpe1, tmpe2, tmpi1, tmpi2, Pxt, Pzt, Nx, Ny, Nz);
	MainProcedure(Cu, Vortex, 3.930, 8.92, 2.733e-2, 3.568, 1.0e+11, 1.5, 1.96, 4.73e+5, 1e-13, V2, a1y, a2y, a1z, a2z, a2x, a1x, b1x, b2x, b1y, b2y, b1z, b2z, e2, e1, F, tmpe0, tmpe1, tmpe2, tmpi1, tmpi2, Pxt, Pzt, Nx, Ny, Nz);
	MainProcedure(Cu, Vortex, 3.930, 8.92, 2.733e-2, 3.568, 1.0e+11, 1.5, 1.96, 4.73e+5, 2e-13, V2, a1y, a2y, a1z, a2z, a2x, a1x, b1x, b2x, b1y, b2y, b1z, b2z, e2, e1, F, tmpe0, tmpe1, tmpe2, tmpi1, tmpi2, Pxt, Pzt, Nx, Ny, Nz);	// продублировать по разным длительностям импульса
	*/
	system("pause");
	return 0;
}
