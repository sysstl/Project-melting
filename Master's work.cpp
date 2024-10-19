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
#include "Integral.h"
#include "Splayn.h"

using namespace std;

double pi = 3.141592654;

enum Metal { Au, Al, Cu, Ni };
enum TypeBeam { Gauss, Vortex };

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

string ConvertNumToString(int i)
{
	string tmp_num_line;

	int length = snprintf(NULL, 0, "%i", i);
	char* str = new char[length + 1];
	snprintf(str, length + 1, "%i", i);
	for (int j = 0; j < length; j++)
	{
		tmp_num_line.push_back(str[j]); // номер столбца
	}

	return tmp_num_line;
}

void MyCreateFile(int number_of_files_, vector<string> filename)
{
	vector<ofstream> file;
	file.reserve(number_of_files_);

	for (int i = 0; i < number_of_files_; i++)
	{
		file.emplace_back(ofstream{ filename[i] });
	}
}

void MyCreateFile_Animate(int number_of_plots, int count_frame_, int start_number_plot, int start_frame)
{
	vector<ofstream> file;
	file.reserve(number_of_plots * count_frame_);
	for (int i = start_number_plot; i < number_of_plots + start_number_plot; ++i)// цикл для создания системных файлов и файлов для данных для Анимации
	{
		string type_plot_Gif;
		string file_name;
		int length = snprintf(NULL, 0, "%i", i);
		char* str = new char[length + 1];
		snprintf(str, length + 1, "%i", i);
		for (int j = 0; j < length; j++)
		{
			type_plot_Gif.push_back(str[j]);
		}

		for (int j = start_frame; j < count_frame_ + start_frame; ++j)
		{
			string number_frame_Gif /*= "GIF"*/;

			int length = snprintf(NULL, 0, "%i", j);
			char* str = new char[length + 1];
			snprintf(str, length + 1, "%i", j);
			for (int k = 0; k < length; k++)
			{
				number_frame_Gif.push_back(str[k]);
			}
											// z_str		// mom_time
			//file_name = "e(xy) z(index) = " + type_plot_Gif + " + " + number_frame_Gif + ".txt";
			file_name = "Data_Gif" + type_plot_Gif + " + " + number_frame_Gif + ".txt";
			file.emplace_back(ofstream{ file_name });
		}
	}
}

void Null_Array(double** Array, int Nx, int Ny)
{
	for (int i = 0; i < Nx; i++)
		for (int j = 0; j < Ny; j++)
			Array[i][j] = 1e-16;
}

void SetDataArray(double** Array, int Nx, int Ny, double a)
{
	for (int i = 0; i < Nx; i++)
		for (int j = 0; j < Ny; j++)
			Array[i][j] = a;
}

void SelectDataForPlot(vector<string> & Initial_Data, vector<string> & Final_Data, double x_left_boundary, double x_right_boundary)
{
	// Initial_Data - вектор названий файлов с данными
	vector<string>::iterator PP_it, PP_it_res;
	int i = 1;
	for (PP_it = Initial_Data.begin(); PP_it != Initial_Data.end(); PP_it++)
	{
		string text;
		text = ConvertNumToString(i);
		string file_name = "AData" + text + ".txt";
		Final_Data.push_back(file_name);
		text.clear();
		i++;
	}

	MyCreateFile(Initial_Data.size(), Final_Data);
	PP_it_res = Final_Data.begin();

	for (PP_it = Initial_Data.begin(); PP_it != Initial_Data.end(); PP_it++)
	{
		ifstream in((*PP_it)); // окрываем файл для чтения

		VecD X;
		VecD Y;
		VecD::iterator it_x, it_y;

		if (in.is_open())
		{
			double x, y;
			while (in >> x >> y)
			{
				/*	new_points.push_back(Point{ x, y });*/
				X.push_back(x);
				Y.push_back(y);
				//cout << x << "   " << y << endl;
			}
		}
		in.close();

		ofstream fout((*PP_it_res));
		it_y = Y.begin();
		for (it_x = X.begin(); it_x != X.end(); it_x++)
		{
			if ((*it_x) >= x_left_boundary && (*it_x) <= x_right_boundary)
			{
				fout << (*it_x) << "   " << (*it_y) << endl;
			}
			it_y++;
		}

		PP_it_res++;
		X.clear();
		Y.clear();
		fout.close();
	}
}

struct P0V
{
	double p0v_s;
	double p0v_sl;
	double p0v_l;
};

struct Point3D // структура необходимая для определение точек, где произошел разрыв
{
	Point3D(int index_x, int index_y, int index_z) : index_x{ index_x }, index_y{ index_y }, index_z{ index_z } {}
	int index_x;
	int index_y;
	int index_z;
};

bool operator ==(const Point3D & p1, const Point3D & p2)
{
	return p1.index_x == p2.index_x && p1.index_y == p2.index_y && p1.index_z == p2.index_z;
}

bool comp_x(const Point3D & pt1, const Point3D & pt2)
{
	return (pt1.index_x < pt2.index_x);
}

bool comp_y(const Point3D & pt1, const Point3D & pt2)
{
	return (pt1.index_y < pt2.index_y);
}

bool comp_z(const Point3D & pt1, const Point3D & pt2)
{
	return (pt1.index_z < pt2.index_z);
}

void MySort_Point3D_y(vector<Point3D> & v)
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
		if (((*it).index_x == (*(it + 1)).index_x) && (it != copy_v.end()))
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

	Integral C_l_fun(My_function, a, b);

	return 9 * n_a[i] * k_b* pow(T_l / T_D[i], 3)* C_l_fun.GaussLegendre_IntervalVariety(5);
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

double GaussBeam(int& i, int& j, int& k, int& Nx, int& Ny, double& dx, double& dy, double& dz, double& dt, int& n, double& beta)
{
	return exp(-dx * dx * (i - Nx / 2) * (i - Nx / 2)) * exp(-dy * dy * (j - Ny / 2) * (j - Ny / 2)) * exp(-dz * k) * dt * n * exp(-beta * dt * n);
}

double VortexBeam(int& i, int& j, int& k, int& Nx, int& Ny, int& Nz, double& dx, double& dy, double& dz, double& dt, int& n, double& beta, int m)
{
	double h11 = 0, h12 = 0, hf1 = 0;

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

	return (pow(sqrt(pow(h11 * h11 + h12 * h12, abs(m))) * cos(m * hf1) * (exp(-dx * dx * (i - Nx / 2) * (i - Nx / 2)) * exp(-dy * dy * (j - Ny / 2) * (j - Ny / 2))), 2) +
		pow(sqrt(pow(h11 * h11 + h12 * h12, abs(m))) * sin(m * hf1) * (exp(-dx * dx * (i - Nx / 2) * (i - Nx / 2)) * exp(-dy * dy * (j - Ny / 2) * (j - Ny / 2))), 2)) * exp(-dz * k) * dt * n * exp(-beta * dt * n);
}

void Calculation00(Metal mt, Parametrs param, Splayn spl_C_e_on_T, Splayn spl_G_e_on_T, TypeBeam tbeam, double*** tmpe0, double*** tmpe1, double*** tmpi0, double*** tmpi1, int& Nx, int& Ny, int& Nz, double& dx, double& dy, double& dz, double& dt, double& A2, int& n, double& beta)
{// Calculation of heat conduction by explicit schem

	double A111 = 0, B111 = 0, A1i111 = 0, CC1111 = 0, CC2111 = 0, FFF = 0;

	for (int i = 1; i < Nx - 1; i++)
	{
		for (int j = 1; j < Ny - 1; j++)
		{
			for (int k = 1; k < Nz - 1; k++)
			{
				A111 = (Dependence_k_e_on_T(mt, param.T00 * tmpe0[i][j][k]) / (100.)) * param.t0 / ((spl_C_e_on_T.GetY(param.T00 * tmpe0[i][j][k]) / (100. * 100. * 100.)) * param.r0 * param.r0);
				B111 = param.kabs * param.P0 * param.t0 * beta / ((spl_C_e_on_T.GetY(param.T00 * tmpe0[i][j][k]) / (100. * 100. * 100.)) * param.T00);
				A1i111 = (Dependence_k_e_on_T(mt, 300.) / (100. * 100.)) * param.t0 / ((Dependence_C_l_on_T(mt, param.T00 * tmpi0[i][j][k]) / (100. * 100. * 100.)) * param.r0 * param.r0);
				if (tbeam == Gauss)
				{
					FFF = GaussBeam(i, j, k, Nx, Ny, dx, dy, dz, dt, n, beta);
				}

				if (tbeam == Vortex)
				{
					FFF = VortexBeam(i, j, k, Nx, Ny, Nz, dx, dy, dz, dt, n, beta, 1);
				}

				tmpe1[i][j][k] = tmpe0[i][j][k] +
					dt * A111 * ((tmpe0[i + 1][j][k] + tmpe0[i - 1][j][k] - 2 * tmpe0[i][j][k]) / (dx * dx) + (tmpe0[i][j + 1][k] + tmpe0[i][j - 1][k] - 2 * tmpe0[i][j][k]) / (dy * dy) + A2 * (tmpe0[i][j][k - 1] + tmpe0[i][j][k + 1] - 2 * tmpe0[i][j][k]) / (dz * dz))
					+ dt * B111 * FFF;

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

void fun2(int& Nx, int& Ny, int& Nz, double& dx, double& dy, double& dz, double& dt, double*** tmpe0, double*** tmpe1, double*** tmpe2, double*** tmpi0, double*** tmpi1, double*** tmpi2, double& A2,/*double*** A1, double*** A1i, double& A2, double*** B1, double*** CC1, double*** CC2,*/ /*double*** F,*/ /*double*** copy_tmpe1, double*** copy_tmpi1,*/ int& begin, int& end, Parametrs & param, vector<Melting> & Melt_metal, Metal mt, /*double*** Ci, double*** gamma,*/ Splayn & spl_G_e_on_T, Splayn spl_C_e_on_T, double beta, int n, Splayn spl_C_l_on_T, TypeBeam tbeam)
{
	double A111 = 0, B111 = 0, A1i111 = 0, CC1111 = 0, CC2111 = 0, FFF = 0;
	for (int k = begin; k < end; k++)
		//for (int i = 1; i < Nx - 1; i++)
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
				if (tbeam == Gauss)
				{
					FFF = GaussBeam(i, j, k, Nx, Ny, dx, dy, dz, dt, n, beta);
				}

				if (tbeam == Vortex)
				{
					FFF = VortexBeam(i, j, k, Nx, Ny, Nz, dx, dy, dz, dt, n, beta, 1);
				}

				tmpe2[i][j][k] = tmpe0[i][j][k] * (1 - 2 * dt * A111 / (dx * dx) - 2 * dt * A111 / (dy * dy) - 2 * dt * A111 * A2 / (dz * dz))
					+ 2 * dt * A111 * ((tmpe1[i + 1][j][k] + tmpe1[i - 1][j][k]) / (dx * dx) + (tmpe1[i][j + 1][k] + tmpe1[i][j - 1][k]) / (dy * dy) + A2 * (tmpe1[i][j][k - 1] + tmpe1[i][j][k + 1]) / (dz * dz))
					+ 2 * dt * B111 * FFF - 2 * dt * CC1111 * (tmpe1[i][j][k] - tmpi1[i][j][k]);

				tmpe2[i][j][k] = tmpe2[i][j][k] / (1 + 2 * dt * A111 / (dx * dx) + 2 * dt * A111 / (dy * dy) + 2 * dt * A111 * A2 / (dz * dz));

				double delta = 1.;
				if ((tmpi1[i][j][k] * param.T00) < (Melt_metal[mt].T_melting - delta))
				{

					tmpi2[i][j][k] = tmpi0[i][j][k] * (1 - 2 * dt * A1i111 / (dx * dx) - 2 * dt * A1i111 / (dy * dy) - 2 * dt * A1i111 * A2 / (dz * dz))
						+ 2 * dt * A1i111 * ((tmpi1[i + 1][j][k] + tmpi1[i - 1][j][k]) / (dx * dx) + (tmpi1[i][j + 1][k] + tmpi1[i][j - 1][k]) / (dy * dy) + A2 * (tmpi1[i][j][k - 1] + tmpi1[i][j][k + 1]) / (dz * dz))
						+ 2 * dt * CC2111 * (tmpe1[i][j][k] - tmpi1[i][j][k]);

					tmpi2[i][j][k] = tmpi2[i][j][k] / (1 + 2 * dt * A1i111 / (dx * dx) + 2 * dt * A1i111 / (dy * dy) + 2 * dt * A1i111 * A2 / (dz * dz));

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

					tmpi2[i][j][k] = tmpi0[i][j][k] * (1 - 2 * dt * A1i111 / (dx * dx) - 2 * dt * A1i111 / (dy * dy) - 2 * dt * A1i111 * A2 / (dz * dz))
						+ 2 * dt * A1i111 * ((tmpi1[i + 1][j][k] + tmpi1[i - 1][j][k]) / (dx * dx) + (tmpi1[i][j + 1][k] + tmpi1[i][j - 1][k]) / (dy * dy) + A2 * (tmpi1[i][j][k - 1] + tmpi1[i][j][k + 1]) / (dz * dz))
						+ 2 * dt * CC2111 * (tmpe1[i][j][k] - tmpi1[i][j][k]);

					tmpi2[i][j][k] = tmpi2[i][j][k] / (1 + 2 * dt * A1i111 / (dx * dx) + 2 * dt * A1i111 / (dy * dy) + 2 * dt * A1i111 * A2 / (dz * dz));

				}

				if ((tmpi1[i][j][k] * param.T00) > (Melt_metal[mt].T_melting + delta))
				{
					A1i111 = (Dependence_k_e_on_T(mt, 300.) / (100. * 100.)) * param.t0 / ((Dependence_C_l_on_T(mt, param.T00 * tmpi1[i][j][k]) / (100. * 100. * 100.)) * param.r0 * param.r0);
					CC2111 = (spl_G_e_on_T.GetY(param.T00 * tmpe1[i][j][k]) / (100. * 100. * 100.)) * param.t0 / ((Dependence_C_l_on_T(mt, param.T00 * tmpi1[i][j][k]) / (100. * 100. * 100.)));

					/*A1i111 = (Dependence_k_e_on_T(mt, 300.) / (100. * 100.)) * param.t0 / ((spl_C_l_on_T.GetY(param.T00 * tmpi1[i][j][k])) * param.r0 * param.r0);
					CC2111 = (spl_G_e_on_T.GetY(param.T00 * tmpe1[i][j][k]) / (100. * 100. * 100.)) * param.t0 / ((spl_C_l_on_T.GetY(param.T00 * tmpi1[i][j][k])));*/

					tmpi2[i][j][k] = tmpi0[i][j][k] * (1 - 2 * dt * A1i111 / (dx * dx) - 2 * dt * A1i111 / (dy * dy) - 2 * dt * A1i111 * A2 / (dz * dz))
						+ 2 * dt * A1i111 * ((tmpi1[i + 1][j][k] + tmpi1[i - 1][j][k]) / (dx * dx) + (tmpi1[i][j + 1][k] + tmpi1[i][j - 1][k]) / (dy * dy) + A2 * (tmpi1[i][j][k - 1] + tmpi1[i][j][k + 1]) / (dz * dz))
						+ 2 * dt * CC2111 * (tmpe1[i][j][k] - tmpi1[i][j][k]);

					tmpi2[i][j][k] = tmpi2[i][j][k] / (1 + 2 * dt * A1i111 / (dx * dx) + 2 * dt * A1i111 / (dy * dy) + 2 * dt * A1i111 * A2 / (dz * dz));
				}
			}
		}
	}
}

void fun21(int& Nx, int& Ny, int& Nz, double*** tmpe2, double*** tmpi1, double*** tmpi2, double& dz, Parametrs & param, Metal mt)
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

void Calculation0(Metal mt, Parametrs & param, Splayn & spl_C_e_on_T, Splayn & spl_G_e_on_T, TypeBeam tbeam, double*** tmpe0, double*** tmpe1, double*** tmpe2, double*** tmpi0, double*** tmpi1, double*** tmpi2, double& A2, double& dx, double& dy, double& dz, double& dt, int& Nx, int& Ny, int& Nz, int& n, double& beta, vector<Point3D> & points_rupture_heat, vector<Interval> & new_interv, vector<Melting> & Melt_metal, Splayn spl_C_l_on_T)
{

	std::thread t26(fun2, ref(Nx), ref(Ny), ref(Nz), ref(dx), ref(dy), ref(dz), ref(dt), tmpe0, tmpe1, tmpe2, tmpi0, tmpi1, tmpi2, ref(A2), ref(new_interv[0].begin), ref(new_interv[0].end), ref(param), ref(Melt_metal), mt, ref(spl_G_e_on_T), spl_C_e_on_T, beta, n, spl_C_l_on_T, tbeam);
	std::thread t27(fun2, ref(Nx), ref(Ny), ref(Nz), ref(dx), ref(dy), ref(dz), ref(dt), tmpe0, tmpe1, tmpe2, tmpi0, tmpi1, tmpi2, ref(A2), ref(new_interv[1].begin), ref(new_interv[1].end), ref(param), ref(Melt_metal), mt, ref(spl_G_e_on_T), spl_C_e_on_T, beta, n, spl_C_l_on_T, tbeam);
	std::thread t28(fun2, ref(Nx), ref(Ny), ref(Nz), ref(dx), ref(dy), ref(dz), ref(dt), tmpe0, tmpe1, tmpe2, tmpi0, tmpi1, tmpi2, ref(A2), ref(new_interv[2].begin), ref(new_interv[2].end), ref(param), ref(Melt_metal), mt, ref(spl_G_e_on_T), spl_C_e_on_T, beta, n, spl_C_l_on_T, tbeam);
	std::thread t29(fun2, ref(Nx), ref(Ny), ref(Nz), ref(dx), ref(dy), ref(dz), ref(dt), tmpe0, tmpe1, tmpe2, tmpi0, tmpi1, tmpi2, ref(A2), ref(new_interv[3].begin), ref(new_interv[3].end), ref(param), ref(Melt_metal), mt, ref(spl_G_e_on_T), spl_C_e_on_T, beta, n, spl_C_l_on_T, tbeam);
	std::thread t30(fun2, ref(Nx), ref(Ny), ref(Nz), ref(dx), ref(dy), ref(dz), ref(dt), tmpe0, tmpe1, tmpe2, tmpi0, tmpi1, tmpi2, ref(A2), ref(new_interv[4].begin), ref(new_interv[4].end), ref(param), ref(Melt_metal), mt, ref(spl_G_e_on_T), spl_C_e_on_T, beta, n, spl_C_l_on_T, tbeam);
	std::thread t31(fun2, ref(Nx), ref(Ny), ref(Nz), ref(dx), ref(dy), ref(dz), ref(dt), tmpe0, tmpe1, tmpe2, tmpi0, tmpi1, tmpi2, ref(A2), ref(new_interv[5].begin), ref(new_interv[5].end), ref(param), ref(Melt_metal), mt, ref(spl_G_e_on_T), spl_C_e_on_T, beta, n, spl_C_l_on_T, tbeam);
	std::thread t32(fun2, ref(Nx), ref(Ny), ref(Nz), ref(dx), ref(dy), ref(dz), ref(dt), tmpe0, tmpe1, tmpe2, tmpi0, tmpi1, tmpi2, ref(A2), ref(new_interv[6].begin), ref(new_interv[6].end), ref(param), ref(Melt_metal), mt, ref(spl_G_e_on_T), spl_C_e_on_T, beta, n, spl_C_l_on_T, tbeam);
	std::thread t33(fun2, ref(Nx), ref(Ny), ref(Nz), ref(dx), ref(dy), ref(dz), ref(dt), tmpe0, tmpe1, tmpe2, tmpi0, tmpi1, tmpi2, ref(A2), ref(new_interv[7].begin), ref(new_interv[7].end), ref(param), ref(Melt_metal), mt, ref(spl_G_e_on_T), spl_C_e_on_T, beta, n, spl_C_l_on_T, tbeam);
	std::thread t34(fun2, ref(Nx), ref(Ny), ref(Nz), ref(dx), ref(dy), ref(dz), ref(dt), tmpe0, tmpe1, tmpe2, tmpi0, tmpi1, tmpi2, ref(A2), ref(new_interv[8].begin), ref(new_interv[8].end), ref(param), ref(Melt_metal), mt, ref(spl_G_e_on_T), spl_C_e_on_T, beta, n, spl_C_l_on_T, tbeam);
	std::thread t35(fun2, ref(Nx), ref(Ny), ref(Nz), ref(dx), ref(dy), ref(dz), ref(dt), tmpe0, tmpe1, tmpe2, tmpi0, tmpi1, tmpi2, ref(A2), ref(new_interv[9].begin), ref(new_interv[9].end), ref(param), ref(Melt_metal), mt, ref(spl_G_e_on_T), spl_C_e_on_T, beta, n, spl_C_l_on_T, tbeam);

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

void QuasiGrid(double& r0, double& kabs, double& z0, double& xy0, vector <double>& Mesh_dz, vector <double>& Mesh_dz_dimension, vector <double>& Mesh_xy_dimension, vector <double>& Mesh_xy, int& Nz_ac, int Nxy_ac)
{
	//vector <double> Mesh_dz, Mesh_dz_dimension;
	vector <double> ::iterator it1;
	Mesh_dz.push_back(0.);
	double NN = Nz_ac;
	int total_grid_z = NN;
	int c = 8; //dfd
	double* s = new double[NN];
	for (int i = 0; i < total_grid_z; i++)
	{
		s[i] = i * 1 / NN;
	}

	cout << endl << endl;

	double* fs = new double[total_grid_z];
	for (int i = 0; i < total_grid_z; i++)
	{
		fs[i] = z0 * kabs * (exp(c * s[i]) - 1) / (exp(c) - 1);
		//cout << fs[i] << endl;
	}

	for (int i = 0; i < total_grid_z - 1; i++)
	{
		Mesh_dz.push_back(fs[i + 1] - fs[i]);
	}

	for (int i = 0; i < total_grid_z; i++)
	{
		fs[i] *= 1e+4 / kabs;
		Mesh_dz_dimension.push_back(fs[i]);
	}

	cout << endl << endl;

	cout << " Bezrazmern  shag " << endl;

	for (it1 = Mesh_dz.begin(); it1 != Mesh_dz.end(); it1++)
	{
		cout << *it1 << "   ";
	}
	cout << endl << endl;

	cout << " Razmern  setka " << endl;

	for (it1 = Mesh_dz_dimension.begin(); it1 != Mesh_dz_dimension.end(); it1++)
	{
		cout << *it1 << "   ";
	}
	cout << endl << endl;

	NN = Nxy_ac;
	c = 8; //dfd
	double* ss = new double[NN];
	for (int i = 0; i < NN; i++)
	{
		ss[i] = i * 1 / NN;
		//	cout << ss[i] << "   ";
	}

	cout << endl << endl;

	double* fss = new double[NN];
	for (int i = 0; i < NN; i++)
	{
		fss[i] = ((xy0 / r0) / 2) + (xy0 / r0 - ((xy0 / r0) / 2)) * (exp(c * ss[i]) - 1) / (exp(c) - 1);
		fss[i] *= 1e+4 * r0;
		//cout << fss[i] << "  ";
	}

	cout << endl << endl;

	double* fssf = new double[NN];
	for (int i = 0; i < NN; i++)
	{
		fssf[i] = ((xy0 / r0) / 2) + (0. - ((xy0 / r0) / 2)) * (exp(c * ss[i]) - 1) / (exp(c) - 1);
		//cout << fssf[i] << "  ";
		fssf[i] *= 1e+4 * r0;
	}

	cout << endl << endl;

	int total_grid_xy = 2 * NN;
	double* res_mesh_xy = new double[total_grid_xy];
	res_mesh_xy[0] = 0;
	int N_mesh_xy = NN; // 25
	cout << res_mesh_xy[0] << "   ";
	for (int i = 1; i < NN; i++)
	{
		res_mesh_xy[i] = fssf[N_mesh_xy - i];
		//cout << res_mesh[i] << "   ";
	}

	for (int i = NN; i < total_grid_xy; i++)
	{
		res_mesh_xy[i] = fss[i - N_mesh_xy];
		//cout << res_mesh[i] << "   ";
	}

	cout << endl << endl;

	//vector <double> Mesh_xy_dimension, Mesh_xy;
	Mesh_xy.push_back(0.);
	for (int i = 0; i < total_grid_xy; i++)
	{
		//	cout << res_mesh_xy[i] << "   "; // Razmern setka
		Mesh_xy_dimension.push_back(res_mesh_xy[i]);
	}

	cout << endl << endl;

	for (int i = 0; i < total_grid_xy - 1; i++)
	{
		Mesh_xy.push_back((res_mesh_xy[i + 1] - res_mesh_xy[i]) / (1e+4 * r0));
		//cout << res_mesh_xy[i+1] - res_mesh_xy[i] << "   ";
	}

	cout << endl << endl;
	cout << Mesh_xy.size() << endl;

	cout << " Bezrazmern  shag " << endl;

	for (it1 = Mesh_xy.begin(); it1 != Mesh_xy.end(); it1++)
	{
		cout << *it1 << "   ";
	}
	cout << endl << endl;

	cout << " Razmern  setka " << endl;

	for (it1 = Mesh_xy_dimension.begin(); it1 != Mesh_xy_dimension.end(); it1++)
	{
		cout << *it1 << "   ";
	}
	cout << endl << endl;
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

void Calculation1(double*** a1x, double*** a2x, double*** b1x, double*** b2x, double*** a1z, double*** b1y, double*** b2y, double*** b1z, double*** b2z, double*** a1y, double*** e1, double& dx, double& dy, double& dz, double& dt, double& CC0, int& Nx, int& Ny, int& Nz, vector<Interval> & new_interv_z)
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
				b2x[0][j][l] = b2x[0][j][l] + dt * CC0 * (e1[0][j][l] - e1[0][j - 1][l]) * ((a1y[1][j][l] - a1y[0][j][l]) * (a1z[0][j][l + 1] - a1z[0][j][l]) - (a1y[0][j][l + 1] - a1y[0][j][l]) * (a1z[0][j][l + 1] - a1z[0][j][l])) / (dx * dy * dz); // было a1z[0][j][l + 1] - a1z[0][j][l])
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
					b2x[i][j][l] = b2x[i][j][l] + dt * CC0 * (e1[i][j][l] - e1[i][j - 1][l]) * ((a1y[i + 1][j][l] - a1y[i][j][l]) * (a1z[i][j][l + 1] - a1z[i][j][l]) - (a1y[i][j][l + 1] - a1y[i][j][l]) * (a1z[i + 1][j][l] - a1z[i][j][l])) / (dx * dy * dz); //	a1z[i + 1][j][l] - a1z[i][j][l] ???????????? а было a1z[i][j][l + 1] - a1z[i][j][l])
					b2x[i][j][l] = b2x[i][j][l] + dt * CC0 * (e1[i][j][l] - e1[i][j][l - 1]) * ((a1y[i][j + 1][l] - a1y[i][j][l]) * (a1z[i + 1][j][l] - a1z[i][j][l]) - (a1y[i + 1][j][l] - a1y[i][j][l]) * (a1z[i][j + 1][l] - a1z[i][j][l])) / (dx * dy * dz); //a1y[i + 1][j][l] - a1y[i][j][l] ??????????? а было  a1y[i][j + 1][l] - a1y[i][j][l]
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

void Calculation10(double*** a1x, double*** a2x, double*** b1x, double*** b2x, double*** a1z, double*** b1y, double*** b2y, double*** b1z, double*** b2z, double*** a1y, double*** e1, double*** V1, double*** V2, double& dx, double& dy, double& dz, double& dt, double& CC0, int& Nx, int& Ny, int& Nz, vector<Interval> & new_interv_z)
{
	// Calculation 1
	// solving the equation of motion

	double a = 1.5; //к-т в искусственной вязкости
	double slag1 = 0;
	double slag2 = 0;
	double q = 0;

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

				// 1-е слагаемое

				if ((b1x[1][j][l] - b1x[0][j][l]) < 0)
				{
					slag1 = pow(b1x[1][j][l] - b1x[0][j][l], 2) / (V2[1][j][l] + V1[0][j][l]);
				}
				else
				{
					slag1 = 0;
				}

				/*if ((b1x[i][j][l] - b1x[i - 1][j][l]) < 0) //???
				{
					slag2 = pow(b1x[i][j][l] - b1x[i - 1][j][l], 2) / (V2[i - 1][j - 1][l - 1] + V1[i - 1][j - 1][l - 1]);
				}
				else
				{
					slag2 = 0;
				}*/

				q = 2 * pow(a, 2) * (slag1 - slag2);
				b2x[0][j][l] = b1x[0][j][l] + dt * CC0 * (e1[1][j][l] - e1[0][j][l]) * (-(a1y[0][j + 1][l] - a1y[0][j][l]) * (a1z[0][j][l + 1] - a1z[0][j][l]) + (a1y[0][j][l + 1] - a1y[0][j][l]) * (a1z[0][j + 1][l] - a1z[0][j][l])) / (dx * dy * dz);



				// 2-е слагаемое

				if ((b1x[0][j + 1][l] - b1x[0][j][l]) < 0)
				{
					slag1 = pow(b1x[0][j + 1][l] - b1x[0][j][l], 2) / (V2[0][j][l] + V1[0][j][l]);
				}
				else
				{
					slag1 = 0;
				}

				if ((b1x[0][j][l] - b1x[0][j - 1][l]) < 0)  // ??
				{
					slag2 = pow(b1x[0][j][l] - b1x[0][j - 1][l], 2) / (V2[0][j - 1][l - 1] + V1[0][j - 1][l - 1]);
				}
				else
				{
					slag2 = 0;
				}

				q = 2 * pow(a, 2) * (slag1 - slag2);
				b2x[0][j][l] = b2x[0][j][l] + dt * CC0 * (e1[0][j][l] - e1[0][j - 1][l]) * ((a1y[1][j][l] - a1y[0][j][l]) * (a1z[0][j][l + 1] - a1z[0][j][l]) - (a1y[0][j][l + 1] - a1y[0][j][l]) * (a1z[1][j][l] - a1z[0][j][l])) / (dx * dy * dz); // было a1z[0][j][l + 1] - a1z[0][j][l])



				// 3-е слагаемое

				if ((b1x[0][j][l + 1] - b1x[0][j][l]) < 0)
				{
					slag1 = pow(b1x[0][j][l + 1] - b1x[0][j][l], 2) / (V2[0][j][l] + V1[0][j][l]);
				}
				else
				{
					slag1 = 0;
				}

				if ((b1x[0][j][l] - b1x[0][j][l - 1]) < 0)//???
				{
					slag2 = pow(b1x[0][j][l] - b1x[0][j][l - 1], 2) / (V2[0][j - 1][l - 1] + V1[0][j - 1][l - 1]);
				}
				else
				{
					slag2 = 0;
				}

				q = 2 * pow(a, 2) * (slag1 - slag2);
				b2x[0][j][l] = b2x[0][j][l] + dt * CC0 * (e1[0][j][l] - e1[0][j][l - 1]) * ((a1y[0][j + 1][l] - a1y[0][j][l]) * (a1z[1][j][l] - a1z[0][j][l]) - (a1y[1][j][l] - a1y[0][j][l]) * (a1z[0][j + 1][l] - a1z[0][j][l])) / (dx * dy * dz); // было a1y[0][j+1][l] - a1y[0][j][l]

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
					// 1-е слагаемое

					if ((b1x[i + 1][j][l] - b1x[i][j][l]) < 0)
					{
						slag1 = pow(b1x[i + 1][j][l] - b1x[i][j][l], 2) / (V2[i][j][l] + V1[i][j][l]);
					}
					else
					{
						slag1 = 0;
					}

					if ((b1x[i][j][l] - b1x[i - 1][j][l]) < 0)
					{
						slag2 = pow(b1x[i][j][l] - b1x[i - 1][j][l], 2) / (V2[i - 1][j - 1][l - 1] + V1[i - 1][j - 1][l - 1]);
					}
					else
					{
						slag2 = 0;
					}

					q = 2 * pow(a, 2) * (slag1 - slag2);
					b2x[i][j][l] = b1x[i][j][l] + dt * CC0 * (e1[i][j][l] - e1[i - 1][j][l] + q) * (-(a1y[i][j + 1][l] - a1y[i][j][l]) * (a1z[i][j][l + 1] - a1z[i][j][l]) + (a1y[i][j][l + 1] - a1y[i][j][l]) * (a1z[i][j + 1][l] - a1z[i][j][l])) / (dx * dy * dz);


					// 2-е слагаемое

					if ((b1x[i][j + 1][l] - b1x[i][j][l]) < 0)
					{
						slag1 = pow(b1x[i][j + 1][l] - b1x[i][j][l], 2) / (V2[i][j][l] + V1[i][j][l]);
					}
					else
					{
						slag1 = 0;
					}

					if ((b1x[i][j][l] - b1x[i][j - 1][l]) < 0)
					{
						slag2 = pow(b1x[i][j][l] - b1x[i][j - 1][l], 2) / (V2[i - 1][j - 1][l - 1] + V1[i - 1][j - 1][l - 1]);
					}
					else
					{
						slag2 = 0;
					}

					q = 2 * pow(a, 2) * (slag1 - slag2);																																					//	a1z[i + 1][j][l] - a1z[i][j][l] ???????????? а было a1z[i][j][l + 1] - a1z[i][j][l])
					b2x[i][j][l] = b2x[i][j][l] + dt * CC0 * (e1[i][j][l] - e1[i][j - 1][l] + q) * ((a1y[i + 1][j][l] - a1y[i][j][l]) * (a1z[i][j][l + 1] - a1z[i][j][l]) - (a1y[i][j][l + 1] - a1y[i][j][l]) * (a1z[i + 1][j][l] - a1z[i][j][l])) / (dx * dy * dz);


					// 3-е слагаемое

					if ((b1x[i][j][l + 1] - b1x[i][j][l]) < 0)
					{
						slag1 = pow(b1x[i][j][l + 1] - b1x[i][j][l], 2) / (V2[i][j][l] + V1[i][j][l]);
					}
					else
					{
						slag1 = 0;
					}

					if ((b1x[i][j][l] - b1x[i][j][l - 1]) < 0)
					{
						slag2 = pow(b1x[i][j][l] - b1x[i][j][l - 1], 2) / (V2[i - 1][j - 1][l - 1] + V1[i - 1][j - 1][l - 1]);
					}
					else
					{
						slag2 = 0;
					}

					q = 2 * pow(a, 2) * (slag1 - slag2);																														//a1y[i + 1][j][l] - a1y[i][j][l] ??????????? а было  a1y[i][j + 1][l] - a1y[i][j][l]
					b2x[i][j][l] = b2x[i][j][l] + dt * CC0 * (e1[i][j][l] - e1[i][j][l - 1] + q) * ((a1y[i][j + 1][l] - a1y[i][j][l]) * (a1z[i + 1][j][l] - a1z[i][j][l]) - (a1y[i + 1][j][l] - a1y[i][j][l]) * (a1z[i][j + 1][l] - a1z[i][j][l])) / (dx * dy * dz);

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

					// 1-е слагаемое

					if ((b1y[i + 1][j][l] - b1y[i][j][l]) < 0)
					{
						slag1 = pow(b1y[i + 1][j][l] - b1y[i][j][l], 2) / (V2[i][j][l] + V1[i][j][l]);
					}
					else
					{
						slag1 = 0;
					}

					if ((b1y[i][j][l] - b1y[i - 1][j][l]) < 0)
					{
						slag2 = pow(b1y[i][j][l] - b1y[i - 1][j][l], 2) / (V2[i - 1][j - 1][l - 1] + V1[i - 1][j - 1][l - 1]);
					}
					else
					{
						slag2 = 0;
					}

					q = 2 * pow(a, 2) * (slag1 - slag2);
					b2y[i][j][l] = b1y[i][j][l] + dt * CC0 * (e1[i][j][l] - e1[i - 1][j][l] + q) * ((a1x[i][j + 1][l] - a1x[i][j][l]) * (a1z[i][j][l + 1] - a1z[i][j][l]) - (a1x[i][j][l + 1] - a1x[i][j][l]) * (a1z[i][j + 1][l] - a1z[i][j][l])) / (dx * dy * dz);


					// 2-е слагаемое

					if ((b1y[i][j + 1][l] - b1y[i][j][l]) < 0)
					{
						slag1 = pow(b1y[i][j + 1][l] - b1y[i][j][l], 2) / (V2[i][j][l] + V1[i][j][l]);
					}
					else
					{
						slag1 = 0;
					}

					if ((b1y[i][j][l] - b1y[i][j - 1][l]) < 0)
					{
						slag2 = pow(b1y[i][j][l] - b1y[i][j - 1][l], 2) / (V2[i - 1][j - 1][l - 1] + V1[i - 1][j - 1][l - 1]);
					}
					else
					{
						slag2 = 0;
					}

					q = 2 * pow(a, 2) * (slag1 - slag2);
					b2y[i][j][l] = b2y[i][j][l] + dt * CC0 * (e1[i][j][l] - e1[i][j - 1][l] + q) * ((a1x[i][j][l + 1] - a1x[i][j][l]) * (a1z[i + 1][j][l] - a1z[i][j][l]) - (a1x[i + 1][j][l] - a1x[i][j][l]) * (a1z[i][j][l + 1] - a1z[i][j][l])) / (dx * dy * dz);


					// 3-е слагаемое

					if ((b1y[i][j][l + 1] - b1y[i][j][l]) < 0)
					{
						slag1 = pow(b1y[i][j][l + 1] - b1y[i][j][l], 2) / (V2[i][j][l] + V1[i][j][l]);
					}
					else
					{
						slag1 = 0;
					}

					if ((b1y[i][j][l] - b1y[i][j][l - 1]) < 0)
					{
						slag2 = pow(b1y[i][j][l] - b1y[i][j][l - 1], 2) / (V2[i - 1][j - 1][l - 1] + V1[i - 1][j - 1][l - 1]);
					}
					else
					{
						slag2 = 0;
					}

					q = 2 * pow(a, 2) * (slag1 - slag2);
					b2y[i][j][l] = b2y[i][j][l] + dt * CC0 * (e1[i][j][l] - e1[i][j][l - 1] + q) * ((a1x[i + 1][j][l] - a1x[i][j][l]) * (a1z[i][j + 1][l] - a1z[i][j][l]) - (a1x[i][j + 1][l] - a1x[i][j][l]) * (a1z[i + 1][j][l] - a1z[i][j][l])) / (dx * dy * dz);


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
					//double 

			/*		if ((b1z[i + 1][j][l] - b1z[i][j][l]) >= 0)
					{
						b2z[i][j][l] = b1z[i][j][l] + dt * (e1[i][j][l] - e1[i - 1][j][l]) * ((a1x[i][j][l + 1] - a1x[i][j][l]) * (a1y[i][j + 1][l] - a1y[i][j][l]) - (a1x[i][j + 1][l] - a1x[i][j][l]) * (a1y[i][j][l + 1] - a1y[i][j][l])) / (dx * dy * dz);
					}

					if ((b1z[i][j + 1][l] - b1z[i][j][l]) >= 0)
					{
						b2z[i][j][l] = b2z[i][j][l] + dt * (e1[i][j][l] - e1[i][j - 1][l]) * ((a1x[i + 1][j][l] - a1x[i][j][l]) * (a1y[i][j][l + 1] - a1y[i][j][l]) - (a1x[i][j][l + 1] - a1x[i][j][l]) * (a1y[i + 1][j][l] - a1y[i][j][l])) / (dx * dy * dz);
					}

					if ((b1z[i][j][l + 1] - b1z[i][j][l]) >= 0)
					{
						b2z[i][j][l] = b2z[i][j][l] + dt * (e1[i][j][l] - e1[i][j][l - 1]) * ((a1x[i][j + 1][l] - a1x[i][j][l]) * (a1y[i + 1][j][l] - a1y[i][j][l]) - (a1x[i + 1][j][l] - a1x[i][j][l]) * (a1y[i][j + 1][l] - a1y[i][j][l])) / (dx * dy * dz);

					}

					if ((b1z[i][j][l] - b1z[i - 1][j][l]) >= 0)
					{
						b2z[i][j][l] = b1z[i][j][l] + dt * (e1[i][j][l] - e1[i - 1][j][l]) * ((a1x[i][j][l + 1] - a1x[i][j][l]) * (a1y[i][j + 1][l] - a1y[i][j][l]) - (a1x[i][j + 1][l] - a1x[i][j][l]) * (a1y[i][j][l + 1] - a1y[i][j][l])) / (dx * dy * dz);
					}

					if ((b1z[i][j][l] - b1z[i][j - 1][l]) >= 0)
					{
						b2z[i][j][l] = b2z[i][j][l] + dt * (e1[i][j][l] - e1[i][j - 1][l]) * ((a1x[i + 1][j][l] - a1x[i][j][l]) * (a1y[i][j][l + 1] - a1y[i][j][l]) - (a1x[i][j][l + 1] - a1x[i][j][l]) * (a1y[i + 1][j][l] - a1y[i][j][l])) / (dx * dy * dz);
					}

					if ((b1z[i][j][l] - b1z[i][j][l - 1]) >= 0)
					{
						b2z[i][j][l] = b2z[i][j][l] + dt * (e1[i][j][l] - e1[i][j][l - 1]) * ((a1x[i][j + 1][l] - a1x[i][j][l]) * (a1y[i + 1][j][l] - a1y[i][j][l]) - (a1x[i + 1][j][l] - a1x[i][j][l]) * (a1y[i][j + 1][l] - a1y[i][j][l])) / (dx * dy * dz);

					}*/

					// 1-е слагаемое

					if ((b1z[i + 1][j][l] - b1z[i][j][l]) < 0)
					{
						slag1 = pow(b1z[i + 1][j][l] - b1z[i][j][l], 2) / (V2[i][j][l] + V1[i][j][l]);
					}
					else
					{
						slag1 = 0;
					}

					if ((b1z[i][j][l] - b1z[i - 1][j][l]) < 0)
					{
						slag2 = pow(b1z[i][j][l] - b1z[i - 1][j][l], 2) / (V2[i - 1][j - 1][l - 1] + V1[i - 1][j - 1][l - 1]);
					}
					else
					{
						slag2 = 0;
					}

					q = 2 * pow(a, 2) * (slag1 - slag2);
					b2z[i][j][l] = b1z[i][j][l] + dt * (e1[i][j][l] - e1[i - 1][j][l] + q) * ((a1x[i][j][l + 1] - a1x[i][j][l]) * (a1y[i][j + 1][l] - a1y[i][j][l]) - (a1x[i][j + 1][l] - a1x[i][j][l]) * (a1y[i][j][l + 1] - a1y[i][j][l])) / (dx * dy * dz);



					// 2-е слагаемое

					if ((b1z[i][j + 1][l] - b1z[i][j][l]) < 0)
					{
						slag1 = pow(b1z[i][j + 1][l] - b1z[i][j][l], 2) / (V2[i][j][l] + V1[i][j][l]);
					}
					else
					{
						slag1 = 0;
					}

					if ((b1z[i][j][l] - b1z[i][j - 1][l]) < 0)
					{
						slag2 = pow(b1z[i][j][l] - b1z[i][j - 1][l], 2) / (V2[i - 1][j - 1][l - 1] + V1[i - 1][j - 1][l - 1]);
					}
					else
					{
						slag2 = 0;
					}

					q = 2 * pow(a, 2) * (slag1 - slag2);
					b2z[i][j][l] = b2z[i][j][l] + dt * (e1[i][j][l] - e1[i][j - 1][l] + q) * ((a1x[i + 1][j][l] - a1x[i][j][l]) * (a1y[i][j][l + 1] - a1y[i][j][l]) - (a1x[i][j][l + 1] - a1x[i][j][l]) * (a1y[i + 1][j][l] - a1y[i][j][l])) / (dx * dy * dz);


					// 3-е слагаемое

					if ((b1z[i][j][l + 1] - b1z[i][j][l]) < 0)
					{
						slag1 = pow(b1z[i][j][l + 1] - b1z[i][j][l], 2) / (V2[i][j][l] + V1[i][j][l]);
					}
					else
					{
						slag1 = 0;
					}

					if ((b1z[i][j][l] - b1z[i][j][l - 1]) < 0)
					{
						slag2 = pow(b1z[i][j][l] - b1z[i][j][l - 1], 2) / (V2[i - 1][j - 1][l - 1] + V1[i - 1][j - 1][l - 1]);
					}
					else
					{
						slag2 = 0;
					}

					q = 2 * pow(a, 2) * (slag1 - slag2);
					b2z[i][j][l] = b2z[i][j][l] + dt * (e1[i][j][l] - e1[i][j][l - 1] + q) * ((a1x[i][j + 1][l] - a1x[i][j][l]) * (a1y[i + 1][j][l] - a1y[i][j][l]) - (a1x[i + 1][j][l] - a1x[i][j][l]) * (a1y[i][j + 1][l] - a1y[i][j][l])) / (dx * dy * dz);

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

void TransformGridFrom_AcToHeat(int& Nx_old, int& Ny_old, int& Nz_old, double& dx_old, double& dy_old, double& dz_old, int& Nx_new, int& Ny_new, int& Nz_new, double& dx_new, double& dy_new, double& dz_new, double*** Initial_Data, double*** Final_Data, double*** Te_old_xy_new_z, double*** Te_old_y_new_xz, Splayn spl)
{
	// инетрпол вдоль Y

	VecD X, Y;

	for (int i = 0; i < Nx_old; i++) // ac
	{
		for (int j = 0; j < Nz_old; j++) // ac
		{
			for (size_t k = 0; k < Ny_old; k++)
			{
				X.push_back(k * dy_old);
				Y.push_back(Initial_Data[i][k][j]);
			}

			spl.InterpolateFast1D(X, Y);

			for (size_t k = 0; k < Ny_new; k++)
			{
				Te_old_y_new_xz[i][k][j] = spl.GetY(k * dy_new);
			}
			X.clear();
			X.erase(X.begin(), X.end());
			Y.clear();
			Y.erase(Y.begin(), Y.end());
		}
	}

	// инетрпол вдоль X

	for (int i = 0; i < Ny_new; i++) // heat
	{
		for (int j = 0; j < Nz_old; j++) // ac
		{
			for (size_t k = 0; k < Nx_old; k++)
			{
				X.push_back(k * dx_old);
				Y.push_back(Te_old_y_new_xz[k][i][j]);
			}

			spl.InterpolateFast1D(X, Y);

			for (size_t k = 0; k < Nx_new; k++)
			{
				Te_old_xy_new_z[k][i][j] = spl.GetY(k * dx_new);
			}
			X.clear();
			X.erase(X.begin(), X.end());
			Y.clear();
			Y.erase(Y.begin(), Y.end());
		}
	}


	// инетрпол вдоль Z

	for (int i = 0; i < Nx_new; i++)
	{
		for (int j = 0; j < Ny_new; j++)
		{
			for (size_t k = 0; k < Nz_old; k++)
			{
				X.push_back(k * dz_old);
				Y.push_back(Te_old_xy_new_z[i][j][k]);
			}

			spl.InterpolateFast1D(X, Y);
			for (int k = 0; k < Nz_new; k++)
			{
				Final_Data[i][j][k] = spl.GetY(k * dz_new);
			}
			X.clear();
			X.erase(X.begin(), X.end());
			Y.clear();
			Y.erase(Y.begin(), Y.end());
		}
	}
}

void TransformGrid(int& Nx_old, int& Ny_old, int& Nz_old, double& dx_old, double& dy_old, double& dz_old, int& Nx_new, int& Ny_new, int& Nz_new, double& dx_new, double& dy_new, double& dz_new, double*** Initial_Data, double*** Final_Data, double*** Te_old_xy_new_z, double*** Te_old_y_new_xz, Splayn spl)
{
	// инетрпол вдоль Z
	VecD X, Y;

	for (size_t i = 0; i < Nx_old; i++)
	{
		for (size_t j = 0; j < Ny_old; j++)
		{
			for (size_t k = 0; k < Nz_old; k++)
			{
				X.push_back(k * dz_old);
				Y.push_back(Initial_Data[i][j][k]);
			}

			spl.InterpolateFast1D(X, Y);

			for (size_t k = 0; k < Nz_new; k++)
			{
				Te_old_xy_new_z[i][j][k] = spl.GetY(k * dz_new);
			}
			X.clear();
			X.erase(X.begin(), X.end());
			Y.clear();
			Y.erase(Y.begin(), Y.end());
		}
	}

	// инетрпол вдоль X

	for (int i = 0; i < Ny_old; i++)
	{
		for (int j = 0; j < Nz_new; j++)
		{
			for (size_t k = 0; k < Nx_old; k++)
			{
				X.push_back(k * dx_old);
				Y.push_back(Te_old_xy_new_z[k][i][j]);
			}

			spl.InterpolateFast1D(X, Y);

			for (int k = 0; k < Nx_new; k++)
			{
				Te_old_y_new_xz[k][i][j] = spl.GetY(k * dx_new);
			}
			X.clear();
			X.erase(X.begin(), X.end());
			Y.clear();
			Y.erase(Y.begin(), Y.end());
		}
	}

	// инетрпол вдоль Y

	for (int i = 0; i < Nx_new; i++)
	{
		for (int j = 0; j < Nz_new; j++)
		{
			for (size_t k = 0; k < Ny_old; k++)
			{
				X.push_back(k * dy_old);
				Y.push_back(Te_old_y_new_xz[i][k][j]);
			}

			spl.InterpolateFast1D(X, Y);
	
			for (int k = 0; k < Ny_new; k++)
			{
				Final_Data[i][k][j] = spl.GetY(k * dy_new);
			}
			X.clear();
			X.erase(X.begin(), X.end());
			Y.clear();
			Y.erase(Y.begin(), Y.end());
		}
	}

	//for (size_t i = 0; i < Nx_old; i++)
	//{
	//	for (size_t j = 0; j < Ny_old; j++)
	//	{
	//		spl.InterpolateFast(1, Initial_Data, i, j, 0, z);
	//		for (size_t k = 0; k < Nz_new; k++)
	//		{
	//			Te_old_xy_new_z[i][j][k] = spl.GetY(k * dz_new);
	//			if (i == Nx_old / 2 && j == Ny_old / 2)
	//			{
	//				spl.InterpolateFast(1, Initial_Data, i, j, 0, z);
	//			//	fout << dz_new * k  * 1e+4 / (5e+5) <<"    "<< 300. * spl.GetY(k * dz_new) << endl;
	//			}
	//		}
	//	}
	//}

	//// инетрпол вдоль X стар

	//for (int i = 0; i < Ny_old; i++)
	//{
	//	for (int j = 0; j < Nz_new; j++)
	//	{
	//		spl.InterpolateFast(1, Te_old_xy_new_z, 0, i, j, x);
	//		for (int k = 0; k < Nx_new; k++)
	//		{
	//			Te_old_y_new_xz[k][i][j] = spl.GetY(k * dx_new);
	//		}
	//	}
	//}


	//// инетрпол вдоль Y стар

	//for (int i = 0; i < Nx_new; i++)
	//{
	//	for (int j = 0; j < Nz_new; j++)
	//	{
	//		spl.InterpolateFast(1, Te_old_y_new_xz, i, 0, j, y);
	//		for (int k = 0; k < Ny_new; k++)
	//		{
	//			Final_Data[i][k][j] = spl.GetY(k * dy_new);
	//		}
	//	}
	//}
}

void Calculation2(double*** a1x, double*** a1y, double*** a1z, double*** a2x, double*** a2y, double*** a2z, double*** b2x, double*** b2y, double*** b2z, double*** V2, double*** e2, double*** tmpe2, double*** tmpi2, double& dx, double& dy, double& dz, double& dx_heat, double& dy_heat, double& dz_heat, double& dt, double& CC0, int& Nx, int& Ny, int& Nz, int& Nx_heat, int& Ny_heat, int& Nz_heat, double& V0, P0V& p0v, vector<Melting> & Melt_metal, Metal mt, Parametrs & param, Splayn & spl_C_e_on_T, Splayn & spl_Te, Splayn & spl_Ti, Splayn & spl_e, double*** e_heat, double*** Te_acoustic, double*** Ti_acoustic, double*** Tei_old_xy_new_z, double*** Tei_old_y_new_xz, int& current_count_frame, vector<Point3D>& points_rupture_ac, vector<Point3D>& points_rupture_heat, vector<Interval>& new_interv_z, Splayn spl_C_l_on_T)
{// Calculation of 
 // Calculatuion 2

	vector<Point3D>::iterator it_ac, it_heat;

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

	TransformGrid(Nx_heat, Ny_heat, Nz_heat, dx_heat, dy_heat, dz_heat, Nx, Ny, Nz, dx, dy, dz, tmpe2, Te_acoustic, Tei_old_xy_new_z, Tei_old_y_new_xz, spl_Te);
	TransformGrid(Nx_heat, Ny_heat, Nz_heat, dx_heat, dy_heat, dz_heat, Nx, Ny, Nz, dx, dy, dz, tmpi2, Ti_acoustic, Tei_old_xy_new_z, Tei_old_y_new_xz, spl_Ti);

	//TransformGrid(Nx_heat, Ny_heat, Nx, Ny, Nz, dx, dy, dz, tmpe2, Te_acoustic, Tei_old_xy_new_z, Tei_old_y_new_xz, spl_Te); // это оставить
	//TransformGrid(Nx_heat, Ny_heat, Nx, Ny, Nz, dx, dy, dz, tmpi2, Ti_acoustic, Tei_old_xy_new_z, Tei_old_y_new_xz, spl_Ti);

	double delta = 1.;

	// способы расчета дваления для случая когда еще не учитывали отрыв
	/*if ((Ti_acoustic[i][j][k] * param.T00) < (Melt_metal[mt].T_melting - delta) )
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
			(Te_acoustic[i][j][k] - 1.);// *//*/

		e2[i][j][k] = (1. - V2[i][j][k]) +
			(Melt_metal[mt].gi * param.T00 * (((((spl_C_l_on_T.GetY(param.T00 * Ti_acoustic[i][j][k])))))) / (pow((Melt_metal[mt].Density + Melt_metal[mt].DensityLiquid) / 2, 2) * pow((Melt_metal[mt].u0 + Melt_metal[mt].u0_Liquid) / 2, 2) * V0 * 1e-6)) *
			(Ti_acoustic[i][j][k] - 1.) / V2[i][j][k] + (Melt_metal[mt].ge * param.T00 * (spl_C_e_on_T.GetY(param.T00 * Te_acoustic[i][j][k]) / (100. * 100. * 100.)) / (((Melt_metal[mt].Density + Melt_metal[mt].DensityLiquid) / 2) * pow((Melt_metal[mt].u0 + Melt_metal[mt].u0_Liquid) / 2, 2) * 1e-6)) *
			(Te_acoustic[i][j][k] - 1.);
	}

	if ((Ti_acoustic[i][j][k] * param.T00) > (Melt_metal[mt].T_melting + delta))
	{
		e2[i][j][k] = (1. - V2[i][j][k]) +
			(Melt_metal[mt].gi * param.T00 * ((Dependence_C_l_on_T(mt, param.T00 * Ti_acoustic[i][j][k]) / (100. * 100. * 100.))) / (Melt_metal[mt].DensityLiquid * Melt_metal[mt].DensityLiquid * pow(Melt_metal[mt].u0_Liquid, 2) * V0 * 1e-6)) *
			(Ti_acoustic[i][j][k] - 1.) / V2[i][j][k] + (Melt_metal[mt].ge * param.T00 * (spl_C_e_on_T.GetY(param.T00 * Te_acoustic[i][j][k]) / (100. * 100. * 100.)) / (Melt_metal[mt].DensityLiquid * pow(Melt_metal[mt].u0_Liquid, 2) * 1e-6)) *
			(Te_acoustic[i][j][k] - 1.);
	}*/

	it_ac = points_rupture_ac.begin();
	ofstream fout_point("Proverka_sovpad.txt");
	for (int i = 0; i < Nx; i++) // Nz_heat, dz_heat - не нужны
	{
		for (int j = 0; j < Ny; j++)// циклы по узлам акустики
		{
			for (int k = 1; k < Nz; k++)//l=1
			{
				if (points_rupture_ac.empty()) // вообще ещё нигде не отровало
				{
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
								(Te_acoustic[i][j][k] - 1.);// */

						e2[i][j][k] = (1. - V2[i][j][k]) +
							(Melt_metal[mt].gi * param.T00 * (((((spl_C_l_on_T.GetY(param.T00 * Ti_acoustic[i][j][k])))))) / (pow((Melt_metal[mt].Density + Melt_metal[mt].DensityLiquid) / 2, 2) * pow((Melt_metal[mt].u0 + Melt_metal[mt].u0_Liquid) / 2, 2) * V0 * 1e-6)) *
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
				else // ! empty
				{
					if ((*it_ac).index_x == i && (*it_ac).index_y == j && (*it_ac).index_z == k)
					{
						fout_point << (*it_ac).index_x << " == " << i << "    " << (*it_ac).index_y << " == " << j << (*it_ac).index_z << " == " << k << endl;
						e2[i][j][k] = 1e-16;
						it_ac++;
					}
					else
					{
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
									(Te_acoustic[i][j][k] - 1.);// */

							e2[i][j][k] = (1. - V2[i][j][k]) +
								(Melt_metal[mt].gi * param.T00 * (((((spl_C_l_on_T.GetY(param.T00 * Ti_acoustic[i][j][k])))))) / (pow((Melt_metal[mt].Density + Melt_metal[mt].DensityLiquid) / 2, 2) * pow((Melt_metal[mt].u0 + Melt_metal[mt].u0_Liquid) / 2, 2) * V0 * 1e-6)) *
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
				/*
				//if ((Ti_acoustic[i][j][k] * param.T00) < (Melt_metal[mt].T_melting - delta) && points_rupture_ac.empty())
				//{
				//	e2[i][j][k] = (1. - V2[i][j][k]) +
				//		(Melt_metal[mt].gi * param.T00 * (Dependence_C_l_on_T(mt, param.T00 * Ti_acoustic[i][j][k]) / (100. * 100. * 100.)) / (Melt_metal[mt].Density * Melt_metal[mt].Density * pow(Melt_metal[mt].u0, 2) * V0 * 1e-6)) *
				//		(Ti_acoustic[i][j][k] - 1.) / V2[i][j][k] + (Melt_metal[mt].ge * param.T00 * (spl_C_e_on_T.GetY(param.T00 * Te_acoustic[i][j][k]) / (100. * 100. * 100.)) / (Melt_metal[mt].Density * pow(Melt_metal[mt].u0, 2) * 1e-6)) *
				//		(Te_acoustic[i][j][k] - 1.);
				//}

				//if ((Ti_acoustic[i][j][k] * param.T00) < (Melt_metal[mt].T_melting - delta) && (!points_rupture_ac.empty())) // точки в веткоре расположены в порядке следования цикла
				//{	// разрыв тв фазы
				//	if ((*it_ac).index_x == i && (*it_ac).index_y == j && (*it_ac).index_z == k)
				//	{
				//		fout_point << (*it_ac).index_x << " == " << i << "    " << (*it_ac).index_y << " == " << j << (*it_ac).index_z << " == " << k << endl;
				//		e2[i][j][k] = 1e-16;
				//		it_ac++;
				//	}
				//	else
				//	{
				//		e2[i][j][k] = (1. - V2[i][j][k]) +
				//			(Melt_metal[mt].gi * param.T00 * (Dependence_C_l_on_T(mt, param.T00 * Ti_acoustic[i][j][k]) / (100. * 100. * 100.)) / (Melt_metal[mt].Density * Melt_metal[mt].Density * pow(Melt_metal[mt].u0, 2) * V0 * 1e-6)) *
				//			(Ti_acoustic[i][j][k] - 1.) / V2[i][j][k] + (Melt_metal[mt].ge * param.T00 * (spl_C_e_on_T.GetY(param.T00 * Te_acoustic[i][j][k]) / (100. * 100. * 100.)) / (Melt_metal[mt].Density * pow(Melt_metal[mt].u0, 2) * 1e-6)) *
				//			(Te_acoustic[i][j][k] - 1.);
				//	}
				//}


				//if ((Ti_acoustic[i][j][k] * param.T00) >= (Melt_metal[mt].T_melting - delta) && (Ti_acoustic[i][j][k] * param.T00) <= (Melt_metal[mt].T_melting + delta))
				//{
				//	//Ci[J/cm3 K] = ((Dependence_C_l_on_T(mt, param.T00 * tmpi2[i][j][k]) / (100. * 100. * 100.) + (Melt_metal[mt].Q_fusion * Melt_metal[mt].Density) / (2 * delta)))
				//	// spl_C_l_on_T.GetY(param.T00 * Ti_acoustic[i][j][k]) 
				//	// spl_Ti.GetY(k * dz) = tmpi2

				//	/*e2[i][j][k] = (1. - V2[i][j][k]) +
				//		(Melt_metal[mt].gi * param.T00 * (((((Dependence_C_l_on_T(mt, param.T00* tmpi2[i][j][k]) / (100. * 100. * 100.) + (Melt_metal[mt].Q_fusion * Melt_metal[mt].Density) / (2 * delta)))))) / (pow((Melt_metal[mt].Density + Melt_metal[mt].DensityLiquid) / 2, 2) * pow((Melt_metal[mt].u0 + Melt_metal[mt].u0_Liquid) / 2, 2) * V0 * 1e-6)) *
				//		(Ti_acoustic[i][j][k] - 1.) / V2[i][j][k] + (Melt_metal[mt].ge * param.T00 * (spl_C_e_on_T.GetY(param.T00 * Te_acoustic[i][j][k]) / (100. * 100. * 100.)) / (((Melt_metal[mt].Density + Melt_metal[mt].DensityLiquid) / 2) * pow((Melt_metal[mt].u0 + Melt_metal[mt].u0_Liquid) / 2, 2) * 1e-6)) *
				//		(Te_acoustic[i][j][k] - 1.);

				//	e2[i][j][k] = (1. - V2[i][j][k]) +
				//		(Melt_metal[mt].gi * param.T00 * (((((spl_C_l_on_T.GetY(param.T00 * Ti_acoustic[i][j][k])))))) / (pow((Melt_metal[mt].Density + Melt_metal[mt].DensityLiquid) / 2, 2) * pow((Melt_metal[mt].u0 + Melt_metal[mt].u0_Liquid) / 2, 2) * V0 * 1e-6)) *
				//		(Ti_acoustic[i][j][k] - 1.) / V2[i][j][k] + (Melt_metal[mt].ge * param.T00 * (spl_C_e_on_T.GetY(param.T00 * Te_acoustic[i][j][k]) / (100. * 100. * 100.)) / (((Melt_metal[mt].Density + Melt_metal[mt].DensityLiquid) / 2) * pow((Melt_metal[mt].u0 + Melt_metal[mt].u0_Liquid) / 2, 2) * 1e-6)) *
				//		(Te_acoustic[i][j][k] - 1.);
				//}

				//if ((Ti_acoustic[i][j][k] * param.T00) > (Melt_metal[mt].T_melting + delta))
				//{
				//	e2[i][j][k] = (1. - V2[i][j][k]) +
				//		(Melt_metal[mt].gi * param.T00 * ((Dependence_C_l_on_T(mt, param.T00 * Ti_acoustic[i][j][k]) / (100. * 100. * 100.))) / (Melt_metal[mt].DensityLiquid * Melt_metal[mt].DensityLiquid * pow(Melt_metal[mt].u0_Liquid, 2) * V0 * 1e-6)) *
				//		(Ti_acoustic[i][j][k] - 1.) / V2[i][j][k] + (Melt_metal[mt].ge * param.T00 * (spl_C_e_on_T.GetY(param.T00 * Te_acoustic[i][j][k]) / (100. * 100. * 100.)) / (Melt_metal[mt].DensityLiquid * pow(Melt_metal[mt].u0_Liquid, 2) * 1e-6)) *
				//		(Te_acoustic[i][j][k] - 1.);
				}*/
			}
		}
	}

	for (int i = 0; i < Nx; i++)// циклы по узлам акустики (определяем, где произошел разрыв металла)
	{
		for (int j = 0; j < Ny; j++)
		{
			for (int k = 0; k < Nz; k++)
			{
				if ( (e2[i][j][k] * p0v.p0v_s) <= -1400)
				{
					e2[i][j][k] = 1e-16;
					Point3D tmp(i, j, k);
					points_rupture_ac.push_back(tmp);
				}
			}
		}
	}

	// Ниже процедур анужна для того, чтобы в сетке тепловой определить где не считать темпер-ру из-за отрыва 
	TransformGridFrom_AcToHeat(Nx, Ny, Nz, dx, dy, dz, Nx_heat, Ny_heat, Nz_heat, dx_heat, dy_heat, dz_heat, e2, e_heat, Tei_old_xy_new_z, Tei_old_y_new_xz, spl_e);
	///TransformGridFrom_AcToHeat(Nx, Nz, Nx_heat, Ny_heat, Nz_heat, dx_heat, dy_heat, dz_heat, e2, e_heat, Tei_old_xy_new_z, Tei_old_y_new_xz, spl_e);

	for (int i = 0; i < Nx_heat; i++)// циклы по узлам акустики (определяем, где произошел разрыв металла)
	{
		for (int j = 0; j < Ny_heat; j++)
		{
			for (int k = 0; k < Nz_heat; k++)
			{
				if ((e_heat[i][j][k] * p0v.p0v_s) <= -1400)
				{
					e_heat[i][j][k] = 1e-16;
					Point3D tmp(i, j, k);
					points_rupture_heat.push_back(tmp);
				}
			}
		}
	}

	if (!points_rupture_ac.empty() /*&& current_count_frame <= 1599*/)
	{
		string current_count_frame_str = ConvertNumToString(current_count_frame);
		ofstream file_animAc("PointsAc" + current_count_frame_str + ".txt");
		file_animAc << "Ac" << endl;
		for (it_ac = points_rupture_ac.begin(); it_ac != points_rupture_ac.end(); it_ac++)
		{
			file_animAc << "(" << (*it_ac).index_x << "," << (*it_ac).index_y << "," << (*it_ac).index_z << ")" << endl;
		}
	}

	if (!points_rupture_heat.empty() /*&& current_count_frame <= 1599*/)
	{
		string current_count_frame_str = ConvertNumToString(current_count_frame);
		ofstream file_animHt("PointsHt" + current_count_frame_str + ".txt");
		file_animHt << "Ht" << endl;
		for (it_heat = points_rupture_heat.begin(); it_heat != points_rupture_heat.end(); it_heat++)
		{
			file_animHt << "(" << (*it_heat).index_x << "," << (*it_heat).index_y << "," << (*it_heat).index_z << ")" << endl;
		}
	}

}

void MainProcedure(Metal mt, TypeBeam tbeam, double kte_kte, double ro0, double CeCe, double CiCi, double gammagamma, double g_e, double g_i, double u00, double tptp, double P00, double*** V00, double*** V1, double*** V2, double*** a1y, double*** a2y, double*** a1z, double*** a2z, double*** a2x, double*** a1x, double*** b1x, double*** b2x, double*** b1y, double*** b2y, double*** b1z, double*** b2z, double*** e2, double*** e1, double*** tmpe0, double*** tmpe1, double*** tmpe2, double*** tmpi0, double*** tmpi1, double*** tmpi2, int Nx_heat, int  Ny_heat, int  Nz_heat, int  Nx_acoustic, int  Ny_acoustic, int  Nz_acoustic, Splayn spl_C_l_on_T, double*** e_heat, double*** Te_acoustic, double*** Ti_acoustic, double*** Tei_old_xy_new_z, double*** Tei_old_y_new_xz, double** massive_melting_tmp1, double** massive_melting_tmp2,  string current_namefile)    // davlenie, gorizontal'naja skorost', vertikal'naja skorost', temperatura
{
	// parametrs of laser and area
	// 20 mkm = 0.002      200 mkm = 0.02 cm 
	double r0 = 1e-2;// 1e-2; // radius of light beam, cm
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
	double p0v_s = Melt_metal[mt].Density * Melt_metal[mt].u0 * Melt_metal[mt].u0 * (1e-6); // 
	cout << " p0v = " << p0v_s << endl;
	double  p0v_sl = ((Melt_metal[mt].Density + Melt_metal[mt].DensityLiquid) / 2) * pow((Melt_metal[mt].u0 + Melt_metal[mt].u0_Liquid) / 2, 2) * (1e-6);
	cout << " p0v 1 = " << p0v_sl << endl;
	double  p0v_l = Melt_metal[mt].DensityLiquid * Melt_metal[mt].u0_Liquid * Melt_metal[mt].u0_Liquid * (1e-6);
	cout << " p0v = " << p0v_l << endl;

	P0V p0v;
	p0v.p0v_s = p0v_s;
	p0v.p0v_sl = p0v_sl;
	p0v.p0v_l = p0v_l;

	double sigma = 1400; // bar Au - предел прочности
	bool melting_split = false; // индикатор того, что (произошел/не произошло) плавление и откол
	bool melting = false;
	bool split = false;

	// откуда вычисляются???????
	//double at1i = gi * T00 * Ci / p0v;// 1.11* 0.132 * 300 * 19.32 / p0v;  // parameter vo vtorom (teplom) slagaemom uravnenija sostojanija  dlja Au
	//double at1e = ge * T00 * Ce / p0v;//1.5 * 300 * 4.56 * 1.2 * 0.01 / p0v;   // parameter vo vtorom (teplom) slagaemom uravnenija sostojanija  dlja Au

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

	//double dx_acoustic_melting = 0.02 / (r0 * 100/*Nx_heat*/);//0.02 / (r0 * 20/*Nx_heat*/);
	//double dy_acoustic_melting = 0.02 / (r0 * 100/*Nx_heat*/);//0.02 / (r0 * 20/*Nx_heat*/);
	//double dz_acoustic_melting = 1e-5 * kabs / 400;//1e-5 * kabs / 40;// Nz_acoustic;

	cout << dx_heat << "  " << dy_heat << "  " << dz_heat << endl;
	cout << dx_acoustic << "  " << dy_acoustic << "  " << dz_acoustic << endl;
	cout << " Razmern shag x y (mkm) heat = " << 1e+4 * r0 * dx_heat << endl;
	cout << " Razmern shag z (mkm) heat = " << 1e+4 * dz_heat / kabs << endl;
	cout << " Razmern shag x y (mkm) acoustic = " << 1e+4 * r0 * dx_acoustic << endl;
	cout << " Razmern shag z (mkm) acoustic = " << 1e+4 * dz_acoustic / kabs;
	cout << endl << endl;
	cout << " Razmern shag x y (nm) heat = " << 1e+4 * r0 * dx_heat * 1000 << endl;
	cout << " Razmern shag z (nm) heat = " << 1e+4 * dz_heat / kabs * 1000 << endl;
	cout << " Razmern shag x y (nm) acoustic = " << 1e+4 * r0 * dx_acoustic * 1000 << endl;
	cout << " Razmern shag z (nm) acoustic = " << 1e+4 * dz_acoustic / kabs * 1000;
	cout << endl << endl;

	double dt = 1. / (t0 * 1e+15);// 2e-4;
	cout << " dt = " << dt << endl;
	int n = 1;
	double tt = 0;
	// 32 переменные double

	Splayn spl_G_e_on_T;// = Calculation_Interpolation("Ge_" + metall + "_new.txt");
	Splayn spl_C_e_on_T;// = Calculation_Interpolation("Ce_" + metall + "_new.txt");
	spl_C_e_on_T.Calculation_Interpolation("Ce_" + metall + "_new.txt");
	spl_G_e_on_T.Calculation_Interpolation("Ge_" + metall + "_new.txt");
	Splayn spl_Te, spl_Ti, spl_e;

	// это для того, чтобы сетку координат один раз задать для интерполяции
	VecD X, Y, Z;
	for (int i = 0; i < Nx_heat; i++)//dz_heat
	{
		X.push_back(i * dx_heat); // если будет неравномерная сетка, то нужно одпраить в методе способ заполнения диагонали
	}

	for (int i = 0; i < Ny_heat; i++)//dz_heat
	{
		Y.push_back(i * dy_heat); // если будет неравномерная сетка, то нужно одпраить в методе способ заполнения диагонали
	}

	for (int i = 0; i < Nz_heat; i++)//dz_heat
	{
		Z.push_back(i * dz_heat); // если будет неравномерная сетка, то нужно одпраить в методе способ заполнения диагонали
	}

	spl_Te.SetInitialData(X, Y, Z, 1);
	spl_Ti.SetInitialData(X, Y, Z, 1);
	Z.clear();
	X.clear();
	Y.clear();

	for (int i = 0; i < Nx_acoustic; i++)//dz_ac
	{
		X.push_back(i * dx_acoustic); // если будет неравномерная сетка, то нужно одпраить в методе способ заполнения диагонали
	}

	for (int i = 0; i < Ny_acoustic; i++)//dz_ac
	{
		Y.push_back(i * dy_acoustic); // если будет неравномерная сетка, то нужно одпраить в методе способ заполнения диагонали
	}

	for (int i = 0; i < Nz_acoustic; i++)//dz_ac
	{
		Z.push_back(i * dz_acoustic); // если будет неравномерная сетка, то нужно одпраить в методе способ заполнения диагонали
	}

	spl_e.SetInitialData(X, Y, Z, 1);
	Z.clear();
	X.clear();
	Y.clear();


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
			massive_melting_tmp2[i][j] = 1e-16;
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
				V1[i][j][k] = V0;
				V00[i][j][k] = V0;
				massive_melting_tmp1[k][i] = 1e-16;
				/*massive_melting_tmp2[i][j][k] = 1e-16;*/
			}
		}
	}
	cout << " hello" << endl;
	system("pause");

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

	vector<Point3D> points_rupture_ac, points_rupture_heat;
	int current_count_frame_begin = 0;

	Calculation00(mt, param, spl_C_e_on_T, spl_G_e_on_T, tbeam, tmpe0, tmpe1, tmpi0, tmpi1, Nx_heat, Ny_heat, Nz_heat, dx_heat, dy_heat, dz_heat, dt, A2, n, beta);
	Calculation1(a1x, a2x, b1x, b2x, a1z, b1y, b2y, b1z, b2z, a1y, e1, dx_acoustic, dy_acoustic, dz_acoustic, dt, CC0, Nx_acoustic, Ny_acoustic, Nz_acoustic, new_interv_acoustic_z);// вязкости еще не будет т.к. поле скоростей берутся из н/у 
	GRU1(e2, Nx_acoustic, Ny_acoustic);
	GRU2(b2x, b2y, b2z, Nx_acoustic, Ny_acoustic, Nz_acoustic);
	Calculation2(a1x, a1y, a1z, a2x, a2y, a2z, b2x, b2y, b2z, V2, e2, tmpe2, tmpi2, dx_acoustic, dy_acoustic, dz_acoustic, dx_heat, dy_heat, dz_heat, dt, CC0, Nx_acoustic, Ny_acoustic, Nz_acoustic, Nx_heat, Ny_heat, Nz_heat, V0, p0v, Melt_metal, mt, param, spl_C_e_on_T, spl_Te, spl_Ti, spl_e, e_heat, Te_acoustic, Ti_acoustic, Tei_old_xy_new_z, Tei_old_y_new_xz, current_count_frame_begin, points_rupture_ac, points_rupture_heat, new_interv_acoustic_z, spl_C_l_on_T);

	tt = n * dt * t0 * 1e+15; // vremja v fs

	int number_plots = 24; //24
	int count_of_lines = 9;
	int npxzt = 1;
	int Ti_P_x_profile = 1;
	int current_number_line_melting = 1;
	int current_time_frame = 800; // номер кадра 
	int total_time_frame = 25000;// 2000; // формальный апарметр теперь (уже)
	vector<double> moments_time; // хранит моменты времени
	vector<int> moments_fix_time_anim;
	//int count_time = 1;
	for (int i = 1; i <= 1800; i++) // сетка времени
	{
		tt = i * dt * t0 * 1e+15; // vremja v fs
		moments_time.push_back(tt);
		//list_namefile.push_back()
	}

	for (int i = current_time_frame; i <= 1600; i++ /*= 10*/) // моменты времени для отрисовки и анимации
	{// аналог n
		moments_fix_time_anim.push_back(i);
	}

	//MyCreateFile(list_namefile.size(), list_namefile);
	vector<double>::iterator it_time, it_last_time; // хранит моменты времени
	vector<int>::iterator it_fix_time_anim;
	it_last_time = moments_time.end() - 1;
	it_fix_time_anim = moments_fix_time_anim.begin();
	cout << endl << endl << " last time = " << (*it_last_time) << endl;

	cout << endl << endl << " first fix time = " << (*(moments_fix_time_anim.begin())) << endl;
	cout << endl << endl << " last fix time = " << (*(moments_fix_time_anim.end() - 1)) << endl;

	cout << endl << endl;

	//GnuPlot plt(number_plots); // объект хранит 5 плотиков (объект на общий случай  без записи данных)
	GnuPlot plt(number_plots, total_time_frame); // с записью в файлы
	system("pause");
	plt.SetParametrs2D(0, 10, 3, "Ti", "x,mkm", "Ti, (K)");
	plt.SetParametrs2D(1, 2, 3, "Te,Ti", "time,fs", "Te,Ti (K)");
	plt.SetParametrsOnPlotColor(2, "Te", "z,mkm", "x,mkm", (1e+4 * dz_heat * Nz_heat / kabs), (1e+4 * r0 * dx_heat * Nx_heat));
	plt.SetParametrsOnPlotColor(3, "Ti", "z,mkm", "x,mkm", (1e+4 * dz_heat * Nz_heat / kabs), (1e+4 * r0 * dx_heat * Nx_heat));
	plt.SetParametrsOnPlotColor(4, "P(x,N/2,z)", "z,mkm", "x,mkm", (1e+4 * dz_acoustic * Nz_acoustic / kabs), (1e+4 * r0 * dx_acoustic * Nx_acoustic));
	plt.SetParametrsOnPlotColor(5, "P(x,y,1)", "y,mkm", "x,mkm", (1e+4 * r0 * dy_acoustic * Ny_acoustic), (1e+4 * r0 * dx_acoustic * Nx_acoustic));

	plt.SetParametrsOnPlotColor(6, "Te", "y,mkm", "x,mkm", (1e+4 * r0 * dy_heat * Ny_heat), (1e+4 * r0 * dx_heat * Nx_heat));
	plt.SetParametrsOnPlotColor(7, "Ti", "y,mkm", "x,mkm", (1e+4 * r0 * dy_heat * Ny_heat), (1e+4 * r0 * dx_heat * Nx_heat));
	plt.SetParametrs2D(8, count_of_lines, 3, "P,bar", "x,mkm", "P,bar");
	plt.SetParametrs2D(9, count_of_lines, 3, "P, bar", "z,mkm", "P, bar");
	plt.SetParametrs2D(10, 1, 3, "Temperature in central zone", "time,fs", "Ti, K"); // копии 1 и 10 графиков
	plt.SetParametrs2D(11, 11, 3, "Temperature T(z)", "z,mkm", "Ti, K");
	plt.SetParametrs2D(12, 11, 3, "Temperature T(x)", "x,mkm", "Ti, K");

	/*plt.SetParametrsOnPlotColor(13, "Melting zone (xz)", "z,mkm", "x,mkm", (1e+4* dz_heat* Nz_heat / kabs), (1e+4* r0* dx_heat* Nx_heat));
	plt.SetParametrsOnPlotColor(14, "Melting zone (xy)", "y,mkm", "x,mkm", (1e+4* r0* dy_heat* Ny_heat), (1e+4* r0* dx_heat* Nx_heat));*/
	plt.SetParametrsOnPlotColor(13, "Melting zone (xz)", "z,mkm", "x,mkm", (1e+4 * dz_acoustic * Nz_acoustic / kabs), (1e+4 * r0 * dx_acoustic * Nx_acoustic));
	plt.SetParametrsOnPlotColor(14, "Melting zone (xy)", "y,mkm", "x,mkm", (1e+4 * r0 * dy_acoustic * Ny_acoustic), (1e+4 * r0 * dx_acoustic * Nx_acoustic));

	plt.SetParametrs2D(15, 2, 3, "Te,Ti", "x,mkm", "Te,Ti (K)");
	plt.SetParametrs2D(16, 2, 3, "Te,Ti", "z,mkm", "Te,Ti (K)");

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	plt.SetParametrs2D(17, 10, 3, "P", "x,mkm", "P (bar)");
	plt.SetParametrs2D(18, 1, 3, "P", "x,mkm", "P (bar)");
	plt.SetParametrs2D(19, 1, 3, "P", "z,mkm", "P (bar)");

	plt.SetParametrsOnPlotColor(20, "Preasure zone (xz)", "z,mkm", "x,mkm", (1e+4 * dz_acoustic * Nz_acoustic / kabs), (1e+4 * r0 * dx_acoustic * Nx_acoustic));
	plt.SetParametrsOnPlotColor(21, "Preasure zone (xy)", "y,mkm", "x,mkm", (1e+4 * r0 * dy_acoustic * Ny_acoustic), (1e+4 * r0 * dx_acoustic * Nx_acoustic));

	plt.SetParametrs2D(22, 10, 3, "P", "x,mkm", "P (bar)");
	plt.SetParametrs2D(23, 10, 3, "Ti", "x,mkm", "Ti, (K)");

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// возможо здесь нужно отобразить в начальный момент времени, но это не обзательно

	plt.SetGridOnPlot3D(8, Nx_acoustic, Ny_acoustic, Nz_acoustic, 1e+4 * dx_acoustic * r0, 1e+4 * dy_acoustic * r0, 1e+4 * dz_acoustic / kabs, count_of_lines, fun_on_x);
	plt.SetGridOnPlot3D(9, Nx_acoustic, Ny_acoustic, Nz_acoustic, 1e+4 * dx_acoustic * r0, 1e+4 * dy_acoustic * r0, 1e+4 * dz_acoustic / kabs, count_of_lines, fun_on_z);
	plt.SetGridOnPlot3D(11, Nx_heat, Ny_heat, Nz_heat, 1e+4 * dx_heat * r0, 1e+4 * dy_heat * r0, 1e+4 * dz_heat / kabs, 11, fun_on_z);
	plt.SetGridOnPlot3D(12, Nx_heat, Ny_heat, Nz_heat, 1e+4 * dx_heat * r0, 1e+4 * dy_heat * r0, 1e+4 * dz_heat / kabs, 11, fun_on_x);
	plt.SetGridOnPlot3D(0, Nx_heat, Ny_heat, Nz_heat, 1e+4 * dx_heat * r0, 1e+4 * dy_heat * r0, 1e+4 * dz_heat / kabs, 100, fun_on_x);
	plt.SetGridOnPlot3D(17, Nx_acoustic, Ny_acoustic, Nz_acoustic, 1e+4 * dx_acoustic * r0, 1e+4 * dy_acoustic * r0, 1e+4 * dz_acoustic / kabs, 100, fun_on_x);
	plt.SetGridOnPlot3D(22, Nx_heat, Ny_heat, Nz_heat, 1e+4 * dx_heat * r0, 1e+4 * dy_heat * r0, 1e+4 * dz_heat / kabs, 100, fun_on_x);
	plt.SetGridOnPlot3D(23, Nx_acoustic, Ny_acoustic, Nz_acoustic, 1e+4 * dx_acoustic * r0, 1e+4 * dy_acoustic * r0, 1e+4 * dz_acoustic / kabs, 100, fun_on_x);

	clock_t t;
	double sum_t = 0;
	t = clock();

	vector<double***> TeTi = { tmpe2 , tmpi2 };//????
	vector<double***> Ti = { tmpi2 };//????
	//vector<double***> P = { e1 };/////
	vector<double***> empty;
	vector<string> TeTiTfs = { "Te","Ti" };
	vector<string> Pfs = { "P" };
	vector<string> TiTfs = { "Ti" };
	vector<string> Legenda_melting_tmp;
	vector<string> Legenda;
	vector<string> Legenda_Ti_x;
	vector<string> Legenda_P_x;
	vector<string> FilesForWritting;  // имена файлов куда запишутся поля темепратур и т.д
	vector<string> FilesForReading;
	vector<int> for_close_and_open = { 0, 2, 3, 4, 5, 6, 7, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 /*, 24, 25, 27, 28, 29, 30, 31, 32, 33, 34 */ };
	vector<int> for_close_and_open_ac = { 0, 1, 2 };
	string namefile;
	string fileoflistnamefiles = "List.txt";
	string depth = "Depth.txt";
	vector<Point3D> points_rupture_for_plot_z;
	vector<Point3D> points_rupture_for_plot_x;
	vector<Point3D>::iterator it;
	vector<Point3D>::iterator it_z, it_x;
	//points_rupture_ac.clear();
	//points_rupture_ac.erase(points_rupture_ac.begin(), points_rupture_ac.end());

	vector<string> list_namefile = { "P1.txt", "P2.txt", "P3.txt", "P4.txt", "P5.txt" };//?????

	string depth2 = current_namefile;// current_namefile;//"Depth2.txt"; // I0 = 0.25e+10
	ofstream fout_depth2(depth2);

	vector<string> Pacteti = { "Pac.txt", "PTi.txt", "PTe.txt", "Pacx.txt", "PTix.txt", "PTex.txt" };
	MyCreateFile(Pacteti.size(), Pacteti);
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
	Pacteti.clear();
	Pacteti = { "Points.txt" };
	MyCreateFile(Pacteti.size(), Pacteti);
	ofstream fout_point("Points.txt");

	int count_depth = 1;
	bool change_time_step = true;
	cout << endl << endl;

	it_time = moments_time.begin();
	do // цикл по времени
	{
		if (points_rupture_ac.empty()){}
		else
		{
			My_unique(points_rupture_ac);
			sort(points_rupture_ac.begin(), points_rupture_ac.end(), comp_x);
			MySort_Point3D_y(points_rupture_ac);
			MySort_Point3D_z(points_rupture_ac);
		}

		t = clock();

		Calculation0(mt, param, spl_C_e_on_T, spl_G_e_on_T, tbeam, tmpe0, tmpe1, tmpe2, tmpi0, tmpi1, tmpi2, A2, dx_heat, dy_heat, dz_heat, dt, Nx_heat, Ny_heat, Nz_heat, n, beta, points_rupture_heat, new_interv, Melt_metal, spl_C_l_on_T);
		if (n == 1)
		{
			Calculation1(a1x, a2x, b1x, b2x, a1z, b1y, b2y, b1z, b2z, a1y, e1, dx_acoustic, dy_acoustic, dz_acoustic, dt, CC0, Nx_acoustic, Ny_acoustic, Nz_acoustic, new_interv_acoustic_z);
		}
		else
		{
			Calculation10(a1x, a2x, b1x, b2x, a1z, b1y, b2y, b1z, b2z, a1y, e1, V00, V1, dx_acoustic, dy_acoustic, dz_acoustic, dt, CC0, Nx_acoustic, Ny_acoustic, Nz_acoustic, new_interv_acoustic_z);
		}

		GRU1(e2, Nx_acoustic, Ny_acoustic);
		GRU2(b2x, b2y, b2z, Nx_acoustic, Ny_acoustic, Nz_acoustic);
		Calculation2(a1x, a1y, a1z, a2x, a2y, a2z, b2x, b2y, b2z, V2, e2, tmpe2, tmpi2, dx_acoustic, dy_acoustic, dz_acoustic, dx_heat, dy_heat, dz_heat, dt, CC0, Nx_acoustic, Ny_acoustic, Nz_acoustic, Nx_heat, Ny_heat, Nz_heat, V0, p0v, Melt_metal, mt, param, spl_C_e_on_T, spl_Te, spl_Ti, spl_e, e_heat, Te_acoustic, Ti_acoustic, Tei_old_xy_new_z, Tei_old_y_new_xz, current_time_frame, points_rupture_ac, points_rupture_heat, new_interv_acoustic_z, spl_C_l_on_T);

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

		plt.SetDataOnPlot3D(1, 0, NULL, T00, Nx_heat / 2, Ny_heat / 2, 1, 2, 0, tt, TeTi, fun_on_t);
		plt.SetDataOnPlot3D(10, 0, NULL, T00, Nx_heat / 2, Ny_heat / 2, 1, 1, 0, tt, Ti, fun_on_t);

		//if (melting)
		// Потом попробуем брать температуру с исходной сетки температуры (в мэйн поменяь сетку массивов)
		// сейчас т.к. стоит сетка температуры акустическая то мб увидим не из-за ни ее ли идет разваливание в зоне плавления?
		if ((Ti_acoustic[Nx_acoustic / 2][Ny_acoustic / 2][1] * param.T00) >= (Melt_metal[mt].T_melting + 1.)) // Зона плавления
		//if ((tmpi1[Nx / 2][Ny / 2][1] * param.T00) >= (Melt_metal[mt].T_melting + 1.)) // Зона плавления
		{
			for (int i = 0; i < Nx_acoustic; i++)
			{
				for (int k = 0; k < Nz_acoustic; k++)
				{
					///if ((tmpi1[i][Ny / 2][k] * param.T00) >= (Melt_metal[mt].T_melting - 1.) && (tmpi1[i][Ny / 2][k] * param.T00) <= (Melt_metal[mt].T_melting + 1.))
					//if ((tmpi1[i][Ny_acoustic / 2][k] * param.T00) >= (Melt_metal[mt].T_melting + 1.))
					if ((Ti_acoustic[i][Ny_acoustic / 2][k] * param.T00) >= (Melt_metal[mt].T_melting + 1.))
					{
						//massive_melting_tmp1[i][Ny_acoustic / 2][k] = tmpi2[i][Ny_acoustic / 2][k];
						massive_melting_tmp1[k][i] = Ti_acoustic[i][Ny_acoustic / 2][k];
					}
				}
			}

			for (int i = 0; i < Nx_acoustic; i++)
			{
				for (int j = 0; j < Ny_acoustic; j++)
				{
					//if ((tmpi1[i][j][1] * param.T00) >= (Melt_metal[mt].T_melting + 1.))
					if ((Ti_acoustic[i][j][1] * param.T00) >= (Melt_metal[mt].T_melting + 1.))
					{
						//massive_melting_tmp2[i][j][1] = tmpi2[i][j][1];
						massive_melting_tmp2[i][j] = Ti_acoustic[i][j][1];
					}
				}
			}
			// 1e+4 * dx_acoustic * r0, 1e+4 * dy_acoustic * r0, 1e+4 * dz_acoustic / kabs,
			plt.SetDataOnPlotColor2D(13, current_time_frame, Nz_acoustic, Nx_acoustic, 1e+4 * dz_acoustic / kabs, 1e+4 * dy_acoustic * r0, massive_melting_tmp1, param.T00);
			plt.SetDataOnPlotColor2D(14, current_time_frame, Nx_acoustic, Ny_acoustic, 1e+4 * dx_acoustic * r0, 1e+4 * dy_acoustic * r0, massive_melting_tmp2, param.T00);

			namefile.clear();
			namefile = "Melting zone (zx) " + metall + "," + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs";
			//plt.ShowDataOnPlotColor(13, namefile, true);
			Sleep(1000);

			namefile.clear();
			namefile = "Melting zone (yx) " + metall + "," + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs";
			//plt.ShowDataOnPlotColor(14, namefile, true);
			Sleep(1000);

			Null_Array(massive_melting_tmp1, Nz_acoustic, Nx_acoustic);
			Null_Array(massive_melting_tmp2, Nx_acoustic, Ny_acoustic);

			plt.Close_and_open_files_for_replot(for_close_and_open);
		}

		// Preasure zone
		if (!points_rupture_ac.empty())/*(e2[Nx_acoustic / 2][Ny_acoustic / 2][1] * p0v) <= -sigma*/
		{
			SetDataArray(massive_melting_tmp1, Nz_acoustic, Nx_acoustic, 2000. / p0v_s);
			SetDataArray(massive_melting_tmp2, Nx_acoustic, Ny_acoustic, 2000. / p0v_s);
			for (int i = 0; i < Nx_acoustic; i++)
			{
				for (int k = 0; k < Nz_acoustic; k++)
				{
					if ((e2[i][Ny_acoustic / 2][k] * p0v_s) == (1e-16 * p0v_s))
					{
						//massive_melting_tmp1[i][Ny_acoustic / 2][k] = tmpi1[i][Ny_acoustic / 2][k];
						massive_melting_tmp1[k][i] = e2[i][Ny_acoustic / 2][k];
					}
				}
			}

			for (int i = 0; i < Nx_acoustic; i++)
			{
				for (int j = 0; j < Ny_acoustic; j++)
				{
					//if ((tmpi1[i][j][1] * param.T00) >= (Melt_metal[mt].T_melting + 1.))
					if ((e2[i][j][1] * p0v_s) == (1e-16 * p0v_s))
					{
						//massive_melting_tmp2[i][j][1] = tmpi1[i][j][1];
						massive_melting_tmp2[i][j] = e2[i][j][1];
					}
				}
			}

			plt.SetDataOnPlotColor2D(20, current_time_frame, Nz_acoustic, Nx_acoustic, 1e+4 * dz_acoustic / kabs, 1e+4 * dy_acoustic * r0, massive_melting_tmp1, p0v_s);
			plt.SetDataOnPlotColor2D(21, current_time_frame, Nx_acoustic, Ny_acoustic, 1e+4 * dx_acoustic * r0, 1e+4 * dy_acoustic * r0, massive_melting_tmp2, p0v_s);
			plt.SetDataOnPlotColor3D(21, current_time_frame, Nx_acoustic, Ny_acoustic, Nz_acoustic, 1e+4* dx_acoustic* r0, 1e+4* dy_acoustic* r0, 1e+4* dz_acoustic / kabs, e2, p0v_s, 1, xy);

			namefile.clear();
			namefile = "Preasure zone (zx) " + metall + "," + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs";
			//plt.ShowDataOnPlotColor(20, namefile, true);
			Sleep(1000);

			namefile.clear();
			namefile = "Preasure zone (yx) " + metall + "," + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs";
			plt.ShowDataOnPlotColor(21, namefile, true);
			Sleep(1000);

			Null_Array(massive_melting_tmp1, Nz_acoustic, Nx_acoustic);
			Null_Array(massive_melting_tmp2, Nx_acoustic, Ny_acoustic);

			plt.Close_and_open_files_for_replot(for_close_and_open);
		}

		// Профили температур ионной решетки в разные моменты
		if ((tmpi1[Nx_heat / 2][Ny_heat / 2][1] * param.T00) >= (Melt_metal[mt].T_melting + 1.)) // Профили температур ионной решетки в разные моменты
		{
			if ((int(ceil(tt)) % 1000 == 0 || int(trunc(tt)) % 1000 == 0) && current_number_line_melting <= 10)
			{
				plt.SetDataOnPlot3D(12, 0, tmpi1, param.T00, NULL, Ny_heat / 2, 1, 12, current_number_line_melting, NULL, empty, fun_on_x);
				plt.SetDataOnPlot3D(11, 0, tmpi1, param.T00, Nx_heat / 2, Ny_heat / 2, NULL, 12, current_number_line_melting, NULL, empty, fun_on_z);
				current_number_line_melting++;
				Legenda_melting_tmp.push_back(ConvertNumToStringdouble(tt) + " fs");
			}
		}

		//!!!!!!!!!!

		if (n == (*it_fix_time_anim)) // блок отрисовки данных и запись данных в файлы
		{//!!!!!!!!!!!!!!!!!!!!!!!!!!	tp * 1e+15
			//plt.Close_and_open_files_for_replot(for_close_and_open);
			namefile.clear();
			namefile = "Te,Ti (t,fms) " + metall + "," + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs, 1";
			plt.ShowDataOnPlot2D(1, 2, TeTiTfs, namefile, true);
			namefile.clear();
			namefile = "Ti (t,fms) " + metall + "," + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs, 1";
			plt.ShowDataOnPlot2D(10, 1, TiTfs, namefile, true);

			plt.SetDataOnPlotColor3D(2, current_time_frame, Nx_heat, Ny_heat, Nz_heat, 1e+4 * dx_heat * r0, 1e+4 * dy_heat * r0, 1e+4 * dz_heat / kabs, tmpe2, param.T00, Ny_heat / 2, xz);
			plt.SetDataOnPlotColor3D(3, current_time_frame, Nx_heat, Ny_heat, Nz_heat, 1e+4 * dx_heat * r0, 1e+4 * dy_heat * r0, 1e+4 * dz_heat / kabs, tmpi2, param.T00, Ny_heat / 2, xz);
			plt.SetDataOnPlotColor3D(4, current_time_frame, Nx_acoustic, Ny_acoustic, Nz_acoustic, 1e+4 * dx_acoustic * r0, 1e+4 * dy_acoustic * r0, 1e+4 * dz_acoustic / kabs, e2, p0v_s, Ny_acoustic / 2, xz);
			plt.SetDataOnPlotColor3D(5, current_time_frame, Nx_acoustic, Ny_acoustic, Nz_acoustic, 1e+4 * dx_acoustic * r0, 1e+4 * dy_acoustic * r0, 1e+4 * dz_acoustic / kabs, e2, p0v_s, 1, xy);

			plt.SetDataOnPlotColor3D(6, current_time_frame, Nx_heat, Ny_heat, Nz_heat, 1e+4 * dx_heat * r0, 1e+4 * dy_heat * r0, 1e+4 * dz_heat / kabs, tmpe2, param.T00, 1, xy);
			plt.SetDataOnPlotColor3D(7, current_time_frame, Nx_heat, Ny_heat, Nz_heat, 1e+4 * dx_heat * r0, 1e+4 * dy_heat * r0, 1e+4 * dz_heat / kabs, tmpi2, param.T00, 1, xy);

			int z_point = 0;// 5;
			Legenda_Ti_x.clear();
			Legenda_P_x.clear();
			int total_count_z = 100;
			for (int Ti_P_x_profile = 1; Ti_P_x_profile <= total_count_z; Ti_P_x_profile++)// цикл по номеру столбца 
			{
				// Ny_heat / 2 = 50   Ny_ac / 2 = 100
				// если меняю общее кол-во линий, то меняем в методе plt.SetGridOnPlot3D
				plt.SetDataOnPlot3D(0, current_time_frame, tmpi2, param.T00, NULL, Ny_heat / 2, z_point, total_count_z, Ti_P_x_profile, NULL, empty, fun_on_x);
				plt.SetDataOnPlot3D(17, current_time_frame, e2, p0v_s, NULL, Ny_acoustic / 2, z_point, total_count_z, Ti_P_x_profile, NULL, empty, fun_on_x);
				plt.SetDataOnPlot3D(22, current_time_frame, e_heat, p0v_s, NULL, Ny_heat / 2, z_point, total_count_z, Ti_P_x_profile, NULL, empty, fun_on_x);
				plt.SetDataOnPlot3D(23, current_time_frame, Ti_acoustic, param.T00, NULL, Ny_acoustic / 2, z_point, total_count_z, Ti_P_x_profile, NULL, empty, fun_on_x);
				Legenda_Ti_x.push_back(ConvertNumToStringdouble(1e+4 * z_point * dz_heat / kabs * 1000) + " nm");
				Legenda_P_x.push_back(ConvertNumToStringdouble(1e+4 * z_point * dz_acoustic / kabs * 1000) + " nm");
				z_point++;// = 5;// номер узла по легенде глубин (по z)
			}

			namefile.clear();
			namefile = "Ti (x,mkm) Legenda z" + metall + "," + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs";
			//plt.ShowDataOnPlot2D(0, 10, Legenda_Ti_x, namefile, true);
			//Sleep(1000);
			namefile.clear();
			namefile = "P (x,mkm) Legenda z" + metall + "," + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs";
			//plt.ShowDataOnPlot2D(17, 10, Legenda_P_x, namefile, true);
			//Sleep(1000);

			if (npxzt <= count_of_lines)
			{
				plt.SetDataOnPlot3D(8, 0, e2, p0v_s, NULL, Ny_acoustic / 2, 1, count_of_lines, npxzt, NULL, empty, fun_on_x);
				plt.SetDataOnPlot3D(9, 0, e2, p0v_s, Nx_acoustic / 2, Ny_acoustic / 2, NULL, count_of_lines, npxzt, NULL, empty, fun_on_z);
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
			namefile = "file P (x,y,1) " + metall + "," + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs";
			plt.ShowDataOnPlotColor(5, namefile, true);
			Sleep(1000);

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

			plt.SetDataOnPlot3D(15, current_time_frame, tmpe2, param.T00, NULL, Ny_heat / 2, 1, 2, 1, NULL, empty, fun_on_x);
			plt.SetDataOnPlot3D(15, current_time_frame, tmpi2, param.T00, NULL, Ny_heat / 2, 1, 2, 2, NULL, empty, fun_on_x);
			plt.SetDataOnPlot3D(16, current_time_frame, tmpe2, param.T00, Nx_heat / 2, Ny_heat / 2, NULL, 2, 1, NULL, empty, fun_on_z);
			plt.SetDataOnPlot3D(16, current_time_frame, tmpi2, param.T00, Nx_heat / 2, Ny_heat / 2, NULL, 2, 2, NULL, empty, fun_on_z);

			plt.SetDataOnPlot3D(18, current_time_frame, e2, p0v_s, NULL, Ny_acoustic / 2, 1, 1, 1, NULL, empty, fun_on_x);
			plt.SetDataOnPlot3D(19, current_time_frame, e2, p0v_s, Nx_acoustic / 2, Ny_acoustic / 2, NULL, 1, 1, NULL, empty, fun_on_z);

			// ЗАПИСЬ ОТДЕЛЬНЫХ СЛАГАЕМЫХ УРАВНЕНИЯ СОСТОЯНИЯ МИ-ГРЮНАЙЗЕНА 
			///spl_Te.InterpolateFast(1, tmpe2, Nx_acoustic / 2, Ny_acoustic / 2);
		///	spl_Ti.InterpolateFast(1, tmpi2, Nx_acoustic / 2, Ny_acoustic / 2);
			// может значения Т перезаписать в матрицу под сетку акустики?
			// если есть изменения в акустике не забывать менять там где идет расчет

			if (!points_rupture_ac.empty())
			{
				for (it = points_rupture_ac.begin(); it != points_rupture_ac.end(); it++)// цикл по точкам где произошел разрыв
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
				if (points_rupture_ac.empty())
				{
					fout_Pac << 1e+4 * dz_acoustic * k / kabs << "   " << p0v_s * (1. - V2[Nx_acoustic / 2][Ny_acoustic / 2][k]) << endl;

					if ((Ti_acoustic[Nx_acoustic / 2][Ny_acoustic / 2][k] * param.T00) < (Melt_metal[mt].T_melting - 1.))
					{
						/*if ((spl_Ti.GetY(k * dz_acoustic) * param.T00) >= 1336.)
						{
							fout_PTi << 1e+4 * dz_acoustic * k / kabs << "   " << p0v * ((Melt_metal[mt].gi * param.T00 * (spl_C_l_on_T.GetY(param.T00 * spl_Ti.GetY(k * dz_acoustic))) / (Melt_metal[mt].Density * Melt_metal[mt].Density * pow(Melt_metal[mt].u0, 2) * V0 * 1e-6)) *
								(spl_Ti.GetY(k * dz_acoustic) - 1.) / V2[Nx_acoustic / 2][Ny_acoustic / 2][k]) << endl;
						}
						else*/
						//{
						fout_PTi << 1e+4 * dz_acoustic * k / kabs << "   " << p0v_s * ((Melt_metal[mt].gi * param.T00 * (Dependence_C_l_on_T(mt, param.T00 * Ti_acoustic[Nx_acoustic / 2][Ny_acoustic / 2][k]) / (100. * 100. * 100.)) / (Melt_metal[mt].Density * Melt_metal[mt].Density * pow(Melt_metal[mt].u0, 2) * V0 * 1e-6)) *
							(Ti_acoustic[Nx_acoustic / 2][Ny_acoustic / 2][k] - 1.) / V2[Nx_acoustic / 2][Ny_acoustic / 2][k]) << endl;
						//}

						fout_PTe << 1e+4 * dz_acoustic * k / kabs << "   " << p0v_s * ((Melt_metal[mt].ge * param.T00 * (spl_C_e_on_T.GetY(param.T00 * Te_acoustic[Nx_acoustic / 2][Ny_acoustic / 2][k]) / (100. * 100. * 100.)) / (Melt_metal[mt].Density * pow(Melt_metal[mt].u0, 2) * 1e-6)) *
							(Te_acoustic[Nx_acoustic / 2][Ny_acoustic / 2][k] - 1.)) << endl;
					}

					if ((Ti_acoustic[Nx_acoustic / 2][Ny_acoustic / 2][k] * param.T00) >= (Melt_metal[mt].T_melting - 1.) && (Ti_acoustic[Nx_acoustic / 2][Ny_acoustic / 2][k] * param.T00) <= (Melt_metal[mt].T_melting + 1.))
					{
						//if ((spl_Ti.GetY(k * dz_acoustic) * param.T00) >= 1336.)
						//{
						// новый способ плавления
						fout_PTi << 1e+4 * dz_acoustic * k / kabs << "   " << p0v_s * ((Melt_metal[mt].gi * param.T00 * ((( /*Dependence_C_l_on_T(mt, param.T00* spl_Ti.GetY(k* dz_acoustic)) / (100. * 100. * 100.) + (Melt_metal[mt].Q_fusion * Melt_metal[mt].Density) / (2 * 1.0)*/  spl_C_l_on_T.GetY(param.T00 * Ti_acoustic[Nx_acoustic / 2][Ny_acoustic / 2][k])))) / (pow((Melt_metal[mt].Density + Melt_metal[mt].DensityLiquid) / 2, 2) * pow((Melt_metal[mt].u0 + Melt_metal[mt].u0_Liquid) / 2, 2) * V0 * 1e-6)) *
							(Ti_acoustic[Nx_acoustic / 2][Ny_acoustic / 2][k] - 1.) / V2[Nx_acoustic / 2][Ny_acoustic / 2][k]) << endl;

						// старый способ плавления
							/*fout_PTi << 1e+4 * dz_acoustic * k / kabs << "   " << p0v * ((Melt_metal[mt].gi * param.T00 * (((Dependence_C_l_on_T(mt, param.T00 * spl_Ti.GetY(k * dz_acoustic)) / (100. * 100. * 100.) + (Melt_metal[mt].Q_fusion * Melt_metal[mt].Density) / (2 * 1.0)))) / (pow((Melt_metal[mt].Density + Melt_metal[mt].DensityLiquid) / 2, 2) * pow((Melt_metal[mt].u0 + Melt_metal[mt].u0_Liquid) / 2, 2) * V0 * 1e-6)) *
								(spl_Ti.GetY(k * dz_acoustic) - 1.) / V2[Nx_acoustic / 2][Ny_acoustic / 2][k]) << endl;*/
								//}

						fout_PTe << 1e+4 * dz_acoustic * k / kabs << "   " << p0v_s * ((Melt_metal[mt].ge * param.T00 * (spl_C_e_on_T.GetY(param.T00 * Te_acoustic[Nx_acoustic / 2][Ny_acoustic / 2][k]) / (100. * 100. * 100.)) / (((Melt_metal[mt].Density + Melt_metal[mt].DensityLiquid) / 2) * pow((Melt_metal[mt].u0 + Melt_metal[mt].u0_Liquid) / 2, 2) * 1e-6)) *
							(Te_acoustic[Nx_acoustic / 2][Ny_acoustic / 2][k] - 1.)) << endl;
					}

					if ((Ti_acoustic[Nx_acoustic / 2][Ny_acoustic / 2][k] * param.T00) > (Melt_metal[mt].T_melting + 1.))
					{
						fout_PTi << 1e+4 * dz_acoustic * k / kabs << "   " << p0v_s * ((Melt_metal[mt].gi * param.T00 * (Dependence_C_l_on_T(mt, param.T00 * Ti_acoustic[Nx_acoustic / 2][Ny_acoustic / 2][k]) / (100. * 100. * 100.)) / (Melt_metal[mt].DensityLiquid * Melt_metal[mt].DensityLiquid * pow(Melt_metal[mt].u0_Liquid, 2) * V0 * 1e-6)) *
							(Ti_acoustic[Nx_acoustic / 2][Ny_acoustic / 2][k] - 1.) / V2[Nx_acoustic / 2][Ny_acoustic / 2][k]) << endl;
						//}

						fout_PTe << 1e+4 * dz_acoustic * k / kabs << "   " << p0v_s * ((Melt_metal[mt].ge * param.T00 * (spl_C_e_on_T.GetY(param.T00 * Te_acoustic[Nx_acoustic / 2][Ny_acoustic / 2][k]) / (100. * 100. * 100.)) / (Melt_metal[mt].DensityLiquid * pow(Melt_metal[mt].u0_Liquid, 2) * 1e-6)) *
							(Te_acoustic[Nx_acoustic / 2][Ny_acoustic / 2][k] - 1.)) << endl;
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
				///spl_Te.InterpolateFast(1, tmpe2, i, Ny_acoustic / 2);
				////spl_Ti.InterpolateFast(1, tmpi2, i, Ny_acoustic / 2);
				//if (points_rupture.empty())
				{
					fout_Pacx << 1e+4 * r0 * dx_acoustic * i << "   " << p0v_s * (1. - V2[i][Ny_acoustic / 2][1]) << endl;

					if ((Ti_acoustic[i][Ny_acoustic / 2][1] * param.T00) < (Melt_metal[mt].T_melting - 1.))
					{
						/*if ((spl_Ti.GetY(1 * dz_acoustic) * param.T00) >= 1336.)
						{
							fout_PTix << 1e+4 * r0 * dx_acoustic * i << "   " << p0v * ((Melt_metal[mt].gi * param.T00 * (spl_C_l_on_T.GetY(spl_Ti.GetY(1 * dz_acoustic) * param.T00)) / (Melt_metal[mt].Density * Melt_metal[mt].Density * pow(Melt_metal[mt].u0, 2) * V0 * 1e-6)) *
								(spl_Ti.GetY(1 * dz_acoustic) - 1.) / V2[i][Ny_acoustic / 2][1]) << endl;
						}
						else
						{*/
						fout_PTix << 1e+4 * r0 * dx_acoustic * i << "   " << p0v_s * ((Melt_metal[mt].gi * param.T00 * (Dependence_C_l_on_T(mt, param.T00 * Ti_acoustic[i][Ny_acoustic / 2][1]) / (100. * 100. * 100.)) / (Melt_metal[mt].Density * Melt_metal[mt].Density * pow(Melt_metal[mt].u0, 2) * V0 * 1e-6)) *
							(Ti_acoustic[i][Ny_acoustic / 2][1] - 1.) / V2[i][Ny_acoustic / 2][1]) << endl;
						//}

						fout_PTex << 1e+4 * r0 * dx_acoustic * i << "   " << p0v_s * ((Melt_metal[mt].ge * param.T00 * (spl_C_e_on_T.GetY(param.T00 * Te_acoustic[i][Ny_acoustic / 2][1]) / (100. * 100. * 100.)) / (Melt_metal[mt].Density * pow(Melt_metal[mt].u0, 2) * 1e-6)) *
							(Te_acoustic[i][Ny_acoustic / 2][1] - 1.)) << endl;
					}

					if ((Ti_acoustic[i][Ny_acoustic / 2][1] * param.T00) >= (Melt_metal[mt].T_melting - 1.) && (Ti_acoustic[i][Ny_acoustic / 2][1] * param.T00) <= (Melt_metal[mt].T_melting + 1.))
					{
						// новый способ плавления
						fout_PTix << 1e+4 * r0 * dx_acoustic * i << "   " << p0v_s * ((Melt_metal[mt].gi * param.T00 * ((spl_C_l_on_T.GetY(Ti_acoustic[i][Ny_acoustic / 2][1] * param.T00))) / (pow((Melt_metal[mt].Density + Melt_metal[mt].DensityLiquid) / 2, 2) * pow((Melt_metal[mt].u0 + Melt_metal[mt].u0_Liquid) / 2, 2) * V0 * 1e-6)) *
							(Ti_acoustic[i][Ny_acoustic / 2][1] - 1.) / V2[i][Ny_acoustic / 2][1]) << endl;

						// старый способ плавления
						/*fout_PTix << 1e+4 * r0 * dx_acoustic * i << "   " << p0v * ((Melt_metal[mt].gi * param.T00 * (((Dependence_C_l_on_T(mt, param.T00 * spl_Ti.GetY(1 * dz_acoustic)) / (100. * 100. * 100.) + (Melt_metal[mt].Q_fusion * Melt_metal[mt].Density) / (2 * 1.0)))) / (pow((Melt_metal[mt].Density + Melt_metal[mt].DensityLiquid) / 2, 2) * pow((Melt_metal[mt].u0 + Melt_metal[mt].u0_Liquid) / 2, 2) * V0 * 1e-6)) *
							(spl_Ti.GetY(1 * dz_acoustic) - 1.) / V2[i][Ny_acoustic / 2][1]) << endl;*/

						fout_PTex << 1e+4 * r0 * dx_acoustic * i << "   " << p0v_s * ((Melt_metal[mt].ge * param.T00 * (spl_C_e_on_T.GetY(param.T00 * Te_acoustic[i][Ny_acoustic / 2][1]) / (100. * 100. * 100.)) / (((Melt_metal[mt].Density + Melt_metal[mt].DensityLiquid) / 2) * pow((Melt_metal[mt].u0 + Melt_metal[mt].u0_Liquid) / 2, 2) * 1e-6)) *
							(Te_acoustic[i][Ny_acoustic / 2][1] - 1.)) << endl;
					}

					if ((Ti_acoustic[i][Ny_acoustic / 2][1] * param.T00) > (Melt_metal[mt].T_melting + 1.))
					{
						fout_PTix << 1e+4 * r0 * dx_acoustic * i << "   " << p0v_s * ((Melt_metal[mt].gi * param.T00 * (Dependence_C_l_on_T(mt, param.T00 * Ti_acoustic[i][Ny_acoustic / 2][1]) / (100. * 100. * 100.)) / (Melt_metal[mt].DensityLiquid * Melt_metal[mt].DensityLiquid * pow(Melt_metal[mt].u0_Liquid, 2) * V0 * 1e-6)) *
							(Ti_acoustic[i][Ny_acoustic / 2][1] - 1.) / V2[i][Ny_acoustic / 2][1]) << endl;

						fout_PTex << 1e+4 * r0 * dx_acoustic * i << "   " << p0v_s * ((Melt_metal[mt].ge * param.T00 * (spl_C_e_on_T.GetY(param.T00 * Te_acoustic[i][Ny_acoustic / 2][1]) / (100. * 100. * 100.)) / (Melt_metal[mt].DensityLiquid * pow(Melt_metal[mt].u0_Liquid, 2) * 1e-6)) *
							(Te_acoustic[i][Ny_acoustic / 2][1] - 1.)) << endl;
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
			//pac.ShowDataOnPlot2D(0, 1, Pfs, namefile, true);
			Sleep(1000);
			namefile.clear();
			namefile = "PTi (z,mkm) " + metall + "," + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs";
			//pac.ShowDataOnPlot2D(1, 1, Pfs, namefile, true);
			Sleep(1000);
			namefile.clear();
			namefile = "PTe (z,mkm) " + metall + "," + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs";
			//pac.ShowDataOnPlot2D(2, 1, Pfs, namefile, true);
			Sleep(1000);
			namefile.clear();
			namefile = "Pacoustic (x,mkm) " + metall + "," + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs";
			//pac.ShowDataOnPlot2D(3, 1, Pfs, namefile, true);
			Sleep(1000);
			namefile.clear();
			namefile = "PTi (x,mkm) " + metall + "," + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs";
			//pac.ShowDataOnPlot2D(4, 1, Pfs, namefile, true);
			Sleep(1000);
			namefile.clear();
			namefile = "PTe (x,mkm) " + metall + "," + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs";
			//pac.ShowDataOnPlot2D(5, 1, Pfs, namefile, true);
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

			current_time_frame++;
			it_fix_time_anim++;
			plt.Close_and_open_files_for_replot(for_close_and_open);
		}
		// next step
		n++;
		cout << " step = " << n << endl;

		// Формирование файлов анимации
		if (current_time_frame == total_time_frame)
		{
			namefile.clear();

			namefile = "Animate Te,Ti (x,mkm) " + metall + "," + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs";
			plt.CreateGifOnPlot2D(15, 2, 3, total_time_frame, moments_fix_time_anim, "Te,Ti", "x,mkm", "Te,Ti (K)", TeTiTfs, namefile, true);
			Sleep(1000);
			namefile.clear();

			namefile = "Animate Te,Ti (z,mkm) " + metall + "," + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs";
			plt.CreateGifOnPlot2D(16, 2, 3, total_time_frame, moments_fix_time_anim, "Te,Ti", "z,mkm", "Te,Ti (K)", TeTiTfs, namefile, true);
			Sleep(1000);
			namefile.clear();

			namefile = "Animate P (x,mkm) " + metall + "," + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs";
			plt.CreateGifOnPlot2D(18, 1, 3, total_time_frame, moments_fix_time_anim, "P", "x,mkm", "P (bar)", Pfs, namefile, true);
			Sleep(1000);
			namefile.clear();

			namefile = "Animate P (z,mkm) " + metall + "," + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs";
			plt.CreateGifOnPlot2D(19, 1, 3, total_time_frame, moments_fix_time_anim, "P", "z,mkm", "P (bar)", Pfs, namefile, true);
			Sleep(1000);
			namefile.clear();




			namefile = "Animate Te(x, N div 2, z) " + metall + "," + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs";
			plt.CreateGifOnPlotColor(2, total_time_frame, moments_fix_time_anim, "Te", "z,mkm", "x,mkm", (1e+4 * dz_heat * Nz_heat / kabs), (1e+4 * r0 * dx_heat * Nx_heat), namefile);
			Sleep(1000);

			namefile.clear();
			namefile = "Animate Ti(x, N div 2, z) " + metall + "," + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs";
			plt.CreateGifOnPlotColor(3, total_time_frame, moments_fix_time_anim, "Ti", "z,mkm", "x,mkm", (1e+4 * dz_heat * Nz_heat / kabs), (1e+4 * r0 * dx_heat * Nx_heat), namefile);
			Sleep(1000);

			namefile.clear();
			namefile = "Animate P (x,N div 2,z) " + metall + "," + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs";
			plt.CreateGifOnPlotColor(4, total_time_frame, moments_fix_time_anim, "P(x,N/2,z)", "z,mkm", "x,mkm", (1e+4 * dz_acoustic * Nz_acoustic / kabs), (1e+4 * r0 * dx_acoustic * Nx_acoustic), namefile);
			Sleep(1000);

			namefile.clear();
			namefile = "Animate P (x,y,1) " + metall + "," + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs";
			plt.CreateGifOnPlotColor(5, total_time_frame, moments_fix_time_anim, "P(x,y,1)", "y,mkm", "x,mkm", (1e+4 * r0 * dy_acoustic * Ny_acoustic), (1e+4 * r0 * dx_acoustic * Nx_acoustic), namefile);
			Sleep(1000);

			namefile.clear();
			namefile = "Animate Te (x,y,1) " + metall + "," + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs";
			plt.CreateGifOnPlotColor(6, total_time_frame, moments_fix_time_anim, "Te", "y,mkm", "x,mkm", (1e+4 * r0 * dy_heat * Ny_heat), (1e+4 * r0 * dx_heat * Nx_heat), namefile);
			Sleep(1000);

			namefile.clear();
			namefile = "Animate Ti (x,Ny,1) " + metall + "," + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs";
			plt.CreateGifOnPlotColor(7, total_time_frame, moments_fix_time_anim, "Ti", "y,mkm", "x,mkm", (1e+4 * r0 * dy_heat * Ny_heat), (1e+4 * r0 * dx_heat * Nx_heat), namefile);
			Sleep(1000);
			namefile.clear();
			current_time_frame++;
		}

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

		V00 = V1;
		V1 = V2;

		e1 = e2;
		it_time++;

	/*	if (n == 5)
		{
			ofstream fout_Test("Te_test.txt");
			for (int i = 0; i < Nx_heat; i++)
			{
				for (int j = 0; j < Ny_heat; j++)
				{
					for (int k = 0; k < Nz_heat; k++)
					{
						fout_Test << tmpe2[i][j][k] * T00 << "   ";
						cout<< tmpe2[i][j][k] * T00 << "   ";
					}
				}
			}
			system("pause");
		}*/

		//if (n == 2500)
		//{
		//	system("pause");// ОСТАНОВИЛСЯ ЗДЕСЬ
		//}

		cout << endl;
	} while (tt <= (*it_last_time) /*8000*//*nn < 1000*/);

	t = clock() - t;

	cout << endl << endl;
	cout << "It took me %d clicks (%f seconds)." << endl <<
		(int)t << "   " << ((double)t) / CLOCKS_PER_SEC;
	cout << endl;

	// могут приодитьс (если вдруг моментов времени меньше заявленных линий на графике, и как результат мы не запишем данные из матрицы в файл )
	if (npxzt < count_of_lines)
	{
		npxzt = count_of_lines;
		plt.SetDataOnPlot3D(8, 0, e1, p0v_s, NULL, Ny_acoustic / 2, 1, count_of_lines, npxzt, NULL, empty, fun_on_x);
		plt.SetDataOnPlot3D(9, 0, e1, p0v_s, Nx_acoustic / 2, Ny_acoustic / 2, NULL, count_of_lines, npxzt, NULL, empty, fun_on_z);
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
	//plt.ShowDataOnPlot2D(8, count_of_lines, Legenda, namefile, true);
	namefile.clear();
	namefile = "file P(z), bar " + metall + ", " + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs";
	//plt.ShowDataOnPlot2D(9, count_of_lines, Legenda, namefile, true);
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

inline int random_ab(int a, int b, std::default_random_engine & g)
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

	//int dxy = 100;
	int Nx_heat = 100;// dxy; // 100 //100; // chislo uzlov po x // 100
	int Ny_heat = 100;// dxy;// 100// 100; // chislo uzlov po y //100
	int Nz_heat = 250;//250; //!!!!!!!!!!!!!!

	int Nx_acoustic = 200;// 200;//400; //dxy; // 100 //100; // chislo uzlov po x // 100
	int Ny_acoustic = 200;// 200;// 400;// dxy;// 100// 100; // chislo uzlov po y //100
	int Nz_acoustic = 400;// 400;// 400;// 10000; // !!!!!!
	// по 400 узлов пока ломануло (посмотреть новый конструткор гнуплот)

	double*** V00, *** V1, *** V2, *** a1y, *** a2y, *** a1z, *** a2z, *** a2x, *** a1x, *** b1x, *** b2x, *** b1y, *** b2y, *** b1z, *** b2z, *** e2, *** e1, *** e_heat, *** tmpe0, *** tmpe1, *** tmpe2, *** tmpi0, *** tmpi1, *** tmpi2;    // davlenie, gorizontal'naja skorost', vertikal'naja skorost', temperatura
	double*** Ti_acoustic, *** Te_acoustic, *** Tei_old_xy_new_z, *** Tei_old_y_new_xz; // массивы для переопрелеления сетки температуры под акустику
	double** massive_melting_tmp1, ** massive_melting_tmp2;

	// 21 массива

	V00 = new double** [Nx_acoustic];
	V1 = new double** [Nx_acoustic];
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
	e_heat = new double** [Nx_heat];
	tmpe0 = new double** [Nx_heat]; // temperature of electrons 
	tmpe1 = new double** [Nx_heat]; // temperature of electrons 
	tmpe2 = new double** [Nx_heat]; // temperature of electrons 
	tmpi0 = new double** [Nx_heat]; // temperature of ions
	tmpi1 = new double** [Nx_heat]; // temperature of ions
	tmpi2 = new double** [Nx_heat]; // temperature of ions
	Ti_acoustic = new double** [Nx_acoustic];
	Te_acoustic = new double** [Nx_acoustic];
	Tei_old_xy_new_z = new double** [Nx_heat];
	Tei_old_y_new_xz = new double** [Nx_acoustic];
	massive_melting_tmp1 = new double* [Nz_acoustic];
	massive_melting_tmp2 = new double* [Nx_acoustic];

	for (int i = 0; i < Nz_acoustic; i++)
	{
		massive_melting_tmp1[i] = new double[Nx_acoustic];
	}

	for (int i = 0; i < Nx_acoustic; i++)
	{
		V00[i] = new double* [Ny_acoustic];
		V1[i] = new double* [Ny_acoustic];
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
		massive_melting_tmp2[i] = new double [Ny_acoustic];
		for (int j = 0; j < Ny_acoustic; j++)
		{
			V00[i][j] = new double[Nz_acoustic];
			V1[i][j] = new double[Nz_acoustic];
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
		e_heat[i] = new double* [Ny_heat];
		for (int j = 0; j < Ny_heat; j++)
		{
			tmpe0[i][j] = new double[Nz_heat];
			tmpe1[i][j] = new double[Nz_heat];
			tmpe2[i][j] = new double[Nz_heat];
			tmpi0[i][j] = new double[Nz_heat];
			tmpi1[i][j] = new double[Nz_heat];
			tmpi2[i][j] = new double[Nz_heat];
			e_heat[i][j] = new double[Nz_heat];
		}
	}

	for (int i = 0; i < Nx_heat; i++)
	{
		Tei_old_xy_new_z[i] = new double* [Ny_heat];
		for (int j = 0; j < Ny_heat; j++)
		{
			Tei_old_xy_new_z[i][j] = new double[Nz_acoustic];
		}
	}

	for (int i = 0; i < Nx_acoustic; i++)
	{
		Tei_old_y_new_xz[i] = new double* [Ny_heat];
		for (int j = 0; j < Ny_heat; j++)
		{
			Tei_old_y_new_xz[i][j] = new double[Nz_acoustic];
		}
	}

	double TTT = 300.;
	/*  отрисовка графиков */

	Splayn spl_G_e_on_T, spl_C_e_on_T; //= Calculation_Interpolation("Ge_Au_new.txt");
	//spl_C_e_on_T.Calculation_Interpolation("Ce_Au_new.txt");
	spl_C_e_on_T.InterpolateFast1D("Ce_Au_new.txt");
	spl_G_e_on_T.Calculation_Interpolation("Ge_Au_new.txt");
	cout << "C_e_on_T = " << spl_C_e_on_T.GetY(TTT) / (100. * 100. * 100.) << "   " << 5.472e-2 << endl;
	cout << "G_e_on_T = " << spl_G_e_on_T.GetY(TTT) / (100. * 100. * 100.) << "   " << 2.5e+10 << endl;
	cout << " kte = " << Dependence_k_e_on_T(Au, TTT) / (100.) << "   " << 3.115 << endl;
	cout << " Cl = " << Dependence_C_l_on_T(Al, TTT) / (100. * 100. * 100.) << "   " << 2.550 << endl;

	ofstream fout1("Ce_interp.txt");
	ofstream fout2("Ge_interp.txt");
	ofstream fout3("ke_depend.txt");
	ofstream fout4("Cl_depend.txt");

	for (double i = 300.; i <= 50250.; i++)
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

	//system("pause");

	vector<string> vvv = { "100 fs.txt", "1 ps.txt", "3 ps.txt", "5 ps.txt", "10 ps.txt" };
	GnuPlot pltt(vvv);
	pltt.SetParametrs2D(0, 5, 3, "h_max(F)", "F, J/cm^2", "h_max, nm");
	vector<string> Melting1 = { "100 fs", "1 ps", "3 ps", "5 ps", "10 ps" };
	namefile = "Dependence h_max on F";
	pltt.ShowDataOnPlot2D(0, 0, 5, Melting1, namefile, true);
	namefile.clear();

	vector<Point3D> v;
	std::default_random_engine g(std::chrono::system_clock::now().time_since_epoch().count());
	for (int i = 0; i < 60; i++)
	{
		srand(static_cast<unsigned int>(time(nullptr)));
		int x = random_ab(47, 53, g);//47 + rand() % 53; // случайное число - формула
		srand(static_cast<unsigned int>(time(nullptr)));
		int y = random_ab(47, 53, g); ;// 47 + rand() % 53; // случайное число - формула
		srand(static_cast<unsigned int>(time(nullptr)));
		int z = random_ab(1, 6, g); ;// 1 + rand() % 6; // случайное число - формула

		Point3D tmp(x, y, z);
		//cout << tmp.index_x << "   " << tmp.index_y << "   " << tmp.index_z << endl;
		v.push_back(tmp);
	}

	if (v.empty())
	{
		cout << " empty " << endl;
	}
	else
	{
		cout << " not empty " << endl;
	}
	cout << v.empty() << endl;
	cout << endl << " size " << v.size() << endl;
	//Point3D tmp();
	vector<Point3D>::iterator iter;
	My_unique(v);
	cout << endl << v.size() << endl;

	cout << endl << endl << endl << endl;
	//system("pause");

	sort(v.begin(), v.end(), comp_x);
	//My_unique(v);
	for (iter = v.begin(); iter != v.end(); iter++)
	{
		//cout << (*iter).index_x << "   " << (*iter).index_y << "   " << (*iter).index_z << endl;
	}
	cout << endl << endl;
	cout << endl << v.size() << endl;

	//system("pause");

	MySort_Point3D_y(v);
	cout << endl;
	for (iter = v.begin(); iter != v.end(); iter++)
	{
		//cout << (*iter).index_x << "   " << (*iter).index_y << "   " << (*iter).index_z << endl;
	}
	cout << endl << v.size() << endl;

	//system("pause");

	MySort_Point3D_z(v);
	cout << endl;
	cout << endl << v.size() << endl;
	for (iter = v.begin(); iter != v.end(); iter++)
	{
		//cout << (*iter).index_x << "   " << (*iter).index_y << "   " << (*iter).index_z << endl;
	}
	cout << endl << v.size() << endl;

	cout << "the end " << endl;
	//system("pause");

	vector<string> Pacteti = { "Data1.txt", "Data2.txt", "Data3.txt", "Data4.txt", "Data5.txt", "Data6.txt" };
	ofstream file1(Pacteti[0]);
	ofstream file2(Pacteti[1]);
	ofstream file3(Pacteti[2]);
	ofstream file4(Pacteti[3]);
	ofstream file5(Pacteti[4]);

	for (int i = 0; i < 50; i++)
	{
		file1 << i * 0.5 << "   " << sin(i * 0.5) << "   " << cos(i * 0.5) << endl;
		file2 << i * 0.5 << "   " << 2 * cos(i * 0.5) << "   " << 2 * cos(i * 0.5) << endl;
		file3 << i * 0.5 << "   " << 2 * sin(i * 0.5) + cos(i * 0.5) << "   " << 5 * cos(i * 0.5) << endl;
		file4 << i * 0.5 << "   " << 2 * sin(i * 0.5) * cos(i * 0.5) << "   " << 10 * cos(i * 0.5) << endl;
		file5 << i * 0.5 << "   " << 2 * sin(i * 0.5) << "   " << 2 * cos(i * 0.5) << endl;
	}

	//vector<string> list_namefile = { "P1.txt", "P2.txt", "P3.txt", "P4.txt", "P5.txt" };
	//MyCreateFile(list_namefile.size(), list_namefile);

	//vector<string> moments_time_string; // хранит моменты времени
	//map<string, string> momtime_n_t;// задает соответствие времени в фс с номером кадра в аниации
	//map<string, string>::iterator it_map;
	//it_map = momtime_n_t.begin();
	//int count_time_frame = 1;
	//for (double i = 100; i <= 100000; i += 100)
	//{
	//	//moments_time_string.push_back(ConvertNumToStringdouble(i) + "fs");
	//	string tmp_first = ConvertNumToStringdouble(i) + "fs";
	//	//(*it_map).second = tmp_second + "fs";
	//	string tmp_second = ConvertNumToString(count_time_frame);
	//	///cout << tmp_first << "   " << tmp_second << endl;
	//	//momtime_n_t[tmp_first] = tmp_second + "fs";
	//	momtime_n_t.insert(make_pair(tmp_first, tmp_second));
	//	cout << tmp_first << "   " << tmp_second << endl;
	//	//cout << (*it_map).first << "  =  " << (*it_map).second << endl;
	//	count_time_frame++;
	//	it_map++;
	//	//list_namefile.push_back()
	//}


	//cout << " map " << endl;
	//for (it_map = momtime_n_t.begin(); it_map != momtime_n_t.end(); it_map++)
	//{
	//	cout << (*it_map).first << " = " << (*it_map).second << endl;
	//}

	//system("pause");
	//GnuPlot Animate2D(Pacteti);
	//Animate2D.CreateGifOnPlot2D(0, 2, 3, 5, "y(x)", "x", "y", Pacteti, "Animate 2.10.2024", true);

	//MyCreateFile(Points_str.size(), Points_str);

	///////////////////////////////////////////////////////////////////////////////////

	//vector<string> P_x = { "Data_Gif18 + 4.txt", "Data_Gif18 + 9.txt", "Data_Gif18 + 19.txt", "Data_Gif18 + 49.txt", "Data_Gif18 + 59.txt",
	//"Data_Gif18 + 69.txt", "Data_Gif18 + 79.txt", "Data_Gif18 + 119.txt", "Data_Gif18 + 159.txt", "Data_Gif18 + 199.txt" };
	//vector<string> P_z_result;

	//string namefile1;
	//GnuPlot plt_P_X(P_x);
	//plt_P_X.SetParametrs2D(0, 10, 3, "P", "x, mkm", "P (bar)");
	//vector<string> plt_P_X_legend = { "50 fs", "100 fs", "200 fs", "500 fs", "600 fs", "700 fs", "800 fs",  "1200 fs", "1600 fs", "2000 fs" };
	//namefile1 = "Anim P(x) in dif mom time";
	//plt_P_X.ShowDataOnPlot2D(0, 0, 10, plt_P_X_legend, namefile1, true);
	//namefile1.clear();

	//vector<string> P_z = { "Data_Gif19 + 4.txt", "Data_Gif19 + 9.txt", "Data_Gif19 + 19.txt", "Data_Gif19 + 29.txt", "Data_Gif19 + 39.txt",
	//"Data_Gif19 + 59.txt", "Data_Gif19 + 79.txt", "Data_Gif19 + 119.txt", "Data_Gif19 + 159.txt", "Data_Gif19 + 199.txt" };
	//SelectDataForPlot(P_z, P_z_result, 0., 0.14);
	//GnuPlot plt_P_Z(P_z_result);
	//plt_P_Z.SetParametrs2D(0, 10, 3, "P", "z, mkm", "P (bar)");
	////"100 fs", "200 fs", "300 fs", "400 fs", "500 fs", "600 fs", "700 fs", "800 fs", "900 fs", "1000 fs"
	//vector<string> plt_P_Z_legend = { "50 fs", "100 fs", "200 fs", "300 fs", "400 fs", "600 fs", "800 fs", "1200 fs", "1600 fs", "2000 fs" };
	//namefile1 = "Anim P(z) in dif mom time";
	//plt_P_Z.ShowDataOnPlot2D(0, 0, 10, plt_P_Z_legend, namefile1, true);
	//namefile1.clear();

	////////////////////////////////////////////////////////////////////////////////

	//system("pause");

	//GnuPlot plt_experim(24,40);
	cout << endl << endl << "Gnuplot" << endl;
	system("pause");
	system("pause");
	//type metal, type beam, kte, Ce,Ci,	   gamma,    ge,   gi, us	tp s	
																								//I0 = 1e+7
	// double*** V2, double*** a1y, double*** a2y, double*** a1z, double*** a2z, double*** a2x, double*** a1x, double*** b1x, double*** b2x, double*** b1y, double*** b2y, double*** b1z, double*** b2z, double*** e2, double*** e1, double*** F, double*** tmpe0, double*** tmpe1, double*** tmpe2, double*** tmpi0, double*** tmpi1, double*** tmpi2, double** Pxt, double** Pzt, int Nx, int Ny, int Nz, Splayn spl_C_l_on_T)
	MainProcedure(Au, Gauss, 3.115, 19.32, 5.472e-2, 2.550, 2.5e+10, 1.5, 1.11, 3.24e+5, 1e-13, 0.15e+12, V00, V1, V2, a1y, a2y, a1z, a2z, a2x, a1x, b1x, b2x, b1y, b2y, b1z, b2z, e2, e1, tmpe0, tmpe1, tmpe2, tmpi0, tmpi1, tmpi2, Nx_heat, Ny_heat, Nz_heat, Nx_acoustic, Ny_acoustic, Nz_acoustic, spl_C_l_on_T, e_heat, Te_acoustic, Ti_acoustic, Tei_old_xy_new_z, Tei_old_y_new_xz, massive_melting_tmp1, massive_melting_tmp2, "Depth1.txt");
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

