#pragma once
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include "Matrix.h"

using namespace std;

enum Plane { xy, xz, yz };
enum Profile { fun_on_x, fun_on_y, fun_on_z, fun_on_t };

class GnuPlot {

	int number_of_plots;

	vector<FILE*> gnuplotPipe;
	vector<ofstream> file;
	vector<ifstream> file_in;
	vector<string> filename;
	vector<Matrix> P_xyzt;

	string GetGPPath();
	string GetGPTerminal();
	size_t GetTermnalWidth(); //������ ������ 

public:
	//  �������� ������� ������� � ������� ������� � �������� ���� (����������: fun � SetDataOnPlotColor3D)!!!!!

	// ����������� ��� �������� �������� � ������� ����� ���������� ������ ���������� ����� � ������� � ������� � ��������
	GnuPlot(int number_of_plots_); // ��� ������ ������� �� ���������� ������ // ++++++++++++++++++
	// ����������� ��� �������� �������� � ������� ����� ���������� ����������������� ������ ������� � ������ ����������
	GnuPlot(int number_of_plots_, vector<string> filename_);

	// ���������� ����� ��� ����������� �������(��) y(x)
	void SetParametrs2D(int Current_number_plot, int number_of_lines, int width_line, string title_plot, string xlabel, string ylabel); // +++ ��������� � ����� �����  // Current_number_plot=0,1,2, ....
	// ���������� ����� ��� ����������� ������ �������� �������� (������)
	void SetParametrsOnPlotColor(int Current_number_plot, string title_plot, string xlabel, string ylabel, long float right_bondary_x, long float top_bondary_y); // Current_number_plot=0,1,2, ....

	// ��� ������� ����� ������� ������ 1 ��� ����� SetDataOnPlot2D (��� 2-� �-��� ���� � ����� � SetDataOnPlot1D(2D, 3D))
	// ���������� ���� 2-�� �-���: �-��� SetDataOnPlot3D ���������� � ����� �� ������� + ��������� ��������� ������ ��� ������ � ����
	void SetGridOnPlot2D(int Current_number_plot, int Nx, int Ny, double dx, double dy, int number_of_lines, Profile prof);
	// ��� ������� ����� ������� ������ 1 ��� ����� SetDataOnPlot3D (��� 2-� �-��� ���� � ����� � SetDataOnPlot1D(2D, 3D))
	// ���������� ���� 2-�� �-���: �-��� SetDataOnPlot3D ���������� � ����� �� ������� + ��������� ��������� ������ ��� ������ � ����
	void SetGridOnPlot3D(int Current_number_plot, int Nx, int Ny, int Nz, double dx, double dy, double dz, int number_of_lines, Profile prof);

	// 1D, 2D, 3D - ����������� ������
   ///(����������� �������� ������) y(x) (Ce �� Te) (�������� ��� ��������� ������� �� ������������ ����� ������� ��� � ������(�� ���������� � �����))
   // ��� ��������� ������� ��� 1D - ������ ���� (U(x), U(t))  
	void SetDataOnPlot1D(int Current_number_plot, int Nx, double dx, double** fun_dimensionless, double parametr_for_dimension, int moment_of_time, double dt, Profile prof);
	// ��� ��������� ������� ������� (y(x)) ��� 2D - ������ ���� (U(x), U(y), U(t))  
	void SetDataOnPlot2D(int Current_number_plot, int Nx, int Ny, double dx, double dy, double** fun_dimensionless, double parametr_for_dimension, int fixed_point_on_axis_x, int fixed_point_on_axis_y, int number_of_lines, int current_number_of_line, double moment_of_time, vector<double**> vec, Profile prof);// ������� ������ ���� �������, ������ ������������� � �������� y(x), � �� ���������
	// ��� ��������� ������� ������� (y(x) ��� ����� t) ��� 3D - ������ ���� (U(x), U(y), U(z),U(t)) 
	void SetDataOnPlot3D(int Current_number_plot, int Nx, int Ny, int Nz, double dx, double dy, double dz, double*** fun_dimensionless, double parametr_for_dimension, int fixed_point_on_axis_x, int fixed_point_on_axis_y, int fixed_point_on_axis_z, int number_of_lines, int current_number_of_line, double moment_of_time, vector<double***> vec, Profile prof); // linetype � ���������� �� ������ � ������ (�� ������ ����� ����)

	// ��� ��������� ���������� ������� (���� ����������� �������� 2D)
	void SetDataOnPlotColor2D(int Current_number_plot, int Nx, int Ny, double dx, double dy, double** fun, double parametr_for_dimension);
	// ��� ��������� ���������� ������� (������� ����������� � ������� (���������) 3-� ������� ������� - ����)
	void SetDataOnPlotColor3D(int Current_number_plot, int Nx, int Ny, int Nz, double dx, double dy, double dz, double*** fun_dimensionless, double parametr_for_dimension, int fixed_point_on_axis, Plane plane);

	// ��������� ������ ������� (���� ���� ����� ����� �� �������) 
	void ShowDataOnPlot2D(int Current_number_plot, int number_of_lines, vector<string> list_name_line, string name_of_file, bool png_);
	// ��������� ������ �������� �������� (������)
	void ShowDataOnPlotColor(int Current_number_plot, string name_of_file, bool png_);
	// ����� ���� ��������� ���� �� ������� ������� (��� ��������� �������� �������� ������ �������)
	void Close_and_open_files_for_replot(vector<int> Current_number_of_plots_); //Current_number_of_plots_ - ������ ������� ���������� ������ ������ ������� ����� ������� � ����� �������
	void Close_all_files_and_plots(int number_of_plots_); // ++++++++++++++++++

};