#pragma once
#include <vector>

typedef std::vector<double> VecD;
typedef std::vector<std::vector<double>> VecVecD;

#define EPS ((double)(pow(10, -14)))

class SLAE;
class Matrixx;

class Matrixx
{
public:
	Matrixx(size_t Rows = 0, size_t Cols = 0);
	Matrixx(VecVecD& values);
	Matrixx(VecD values, bool gorizontal = true);
	size_t GetRows();
	size_t GetCols();
	//double GetDeterminant();
	Matrixx GetTransposed();
	Matrixx GetReverse();
	void SwitchRows(size_t Which, size_t Where);
	size_t GetIndexOfBiggestValue(size_t Row, size_t Col);
	double GetNorma();

	//---------------------------//

	Matrixx operator-(const Matrixx& other);
	Matrixx operator*(const Matrixx& other);
	VecD& operator[](size_t Index);

private:

	size_t MakeItTriangle();
	//---------------------------//
	VecVecD _Values;
	size_t _Rows;
	size_t _Cols;
};

class SLAE
{
public:

	SLAE(size_t Rows, size_t Cols);
	SLAE(VecVecD& Values, VecD& FreeCol);
	SLAE(Matrixx& Matrix, VecD& FreeCol);
	SLAE(VecD& Matrix, VecD& FreeCol);
	VecD GetSolution();
	VecD SolveProgon();
	VecD SolveProgonLinearInterpolate();
	Matrixx GetMatrix();
	VecD GetFreeCol();

private:

	void SwitchRows(size_t Which, size_t Where);
	void MakeItTriangle();

	Matrixx _Matrix;
	VecD _Diagonal;
	VecD _FreeCol;
};
