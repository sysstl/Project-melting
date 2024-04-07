#pragma once
#include <vector>

typedef std::vector<double> VecD;
typedef std::vector<std::vector<double>> VecVecD;

#define EPS ((double)(pow(10, -14)))

class SLAE;

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
