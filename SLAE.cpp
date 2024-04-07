#include "SLAE.h"

SLAE::SLAE(size_t Rows, size_t Cols) : _Matrix(Rows, Cols), _FreeCol(VecD(Cols, 0.0))
{}

SLAE::SLAE(VecVecD& Values, VecD& FreeCol) : _Matrix(Values), _FreeCol(FreeCol)
{}

SLAE::SLAE(Matrixx& Matrix, VecD& FreeCol) : _Matrix(Matrix), _FreeCol(FreeCol)
{}

SLAE::SLAE(VecD& Diagonal, VecD& FreeCol) : _Diagonal(Diagonal), _FreeCol(FreeCol)
{}

VecD SLAE::GetSolution()
{
	SLAE temp = *this;
	size_t size = temp._Matrix.GetRows(), size2 = temp._Matrix.GetCols();
	if ((size != size2) || size == 0 || size2 == 0) {
		throw std::exception("Imposible to solve this slae");
	}

	VecD result(size, 0.0);
	temp.MakeItTriangle();
	double Solution;
	for (int64_t i = size - 1; i >= 0; --i)
	{
		Solution = temp._FreeCol[i];
		for (size_t j = size - 1; j > i; --j)
		{
			Solution -= result[j] * temp._Matrix[i][j];
		}
		Solution /= temp._Matrix[i][i];
		if (abs(Solution) < EPS) {
			result[i] = 0;
		}
		else {
			result[i] = Solution;
		}
	}
	return result;
}

VecD SLAE::SolveProgon()
{
	size_t Size = this->_Matrix.GetRows();
	VecD alpha(Size, 0), beta(Size, 0);

	alpha[1] = (-this->_Matrix[0][1]) / this->_Matrix[0][0];
	beta[1] = this->_FreeCol[0] / this->_Matrix[0][0];

	for (size_t i = 2; i < Size; ++i)
	{
		alpha[i] = (-this->_Matrix[i - 1][i]) / (this->_Matrix[i - 1][i - 2] * alpha[i - 1] + this->_Matrix[i - 1][i - 1]);
		beta[i] = (this->_FreeCol[i - 1] - this->_Matrix[i - 1][i - 2] * beta[i - 1]) / (this->_Matrix[i - 1][i - 2] * alpha[i - 1] + this->_Matrix[i - 1][i - 1]);
	}

	VecD Result(Size, 0);

	Result[Size - 1] = (this->_FreeCol[Size - 1] - this->_Matrix[Size - 1][Size - 2] * beta[Size - 1]) / (this->_Matrix[Size - 1][Size - 1] + this->_Matrix[Size - 1][Size - 2] * alpha[Size - 1]);

	for (size_t i = Size - 1; i > 0; --i) {
		Result[i - 1] = alpha[i] * Result[i] + beta[i];
	}

	return Result;
}

VecD SLAE::SolveProgonLinearInterpolate()
{
	for (size_t i = 1; i < this->_Diagonal.size() - 2; i += 2) {
		this->_FreeCol[i] = (this->_FreeCol[i] - this->_FreeCol[i - 1]) / (this->_Diagonal[i]);
	}
	this->_FreeCol[this->_Diagonal.size() - 1] = (this->_FreeCol[this->_Diagonal.size() - 1] - this->_FreeCol[this->_Diagonal.size() - 2]) / (this->_Diagonal[this->_Diagonal.size() - 1]);
	return this->_FreeCol;
}

Matrixx SLAE::GetMatrix()
{
	return this->_Matrix;
}

VecD SLAE::GetFreeCol()
{
	return this->_FreeCol;
}

void SLAE::SwitchRows(size_t Which, size_t Where)
{
	this->_Matrix.SwitchRows(Which, Where);
	double temp = this->_FreeCol[Where];
	this->_FreeCol[Where] = this->_FreeCol[Which];
	this->_FreeCol[Which] = temp;
}

void SLAE::MakeItTriangle()
{
	size_t Rows = this->_FreeCol.size();
	double ValueToMulty;
	for (size_t i = 0; i < Rows; ++i)
	{
		if (abs(this->_Matrix[i][i]) <= EPS) {
			this->SwitchRows(i, this->_Matrix.GetIndexOfBiggestValue(i, i));
			if (abs(this->_Matrix[i][i]) <= EPS) {
				throw std::exception("Determinant = 0");
			}
		}
		for (size_t j = i + 1; j < Rows; ++j)
		{
			ValueToMulty = this->_Matrix[j][i] / this->_Matrix[i][i];
			for (size_t k = i; k < Rows; ++k) {
				this->_Matrix[j][k] -= (ValueToMulty * this->_Matrix[i][k]);
			}
			this->_FreeCol[j] -= (ValueToMulty * this->_FreeCol[i]);
		}
	}
}

Matrixx::Matrixx(size_t Rows, size_t Cols) : _Rows(Rows), _Cols(Cols)
{
	for (size_t i = 0; i < Rows; ++i) {
		this->_Values.push_back(VecD(Cols, 0.0));
	}
}

Matrixx::Matrixx(VecVecD& values) : _Values(values)
{
	this->_Rows = values.size();
	if (this->_Rows != 0) {
		this->_Cols = values[0].size();
	}
	else {
		this->_Cols = 0;
	}
}

Matrixx::Matrixx(VecD values, bool gorizontal)
{
	if (gorizontal)
	{
		this->_Rows = 1;
		this->_Cols = values.size();
		this->_Values.push_back(values);
	}
	else
	{
		this->_Cols = 1;
		this->_Rows = values.size();
		for (size_t i = 0; i < this->_Rows; i++) {
			this->_Values.push_back(VecD(1, values[i]));
		}
	}
}

size_t Matrixx::GetRows()
{
	return this->_Rows;
}

size_t Matrixx::GetCols()
{
	return this->_Cols;
}

//double Matrix::GetDeterminant()
//{
//	size_t Rows = this->_Rows, Cols = this->_Cols, NumberOfChanges = 0;
//	double ValueToMulty;
//	Matrix temp(*this);
//	if (Rows != Cols) {
//		throw std::exception("Imposible to calc determinant");
//	}
//	try
//	{
//		NumberOfChanges = temp.MakeItTriangle();
//		if (abs(temp._Values[Rows - 1][Cols - 1]) <= EPS)
//			return 0;
//		double Result = 1;
//		for (int i = 0; i < Rows; ++i)
//			Result *= temp._Values[i][i];
//		return (Result * pow(-1, NumberOfChanges));
//	}
//	catch (std::exception& exept)
//	{
//		return 0;
//	}
//}

Matrixx Matrixx::GetTransposed()
{
	Matrixx Result(this->_Cols, this->_Rows);
	for (size_t i = 0; i < this->_Rows; ++i) {
		for (size_t j = 0; j < this->_Cols; ++j) {
			Result._Values[j][i] = this->_Values[i][j];
		}
	}
	return Result;
}

Matrixx Matrixx::GetReverse()
{
	int size = this->_Values.size();

	VecVecD mat(size, VecD(2 * size));
	for (size_t i = 0; i < size; i++) {
		for (size_t j = 0; j < size; j++) {
			mat[i][j] = this->_Values[i][j];
		}
		mat[i][i + size] = 1.0;
	}

	for (size_t i = 0; i < size; i++) {
		double Divisor = mat[i][i];
		for (size_t j = 0; j < 2 * size; j++) {
			mat[i][j] /= Divisor;
		}
		for (size_t k = 0; k < size; k++) {
			if (k != i) {
				double Multiplier = mat[k][i];
				for (size_t j = 0; j < 2 * size; j++) {
					mat[k][j] -= Multiplier * mat[i][j];
				}
			}
		}
	}
	VecVecD Reveresed(size, std::vector<double>(size));
	for (size_t i = 0; i < size; i++) {
		for (size_t j = 0; j < size; j++) {
			Reveresed[i][j] = mat[i][j + size];
		}
	}
	return Matrixx(Reveresed);
}

void Matrixx::SwitchRows(size_t Which, size_t Where)
{
	if (Where >= this->_Values.size() || Where >= this->_Values.size()) {
		throw "Invalid Index";
	}
	VecD temp = this->_Values[Where];
	this->_Values[Where] = this->_Values[Which];
	this->_Values[Which] = temp;
}

size_t Matrixx::GetIndexOfBiggestValue(size_t Row, size_t Col)
{
	double BiggestValue = 0;
	size_t IndexOfRow = Row;
	for (size_t i = Row + 1; i < this->_Rows; ++i)
	{
		if (abs(this->_Values[i][Col]) > BiggestValue)
		{
			BiggestValue = abs(this->_Values[i][Col]);
			IndexOfRow = i;
		}
	}
	return IndexOfRow;
}

double Matrixx::GetNorma()
{
	double result = 0;
	for (size_t i = 0; i < this->_Rows; i++) {
		for (size_t j = 0; j < this->_Cols; j++) {
			result += (this->_Values[i][j] * this->_Values[i][j]);
		}
	}
	return sqrt(result);
}

Matrixx Matrixx::operator-(const Matrixx& other)
{
	if (this->_Rows != other._Rows || this->_Cols != other._Cols) {
		throw std::exception("Imposible to subtract these matrices");
	}

	Matrixx Result(this->_Rows, this->_Cols);
	for (size_t i = 0; i < this->_Rows; ++i) {
		for (size_t j = 0; j < this->_Cols; ++j) {
			Result._Values[i][j] = this->_Values[i][j] - other._Values[i][j];
		}
	}
	return Result;
}

Matrixx Matrixx::operator*(const Matrixx& other)
{
	if (this->_Cols != other._Rows) {
		throw std::exception("Imposible to multiply these matrices");
	}

	Matrixx Result(this->_Rows, other._Cols);
	for (size_t i = 0; i < Result._Rows; ++i) {
		for (size_t j = 0; j < Result._Cols; ++j) {
			for (size_t k = 0; k < this->_Cols; ++k) {
				Result._Values[i][j] += this->_Values[i][k] * other._Values[k][j];
			}
		}
	}
	return Result;
}

VecD& Matrixx::operator[](size_t Index)
{
	if (Index >= this->_Rows) {
		throw std::exception("Invalid index");
	}
	return this->_Values[Index];
}

size_t Matrixx::MakeItTriangle()
{
	size_t Rows = this->_Rows, Cols = this->_Cols, NumberOfChanges = 0;
	double ValueToMulty;
	for (size_t i = 0; i < Rows; ++i)
	{
		if (abs(this->_Values[i][i]) <= EPS)
		{
			this->SwitchRows(i, this->GetIndexOfBiggestValue(i, i));
			if (abs(this->_Values[i][i]) <= EPS) {
				throw std::exception("Det = 0");
			}
			NumberOfChanges++;
		}
		for (size_t j = i + 1; j < Rows; ++j)
		{
			ValueToMulty = this->_Values[j][i] / this->_Values[i][i];
			for (int k = i; k < Cols; k++) {
				this->_Values[j][k] -= (ValueToMulty * this->_Values[i][k]);
			}
		}
	}
	return NumberOfChanges;
}