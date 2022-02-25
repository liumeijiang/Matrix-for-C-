#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <cassert>
#include<math.h>

using namespace std;

typedef unsigned int uint;
#define THRESHOLD 1e-6

class Matrix
{
public:
	Matrix();
	Matrix(uint row, uint col);
	Matrix(const Matrix &t);
	Matrix(double *a,double *b);
	~Matrix();

	static Matrix ones(uint row, uint col);
	static Matrix zeros(uint row, uint col);
	static Matrix unit(uint n);
	static uint rank(Matrix t,Matrix **out);

	friend ostream& operator<<(ostream &cout, Matrix &t);
	double* operator[](uint i);
	Matrix operator+(const Matrix &t);
	Matrix operator*(const Matrix &t);
	Matrix &operator=(const Matrix &t);
	Matrix &operator=(double *a);
	Matrix operator*(double n);

	uint getCol() { return col; }
	uint getRow() {return row;}
	uint getSize() {return size;}

	double det();
	Matrix &resize(uint row, uint col);
	Matrix T();
	Matrix adjoint();
	Matrix inv();
	Matrix &swap(uint a, uint b, uint key);
	Matrix &insert(Matrix &t,uint clip,uint key);
	Matrix &deleteClip(uint clip, uint key);
	Matrix getClip(uint clip, uint key);
	Matrix &clipAdd(Matrix t, uint clip, uint key);
	Matrix &clipMulti(double n, uint clip, uint key);

private:
	double* ptr;
	uint row, col;
	uint size;
	static Matrix complementMinor(Matrix &t, uint row, uint col);
	
};

ostream& operator<<(ostream& cout, Matrix& t)
{
	for (uint i = 0; i < t.getRow(); i++)
	{
		for (uint j = 0; j < t.getCol(); j++)
		{
			cout << t[i][j] << "\t";
		}
		cout << "\n";
	}
	return cout;
}

Matrix::Matrix(uint row, uint col) :row(row), col(col)
{
	this->size = col * row;
	this->ptr = new double[size];
}

Matrix Matrix::operator+(const Matrix& t)
{
	Matrix* temp = new Matrix(row, col);
	for (uint j = 0; j < this->size; j++)
	{
		*(temp->ptr + j) = *(this->ptr + j) + *(t.ptr + j);
	}
	return *temp;
}
Matrix& Matrix::operator=(const Matrix& t)
{
	if (this != &t)
	{

		this->row = t.row;
		this->col = t.col;
		this->size = col * row;
		if (ptr != NULL) delete[] ptr;
		this->ptr = new double[size];
		for (uint i = 0; i < t.size; i++)
		{
			*(this->ptr + i) = *(t.ptr + i);
		}
	}
	return *this;
}
Matrix::Matrix(double* a, double* b)//从数组构造（a,a+length），返回行矩阵
{
	assert(b > a);
	uint len = ((uint)b - (uint)a) / 8;
	this->row = 1;
	this->col = len;
	this->size = len;
	this->ptr = new double[len];
	for (uint i = 0; i < len; i++)
	{
		*(this->ptr + i) = *(a + i);
	}
}
Matrix Matrix::operator*(const Matrix& t)
{
	assert(col == t.row);
	Matrix* temp = new Matrix(row, t.col);
	double tt;
	for (uint i = 0; i < temp->row; i++)
	{
		for (uint j = 0; j < temp->col; j++)
		{
			tt = 0;
			for (uint k = 0; k < col; k++)
			{
				tt = tt + *(this->ptr + i * col + k) * *(t.ptr + k * t.col + j);
			}
			*(temp->ptr + j + i * temp->col) = tt;
		}
	}
	return *temp;
}

Matrix Matrix::operator*(double n)
{
	Matrix* temp = new Matrix(row, col);
	for (uint j = 0; j < this->size; j++)
	{
		*(temp->ptr + j) = *(this->ptr + j) * n;
	}
	return *temp;
}
double* Matrix::operator[](uint i) {

	return this->ptr + i * col;
}
Matrix::Matrix(const Matrix& t)
{

	this->row = t.row;
	this->col = t.col;
	this->size = col * row;
	this->ptr = new double[size];
	for (uint i = 0; i < t.size; i++)
	{
		*(this->ptr + i) = *(t.ptr + i);
	}

}

Matrix::~Matrix()
{
	//cout << "Matrix::~Matrix()" << this->size<<endl;
	if (ptr != NULL)
	{
		delete[] ptr;
		ptr = NULL;
	}
}


Matrix Matrix::ones(uint row, uint col)//返回全1矩阵
{
	Matrix* t = new Matrix(row, col);
	for (uint i = 0; i < t->size; i++)
	{
		*(t->ptr + i) = 1;
	}
	return *t;
}
Matrix Matrix::zeros(uint row, uint col)//返回全零矩阵
{
	Matrix* t = new Matrix(row, col);
	for (uint i = 0; i < t->size; i++)
	{
		*(t->ptr + i) = 0;
	}
	return *t;
}
Matrix Matrix::unit(uint n)//返回单位矩阵
{
	Matrix t = Matrix::zeros(n, n);
	for (uint i = 0; i < n; i++)
	{
		for (uint j = 0; j < n; j++)
		{
			if (i == j)
			{
				t[i][j] = 1;
			}
		}
	}
	return t;
}
double Matrix::det() //返还矩阵的行列式
{
	assert(col == row && col > 0);
	if (row == 1)
	{
		return *(this->ptr);
	}
	if (row == 2) {
		return (*this)[0][0] * (*this)[1][1] - (*this)[0][1] * (*this)[1][0];
	}
	double output = 0;
	for (int j = 0; j < col; j++)
	{
		Matrix t = Matrix::complementMinor(*this, 0, j);
		output = output + pow(-1, j) * t.det() * (*this)[0][j];
	}
	return output;
}
Matrix Matrix::complementMinor(Matrix& t, uint row, uint col) //返还余子式矩阵
{
	Matrix* output = new Matrix(t.row - 1, t.col - 1);

	uint n = 0;
	for (uint i = 0; i < t.row; i++)
	{
		for (uint j = 0; j < t.col; j++)
		{
			if (i != row && j != col)
			{
				*(output->ptr + n++) = t[i][j];
			}
		}
	}
	return *output;
}
Matrix& Matrix::resize(uint row, uint col)
{
	assert(row * col == this->size);
	this->row = row;
	this->col = col;
	return *this;
}
Matrix& Matrix::swap(uint a, uint b, uint key)//交换矩阵两行（key =0）或两列（key=1）
{
	if (a == b)
	{
		return *this;
	}
	double temp = 0;
	switch (key)
	{
	case 0:
		assert(a < row&& b < row);
		for (uint i = 0; i < col; i++)
		{
			temp = (*this)[a][i];
			(*this)[a][i] = (*this)[b][i];
			(*this)[b][i] = temp;
		}
		return *this;
	case 1:
		assert(a < col&& b < col);
		for (uint i = 0; i < row; i++)
		{
			temp = (*this)[i][a];
			(*this)[i][a] = (*this)[i][b];
			(*this)[i][b] = temp;
		}
		return *this;
	default:
		assert(key < 2);
	}
}
Matrix& Matrix::insert(Matrix& t, uint clip, uint key)//在clip前插入一个矩阵
{
	double* temp = new double[size + t.size];
	switch (key)
	{
	case 0:
		assert(col == t.col);
		if (clip > row) clip = row;

		for (uint i = 0; i < size + t.size; i++)
		{
			if (i < clip * col) *(temp + i) = *(this->ptr + i);
			else if (i < clip * col + t.size) *(temp + i) = *(t.ptr + i - clip * col);
			else *(temp + i) = *(this->ptr + i - clip * col - t.size);
		}
		row = row + t.row;
	case 1:
		assert(row == t.row);
		if (clip > col) clip = col;
		for (uint i = 0; i < row; i++)
		{
			for (uint j = 0; j < col + t.col; j++)
			{
				if (j < clip) *(temp + i * (col + t.col) + j) = (*this)[i][j];
				else if (j < clip + t.col) *(temp + i * (col + t.col) + j) = t[i][j - clip];
				else  *(temp + i * (col + t.col) + j) = (*this)[i][j - clip - t.col];
			}
		}
		col = col + t.col;
	default:
		assert(key < 2);
	}
	size = size + t.size;
	delete[] ptr;
	ptr = temp;
	return *this;
}
Matrix Matrix::T()//返回转置矩阵
{
	Matrix* temp = new Matrix(row, col);

	for (uint i = 0;i < row;i++)
	{
		for (uint j = 0; j < col; j++)
		{
			*(temp->ptr + j * row + i) = (*this)[i][j];
		}
	}
	return *temp;
}
Matrix Matrix::adjoint()//返回伴随矩阵
{
	Matrix* output = new Matrix(this->row, this->col);
	for (uint i = 0; i < row; i++)
	{
		for (uint j = 0; j < col; j++)
		{
			Matrix t = Matrix::complementMinor(*this, i, j);
			(*output)[j][i] = t.det() * pow(-1, i + j);
		}
	}
	return *output;
}
Matrix Matrix::inv()//返回逆矩阵
{
	double det = this->det();
	assert(det != 0);
	return this->adjoint() * (1 / det);
}
Matrix& Matrix::deleteClip(uint clip, uint key)
{
	uint n = 0;
	switch (key)
	{
	case 0:
		assert(clip < row);
		this->row = row - 1;
		for (uint i = clip; i < row; i++)
		{
			for (uint j = 0; j < col; j++)
			{
				(*this)[i][j] = (*this)[i + 1][j];
			}
		}
		size = row * col;
		return *this;
	case 1:
		assert(clip < col);
		for (uint i = 0; i < row; i++)
		{
			for (uint j = 0; j < col; j++)
			{
				int shift = 0;
				if (i * j >= (col - 1) * row)
				{
					break;
				}
				if (j < clip)
				{
					if ((n + i) % col == clip) shift = 1;
					*(this->ptr + n) = *(this->ptr + n + i + shift);
					n++;
				}
				else {
					if ((n + i + 1) % col == clip) shift = 1;
					*(this->ptr + n) = *(this->ptr + n + i + 1 + shift);
					n++;
				}

			}
		}

		col = col - 1;
		size = col * row;
		return *this;
	default:
		assert(key < 2);
	}
}
Matrix Matrix::getClip(uint clip, uint key)
{
	Matrix* t = NULL;
	switch (key)
	{
	case 0:
		assert(clip < row);
		t = new Matrix(1, col);
		for (uint i = 0; i < col; i++)
		{
			*(t->ptr + i) = (*this)[clip][i];
		}
		break;
	case 1:
		assert(clip < col);
		t = new Matrix(row, 1);
		for (uint i = 0; i < row; i++)
		{
			*(t->ptr + i) = (*this)[i][clip];
		}
		break;
	default:
		assert(key < 2);
		break;
	}
	return *t;
}
Matrix& Matrix::clipAdd(Matrix t, uint clip, uint key = 0)
{
	switch (key)
	{
	case 0:
		assert(clip < row&& t.row == 1);
		for (uint i = 0; i < col; i++)
		{
			(*this)[clip][i] += t[0][i];
		}
		break;
	case 1:
		assert(clip < col&& t.col == 1);
		for (uint i = 0; i < row; i++)
		{
			(*this)[i][clip] += t[i][0];
		}
		break;
	default:
		assert(key < 2);
		break;
	}
	return *this;
}
Matrix& Matrix::clipMulti(double n, uint clip, uint key = 0)
{
	switch (key)
	{
	case 0:
		assert(clip < row);
		for (uint i = 0; i < col; i++)
		{
			(*this)[clip][i] *= n;
		}
		break;
	case 1:
		assert(clip < col);
		for (uint i = 0; i < row; i++)
		{
			(*this)[i][clip] *= n;
		}
		break;
	default:
		assert(key < 2);
		break;
	}
	return *this;
}
uint Matrix::rank(Matrix t, Matrix** out = NULL)
{
	if (t.row > t.col)
	{
		t = t.T();
	}
	uint shift = 0;
	for (uint i = 0; i < t.col; i++)
	{
		double max = 0;
		uint n = i;
		for (uint j = i - shift; j < t.row; j++)
		{
			if (abs(t[j][i]) > abs(max))
			{
				max = t[j][i];
				n = j;
			}
		}
		t.swap(i - shift, n, 0);
		if (abs(max) < THRESHOLD)
		{
			shift++;
			continue;
		}
		t.clipMulti(1 / t[i - shift][i], i - shift);
		for (uint j = i + 1 - shift; j < t.row; j++)
		{
			t.clipAdd(t.getClip(i - shift, 0) * -t[j][i], j);
		}
		cout << t << endl;
	}
	*out = new Matrix(t);
	return t.col - shift;
}
