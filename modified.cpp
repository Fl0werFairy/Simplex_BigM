#include <fstream>
#include <iomanip>
#include <algorithm>
#include "modified.h"
using namespace std;
namespace trial{
	Row::Row():size(0),data(nullptr) {}
	Row::Row(int s, double *d): size(s),data(d) {}

	double& Row::operator[](int index) {
		return data[index];
	}
	double Row::operator[](int index) const {
		return data[index];
	}
	void Row::addRow(Row &another, double ratio) {
		for (int i = 0; i < size; i++)
			data[i] += another.data[i] * ratio;
	}
	void Row::multi(double ratio) {
		for (int i = 0; i < size; i++)
			data[i] *= ratio;
	}

	Matrix::Matrix():row(0),column(0),data(nullptr){}
	Matrix::Matrix(int _row, int _column) : row(_row), column(_column) {
		data = new double[_row * _column];
		fill(data,data+_row*_column,0);
	}
	Matrix::Matrix(int _row, int _column, double constant) : row(_row), column(_column) {
		data = new double[_row * _column];
		fill(data,data+_row*_column,constant);
	}

	Matrix::Matrix(const Matrix &mat) {
		row = mat.row;
		column = mat.column;
		data = new double[row * column];
		copy(mat.data,mat.data+row*column,data);
	}

	Matrix Matrix::getColumn(int index) const {
		Matrix result(row, 1);
		for (int i = 0; i < row; i++)
			result.data[i] = data[i * column + index];
		return result;
	}

	Matrix Matrix::getColumns(int beginIndex, int endIndex) const {
		Matrix result(row, endIndex - beginIndex);
		int count = 0;
		for (int i = 0; i < row; i++)
			for (int j = beginIndex; j < endIndex; j++)
				result.data[count++] = data[i * column + j];
		return result;
	}

	void Matrix::appendColumn(Matrix &mat) {
		int count = 0;
		int total = row * (column + mat.column);
		auto *tempData = new double[total];
		for (int i = 0; i < row; i++) {
			for (int j = 0; j < column; j++)
				tempData[count++] = data[i * column + j];
			for (int j = 0; j < mat.column; j++)
				tempData[count++] = mat.data[i * mat.column + j];
		}
		delete[] data;
		column = column + mat.column;
		data = tempData;
	}

	Matrix Matrix::operator+(const Matrix &mat) const {
		Matrix result(row, column);
		for (int i = 0; i < row * column; i++)
			result.data[i] = data[i] + mat.data[i];
		return result;
	}
	Matrix Matrix::operator-(const Matrix &mat) const {
		Matrix result(row, column);
		for (int i = 0; i < row * column; i++)
			result.data[i] = data[i] - mat.data[i];
		return result;
	}
	Matrix Matrix::operator*(double scale) const {
		Matrix result(row, column);
		for (int i = 0; i < row * column; i++)
			result.data[i] = scale * data[i];
		return result;
	}
	Matrix operator*(double scale, const Matrix &mat) {
		int row = mat.row;
		int column = mat.column;
		Matrix result(row, column);
		for (int i = 0; i < row * column; i++)
			result.data[i] = scale * mat.data[i];
		return result;
	}
	Matrix Matrix::operator*(const Matrix &mat) const {
		Matrix result(row, mat.column);
		for (int i = 0; i < row; i++)
			for (int j = 0; j < mat.column; j++)
				for (int k = 0; k < column; k++)
					result.data[i * row + j] += data[i * row + k] * mat.data[k * row + j];
		return result;
	}
	Row Matrix::operator[](int index) const {
		return {column, data + index * column};
	}
	istream &operator>>(istream &is, Matrix &mat) {
		int row = mat.row;
		int column = mat.column;
		for (int i = 0; i < row * column; i++)
			is >> mat.data[i];
		return is;
	}
	ostream &operator<<(ostream &os, const Matrix &mat) {
		int row = mat.row;
		int column = mat.column;
		for (int i = 0; i < row * column; i++)
			os << mat.data[i] << " \n"[i % column == column - 1];
		return os;
	}

	Solver::Solver(int n, int m, const Matrix& c, const Matrix& a, const Matrix& b, const Matrix& d, const Matrix& e) :
			n(n), m(m), c(c), a(a), b(b), d(d), e(e) {}

	void Solver::relax() {
		Matrix zero = Matrix(1, 1);
		Matrix one = Matrix(1, 1, 1);
		for (int i = 0; i < m; i++) {
			if (equal(d[i][0], 1)) {
				n++;
				c.appendColumn(zero);
				e.appendColumn(one);
				Matrix newColumn(m, 1);
				newColumn[i][0] = -1;
				a.appendColumn(newColumn);
				d[i][0] = 0;
			}
			else if (equal(d[i][0], -1)) {
				n++;
				c.appendColumn(zero);
				e.appendColumn(one);
				Matrix newColumn(m, 1);
				newColumn[i][0] = 1;
				a.appendColumn(newColumn);
				d[i][0] = 0;
			}
		}
	}
	void Solver::normalize() {
		Matrix one = Matrix(1, 1, 1);
		relax();
		for (int i = 0; i < m; i++) {
			if (b[i][0] >= 0) continue;
			b[i][0] *= -1;
			for (int j = 0; j < n; j++)
				a[i][j] *= -1;
		}
		for (int j = 0; j < n; j++) {
			if (equal(e[0][j], -1)) {
				c[0][j] *= -1;
				e[0][j] = 1;
				for (int i = 0; i < m; i++) {
					a[i][j] *= -1;
				}
				negative.push_back(j);
			}
			else if (equal(e[0][j], 0)) {
				n++;
				c.appendColumn(one);
				c[0][n - 1] = c[0][j] * -1;
				e[0][j] = 1;
				e.appendColumn(one);
				Matrix newColumn = a.getColumn(j) * -1;
				a.appendColumn(newColumn);
				noConstraints.emplace_back(j, n - 1);
			}
		}
	}

	void Solver::recover(Matrix &x) {
		for (auto iter : negative) {
			int p = negative[iter];
			x[0][p] *= -1;
		}
		for (auto [p,q] :noConstraints)
			x[0][p] = x[0][p] - x[0][q];
	}


	BigMSolver::BigMSolver(int n, int m, Matrix c, Matrix a, Matrix b, Matrix d,
	                       Matrix e) :
			Solver(n, m, c, a, b, d, e) {}

	void BigMSolver::pivot(int inIndex, int outIndex) {

		baseIndex[outIndex] = inIndex;
		Row pivot = a[outIndex];
		Row pivotB = b[outIndex];
		// normalize to one
		double toOneRatio = 1 / a[outIndex][inIndex];
		pivot.multi(toOneRatio);
		pivotB.multi(toOneRatio);

		for (int i = 0; i < m; i++) {
			if (i != outIndex) {
				Row another = a[i];
				Row anotherB = b[i];
				double ratio = -another[inIndex];
				another.addRow(pivot, ratio);
				anotherB.addRow(pivotB, ratio);
			}

		}
		// transform checked number to zero
		Row cRow = c[0];
		double ratio = -cRow[inIndex];
		cRow.addRow(pivot, ratio);
		value += pivotB[0] * ratio;
	}

	void BigMSolver::solve(int &k, double &y, Matrix &x) {
		nonManualVariableCount = n;
		// add artificial variables
		Matrix zero = Matrix(m, 1);
		Matrix one = Matrix(1, 1, 1);
		Matrix matrixM = Matrix(1, 1, M);
		for (int i = 0; i < m; i++) {
			a.appendColumn(zero);
			a[i][n] = 1;
			c.appendColumn(matrixM);
			e.appendColumn(one);
			baseIndex.push_back(n);
			n++;
		}
		// to max format
		Row checkedArray = c[0];
		for (int i = 0; i < n; i++)
			checkedArray[i] = -checkedArray[i];

		// normalize checked number of manual variables
		for (int i = 0; i < m; i++) {
			auto rowi = a[i];
			checkedArray.addRow(rowi, M);
			double delta = M * b[i][0];
			value += delta;
		}

		while (true) {
			// get max
			int inIndex = -1;
			double maxCheckedNum = 0;
			for (int i = 0; i < n; i++) {
				double num = checkedArray[i];
				if (GNE(num, 0)) {
					// judge unbound
					bool flag = true;
					for (int j = 0; j < m; j++) {
						double aji = a[j][i];
						if (GNE(aji, 0)) {
							flag = false;
							break;
						}
					}
					if (flag)k = 0;
						// All numbers in column i are <= 0: unbound!
					// get max index
					if (!flag && GNE(num, maxCheckedNum)) {
						inIndex = i;
						maxCheckedNum = num;
					}
				}
			}
			if (inIndex == -1) {
				// all of checked numbers are <= 0
				// judge manual variables
				for (int i = 0; i < m; i++) {
					if (baseIndex[i] >= nonManualVariableCount && b[i][0] != 0) {
						// no feasible solution!
						k = -1;
						return;
					}
				}
				if(k==0)return;
				// output answer
				k = 1;
				for (int i = 0; i < m; i++) {
					int base = baseIndex[i];
					x[0][base] = b[i][0];
				}
				y = value;
				return;
			}
			int outIndex = -1;
			int minOutIndex = 0x7fffffff;
			double minRatio = 0;
			for (int i = 0; i < m; i++) {
				double aik = a[i][inIndex];
				if (GEQ(0, aik)) continue;
				double ratio = b[i][0] / aik;
				if (outIndex == -1 || GNE(minRatio, ratio) || (equal(minRatio, ratio) && baseIndex[i] < minOutIndex)) {
					minRatio = ratio;
					outIndex = i;
					minOutIndex = baseIndex[i];
				}
			}
			// outIndex == -1 is impossible, due to this case is unbound!
			pivot(inIndex, outIndex);
		}
	}

	pair<int,vector<double>> work(ifstream& fin)
	{
		int n, m;
		fin >> m >> n;
		trial::Matrix c(1, n);
		fin>>c;
		c=c*-1;

		trial::Matrix mat(m, n + 1);
		fin >> mat;

		trial::Matrix a = mat.getColumns(0, n);
		trial::Matrix b = mat.getColumn(n);
		trial::Matrix d(m,1,-1);
		trial::Matrix e(1, n,1);

		int k;
		double y;
		trial::Matrix x(1, n + 2 * m);
		trial::BigMSolver bigMSolver(n, m, c, a, b, d, e);
		bigMSolver.normalize();
		bigMSolver.solve(k, y, x);
		bigMSolver.recover(x);

		vector<double> ans{-y};
		for (int i = 0; i < n; ++i)ans.emplace_back(x[0][i]);

		cout << fixed << setprecision(6);
		if(k==-1)cout<<"Infeasible\n";
		else if(k==0)cout<<"Unbounded\n";
		else if (k == 1) {
			cout << ans[0] << endl;
			for (int i = 1; i <= n; i++) {
				cout << ans[i] << " ";
			}
			cout << endl;
		}
		return {k, ans};
	}
}
