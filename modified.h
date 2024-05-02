#ifndef ATEST_TRIAL_H
#define ATEST_TRIAL_H
#include <vector>
#include <cmath>
#include <iostream>
namespace trial{
	class Row {
	public:
		Row();
		Row(int size, double *data);
		double &operator[](int index);
		double operator[](int index) const;
		void addRow(Row &another, double ratio);
		void multi(double ratio);

	private:
		int size;
		double *data;
	};

#define GNE(x, y) (x > y && !equal(x, y))
#define GEQ(x, y) (x > y || equal(x, y))

	const double EPS = 1e-8;
	inline bool equal(double a, double b) {
		return fabs(a - b) < EPS;
	}
	class Matrix {
	public:
		Matrix();
		Matrix(int row, int column);
		Matrix(int row, int column, double constant);
		Matrix(const Matrix &mat);
		~Matrix()=default;

		[[nodiscard]] Matrix getColumn(int index) const;
		[[nodiscard]] Matrix getColumns(int beginIndex, int endIndex) const;

		void appendColumn(Matrix &mat);
		Matrix operator+(const Matrix &mat) const;
		Matrix operator-(const Matrix &mat) const;
		Matrix operator*(double scale) const;
		friend Matrix operator*(double scale, const Matrix &mat);
		Matrix operator*(const Matrix &mat) const;
		Row operator[](int index) const;

		friend std::istream &operator>>(std::istream &is, Matrix &mat);
		friend std::ostream &operator<<(std::ostream &os, const Matrix &mat);

//private:
		int row;
		int column;
		double *data;
	};
	class Solver {
	public:
		Solver(int n, int m, const Matrix& c, const Matrix& a, const Matrix& b, const Matrix& d, const Matrix& e);
		~Solver()=default;

		void relax();
		void normalize();
		virtual void solve(int &k, double &y, Matrix &x) = 0;
		void recover(Matrix &x);

	protected:
		int n, m;
		Matrix c;
		Matrix a, b, d;
		Matrix e;
		std::vector<int> negative;  // xi -> -xi, if xi <= 0
		std::vector<std::pair<int, int> > noConstraints;  // xi -> xi - xj, if no constraints on xi
	};
#define M 1e3

	class BigMSolver : public Solver {
	private:
		double value = 0;
		int nonManualVariableCount = 0;
		std::vector<int> baseIndex;

		void pivot(int inIndex, int outIndex);

	public:
		BigMSolver(int n, int m, Matrix c, Matrix a, Matrix b, Matrix d, Matrix e);

		~BigMSolver()=default;

		void solve(int &k, double &y, Matrix &x) override;
	};
	std::pair<int,std::vector<double>> work(std::ifstream&);
}
#endif
