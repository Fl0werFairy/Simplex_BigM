#include<iostream>
#include<iomanip>
#include<vector>
#include<cmath>
#include<utility>
using namespace std;
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

#define GRE(x, y) (x > y && !equal(x, y))
#define GEQ(x, y) (x > y || equal(x, y))
#define SIMPLEX_BIG_M


const double EPS = 1e-5;

bool equal(double a, double b) {
	return fabs(a - b) < EPS;
}


Row::Row() {
	size = 0;
	data = nullptr;
}

Row::Row(int size, double *data): size(size) {
	this->data = data;
}

double& Row::operator[](int index) {
	if (index < 0 || index >= size) {
		cerr << "Error!" << endl;
		exit(0);
	}
	return data[index];
}

double Row::operator[](int index) const {
	if (index < 0 || index >= size) {
		cerr << "Error!" << endl;
		exit(0);
	}
	return data[index];
}

void Row::addRow(Row &another, double ratio) {
	if (another.size != size) {
		cerr << "Error: Row size " << another.size << " " << size;
		exit(-1);
	}
	for (int i = 0; i < size; i++) {
		data[i] += another.data[i] * ratio;
	}
}

void Row::multi(double ratio) {
	for (int i = 0; i < size; i++) {
		data[i] *= ratio;
	}
}


class Matrix {
public:
	Matrix();

	Matrix(int row, int column);

	Matrix(int row, int column, double constant);

	Matrix(int row, int column, const double *data);

	Matrix(const Matrix &mat);

	~Matrix();



	Matrix getColumn(int index) const;


	Matrix getColumns(int beginIndex, int endIndex) const;


	void appendColumn(Matrix &mat);

	Matrix operator+(const Matrix &mat) const;
	Matrix operator-(const Matrix &mat) const;
	Matrix operator*(double scale) const;
	friend Matrix operator*(double scale, const Matrix &mat);
	Matrix operator*(const Matrix &mat) const;

	Row operator[](int index) const;

	friend istream &operator>>(istream &is, const Matrix &mat);

	friend ostream &operator<<(ostream &os, const Matrix &mat);

//private:
	int row;
	int column;
	double *data;
};
Matrix::Matrix() {
	row = 0;
	column = 0;
	data = nullptr;
}

Matrix::Matrix(int row, int column) : row(row), column(column) {
	data = new double[row * column];
	for (int i = 0; i < row * column; i++)
		data[i] = 0;
}

Matrix::Matrix(int row, int column, double constant) : row(row), column(column) {
	data = new double[row * column];
	for (int i = 0; i < row * column; i++)
		data[i] = constant;
}

Matrix::Matrix(int row, int column, const double *data) : row(row), column(column) {
	this->data = new double[row * column];
	for (int i = 0; i < row * column; i++)
		this->data[i] = data[i];
}

Matrix::Matrix(const Matrix &mat) {
	row = mat.row;
	column = mat.column;
	data = new double[row * column];
	for (int i = 0; i < row * column; i++)
		data[i] = mat.data[i];
}

Matrix::~Matrix() {
	if (data != nullptr) {
		delete[] data;
		data = nullptr;
	}
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
	if (row != mat.row) {
		cerr << "Error!" << endl;
		exit(0);
	}
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
	if (row != mat.row || column != mat.column) {
		cerr << "Error!" << endl;
		exit(0);
	}
	Matrix result(row, column);
	for (int i = 0; i < row * column; i++)
		result.data[i] = data[i] + mat.data[i];
	return result;
}

Matrix Matrix::operator-(const Matrix &mat) const {
	if (row != mat.row || column != mat.column) {
		cerr << "Error!" << endl;
		exit(0);
	}
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
	if (column != mat.row) {
		cerr << "Error!" << endl;
		exit(0);
	}
	Matrix result(row, mat.column);
	for (int i = 0; i < row; i++)
		for (int j = 0; j < mat.column; j++)
			for (int k = 0; k < column; k++)
				result.data[i * row + j] += data[i * row + k] * mat.data[k * row + j];
	return result;
}

Row Matrix::operator[](int index) const {
	if (index < 0 || index >= row) {
		cerr << "Index out of range: " << index << endl;
		exit(0);
	}
	return {column, data + index * column};
}
istream &operator>>(istream &is, const Matrix &mat) {
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


class Solver {
public:
	Solver(int n, int m, Matrix c, Matrix a, Matrix b, Matrix d, Matrix e);
	~Solver();

	void relax();
	void normalize();
	virtual void solve(int &k, double &y, Matrix &x) = 0;
	void recover(Matrix &x);


protected:
	int n, m;
	Matrix c;
	Matrix a, b, d;
	Matrix e;
	vector<int> negative;  // xi -> -xi, if xi <= 0
	vector<pair<int, int> > noConstraints;  // xi -> xi - xj, if no constraints on xi
};
Solver::Solver(int n, int m, Matrix c, Matrix a, Matrix b, Matrix d, Matrix e) :
		n(n), m(m), c(c), a(a), b(b), d(d), e(e) {
	negative.clear();
	noConstraints.clear();
}

Solver::~Solver() = default;

void Solver::relax() {
	Matrix zero = Matrix(1, 1);
	Matrix one = Matrix(1, 1, 1);
	for (int i = 0; i < m; i++) {
		if (equal(d[i][0], 0)) continue;
		if (equal(d[i][0], 1)) {
			n++;
			c.appendColumn(zero);
			e.appendColumn(one);
			Matrix newColumn(m, 1);
			newColumn[i][0] = -1;
			a.appendColumn(newColumn);
			d[i][0] = 0;
			continue;
		}
		if (equal(d[i][0], -1)) {
			n++;
			c.appendColumn(zero);
			e.appendColumn(one);
			Matrix newColumn(m, 1);
			newColumn[i][0] = 1;
			a.appendColumn(newColumn);
			d[i][0] = 0;
			continue;
		}
	}
}

void Solver::normalize() {
	Matrix one = Matrix(1, 1, 1);
	relax();
	for (int i = 0; i < m; i++) {
		if (b[i][0] >= 0) continue;
		b[i][0] *= -1;
		for (int j = 0; j < n; j++) {
			a[i][j] *= -1;
		}
	}
	for (int j = 0; j < n; j++) {
		if (equal(e[0][j], 1)) continue;
		if (equal(e[0][j], -1)) {
			c[0][j] *= -1;
			e[0][j] = 1;
			for (int i = 0; i < m; i++) {
				a[i][j] *= -1;
			}
			negative.push_back(j);
			continue;
		}
		if (equal(e[0][j], 0)) {
			n++;
			c.appendColumn(one);
			c[0][n - 1] = c[0][j] * -1;
			e[0][j] = 1;
			e.appendColumn(one);
			Matrix newColumn = a.getColumn(j) * -1;
			a.appendColumn(newColumn);
			noConstraints.emplace_back(j, n - 1);
			continue;
		}
	}
}

void Solver::recover(Matrix &x) {
	for (auto iter : negative) {
		int p = negative[iter];
		x[0][p] *= -1;
	}
	for (auto iter = noConstraints.begin(); iter != noConstraints.end(); iter++) {
		int p = iter->first;
		int q = iter->second;
		x[0][p] = x[0][p] - x[0][q];
	}
}

#define Maxm 2500

class BigMSolver : public Solver {
private:
	double value = 0;
	int nonManualVariableCount = 0;
	vector<int> baseIndex;

	void exchange(int inIndex, int outIndex);

public:
	BigMSolver(int n, int m, Matrix c, Matrix a, Matrix b, Matrix d, Matrix e);

	~BigMSolver();

	void solve(int &k, double &y, Matrix &x) override;
};
BigMSolver::BigMSolver(int n, int m, Matrix c, Matrix a, Matrix b, Matrix d,
                       Matrix e) :
		Solver(n, m, c, a, b, d, e) {}

BigMSolver::~BigMSolver() = default;

void BigMSolver::exchange(int inIndex, int outIndex) {

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
	Matrix matrixM = Matrix(1, 1, Maxm);
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
	for (int i = 0; i < n; i++) {
		checkedArray[i] = -checkedArray[i];
	}
	// normalize checked number of manual variables
	for (int i = 0; i < m; i++) {
		auto rowi = a[i];
		checkedArray.addRow(rowi, Maxm);
		double delta = Maxm * b[i][0];
		value += delta;
	}

	while (true) {
		// get max
		int inIndex = -1;
		double maxCheckedNum = 0;
		for (int i = 0; i < n; i++) {
			double num = checkedArray[i];
			if (GRE(num, 0)) {
				// judge unbound
				bool flag = true;
				for (int j = 0; j < m; j++) {
					double aji = a[j][i];
					if (GRE(aji, 0)) {
						flag = false;
						break;
					}
				}
				if (flag) {
					// All numbers in column i are <= 0: unbound!
					k = 0;
					// return;
				}
				// get max index
				if (!flag&&GRE(num, maxCheckedNum)) {
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
			if (outIndex == -1 || GRE(minRatio, ratio) || (equal(minRatio, ratio) && baseIndex[i] < minOutIndex)) {
				minRatio = ratio;
				outIndex = i;
				minOutIndex = baseIndex[i];
			}
		}
		// outIndex == -1 is impossible, due to this case is unbound!
		exchange(inIndex, outIndex);


	}
}

int main() {
    int N,M,K;
    cin>>N>>M;
    cin>>K;
    int S,T;
    S=1;T=N;

    Matrix Flow(N+M,M+1),con(N+M,1);
    Matrix sign(N+M,1);
    Matrix cost(1,M);
    Matrix pm(1,M,1);
    for(int i=0;i<M;++i)
    {
        static int u,v,c,p;
        cin>>u>>v>>c>>p;
        Flow.data[i*(M+1)+i]=1;
        con.data[i]=Flow.data[i*(M+1)+M]=c;
        sign.data[i]=-1;
        Flow.data[(u-1+M)*(M+1)+i]=1;
        Flow.data[(v-1+M)*(M+1)+i]=-1;
        cost.data[i]=p;
    }

    con.data[M+S-1]=Flow.data[(M+S-1)*(M+1)+M]=K;
    sign.data[M+S-1]=1;
    con.data[M+T-1]=Flow.data[(M+T-1)*(M+1)+M]=-K;
    sign.data[M+T-1]=-1;

    BigMSolver bs(M,N+M,cost,Flow.getColumns(0,M),con,sign,pm);
    int flag;
    double ans;
    Matrix x(1,M+2*(N+M)-1);
    bs.normalize();
    bs.solve(flag,ans,x);
    bs.recover(x);
    if(flag==-1)cout<<"Infeasible";
    else if(flag==1)
    {
        cout<<fixed<<setprecision(3)<<ans;
    }
	// int n, m;
	// cin >> m >> n;

	// Matrix c(1, n);
    // cin>>c;


	// Matrix mat(m, n + 1);
	// cin >> mat;

	// Matrix a = mat.getColumns(0, n);
	// Matrix b = mat.getColumn(n);
	// Matrix d(m,1,-1);

	// Matrix e(1, n,1);
	// int k;
	// double y;
	// Matrix x(1, n + 2 * m);

	// BigMSolver bigMSolver(n, m, c, a, b, d, e);

	// bigMSolver.normalize();

	// bigMSolver.solve(k, y, x);
	// bigMSolver.recover(x);

	// cout << fixed << setprecision(3);

	// if(k==-1)cout<<"Infeasible";
	// else if(k==0)cout<<"Unbounded";
	// else if (k == 1) {
	// 	cout << y << endl;
	// 	for (int i = 0; i < n; i++) {
	// 		cout << x[0][i] << " ";
	// 	}
	// 	cout << endl;
	// }

}