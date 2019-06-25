#include <bits/stdc++.h>

using namespace std;

typedef long long ll;

template<typename T> struct matrix;

template <typename T> ostream &operator<<(ostream &out, matrix<T> mat){
	int n = mat.n;
	int m = mat.m;
	out << "{" << endl;
	for (int i = 0; i < n; i++){
		cout << "\t{";
		for (int j = 0; j < m; j++){
			if (j != 0)
				out << ", ";
			out << mat(i, j);
		}
		out << "}";
		if (i+1 != n)
			out << ",";
		out << endl;
	}
	out << "}";
	return out;
}

template<typename T> struct matrix {
	vector<vector<T>> in;
	int n, m;
	bool is_transposed;
	matrix(){}
	matrix(int n, int m, int op = 0):n(n), m(m), in(n, vector<T>(m, 0)){
		if (op == -1 && n == m)
			for (int i = 0; i < n; i++)
				in[i][i] = 1;
		is_transposed = false;
	}
	matrix(initializer_list<initializer_list<T>> c):
		n(c.size()), m((*c.begin()).size()){
			in = vector<vector<T>>(n, vector<T>(m, 0));
			int i, j;
			i = 0;
			for (auto &it : c){
				j = 0;
				for (auto &jt : it){
					in[i][j] = jt;
					j++;
				}
				i++;
			}
			is_transposed = false;
		}
	matrix(vector<T> v):
		matrix(v.size(), v.size(), 0){
			for (int i = 0; i < v.size(); i++)
				in[i][i] = v[i];
		}
	T &operator()(int i){
		return in[i][0];
	}
	T &operator()(int i, int j){
		if (is_transposed)
			swap(i, j);
		return in[i][j];	
	}
	void operator*=(T r){
		matrix<T> &l = *this;
		for (int i = 0; i < n; i++)
			for (int j = 0; j < m; j++)
				l(i, j) *= r;	
	}
	void operator+=(matrix<T> r){
		if (r.n != n || r.m != m){
			cout << *this << endl;
			cout << r << endl;
			throw logic_error("invalid dimensions while summing matrices!");
		}

		matrix<T> &l = *this;
		for (int i = 0; i < n; i++)
			for (int j = 0; j < m; j++)
				l(i, j) += r(i, j);
	}
	matrix<T> operator+(matrix<T> r){
		matrix<T> l = *this;
		l += r;
		return l;
	}
	matrix<T> operator*(T sc){
		matrix<T> l = *this;
		for (int i = 0; i < n; i++)
			for (int j = 0; j < m; j++)
				l(i, j) *= sc;
		return l;
	}
	matrix<T> operator*(matrix<T> r){
		matrix<T> &l = *this;
		int row = l.n;
		int col = r.m;
		int K = l.m;
		if (l.m != r.n)
			throw logic_error("invalid dimensions while multiplying matrices");
		matrix<T> m(row, col, 0);
		for (int i = 0; i < row; i++)
			for (int j = 0; j < col; j++)
				for (int k = 0; k < K; k++)
					m(i, j) = (m(i, j) + l(i, k)*r(k, j) );
		return m;
	}
	matrix<T> operator-(matrix<T> r){
		matrix<T> l = (*this);
		r *= -1;
		l += r;
		return l;
	}
	matrix<T> operator^(ll e){
		matrix<T> &m = (*this);
		if (e == 0) return matrix(m.n, m.n, 1);
		if (e == 1) return m;
		if (e == 2) return m*m;
		auto m_ = m^(e/2); m_ = m_*m_;
		if (e%2 == 1) m_ = m_ * m;
		return m_;
	}
	T trace(){
		if (n != m) throw logic_error("matrix is not diagonal");
		matrix<T> m = (*this);
		T t = 0;
		for (int i = 0; i < n; i++)
			t += m(i, i);
		return t;
	}
	matrix<T> trans(){
		matrix<T> t = (*this);
		swap(t.n, t.m);
		t.is_transposed = true;
		return t;
	}
	matrix<T> inverse(){
		matrix<T> l(n, n, -1);
		matrix<T> r = *this;
		for (int i = 0; i < n; i++){
			if (abs(r(i, i)) < EPS) continue;
			double inv = 1.0/r(i, i);
			for (int j = 0; j < n; j++){
				l(i, j) = l(i, j) * inv;
				r(i, j) = r(i, j) * inv;
			}
			for (int k = 0; k < n; k++){
				if (k == i) continue;
				if (abs(r(k, i)) < EPS) continue;
				double fat = r(k, i);
				for (int j = 0; j < n; j++){
					l(k, j) = l(k, j) - fat*l(i, j);
					r(k, j) = r(k, j) - fat*r(i, j);
				}
			}
		}
		return l;
	}


	/*vec only*/

	matrix(int n):
		matrix(n, 1){
			for (int i = 0; i < n; i++)
				in[i][0] = 0;
		}
	double norm(){
		matrix<double> v = *this;
		if (n == 2) return hypot(v(0), v(1));
		return sqrt((v.trans()*v).trace());
	}
	matrix<double> unitary(){
		matrix<double> ans = *this;
		return ans*(1/norm());
	}
	matrix<int> dimensions(){
		return {{n}, {m}};
	}
};


template <typename T>
using vec = matrix<T>;

double fabs(vec<double> v){
	return v.norm();
}
