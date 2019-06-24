#include <bits/stdc++.h>

const double EPS = 1e-12;
const double STOPPING_EPS = 1e-12;
const double SQ2 = sqrt(2.0);

#include "matrix.cpp"

using namespace std;

ofstream xout;
ofstream yout;

const int N = 2;

double f0(vec<double> x){
	double ans = 0.0;
	for (int i = 0; i < N-1; i++){
		ans += (100*(x(i+1) - x(i)*x(i))*(x(i+1) - x(i)*x(i))
		    +   (1 - x(i))*(1 - x(i)));
	}
	return ans;
}

double f1(vec<double> x){
	double ans = 0;
	for (int i = 0; i < N; i++)
		ans += x(i)*x(i)*x(i)*x(i) - 16*x(i)*x(i) + 5*x(i);
	return ans;
}

double f2(vec<double> x){
	double l = x(0)*x(0) + x(1) - 11;
	double r = x(0) + x(1)*x(1) - 7;
	return l*l + r*r;
}

vector<vec<double>> canonical_dir(int n){
	vector<vec<double>> ans;
	for (int i = 0; i < n; i++){
		vec<double> v(n);
		v(i) = 1.0;
		ans.push_back(v);
	}
	return ans;
}

//vec<double> gradient(vec<double> x, double (f)(vec<double>)){
//	int n = x.n;
//	if (f == f0){
//		vec<double>	g(n);
//		g(0) = 200*(x(1) - x(0)*x(0))*(-2*x(0)) + 2*(1 - x(0))*-1;
//		for (int i = 1; i < n-1; i++)
//			g(i) = 200*(x(i) - x(i-1)*x(i-1)) + 200*(x(i+1) - x(i)*x(i))*(-2*x(i)) + 2*(1 - x(i))*-1;
//		g(n-1) = 200*(x(n-1) - x(n-2)*x(n-2));
//		return g;
//	}
//	throw logic_error("not implemented yet");
//}

template<typename X> X gradient(X l, double (f)(X), int n){
//	return gradient(l, f);
	X dir = l*0;
	double alpha = f(l);
	for (X &b : canonical_dir(n)){
		X x_ = l + (b*EPS);
		double y_ = f(x_);
		double beta = (1.0/EPS)*(f(x_) - alpha);
		dir = dir + (b*beta);
	}
	return dir;
}


template <typename X> pair<X, double> goldensection(X l, X r, double (f)(X)){
	static double ratioB = 0.618033;
	static double ratioA = 1.000-ratioB;
	X delta = r-l;
	X a = l + delta*ratioA;
	double ya = f(a);
	X b = l + delta*ratioB;
	double yb = f(b);
	char aux;
	X ans = a;
	double f_ans = ya;
	while (fabs(delta) > STOPPING_EPS){
		if (ya < yb){
			r = b;
			delta = r-l;
			b = a;
			yb = ya;
			a = l + delta*ratioA;
			ya = f(a);
			if (ya < f_ans){
				ans = a;
				f_ans = ya;
			}
		}
		else{
			l = a;
			delta = r-l;
			a = b;
			ya = yb;
			b = l + delta*ratioB;
			yb = f(b);
			if (yb < f_ans){
				ans = b;
				f_ans = yb;
			}
		}
	}
	return {ans, f_ans};
}

template <typename X> X right_limit(X l, X dir, double (f)(X) ){
	auto tmpL = l;
	dir = dir.unitary()*EPS;
	double f_x = f(l);
	while (true){
		dir *= 2;
		double f_next_x = f(l + dir);
		if (f_next_x > f_x){
			auto ans = l+dir;
			//if (f(tmpL) < f_next_x) return tmpL;
			return l+dir;
		}
		f_x = f_next_x;
	}
}
/*
https://en.wikipedia.org/wiki/Davidon%E2%80%93Fletcher%E2%80%93Powell_formula
*/
vec<double> DFP(vec<double> l, double (f)(vec<double>)){
	int n = l.n;
	vec<double> grad = gradient(l, f, n);
	matrix<double> I(n, n, -1);
	matrix<double> B_ = I;
	while (true){
		if ((grad.trans()*grad).trace() < STOPPING_EPS) return l;
		vec<double> p = (B_*grad)*-1;
		vec<double>	r = right_limit(l, p, f);
		auto line_search = goldensection(l, r, f);
		vec<double> min_point = line_search.first;
		if ((min_point-l).norm() < STOPPING_EPS) return l;
		if (fabs(f(min_point) - f(l)) < STOPPING_EPS) return l;


		vec<double> s = min_point-l;

		l = min_point;

		vec<double> next_grad = gradient(l, f, n);
		vec<double> y = next_grad - grad;
		grad = next_grad;
		
		double aV = -1.0/(y.trans()*B_*y).trace();
		auto V = (B_*y)*((B_*y).trans())*aV;
		double aW = 1.0/(s.trans()*y).trace();
		auto W = (s*s.trans())*aW;
		B_ += V;
		B_ += W;
	}
}

/*
https://en.wikipedia.org/wiki/Broyden%E2%80%93Fletcher%E2%80%93Goldfarb%E2%80%93Shanno_algorithm
*/
vec<double> BFGS(vec<double> l, double (f)(vec<double>)){
	int n = l.n;
	vec<double> grad = gradient(l, f, n);
	matrix<double> I(n, n, -1);
	matrix<double> B_ = I;
	while (true){
		if ((grad.trans()*grad).trace() < STOPPING_EPS) return l;
		vec<double> p = (B_*grad)*-1;
		vec<double>	r = right_limit(l, p, f);
		auto line_search = goldensection(l, r, f);
		vec<double> min_point = line_search.first;
		if ((min_point-l).norm() < STOPPING_EPS) return l;
		if (fabs(f(min_point) - f(l)) < STOPPING_EPS) return l;


		vec<double> s = min_point-l;

		l = min_point;

		vec<double> next_grad = gradient(l, f, n);
		vec<double> y = next_grad - grad;
		grad = next_grad;
		
		double t = (s.trans()*y).trace();
		double aV = 1.0/(t*t);
		auto V = (s*s.trans())*aV*(s.trans()*y + y.trans()*B_*y).trace();
		double aW = -1.0/(t);
		auto W = (B_*y*s.trans() + s*y.trans()*B_)*aW;

		B_ += V;
		B_ += W;
	}
}


/*
https://en.wikipedia.org/wiki/Ellipsoid_method
*/
vec<double> ellipsoid(vec<double> x, double (f)(vec<double>)){
	int n = x.n;
	matrix<double> aux(n, n), P;
	for (int i = 0; i < n; i++)
		aux(i, i) = 9999;
	double sc = (1.0*n*n)/(n*n-1);
	P = aux;
	while (true){
		if (P.trace() < STOPPING_EPS)
			return x;
		vec<double> g = gradient(x, f, n);
		double den = sqrt((g.trans()*P*g).trace());
		if (den < STOPPING_EPS)
			return x;
		
		vec<double> g_ = g*(1.0/den);	

		x += (P*g_)*(-1.0/(n+1));

		P += (P*(g_*(g_.trans()))*P)*(-2.0/(n+1));

		P *= (1.0*n*n)/(n*n-1);
	}
}


int main(){
	cout << fixed << setprecision(8);
	mt19937 rng((int) chrono::steady_clock::now().time_since_epoch().count());
	uniform_real_distribution<double> distribution(-50.0, 50.0);


	auto f = f2; //funcao otimizada
	auto method = BFGS; //método de otimizacao
	vec<double> i(N);
	for (int j = 0; j < N; j++)
		i(j) = distribution(rng);
	i(0) = abs(i(0));
	auto ans = method(i, f);
	cout << "starting point " << i << endl;
	cout << "minimum found at " << ans << endl;
	cout << "f(x) = " << f(ans) << endl;
	return 0;
}
