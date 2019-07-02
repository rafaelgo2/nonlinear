#include <bits/stdc++.h>

const double EPS = 1e-12;
const double STOPPING_EPS = 1e-12;
const double HESSIAN_EPS = 1e-4;
const double SQ2 = sqrt(2.0);

#include "matrix.cpp"

using namespace std;

int it = 0;

void log(ofstream &logger, vec<double> v, double (f)(vec<double>)){
	logger << fixed << setprecision(5);
	logger << f(v) << '\t';
	logger << endl;
	it++;
}

int N;

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

double f3(vec<double> x){
	return x(0)*x(0) + x(1)*x(1);
}

double f4(vec<double> x){
	return 0.26*(x(0)*x(0) + x(1)*x(1)) - 0.48*x(0)*x(1);
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

vec<double> gradient(vec<double> l, double (f)(vec<double>)){
	//	return gradient(l, f);
	vec<double> dir = l*0;
	double alpha = f(l);
	int n = l.n;
	for (vec<double> &b : canonical_dir(n)){
		vec<double> x_ = l + (b*EPS);
		double y_ = f(x_);
		double beta = (1.0/EPS)*(f(x_) - alpha);
		dir = dir + (b*beta);
	}
	return dir;
}

matrix<double> hessian(vec<double> x, double (f)(vec<double>)){
	matrix<double> h;
	if (false && f == f2){
		h = {
			{4*(x(0)*x(0) + x(1) - 11) + 8*x(0)*x(0) + 2, 4*x(0) + 4*x(1)},
			{4*x(0) + 4*x(1), 4*(x(0) + x(1)*x(1) - 7) + 8*x(1)*x(1)+2}
		};
	}
	int n = x.n;
	auto D = canonical_dir(n);
	h = matrix<double>(n, n);

	double EPS = HESSIAN_EPS;
	for (int i = 0; i < n; i++){
		for (int j = i; j < n; j++){
			double ans = 0.0;
			ans += f(x + D[i]*EPS + D[j]*EPS);
			ans -= f(x + D[i]*EPS);
			ans -= f(x + D[j]*EPS);
			ans += f(x);
			h(j, i) = h(i, j) = ans/(EPS*EPS);
		}
	}
	return h;
}

pair<vec<double>, double> armijo(vec<double> l, vec<double> r, double (f)(vec<double>)){
	auto d = r-l;
	auto alpha = d.norm();
	d = d.unitary();


	double f_l = f(l);

	double beta = 1e-4;
	double sigma = 0.5;

	double f_g = (gradient(l, f).trans()*d).trace();


	while (true){
		if (f(l + d*alpha) <= f_l + beta*alpha*f_g) break;
		alpha *= sigma;
	}

	return {l + d*alpha, f(l + d*alpha)};
}


pair<vec<double>, double> goldensection(vec<double> l, vec<double> r, double (f)(vec<double>)){
	static double ratioB = 0.618033;
	static double ratioA = 1.000-ratioB;
	vec<double> delta = r-l;
	vec<double> a = l + delta*ratioA;
	double ya = f(a);
	vec<double> b = l + delta*ratioB;
	double yb = f(b);
	char aux;
	vec<double> ans = a;
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

auto line_search_method = goldensection;

template <typename X> X right_limit(X l, X dir, double (f)(X) ){
	auto tmpL = l;
	dir = dir.unitary()*EPS;
	double f_x = f(l);
	while (true){
		dir *= 1.6;
		double f_next_x = f(l + dir);
		if (f_next_x > f_x){
			auto ans = l+dir;
			//if (f(tmpL) < f_next_x) return tmpL;
			return l+dir;
		}
		f_x = f_next_x;
	}
}

vec<double> gradient_method(vec<double> l, double (f)(vec<double>)){
	ofstream logger("logs/gradient_method.txt");
	it = 0;
	auto start = chrono::steady_clock::now();
	auto leave = [&]{
		auto end = chrono::steady_clock::now();
		logger
			<< "Elapsed time in milliseconds : " 
			<< chrono::duration_cast<chrono::milliseconds>(end - start).count()
			<< " ms" << endl;
		logger << it << endl;
		return l;
	};
	int n = l.n;
	vec<double> grad = gradient(l, f);
	while (true){
		log(logger, l, f);
		if ((grad.trans()*grad).trace() < STOPPING_EPS) return leave();
		vec<double> p = grad*-1;
		vec<double>	r = right_limit(l, p, f);
		auto line_search = line_search_method(l, r, f);
		vec<double> min_point = line_search.first;
		if ((min_point-l).norm() < STOPPING_EPS) return leave();
		if (fabs(f(min_point) - f(l)) < STOPPING_EPS) return leave();
		l = min_point;
		grad = gradient(l, f);
	}
}

vec<double> newton_method(vec<double> l, double (f)(vec<double>)){
	ofstream logger("logs/newton_method.txt");
	it = 0;
	auto start = chrono::steady_clock::now();
	auto leave = [&]{
		auto end = chrono::steady_clock::now();
		logger
			<< "Elapsed time in milliseconds : " 
			<< chrono::duration_cast<chrono::milliseconds>(end - start).count()
			<< " ms" << endl;
		logger << it << endl;
		return l;
	};
	int n = l.n;
	vec<double> grad = gradient(l, f);
	while (true){
		log(logger, l, f);
		vec<double> p = hessian(l, f).inverse()*grad*-1;
		vec<double> min_point = l + p;
		if ((min_point-l).norm() < STOPPING_EPS) return leave();
		if (fabs(f(min_point) - f(l)) < STOPPING_EPS) return leave();
		l = min_point;
		grad = gradient(l, f);
	}
}

vec<double> newton_method_line_search(vec<double> l, double (f)(vec<double>)){
	ofstream logger("logs/newton_method_line_search.txt"); 
	it = 0;
	auto start = chrono::steady_clock::now();
	auto leave = [&]{
		auto end = chrono::steady_clock::now();
		logger
			<< "Elapsed time in milliseconds : " 
			<< chrono::duration_cast<chrono::milliseconds>(end - start).count()
			<< " ms" << endl;
		logger << it << endl;
		return l;
	};
	int n = l.n;
	vec<double> grad = gradient(l, f);
	while (true){
		log(logger, l, f);
		vec<double> p = hessian(l, f).inverse()*grad*-1;

		vec<double>	r = right_limit(l, p, f);
		auto line_search = line_search_method(l, r, f);
		vec<double> min_point = line_search.first;
		if ((min_point-l).norm() < STOPPING_EPS) return leave();
		if (fabs(f(min_point) - f(l)) < STOPPING_EPS) return leave();
		l = min_point;
		grad = gradient(l, f);
	}
}

/*
https://en.wikipedia.org/wiki/Davidon%E2%80%93Fletcher%E2%80%93Powell_formula
 */
vec<double> DFP(vec<double> l, double (f)(vec<double>)){
	cout << "DFP" << endl;
	ofstream logger("logs/DFP.txt");
	it = 0;
	auto start = chrono::steady_clock::now();
	auto leave = [&]{
		auto end = chrono::steady_clock::now();
		logger
			<< "Elapsed time in milliseconds : " 
			<< chrono::duration_cast<chrono::milliseconds>(end - start).count()
			<< " ms" << endl;
		logger << it << endl;
		return l;
	};
	int n = l.n;
	vec<double> grad = gradient(l, f);
	matrix<double> I(n, n, -1);
	matrix<double> B_ = I;
	int IT_MAX = 1e8;
	while (IT_MAX--){
		log(logger, l, f);
		if ((grad.trans()*grad).trace() < STOPPING_EPS) return leave();
		vec<double> p = (B_*grad)*-1;
		vec<double>	r = right_limit(l, p, f);
		auto line_search = line_search_method(l, r, f);
		vec<double> min_point = line_search.first;
		if ((min_point-l).norm() < STOPPING_EPS) return leave();
		if (fabs(f(min_point) - f(l)) < STOPPING_EPS) return leave();
		//cout << fabs(f(min_point) - f(l)) << endl;
		vec<double> s = min_point-l;

		l = min_point;

		vec<double> next_grad = gradient(l, f);
		vec<double> y = next_grad - grad;
		grad = next_grad;

		double aV = -1.0/(y.trans()*B_*y).trace();
		auto V = (B_*y)*((B_*y).trans())*aV;
		double aW = 1.0/(s.trans()*y).trace();
		auto W = (s*s.trans())*aW;
		B_ += V;
		B_ += W;
	}
	return leave();
}

/*
https://en.wikipedia.org/wiki/Broyden%E2%80%93Fletcher%E2%80%93Goldfarb%E2%80%93Shanno_algorithm
 */
vec<double> BFGS(vec<double> l, double (f)(vec<double>)){
	ofstream logger("logs/BFGS.txt");
	it = 0;
	auto start = chrono::steady_clock::now();
	auto leave = [&]{
		auto end = chrono::steady_clock::now();
		logger
			<< "Elapsed time in milliseconds : " 
			<< chrono::duration_cast<chrono::milliseconds>(end - start).count()
			<< " ms" << endl;
		logger << it << endl;
		return l;
	};
	int n = l.n;
	vec<double> grad = gradient(l, f);
	matrix<double> I(n, n, -1);
	matrix<double> B_ = I;
	while (true){
		log(logger, l, f);
		if ((grad.trans()*grad).trace() < STOPPING_EPS) return leave();
		vec<double> p = (B_*grad)*-1;
		vec<double>	r = right_limit(l, p, f);
		auto line_search = line_search_method(l, r, f);
		vec<double> min_point = line_search.first;
		if ((min_point-l).norm() < STOPPING_EPS) return leave();
		if (fabs(f(min_point) - f(l)) < STOPPING_EPS) return leave();


		vec<double> s = min_point-l;

		l = min_point;

		vec<double> next_grad = gradient(l, f);
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

vec<double> conjugate_gradient_method_sq(vec<double> l, double(f)(vec<double>)){
	ofstream logger("logs/conjugate_gradient_method_sq.txt");
	it = 0;
	auto start = chrono::steady_clock::now();
	auto leave = [&]{
		auto end = chrono::steady_clock::now();
		logger
			<< "Elapsed time in milliseconds : " 
			<< chrono::duration_cast<chrono::milliseconds>(end - start).count()
			<< " ms" << endl;
		logger << it << endl;
		return l;
	};
	auto grad = gradient(l, f);
	auto d = grad*-1;
	while (true){
		log(logger, l, f);
		if ((grad.trans()*grad).trace() < STOPPING_EPS) return leave();
		auto B = hessian(l, f);
		double alpha = -(grad.trans()*d).trace()/(d.trans()*(B*d)).trace();
		vec<double> min_point = l + d*alpha;
		if ((min_point-l).norm() < STOPPING_EPS) return leave();
		if (fabs(f(min_point) - f(l)) < STOPPING_EPS) return leave();
		auto r = min_point;
		auto next_grad = gradient(r, f);
		double beta = (next_grad.trans()*(B*d)).trace()/(d.trans()*(B*d)).trace();

		d = next_grad*-1 + d*beta;
		grad = next_grad;
		l = r;
	}
}

vec<double> conjugate_gradient_method_uni(vec<double> l, double(f)(vec<double>)){
	ofstream logger("logs/conjugate_gradient_method_uni.txt");  
	it = 0;
	auto start = chrono::steady_clock::now();
	auto leave = [&]{
		auto end = chrono::steady_clock::now();
		logger
			<< "Elapsed time in milliseconds : " 
			<< chrono::duration_cast<chrono::milliseconds>(end - start).count()
			<< " ms" << endl;
		logger << it << endl;
		return l;
	};
	auto grad = gradient(l, f);
	auto d = grad*-1;
	while (true){
		log(logger, l, f);
		if ((grad.trans()*grad).trace() < STOPPING_EPS) return leave();
		vec<double>	r = right_limit(l, d, f);
		auto line_search = line_search_method(l, r, f);
		vec<double> min_point = line_search.first;
		if ((min_point-l).norm() < STOPPING_EPS) return leave();
		if (fabs(f(min_point) - f(l)) < STOPPING_EPS) return leave();
		l = min_point;
		auto next_grad = gradient(l, f);
		double beta = (next_grad.trans()*d).trace()/(d.trans()*d).trace();

		d = next_grad*-1 + d*beta;
		grad = next_grad;
	}
}

/*
https://en.wikipedia.org/wiki/Ellipsoid_method
 */
vec<double> ellipsoid(vec<double> x, double (f)(vec<double>)){
	ofstream logger("logs/ellipsoid.txt");
	it = 0;
	auto start = chrono::steady_clock::now();
	auto leave = [&]{
		auto end = chrono::steady_clock::now();
		logger
			<< "Elapsed time in milliseconds : " 
			<< chrono::duration_cast<chrono::milliseconds>(end - start).count()
			<< " ms" << endl;
		logger << it << endl;
		return x;
	};
	int n = x.n;
	matrix<double> aux(n, n), P;
	for (int i = 0; i < n; i++)
		aux(i, i) = 9999;
	double sc = (1.0*n*n)/(n*n-1);
	P = aux;
	while (true){
		log(logger, x, f);
		if (P.trace() < STOPPING_EPS) return leave();
		vec<double> g = gradient(x, f);
		double den = sqrt((g.trans()*P*g).trace());
		if (den < STOPPING_EPS) return leave();

		vec<double> g_ = g*(1.0/den);	

		x += (P*g_)*(-1.0/(n+1));

		P += (P*(g_*(g_.trans()))*P)*(-2.0/(n+1));

		P *= (1.0*n*n)/(n*n-1);
	}
}


int main(){
	cout << fixed << setprecision(8);
	mt19937 rng((int) chrono::steady_clock::now().time_since_epoch().count());
	double lim = 5.0;
	uniform_real_distribution<double> distribution(-lim, lim);


	auto f = f1; //funcao otimizada
	N = 20; //dimensao do vetor do dominio
	line_search_method = goldensection; //busca linear

	vec<double> i(N);
	for (int j = 0; j < N; j++){
		i(j) = distribution(rng);
	}
	cout << "starting point " << i << endl;
	for (auto method : {gradient_method, conjugate_gradient_method_sq, conjugate_gradient_method_uni, newton_method, newton_method_line_search, DFP, BFGS, ellipsoid}){
		auto ans = method(i, f);
		cout << "minimum found at " << ans << endl;
		cout << "f(x) = " << f(ans) << endl;
	}
	return 0;
}
