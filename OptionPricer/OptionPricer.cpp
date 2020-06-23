#include<iostream>
#include<vector>
#include<math.h>
#include<iterator>
#include<algorithm>
#include<boost/math/distributions/normal.hpp>
using boost::math::normal;


template <class T>
class EuropeanOption {
private:
	T S_, K_, r_, sig_, T_;

public:
	EuropeanOption(const T &S, const T &K, const T &r, const T &sig, const T &T) {
		std::cout << "Running my constructor..." << std::endl;
		S_ = S;
		K_ = K;
		r_ = r;
		sig_ = sig;
		T_ = T;

	}
	~EuropeanOption() {
		std::cout << "Running my destructor..." << std::endl;
	}
	void CallPricer_BSE();
	void CallPricer_BM();

};

template <class T>
void EuropeanOption<T>::CallPricer_BM() {
	std::cout << "Running Binomial Model call pricer..." << std::endl;

	std::vector<std::vector<float>> stock_matrix;
	std::vector<std::vector<float>> price_matrix;
	std::vector<float> temp_vec;

	int gridcount = 50;
	float u, q;
	u = exp(sig_*pow(T_ / gridcount, 0.5));
	q = (exp(r_*T_ / gridcount) - (1 / u)) / (u - 1 / u);

	for (int ia = 0;ia <= gridcount;ia++) {
		std::vector<float>(gridcount - ia + 1).swap(temp_vec);
		stock_matrix.push_back(temp_vec);
		price_matrix.push_back(temp_vec);
		for (int ib = 0;ib < price_matrix[ia].size();ib++) {
			stock_matrix[ia][ib] = S_ * pow(u, price_matrix[ia].size() - ib) / pow(u, ib);
			std::cout << "stock price [" << ia << "][" << ib << "]: " << stock_matrix[ia][ib] << std::endl;
			if (ia == 0) {
				price_matrix[ia][ib] = std::max(stock_matrix[ia][ib] - K_, 0.0f);
				std::cout << "call price [" << ia << "][" << ib << "]: " << price_matrix[ia][ib] << std::endl;
			}
			else {
				price_matrix[ia][ib] = price_matrix[ia - 1][ib] * q + price_matrix[ia - 1][ib + 1] * (1 - q);
				std::cout << "call price: [" << ia << "][" << ib << "]: " << price_matrix[ia][ib] << std::endl;
			}
		}
	}
	std::cout << "FINAL call price : " << price_matrix[gridcount][0] << std::endl;
}

template <class T>
void EuropeanOption<T>::CallPricer_BSE() {
	std::cout << "Running BSE call pricer..." << std::endl;
	normal sn(0, 1);

	float d1 = (log(S_ / K_) + (r_ + 0.5*pow(sig_, 2)*T_)) / (sig_*pow(T_, 0.5));
	float d2 = d1 - sig_ * pow(T_, 0.5);
	float  thisprice;
	std::vector<float> result;

	std::cout << "S: " << S_ << std::endl;
	std::cout << "K: " << K_ << std::endl;
	std::cout << "r: " << r_ << std::endl;
	std::cout << "sig: " << sig_ << std::endl;
	std::cout << "T: " << T_ << std::endl;
	std::cout << "d1: " << d1 << std::endl;
	std::cout << "d2: " << d2 << std::endl;
	std::cout << "cdf(d1): " << cdf(sn, d1) << std::endl;
	std::cout << "cdf(d2): " << cdf(sn, d2) << std::endl;

	thisprice = S_ * cdf(sn, d1) - K_ * exp(-r_ * T_)*cdf(sn, d2);
	std::cout << "BSE call price: " << thisprice << std::endl;
	result.push_back(thisprice);

}

int main() {

	//std::vector<double> K{ 75, 95, 100, 105, 125 };

	EuropeanOption<float> d1(100, 95, 0.02, 0.2, 1);

	d1.CallPricer_BM();
	d1.CallPricer_BSE();

	//vector<double> call = BSE(S, K, r, sig, T);
	//for (auto it = call.begin(); it != call.end(); it++) {
	//	cout << "call price " << distance(call.begin(), it) << ": " << *it << endl;
	//}

	return 0;
}

