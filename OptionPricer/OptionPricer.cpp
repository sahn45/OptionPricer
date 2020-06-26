#include<iostream>
#include<vector>
#include<math.h>
#include<iterator>
#include<algorithm>
#include<boost/math/distributions/normal.hpp>
#include<Eigen/Dense>
#include<Eigen/Sparse>
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
	void CallPricer_FD();

};


template <class T>
void EuropeanOption<T>::CallPricer_FD() {
	std::cout << "Running finite difference call pricer..." << std::endl;

	int timegridcnt = 10;
	int stockgridcnt = 10;
	int upperMultiple = 3;

	float dt = T_ / timegridcnt;
	float ds = S_ * upperMultiple / stockgridcnt;
	float upperPayoff = max(S_*upperMultiple - K_,0), lowerPayoff = 0;

	Eigen::MatrixXd A(stockgridcnt, stockgridcnt);
	Eigen::MatrixXd K(stockgridcnt, 1);
	Eigen::MatrixXd F(stockgridcnt, timegridcnt);
	Eigen::MatrixXd FSolve(stockgridcnt, 1);
	Eigen::Array3d ABlock;

	//This sets up the A matrix to solve the FD equation
	for (int i=0;i<A.rows();i++) {

		ABlock[0] = 0.5*dt*(pow(sig_, 2)*pow(i+1, 2) - r_ * (i+1));
		ABlock[1] = 1 - dt * (pow(sig_, 2)*pow(i+1, 2) + r_);
		ABlock[2] = 0.5*dt*(pow(sig_, 2)*pow(i+1, 2) + r_ * (i+1));

		std::cout << "aBlock: " << std::endl << ABlock << std::endl;
		std::cout << "Head<2> aBlock: " << std::endl << ABlock.head<2>() << std::endl;
		std::cout << "Tail<2> aBlock: " << std::endl << ABlock.tail<2>() << std::endl;

		if (i == 0) {
			A.block(0, 0, 1, 2) = ABlock.tail<2>().transpose();
			std::cout << "i " << i << std::endl << "A: " << std::endl << A << std::endl;
		}
		else if (i == A.rows()-1) {
			A.block(i - 1, i - 2, 1, 2) = ABlock.head<2>().transpose();
			std::cout << "i " << i << std::endl << "A: " << std::endl << A << std::endl;
		}
		else {

			A.block(i - 1, i - 2, 1, 3) = ABlock.transpose();

			std::cout << "i " << i << std::endl << "A: " << std::endl << A << std::endl;

		}
	}

	for (int i=0;i<timegridcnt;i++) {
		if (i == 0) {
			F.col(i) = 0;
		}
		else {
			FSolve = A * F.col(i - 1);
			F.block(0, i, 1, 1) = F.block(0, i, 1, 1) + 0.5*dt*(pow(sig_, 2) - r_)*upperPayoff;
			F.block(F.rows()-1, i, 1, 1) = F.block(F.rows()-1, i, 1, 1) + 0.5*dt*(pow(sig_, 2)*pow((F.rows()-1), 2) + r_ * (F.rows()-1))*lowerPayoff;
		}

	}




}

template <class T>
void EuropeanOption<T>::CallPricer_BM() {
	std::cout << "Running Cox-Ross-Rubinstein Binomial Model call pricer..." << std::endl;

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


	Eigen::MatrixXd AMat(5, 5);
	Eigen::Array3d aBlock;

	float sig = 0.2;
	float r = 0.02;
	float dt = 0.1;

	for (int ia = 1;ia <= AMat.rows();ia++) {
		
		aBlock[0] = 0.5*dt*(pow(sig, 2)*pow(ia, 2) - r * ia);
		aBlock[1] = 1 - dt * (pow(sig, 2)*pow(ia, 2) + r);
		aBlock[2] = 0.5*dt*(pow(sig, 2)*pow(ia, 2) + r * ia);

		std::cout << "aBlock: " << std::endl << aBlock << std::endl;
		std::cout << "Head<2> aBlock: " << std::endl << aBlock.head<2>() << std::endl;
		std::cout << "Tail<2> aBlock: " << std::endl << aBlock.tail<2>() << std::endl;

		if (ia == 1) {
			AMat.block(0,0,1,2) = aBlock.tail<2>().transpose();
			std::cout << "ia " << ia << std::endl << "AMat: " << std::endl << AMat << std::endl;
		}
		else if (ia == AMat.rows()) {
			AMat.block(ia-1, ia-2,1,2) = aBlock.head<2>().transpose();
			std::cout << "ia " << ia << std::endl << "AMat: " << std::endl << AMat << std::endl;
		}
		else {

			AMat.block(ia-1, ia-2, 1, 3) = aBlock.transpose();

			std::cout << "ia " << ia << std::endl << "AMat: " << std::endl << AMat << std::endl;

		}
	}


	//****************************
	Eigen::MatrixXd d1(3, 3);
	Eigen::MatrixXd d2(3, 3);
	Eigen::MatrixXd d3(5, 5);
	Eigen::MatrixXd d4(4, 4);

	d4.setIdentity();

	d1 << 1, 2, 3,
		4, 5, 6,
		7, 8, 9;

	d2 << 1, 1, 1,
		1, 1, 1,	
		1, 1, 1;

	std::cout << d1.col(0) << std::endl;

	d3.bottomLeftCorner(4, 4) = d4;

	std::cout << d3 << std::endl;

	std::cout << d1.cwiseProduct(d1) << std::endl;

	std::cout << d1 << std::endl;

	std::cout << d1 + d2 << std::endl;

	std::cout << d1.inverse() << std::endl;

	std::cout << d1.adjoint() << std::endl;

	std::cout << (d1 + d2).inverse() << std::endl;

	std::cout << d1.size() << std::endl;

	EuropeanOption<float> t1(100, 95, 0.02, 0.2, 1);

	t1.CallPricer_BM();
	t1.CallPricer_BSE();

	//vector<double> call = BSE(S, K, r, sig, T);
	//for (auto it = call.begin(); it != call.end(); it++) {
	//	cout << "call price " << distance(call.begin(), it) << ": " << *it << endl;
	//}

	return 0;
}

