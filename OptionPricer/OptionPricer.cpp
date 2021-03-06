#include<iostream>
#include<vector>
#include<math.h>
#include<iterator>
#include<algorithm>
#include<boost/math/distributions/normal.hpp>
#include<Eigen/Dense>
#include<Eigen/Sparse>
#include<string>
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
	void CallPricer_BSE(const char *callput);
	void CallPricer_BM(const char *callput, const int &gridcount);
	void CallPricer_FD(const char *callput, const char *fdType, const int &timegridcnt, const int &stockgridcnt, const int &upperMultiple);

};


template <class T>
void EuropeanOption<T>::CallPricer_FD(const char *callput, const char *fdType, const int &timegridcnt, const int &stockgridcnt, const int &upperMultiple) {

	std::cout << "Running '" << fdType << "' finite difference call pricer..." << std::endl;

	//int timegridcnt = 10;
	//int stockgridcnt = 10;
	//int upperMultiple = 3;

	float dt = T_ / timegridcnt;
	float ds = S_ * upperMultiple / stockgridcnt;
	float upperPayoff = std::max(S_*upperMultiple - K_, 0.0f), lowerPayoff = 0;


	Eigen::MatrixXd A(stockgridcnt, stockgridcnt);
	A = Eigen::MatrixXd::Zero(stockgridcnt, stockgridcnt);
	Eigen::MatrixXd K(stockgridcnt, 1);
	Eigen::MatrixXd F(stockgridcnt, timegridcnt);
	Eigen::Array3d ABlock;

	if (strcmp(fdType, "explicit") == 0) {
		//This sets up the tridiagonal A matrix to solve the FD equation under Explicit FD method
		for (int i = 0;i < stockgridcnt;i++) {

			ABlock[0] = 0.5*dt*(pow(sig_, 2)*pow(i + 1, 2) - r_ * (i + 1));
			ABlock[1] = 1 - dt * (pow(sig_, 2)*pow(i + 1, 2) + r_);
			ABlock[2] = 0.5*dt*(pow(sig_, 2)*pow(i + 1, 2) + r_ * (i + 1));

			std::cout << "aBlock: " << std::endl << ABlock << std::endl;
			std::cout << "Head<2> aBlock: " << std::endl << ABlock.head<2>() << std::endl;
			std::cout << "Tail<2> aBlock: " << std::endl << ABlock.tail<2>() << std::endl;

			if (i == 0) {
				A.block(0, 0, 1, 2) = ABlock.tail<2>().transpose();
				std::cout << "i " << i << std::endl << "A: " << std::endl << A << std::endl;
			}
			else if (i == A.rows() - 1) {
				A.block(i, i - 1, 1, 2) = ABlock.head<2>().transpose();
				std::cout << "i " << i << std::endl << "A: " << std::endl << A << std::endl;
			}
			else {

				A.block(i, i - 1, 1, 3) = ABlock.transpose();
				std::cout << "i " << i << std::endl << "A: " << std::endl << A << std::endl;

			}
		}

		for (int i = timegridcnt - 1;i >= 0;i--) {
			//First pass sets up payout vector at expiration of option
			if (i == timegridcnt - 1) {
				for (int ia = 0;ia < stockgridcnt;ia++) {
					if (strcmp(callput, "call") == 0) {
						F.col(i)[ia] = std::max((ia + 1)*ds - K_, 0.0f);
					}
					else if (strcmp(callput, "put") == 0) {
						F.col(i)[ia] = std::max(K_ - (ia + 1)*ds, 0.0f);
					}
				}
				std::cout << "F(" << i << "): " << std::endl << F << std::endl;
			}
			//Iterates through time points using Explicit FD formula: F(i-1) = A*F(i) + K(i)
			else {
				std::cout << "F.col(" << i + 1 << "): " << std::endl << F.col(i + 1) << std::endl;
				std::cout << "A: " << std::endl << A << std::endl;
				F.col(i) = A * F.col(i + 1);
				std::cout << "F.col(" << i << "): " << std::endl << F.col(i) << std::endl;
				F(0, i) = F(0, i) + 0.5*dt*(pow(sig_, 2) - r_)*lowerPayoff;
				F(stockgridcnt - 1, i) = F(stockgridcnt - 1, i) + 0.5*dt*(pow(sig_, 2)*pow((stockgridcnt - 1), 2) + r_ * (stockgridcnt - 1))*upperPayoff;
				std::cout << "F(" << i << "): " << std::endl << F << std::endl;
			}

		}
	}
	else if (strcmp(fdType, "implicit") == 0) {

		Eigen::MatrixXd tempM(stockgridcnt, 1);

		//This sets up the tridiagonal A matrix to solve the FD equation under Implicit FD method
		for (int i = 0;i < stockgridcnt;i++) {

			ABlock[0] = 0.5*dt*(-pow(sig_, 2)*pow(i + 1, 2) + r_ * (i + 1));
			ABlock[1] = 1 + dt * (pow(sig_, 2)*pow(i + 1, 2) + r_);
			ABlock[2] = 0.5*dt*(-pow(sig_, 2)*pow(i + 1, 2) - r_ * (i + 1));

			std::cout << "aBlock: " << std::endl << ABlock << std::endl;
			std::cout << "Head<2> aBlock: " << std::endl << ABlock.head<2>() << std::endl;
			std::cout << "Tail<2> aBlock: " << std::endl << ABlock.tail<2>() << std::endl;

			if (i == 0) {
				A.block(0, 0, 1, 2) = ABlock.tail<2>().transpose();
				std::cout << "i " << i << std::endl << "A: " << std::endl << A << std::endl;
			}
			else if (i == A.rows() - 1) {
				A.block(i, i - 1, 1, 2) = ABlock.head<2>().transpose();
				std::cout << "i " << i << std::endl << "A: " << std::endl << A << std::endl;
			}
			else {

				A.block(i, i - 1, 1, 3) = ABlock.transpose();
				std::cout << "i " << i << std::endl << "A: " << std::endl << A << std::endl;

			}
		}

		for (int i = timegridcnt - 1; i >= 0; i--) {
			//First pass sets up payout vector at expiration of option
			if (i == timegridcnt - 1) {
				for (int ia = 0;ia < stockgridcnt;ia++) {
					if (strcmp(callput, "call") == 0) {
						F.col(i)[ia] = std::max((ia + 1)*ds - K_, 0.0f);
					}
					else if (strcmp(callput, "put") == 0) {
						F.col(i)[ia] = std::max(K_ - (ia + 1)*ds, 0.0f);
					}
				}
				std::cout << "F(" << i << "): " << std::endl << F << std::endl;
			}
			//Iterates through time points using Implicit FD formula: F(i) = A^(-1)*[F(i+1) + K(i)]
			else {
				std::cout << "F.col(" << i + 1 << "): " << std::endl << F.col(i + 1) << std::endl;
				std::cout << "A: " << std::endl << A << std::endl;
				tempM = F.col(i + 1);
				tempM(0) = F(0, i + 1) - 0.5*dt*(-pow(sig_, 2) + r_)*lowerPayoff;
				tempM(stockgridcnt - 1) = F(stockgridcnt - 1, i + 1) - 0.5*dt*(-pow(sig_, 2)*pow((stockgridcnt - 1), 2) - r_ * (stockgridcnt - 1))*upperPayoff;
				F.col(i) = A.inverse()*tempM;
				std::cout << "F(" << i << "): " << std::endl << F << std::endl;

			}

		}
	}
	else if (strcmp(fdType, "crank-nicolson") == 0) {

		Eigen::MatrixXd tempM(stockgridcnt, 1);
		Eigen::MatrixXd B(stockgridcnt, stockgridcnt);
		B = Eigen::MatrixXd::Zero(stockgridcnt, stockgridcnt);
		Eigen::Array3d BBlock;

		std::cout << B << std::endl;

		//This sets up the tridiagonal A matrix to solve the FD equation under Crank-Nicolson method
		for (int i = 0;i < A.rows();i++) {

			ABlock[0] = -0.25*dt*(pow(sig_, 2)*pow(i + 1, 2) - r_ * (i + 1));
			ABlock[1] = 1 + 0.5*dt * (pow(sig_, 2)*pow(i + 1, 2) + r_);
			ABlock[2] = -0.25*dt*(pow(sig_, 2)*pow(i + 1, 2) + r_ * (i + 1));

			std::cout << "aBlock: " << std::endl << ABlock << std::endl;
			std::cout << "Head<2> aBlock: " << std::endl << ABlock.head<2>() << std::endl;
			std::cout << "Tail<2> aBlock: " << std::endl << ABlock.tail<2>() << std::endl;

			BBlock[0] = 0.25*dt*(pow(sig_, 2)*pow(i + 1, 2) - r_ * (i + 1));
			BBlock[1] = 1 - 0.5*dt * (pow(sig_, 2)*pow(i + 1, 2) + r_);
			BBlock[2] = 0.25*dt*(pow(sig_, 2)*pow(i + 1, 2) + r_ * (i + 1));

			std::cout << "bBlock: " << std::endl << BBlock << std::endl;
			std::cout << "Head<2> bBlock: " << std::endl << BBlock.head<2>() << std::endl;
			std::cout << "Tail<2> bBlock: " << std::endl << BBlock.tail<2>() << std::endl;

			if (i == 0) {
				A.block(0, 0, 1, 2) = ABlock.tail<2>().transpose();
				std::cout << "i " << i << std::endl << "A: " << std::endl << A << std::endl;

				B.block(0, 0, 1, 2) = BBlock.tail<2>().transpose();
				std::cout << "i " << i << std::endl << "B: " << std::endl << B << std::endl;
			}
			else if (i == stockgridcnt - 1) {
				A.block(i, i - 1, 1, 2) = ABlock.head<2>().transpose();
				std::cout << "i " << i << std::endl << "A: " << std::endl << A << std::endl;

				B.block(i, i - 1, 1, 2) = BBlock.head<2>().transpose();
				std::cout << "i " << i << std::endl << "B: " << std::endl << B << std::endl;
			}
			else {

				A.block(i, i - 1, 1, 3) = ABlock.transpose();
				std::cout << "i " << i << std::endl << "A: " << std::endl << A << std::endl;

				B.block(i, i - 1, 1, 3) = BBlock.transpose();
				std::cout << "i " << i << std::endl << "B: " << std::endl << B << std::endl;

			}
		}

		for (int i = timegridcnt - 1;i >= 0;i--) {
			//First pass sets up payout vector at expiration of option
			if (i == timegridcnt - 1) {
				for (int ia = 0;ia < stockgridcnt;ia++) {
					if (strcmp(callput, "call") == 0) {
						F.col(i)[ia] = std::max((ia + 1)*ds - K_, 0.0f);
					}
					else if (strcmp(callput, "put") == 0) {
						F.col(i)[ia] = std::max(K_ - (ia + 1)*ds, 0.0f);
					}
				}
				std::cout << "F(" << i << "): " << std::endl << F << std::endl;
			}
			//Iterates through time points using Implicit FD formula: F(i) = A^(-1)*[F(i+1) + K(i)]
			else {
				std::cout << "F.col(" << i + 1 << "): " << std::endl << F.col(i + 1) << std::endl;
				std::cout << "A: " << std::endl << A << std::endl;
				std::cout << "B: " << std::endl << B << std::endl;
				tempM = B * F.col(i + 1);
				//std::cout << "tempM: " << std::endl << tempM << std::endl;
				tempM(0) = F(0, i + 1) + 2 * 0.25*dt*(pow(sig_, 2) - r_)*lowerPayoff;
				//std::cout << "tempM: " << std::endl << tempM << std::endl;
				tempM(stockgridcnt - 1) = F(stockgridcnt - 1, i + 1) + 2 * 0.25*dt*(pow(sig_, 2)*pow((stockgridcnt - 1), 2) + r_ * (stockgridcnt - 1))*upperPayoff;
				//std::cout << "tempM: " << std::endl << tempM << std::endl;
				F.col(i) = A.inverse()*tempM;
				std::cout << "F(" << i << "): " << std::endl << F << std::endl;

			}

		}

	}
}


template <class T>
void EuropeanOption<T>::CallPricer_BM(const char *callput, const int &gridcount) {
	std::cout << "Running Cox-Ross-Rubinstein Binomial Model '" << callput << "' pricer..." << std::endl;

	std::vector<std::vector<float>> stock_matrix;
	std::vector<std::vector<float>> price_matrix;
	std::vector<float> temp_vec;

	//int gridcount = 50;
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
				if (strcmp(callput, "call") == 0) {
					price_matrix[ia][ib] = std::max(stock_matrix[ia][ib] - K_, 0.0f);
				}
				else if (strcmp(callput, "put") == 0) {
					price_matrix[ia][ib] = std::max(K_ - stock_matrix[ia][ib], 0.0f);
				}
				std::cout << "'" << callput << "' price [" << ia << "][" << ib << "]: " << price_matrix[ia][ib] << std::endl;
			}
			else {
				price_matrix[ia][ib] = price_matrix[ia - 1][ib] * q + price_matrix[ia - 1][ib + 1] * (1 - q);
				std::cout << "'" << callput << "' price: [" << ia << "][" << ib << "]: " << price_matrix[ia][ib] << std::endl;
			}
		}
	}
	std::cout << "FINAL call price : " << price_matrix[gridcount][0] << std::endl;
}

template <class T>
void EuropeanOption<T>::CallPricer_BSE(const char *callput) {
	std::cout << "Running BSE '" << callput << "' pricer..." << std::endl;
	normal sn(0, 1);

	float d1 = (log(S_ / K_) + (r_ + 0.5*pow(sig_, 2)*T_)) / (sig_*pow(T_, 0.5));
	float d2 = d1 - sig_ * pow(T_, 0.5);
	float  thisprice;

	std::cout << "S: " << S_ << std::endl;
	std::cout << "K: " << K_ << std::endl;
	std::cout << "r: " << r_ << std::endl;
	std::cout << "sig: " << sig_ << std::endl;
	std::cout << "T: " << T_ << std::endl;
	std::cout << "d1: " << d1 << std::endl;
	std::cout << "d2: " << d2 << std::endl;
	std::cout << "cdf(d1): " << cdf(sn, d1) << std::endl;
	std::cout << "cdf(d2): " << cdf(sn, d2) << std::endl;

	if (strcmp(callput, "call")) {
		thisprice = S_ * cdf(sn, d1) - K_ * exp(-r_ * T_)*cdf(sn, d2);
	}
	else if (strcmp(callput, "put")) {
		thisprice = K_ * exp(-r_ * T_)*cdf(sn, -d2) - S_ * cdf(sn, -d1);
	}
	std::cout << "BSE '" << callput << "' price: " << thisprice << std::endl;

}



