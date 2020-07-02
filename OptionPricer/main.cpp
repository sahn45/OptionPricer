#include<iostream>
#include<vector>
#include<math.h>
#include<iterator>
#include<algorithm>
#include<boost/math/distributions/normal.hpp>
#include<Eigen/Dense>
#include<Eigen/Sparse>
#include<string>
#include<C:\Users\sahn4\source\repos\OptionPricer\OptionPricer\OptionPricer.cpp>
#include<C:\Users\sahn4\source\repos\OptionPricer\OptionPricer\LiborMarketModel.cpp>
#include<C:\Users\sahn4\source\repos\OptionPricer\OptionPricer\SABRModel.cpp>

using boost::math::normal;

int main() {


	Eigen::ArrayXd capfloorSig(5);
	capfloorSig << 0.2, 0.22, 0.23, 0.25, 0.26;
	Eigen::MatrixXd swaptionSig(5, 5);
	swaptionSig << 0.2, 0.2, 0.2, 0.2, 0.2,
		0.2, 0.2, 0.2, 0.2, 0.2,
		0.2, 0.2, 0.2, 0.2, 0.2,
		0.2, 0.2, 0.2, 0.2, 0.2,
		0.2, 0.2, 0.2, 0.2, 0.2;
	Eigen::MatrixXd corr(5, 5);
	corr << 1, 0.98, 0.98, 0.98, 0.98,
		0.98, 1, 0.98, 0.98, 0.98,
		0.98, 0.98, 1, 0.98, 0.98,
		0.98, 0.98, 0.98, 1, 0.98,
		0.98, 0.98, 0.98, 0.98, 1;
	Eigen::ArrayXd term(5);
	term << 1, 2, 3, 4, 5;

	LMM<float> l1(capfloorSig, swaptionSig, corr, term);
	l1.LMM_calibrate();
	//l1.LMM_runMonteCarlo();


	/*
	EuropeanOption<float> t1(100, 95, 0.02, 0.2, 1);

	t1.CallPricer_FD("call", "explicit", 10, 10, 3);
	t1.CallPricer_FD("put", "explicit", 10, 10, 3);
	t1.CallPricer_FD("call", "implicit", 10, 10, 3);
	t1.CallPricer_FD("put", "implicit", 10, 10, 3);
	t1.CallPricer_FD("call", "crank-nicolson", 10, 10, 3);
	t1.CallPricer_FD("put", "crank-nicolson", 10, 10, 3);
	t1.CallPricer_BM("call", 50);
	t1.CallPricer_BM("put", 50);
	t1.CallPricer_BSE("call");
	t1.CallPricer_BSE("put");


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
			AMat.block(0, 0, 1, 2) = aBlock.tail<2>().transpose();
			std::cout << "ia " << ia << std::endl << "AMat: " << std::endl << AMat << std::endl;
		}
		else if (ia == AMat.rows()) {
			AMat.block(ia - 1, ia - 2, 1, 2) = aBlock.head<2>().transpose();
			std::cout << "ia " << ia << std::endl << "AMat: " << std::endl << AMat << std::endl;
		}
		else {

			AMat.block(ia - 1, ia - 2, 1, 3) = aBlock.transpose();

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


	//vector<double> call = BSE(S, K, r, sig, T);
	//for (auto it = call.begin(); it != call.end(); it++) {
	//	cout << "call price " << distance(call.begin(), it) << ": " << *it << endl;
	//}

	*/


	return 0;
}

