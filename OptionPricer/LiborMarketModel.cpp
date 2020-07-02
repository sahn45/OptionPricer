/*
Notes:
http://lesniewski.us/papers/presentations/RiskCourse_LMM.pdf

*/

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
class LMM {
private:
	Eigen::ArrayXd sigCapFloor_;
	Eigen::MatrixXd sigSwaption_;
	Eigen::MatrixXd corr_;
	Eigen::ArrayXd term_;

public:
	LMM(const Eigen::ArrayXd &sigCapFloor, const Eigen::MatrixXd &sigSwaption, const Eigen::MatrixXd &corr, const Eigen::ArrayXd &term) {
		std::cout << "Running constructor..." << std::endl;

		sigCapFloor_ = sigCapFloor;
		sigSwaption_ = sigSwaption;
		corr_ = corr;
		term_ = term;

	}
	~LMM() {
		std::cout << "Running destructor..." << std::endl;
	}

	void LMM_calibrate();
	void LMM_runMonteCarlo();
	
};

template <class T>
void LMM<T>::LMM_calibrate() {
	
	Eigen::MatrixXd kappaLamda(term_.size(), term_.size());
	//kappaLamda.Ones();
	//Eigen::ArrayXd tau(term_.size());
	float zeta;

	//tau(0) = term_(0);
	//tau.tail(term_.size()-1) = term_.tail(term_.size() - 1) - term_.head(term_.size() - 1);
	//sigCapFloor_*sigCapFloor_;//*kappaLamda*tau;*exp(kappaLamda*tau)*tau;


}

