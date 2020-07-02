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
class SABRModel {
private:
	T alpha_, beta_, nu_, rho_;

public:
	SABRModel(const std::vector<T> mktQuotes) {
		std::cout << "Running constructor..." << std::endl;


	}


};