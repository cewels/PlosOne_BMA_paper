#ifndef Dist_H
#define Dist_H

#include <boost/random/mersenne_twister.hpp>
#include <boost/generator_iterator.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/binomial_distribution.hpp>
#include <boost/random/gamma_distribution.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/bernoulli_distribution.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/poisson_distribution.hpp>
#include <boost/random/random_device.hpp>
#include <boost/random/lagged_fibonacci.hpp>
#include <boost/random/lognormal_distribution.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/beta_distribution.hpp>

class Context
{
	
public:
	typedef boost::lagged_fibonacci607  RandomGeneratorType;
		
	int uniform_int_dis(int lowerLimit, int upperLimit, RandomGeneratorType& rg);
	double uniform_real_dis(double lowerLimit, double upperLimit, RandomGeneratorType& rg);
	int bernoullie_dis(double p, RandomGeneratorType& rg);
	int binomial_dis(int n, double p, RandomGeneratorType& rg);
	long long binomial_dis_long(long long n, double p, RandomGeneratorType& rg);
	double normal_dis(double mean, double S, RandomGeneratorType& rg);
	double gamma_dis(double alpha, RandomGeneratorType& rg);
	double poisson_dis(int lambda, RandomGeneratorType& rg);
	double beta_dis(double a, double b, RandomGeneratorType& rg);
	//double lognormal_dis(double mean, double sigma, RandomGeneratorType& rg);
};
#endif