#include "Dist.h"
#include <omp.h>
#include <ctime>
#include <thread>
#include <math.h>


int Context::uniform_int_dis(int lowerLimit, int upperLimit, RandomGeneratorType& rg) {
	boost::uniform_int<> distribution(lowerLimit, upperLimit);
	boost::variate_generator< RandomGeneratorType&, boost::uniform_int<> >
		uni_int_s(rg, distribution);
	distribution.reset();
	return uni_int_s();
}


double Context::uniform_real_dis(double lowerLimit, double upperLimit, RandomGeneratorType& rg) {
	boost::uniform_real<> distribution(lowerLimit, upperLimit);
	boost::variate_generator< RandomGeneratorType&, boost::uniform_real<> >	uni_real_s(rg, distribution);
	distribution.reset();
	return uni_real_s();
}
int Context::bernoullie_dis(double p, RandomGeneratorType& rg) {
	boost::bernoulli_distribution<> distribution(p);
	boost::variate_generator< RandomGeneratorType&, boost::bernoulli_distribution<> >
		ber_s(rg, distribution);
	distribution.reset();
	return ber_s();
}

int Context::binomial_dis(int n, double p, RandomGeneratorType& rg) {
	boost::binomial_distribution<> distribution(n, p);
	boost::variate_generator< RandomGeneratorType&, boost::binomial_distribution<> >
		bino_s(rg, distribution);
	distribution.reset();
	return bino_s();
}

long long Context::binomial_dis_long(long long n, double p, RandomGeneratorType& rg) {
	boost::binomial_distribution<> distribution(n, p);
	boost::variate_generator< RandomGeneratorType&, boost::binomial_distribution<> >
		bino_s(rg, distribution);
	distribution.reset();
	return bino_s();
}

double Context::normal_dis(double mean, double S, RandomGeneratorType& rg) {
	boost::normal_distribution<> distribution(mean, S);
	boost::variate_generator< RandomGeneratorType&, boost::normal_distribution<> >
		norm_s(rg, distribution);
	distribution.reset();
	return norm_s();
}

double Context::gamma_dis(double alpha, RandomGeneratorType& rg) {
	boost::gamma_distribution<> distribution(alpha);
	boost::variate_generator< RandomGeneratorType&, boost::gamma_distribution<> >
		gamma_s(rg, distribution);
	distribution.reset();
	return gamma_s();
}

double Context::poisson_dis(int lambda, RandomGeneratorType& rg) {
	boost::poisson_distribution<> distribution(lambda);
	boost::variate_generator< RandomGeneratorType&, boost::poisson_distribution<> >
		pois_s(rg, distribution);
	distribution.reset();
	return pois_s();
}

double Context::beta_dis(double a, double b, RandomGeneratorType& rg) {
	double x = gamma_dis(a, rg);
	double y = gamma_dis(b, rg);
	return x / (x + y);
}

/*double Context::lognormal_dis(double mean, double sigma, RandomGeneratorType& rg) {
	boost::lognormal_distribution<> distribution(mean, sigma);
	boost::variate_generator< RandomGeneratorType&, boost::lognormal_distribution<> >
		log_s(rg, distribution);
	distribution.reset();
	return log_s();
}*/
