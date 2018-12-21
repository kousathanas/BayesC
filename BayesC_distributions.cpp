#include <math.h>
#include <random>
#include <boost/random.hpp>
#include <boost/math/distributions.hpp>
boost::random::mt19937 gen(time(0));

//distributions
double runif(double lower, double higher)
{
	boost::random::uniform_real_distribution<> dist(lower, higher);
	return dist(gen);
}

double rnorm(double mean, double sd)
{
	boost::random::normal_distribution<> dist(mean, sd);
	return dist(gen);
}


double rbeta(double alpha, double beta)
{

	boost::math::beta_distribution<> dist(alpha, beta);
	double q = quantile(dist, runif(0,1));

	return(q);
}

double rinvchisq(double df, double scale)
{

	boost::math::inverse_chi_squared_distribution<> dist(df, scale);
	double q = quantile(dist, runif(0,1));

	return(q);
}
int rbernoulli(double p)
{
	std::bernoulli_distribution dist(p);
	return dist(gen);
}


int rbinom(int k, double p)
{

	std::binomial_distribution<> dist(k,p);
	return dist(gen);
}
