/*
 * BayesC_distributions.h
 *
 *  Created on: Aug 15, 2018
 *      Author: Athanasios Kousathanas
 */

#ifndef BAYESC_DISTRIBUTIONS_H_
#define BAYESC_DISTRIBUTIONS_H_

//distributions
double runif(double lower, double higher);
double rnorm(double mean, double sd);
double rbeta(double alpha, double beta);
double rinvchisq(double df, double scale);
int rbernoulli(double p);
int rbinom(int k, double p);



#endif /* BAYESC_DISTRIBUTIONS_H_ */
