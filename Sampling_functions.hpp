/*
 * Sampling_functions.hpp
 *
 *  Created on: Dec 21, 2018
 *      Author: thanasis
 */
#include <Eigen/Core>
using namespace Eigen;

#ifndef SAMPLING_FUNCTIONS_HPP_
#define SAMPLING_FUNCTIONS_HPP_

double sample_mu(int N, double Esigma2,const VectorXd& Y,const MatrixXd& X,const VectorXd& beta);
double sample_psi2_chisq(const VectorXd& beta,int NZ,double v0B,double s0B);
double sample_sigma_chisq(int N,const VectorXd& epsilon,double v0E,double s0E);
double sample_w(int M,int NZ);

#endif /* SAMPLING_FUNCTIONS_HPP_ */
