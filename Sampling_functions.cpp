#include "BayesC_distributions.h"
#include <Eigen/Core>
using namespace Eigen;

//sampling functions
double sample_mu(int N, double Esigma2,const VectorXd& Y,const MatrixXd& X,const VectorXd& beta)
{
	double mean=((Y-X*beta).sum())/N;
	double sd=sqrt(Esigma2/N);
	double mu=rnorm(mean,sd);
	return(mu);
}

//sample variance of beta
double sample_psi2_chisq(const VectorXd& beta,int NZ,double v0B,double s0B){
	double df=v0B+NZ;
	double scale=(beta.squaredNorm()*NZ+v0B*s0B)/(v0B+NZ);
	double psi2=rinvchisq(df, scale);
	return(psi2);
}

//sample error variance of Y
double sample_sigma_chisq(int N,const VectorXd& epsilon,double v0E,double s0E){
	double sigma2=rinvchisq(v0E+N, (epsilon.squaredNorm()+v0E*s0E)/(v0E+N));
	return(sigma2);
}

//sample mixture weight
double sample_w(int M,int NZ){
	double w=rbeta(1+NZ,1+(M-NZ));
	return(w);
}
