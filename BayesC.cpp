#include <iostream>
#include <fstream>
#include <stdio.h>
#include <math.h>
#include <omp.h>

#include <random>
#include <map>
#include <string>
#include <iomanip>

#include <unistd.h>
#include <string>
#include <algorithm>
#include <random>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <boost/random.hpp>
#include <boost/math/distributions.hpp>
#include <boost/math/distributions/inverse_chi_squared.hpp>


using namespace std;
using namespace Eigen;

using Eigen::Matrix;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::LLT;
using Eigen::Lower;
using Eigen::Map;
using Eigen::Upper;
typedef Map<MatrixXd> MapMatd;

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
	//cout<<df<<"\t"<<scale<<endl;
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

int main(int argc, char *argv[])
{
	int M=25000;
	int N=5000;
	int MT=20000;
	int i,j,k,l,m=0;
	double sigmaY_true=1;
	double sigmab_true=1;
	//Genotype matrix
	MatrixXd X(N,M);
	for (i=0;i<N;i++){
		for (j=0;j<M;j++){
			X(i,j)=rnorm(0,1);
		}
	}
	//beta coefficients
	VectorXd beta_true(M);
	for (i=0;i<M;i++){
		beta_true[i]=rnorm(0,sigmab_true);
	}
	//set some betas to zero
	for (i=0;i<MT;i++){
		beta_true[i]=0.0;

	}

	//error
	VectorXd error(N);
	for (i=0;i<N;i++){
		error[i]=rnorm(0,sigmaY_true);
	}

	//construct phenotypes
	VectorXd Y(N);
	Y=X*beta_true;
	Y+=error;

	//normalize
	RowVectorXd mean = X.colwise().mean();
	RowVectorXd sd = (X.rowwise() - mean).array().square().colwise().mean();
	X = (X.rowwise() - mean).array().rowwise() / sd.array();
	//cout<<"OK"<<endl;

	int iter=5000;
	int burnin=2500;

	double Emu=0;
	VectorXd vEmu(N);
	vEmu.setZero();
	VectorXd Ebeta(M);
	Ebeta.setZero();
	VectorXd ny(M);
	ny.setZero();
	double Ew=0.5;
	//residual error
	VectorXd epsilon(N);

	epsilon=Y-X*Ebeta-vEmu;

	std::vector<int> markerI;
	for (int i=0; i<M; ++i) {
		markerI.push_back(i);
	}
	int marker=0;

	//non-zero variable NZ
	int NZ=0;

	double Esigma2=epsilon.squaredNorm()/(N*0.5);
	double Epsi2=rbeta(1,1);

	//Standard parameterization of hyperpriors for variances
	double v0E,s0E,v0B,s0B=0.001;


	 /* Alternative parameterization of hyperpriors for variances
				double v0E,v0B=4;
				double s0B=((v0B-2)/v0B)*Epsi2;
				double s0E=((v0E-2)/v0E)*Esigma2;
*/

	//pre-computed elements for calculations
	VectorXd el1(M);
	VectorXd el2(M);
	for (int i=0; i<M; ++i) {
	el1[i]=X.col(i).transpose()*X.col(i);
	el2[i]=X.col(i).transpose();
	}
	//begin GIBBS sampling iterations

	for (i=0;i<iter;i++){

		Emu=sample_mu(N,Esigma2,Y,X,Ebeta);

		//sample effects and probabilities jointly
		std::random_shuffle(markerI.begin(), markerI.end());
		for (j=0;j<M;j++){
			marker=markerI[j];

			epsilon=epsilon+X.col(j)*Ebeta[j];

			double Cj=el1[j]+Esigma2/Epsi2;
			double rj=el2[i]*epsilon;

			epsilon=epsilon-X.col(j)*Ebeta[j];

			double ratio=(((exp(-(pow(rj,2))/(2*Cj*Esigma2))*sqrt((Epsi2*Cj)/Esigma2))));
			ratio=Ew/(Ew+ratio*(1-Ew));
			ny[j]=rbernoulli(ratio);

			if (ny[j]==0){
				Ebeta[j]=0;
			}
			if (ny[j]==1){
				Ebeta[j]=rnorm(rj/Cj,Esigma2/Cj);
			}

		}
		NZ=ny.sum();
		//cout<<NZ<<endl;

		Ew=sample_w(M,NZ);
		epsilon=Y-X*Ebeta;

		Epsi2=sample_psi2_chisq(Ebeta,NZ,v0B,s0B);
		Esigma2=sample_sigma_chisq(N,epsilon,v0E,s0E);

	}


	//write to files
	ofstream myfile1;
	myfile1.open ("Y.txt");
	myfile1 << Y << '\n';
	myfile1.close();
/*
	ofstream myfile2;
	myfile2.open ("X.txt");
	myfile2 << X << '\n';
	myfile2.close();
*/
	ofstream myfile3;
	myfile3.open ("beta.txt");
	myfile3 << beta_true << '\n';
	myfile3.close();

	ofstream myfile4;
	myfile4.open ("Ebeta.txt");
	myfile4 << Ebeta << '\n';
	myfile4.close();

	return 0;
}
