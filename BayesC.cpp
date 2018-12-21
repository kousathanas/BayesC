#include <iostream>
#include <fstream>
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <chrono>

#include <random>
#include <map>
#include <string>
#include <iomanip>

#include <unistd.h>
#include <string>
#include <algorithm>
#include <random>

#include <Eigen/Core>

#include <boost/program_options.hpp>
#include <iterator>

#include "BayesC_distributions.h"
#include "Sampling_functions.hpp"

using namespace std;
using namespace Eigen;
namespace po = boost::program_options;

void ReadFromFile(std::vector<double> &x, const std::string &file_name)
{
	std::ifstream read_file(file_name);
	assert(read_file.is_open());

	std::copy(std::istream_iterator<double>(read_file), std::istream_iterator<double>(),
			std::back_inserter(x));

	read_file.close();
}


int main(int argc, char *argv[])
{

	po::options_description desc("Options");
	desc.add_options()
		("M", po::value<int>()->required(), "No. of simulated markers")
		("N", po::value<int>()->required(), "No. of simulated individuals")
		("iter", po::value<int>()->default_value(5000), "No. of Gibbs iterations")
		("pNZ", po::value<double>()->default_value(0.5), "Proportion nonzero (simulations)")
		("h2", po::value<double>()->default_value(0.6), "Heritability (simulations)")
		("input", po::value<std::string>()->default_value("none"),"Input filename")
		("out", po::value<std::string>()->default_value("BayesC_out"),"Output filename")
	;

	//clock starts
	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

	//map variables
	po::variables_map vm;
	po::store(po::parse_command_line(argc,argv,desc),vm);
	po::notify(vm);

	int M=vm["M"].as<int>();
	int N=vm["N"].as<int>();
	int iter=vm["iter"].as<int>();

	string input=vm["input"].as<string>();
	string output=vm["out"].as<string>();

	MatrixXd X(N,M);
	VectorXd Y(N);

	//beta coefficients
	VectorXd beta_true(M);
	beta_true.setZero();

	int i,j,k,l,m=0;

	//Was an input matrix given?

	if (input!="none"){ //Either read input tables for X and Y
		ifstream f1(input+".X");
		//f1 >> m >> n;
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < M; j++)
			{
				f1 >> X(i,j);
				//cout<<X(i,j)<<endl;
			}
		}
		f1.close();
		cout<<"finished reading matrix X!"<<endl;

		std::vector<double> Y_in;
		ReadFromFile(Y_in, input+".Y");
		double* ptr_Y = &Y_in[0];
		Eigen::Map<Eigen::VectorXd> Y1(ptr_Y, Y_in.size());
		cout<<"finished reading vector Y!"<<endl;
		if(Y_in.size()!=N){cout<<"input Y vector size doesnt much the size indicated in the command line"<<endl;return 0;}
		Y=Y1;
		Y_in.clear();


	}else //or simulate
	{
		double h2=vm["h2"].as<double>();
		double pNZ=vm["pNZ"].as<double>();
		double sigmaY_true=1;

		int MT=pNZ*M;

		//Fill Genotype matrix
		for (i=0;i<N;i++){
			for (j=0;j<M;j++){
				X(i,j)=rnorm(0,1);
			}
		}
		for (i=0;i<MT;i++){
			beta_true[i]=rnorm(0,sqrt(h2/MT));
		}

		//error
		VectorXd error(N);
		for (i=0;i<N;i++){
			error[i]=rnorm(0,sqrt(1-h2));
		}

		//construct phenotypes
		Y=X*beta_true;
		Y+=error;
	}
	//standardize matrix X
	RowVectorXd mean = X.colwise().mean();
	RowVectorXd sd = ((X.rowwise() - mean).array().square().colwise().sum() / (X.rows() - 1)).sqrt();
	X = (X.rowwise() - mean).array().rowwise() / sd.array();

	//standardize vector Y
    Y = (Y.array() - Y.array().mean());
    Y /= sqrt(Y.squaredNorm() / (double(N - 1)));

    //Initialize variables
	double Emu=0;
	VectorXd vEmu(N);
	vEmu.setOnes();

	VectorXd Ebeta(M);
	Ebeta.setZero();
	VectorXd ny(M);
	ny.setZero();
	double Ew=0.5;
	//residual error
	VectorXd epsilon(N);

	epsilon=Y-X*Ebeta-vEmu*Emu;

	std::vector<int> markerI;
	for (int i=0; i<M; ++i) {
		markerI.push_back(i);
	}
	int marker=0;


	int NZ=0;

	double Esigma2=epsilon.squaredNorm()/(N*0.5);
	double Epsi2=rbeta(1,1);

	//Standard parameterization of hyperpriors for variances
	//double v0E=0.001,s0E=0.001,v0B=0.001,s0B=0.001;


	// Alternative parameterization of hyperpriors for variances
	double v0E=4,v0B=4;
	double s0B=((v0B-2)/v0B)*Epsi2;
	double s0E=((v0E-2)/v0E)*Esigma2;


	//pre-computed elements for calculations
	VectorXd el1(M);
	for (int i=0; i<M; ++i) {
		el1[i]=X.col(i).transpose()*X.col(i);
	}

	//open files for writing
	std::ofstream ofs;
	ofs.open(output+"_estimates.txt");
	for (int i=0; i<M; ++i) {
		ofs << "beta_" <<i<< ' ';
	}
	for (int i=0; i<M; ++i) {
		ofs << "incl_" <<i<< ' ';
	}
	ofs << "Ew" << " ";
	ofs << "Epsi2" << " ";
	ofs << "Esigma2" << " ";
	ofs << "\n";
	ofs.close();

	std::chrono::steady_clock::time_point end1= std::chrono::steady_clock::now();
	std::cout << "Time taken for Reading/generating data = " << std::chrono::duration_cast<std::chrono::nanoseconds> (end1 - begin).count()*1e-9 <<" seconds"<<std::endl;

	//begin GIBBS sampling iterations

	ofs.open (output+"_estimates.txt", std::ios_base::app);
	for (i=0;i<iter;i++){

		Emu=sample_mu(N,Esigma2,Y,X,Ebeta);

		//sample effects and probabilities jointly
		std::random_shuffle(markerI.begin(), markerI.end());

		for (j=0;j<M;j++){
			marker=markerI[j];

			epsilon=epsilon+X.col(marker)*Ebeta[marker];

			double Cj=el1[marker]+Esigma2/Epsi2;
			double rj=X.col(marker).transpose()*epsilon;

			double ratio=(((exp(-(pow(rj,2))/(2*Cj*Esigma2))*sqrt((Epsi2*Cj)/Esigma2))));
			ratio=Ew/(Ew+ratio*(1-Ew));
			ny[marker]=rbernoulli(ratio);

			if (ny[marker]==0){
				Ebeta[marker]=0;
			}
			else if (ny[marker]==1){
				Ebeta[marker]=rnorm(rj/Cj,Esigma2/Cj);
			}

			epsilon=epsilon-X.col(marker)*Ebeta[marker];

		}
		for (j=0;j<M;j++){
			ofs << Ebeta[j] << " ";
		}
		for (j=0;j<M;j++){
			ofs << ny[j] << " ";
		}
		NZ=ny.sum();
		//cout<<NZ<<endl;

		Ew=sample_w(M,NZ);
		epsilon=Y-X*Ebeta-vEmu*Emu;

		Epsi2=sample_psi2_chisq(Ebeta,NZ,v0B,s0B);
		Esigma2=sample_sigma_chisq(N,epsilon,v0E,s0E);

		ofs << Ew << " ";
		ofs << Epsi2 << " ";
		ofs << Esigma2 << " ";
		ofs << "\n";

	}
	ofs.close();
//write out simulated data
if (input=="none"){
	//write to files
	ofstream myfile1;
	myfile1.open (output+"_simulated_Y.txt");
	for (i=0;i<N;i++){
		myfile1 << Y[i] << ' ';
	}
	myfile1 << endl;
	myfile1.close();
/*
	ofstream myfile2;
	myfile2.open (output+"_simulated_X.txt");
	for (i=0;i<N;i++){
		for (j=0;j<M;j++){
			myfile2<<X(i,j)<< ' ';
		}
		myfile2<<endl;
	}
	myfile2.close();
*/
	ofstream myfile3;
	myfile3.open (output+"_simulated_betatrue.txt");
	myfile3 << beta_true << ' ';
	myfile3.close();
}


std::chrono::steady_clock::time_point end2= std::chrono::steady_clock::now();
std::cout << "Time taken for full analysis = " << std::chrono::duration_cast<std::chrono::nanoseconds> (end2 - begin).count()*1e-9 <<" seconds"<<std::endl;

	return 0;
}
