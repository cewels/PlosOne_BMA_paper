#include "Eigen/Dense"
#include "Eigen/LU"
# include <cstdlib>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>
#include <omp.h>
#include "Dist.h"
#include <iterator>
using namespace Eigen;
using namespace std;
#include <boost/random/lagged_fibonacci.hpp>

typedef boost::random::lagged_fibonacci607 base_generator_type;

struct symptom{
	int Fam_ID;
	int ID;
	int p_id;
	int m_id;
	int sex;
	int t_pheno;
	int sym[12];
};

void dirichlet_sample(int n, double a[], double *(&out), base_generator_type& eng);
void multi_boost(int n, double p[], int ncat, int *(&ix), base_generator_type& eng);
double ln_binomial(double lambda, int y);
void inv_Hessian_LCA(double *p_star, double **lambda, int no_ind, int Nsym, int n_k, struct symptom *dat);
double ln_lca(struct symptom *dat, double *Theta, int size_of_Theta, int no_ind, int Nsym, int nk);
double DIC3_LC(double mean_log_dev, double **p_star, double ***lambda, int NUM_KEEP, int K, int Nsym, int no_ind, struct symptom *dat);
double DIC3_GoM(double mean_log_dev, double ***gik, double ***gamma, int NUM_KEEP, int K, int Nsym, int no_ind, struct symptom *dat);
