#include "otherDist.h"

void dirichlet_sample(int n, double a[], double *(&out), base_generator_type& eng)

//****************************************************************************80
//
//  Purpose:
//
//    DIRICHLET_SAMPLE samples the Dirichlet PDF-based on gamma

//
//  Author:
//
//    me
//
//  Reference:
//
{
	int i;
	double sum;
	Context dist;

	Map<RowVectorXd> v(a, n);
	VectorXd v1(n);

	for (i = 0; i < n; i++)
	{
		v1(i) = dist.gamma_dis(v(i), eng);

	}

	sum = v1.sum();
	for (i = 0; i < n; i++)
	{
		out[i] = v1(i) / sum;
	}
	return;
}
void multi_boost(int n, double p[], int ncat, int *(&ix), base_generator_type& eng)

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_MULTINOMIAL_SAMPLE generates a multinomial random deviate.
//
// modified, using boost
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 May 2013
//
//  Author:
//
//    Original FORTRAN77 version by Barry Brown, James Lovato.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Luc Devroye,
//    Non-Uniform Random Variate Generation,
//    Springer, 1986,
//    ISBN: 0387963057,
//    LC: QA274.D48.
//
//  Parameters:
//
//    Input, int N, the number of events, which will be
//    classified into one of the NCAT categories.
//
//    Input, double P[NCAT-1].  P(I) is the probability that an event
//    will be classified into category I.  Thus, each P(I) must be between 
//    0.0 and 1.0.  Only the first NCAT-1 values of P must be defined since 
//    P(NCAT) would be 1.0 minus the sum of the first NCAT-1 P's.
//
//    Input, int NCAT, the number of categories.
//
//    Output, int I4VEC_MULTINOMIAL_SAMPLE[NCAT], a random observation from 
//    the multinomial distribution.  All IX(i) will be nonnegative and their 
//    sum will be N.
//
{
	int i;
	int icat;
	int ntot;
	double prob;
	double ptot;
	Context dist;

	if (n < 0)
	{
		cerr << "\n";
		cerr << "I4VEC_MULTINOMIAL_SAMPLE - Fatal error!\n";
		cerr << "  N < 0\n";
		exit(1);
	}

	if (ncat <= 1)
	{
		cerr << "\n";
		cerr << "I4VEC_MULTINOMIAL_SAMPLE - Fatal error!\n";
		cerr << "  NCAT <= 1\n";
		exit(1);
	}

	for (i = 0; i < ncat - 1; i++)
	{
		if (p[i] < 0.0)
		{
			cerr << "\n";
			cerr << "I4VEC_MULTINOMIAL_SAMPLE - Fatal error!\n";
			cerr << "  Some P(i) < 0.\n";
			exit(1);
		}

		if (1.0 < p[i])
		{
			cerr << "\n";
			cerr << "I4VEC_MULTINOMIAL_SAMPLE - Fatal error!\n";
			cerr << "  Some 1 < P(i).\n";
			exit(1);
		}
	}

	ptot = 0.0;
	for (i = 0; i < ncat - 1; i++)
	{
		ptot = ptot + p[i];
	}

	if (1.0 < ptot)
	{
		cerr << "\n";
		cerr << "I4VEC_MULTINOMIAL_SAMPLE - Fatal error!\n";
		cerr << "  1.0 < Sum of P().\n";
		exit(1);
	}
	//
	//  Initialize variables.
	//
	ntot = n;
	ptot = 1.0;

	for (i = 0; i < ncat; i++)
	{
		ix[i] = 0;
	}
	//
	//  Generate the observation.
	//


	for (icat = 0; icat < ncat - 1; icat++)
	{
		prob = p[icat] / ptot;
		ix[icat] = dist.binomial_dis(ntot, prob, eng);
		ntot = ntot - ix[icat];
		if (ntot <= 0)
		{
			return;
		}
		ptot = ptot - p[icat];
	}

	ix[ncat - 1] = ntot;

	return;
}
double ln_binomial(double lambda, int y){

	if (lambda > 1){
		cout << "Lambda is greater than 1" << endl;
		exit(1);
	}
	else{
		if (y == 0){
			return log(1 - lambda);
		}
		else{
			return log(lambda);
		}
	}
}

double ln_lca(struct symptom *dat, double *Theta, int size_of_Theta, int no_ind, int Nsym, int nk){
	// work out the loglike of lca
	int i, j, k;
	double out = 0;
	double temp;
	double *p = new double[nk];// store p_k
	double **lambda;
	lambda = new double *[nk];

	for (k = 0; k < nk; k++){
		lambda[k] = new double[Nsym];
	}

	k = 0;
	j = 0;

	for (i = 0; i < size_of_Theta; i++){
		if (i < nk){
			p[i] = Theta[i];
		}
		else{// lambda space
			lambda[k][j] = Theta[i];
			j++;
			if (j > Nsym){
				j = 0;
				k++;
			}
		}
	}

	for (i = 0; i < no_ind; i++){
		for (j = 0; j < Nsym; j++){
			temp = 0;
			for (k = 0; k < nk; k++){
				if (dat[i].sym[j] == 1){
					temp += p[k] * lambda[k][j];
				}
				else{
					temp += p[k] * (1 - lambda[k][j]);
				}
			}
			out += log(temp);
		}
	}

	return(out);
}
void inv_Hessian_LCA(double *p_star, double **lambda, int no_ind, int Nsym, int n_k, struct symptom *dat){

	// aim: work out Hessian for LCA
	// need to find function to find inverse
	int npara = n_k + n_k*Nsym;
	MatrixXd out(npara, npara);
	int i, j, k, t = 0;

	double tol_1 = 1e-4;
	double tol_2 = 2e-4;
	double * Theta = new double[n_k + n_k*Nsym];

	for (i = 0; i < n_k; i++){
		Theta[t] = p_star[i];
		t++;
	}
	for (k = 0; k < n_k; k++){
		for (j = 0; j < Nsym; j++){
			Theta[t] = lambda[k][j];
			t++;
		}
	}
	if (t != n_k + n_k*Nsym){
		cout << "something fishy with t" << endl;
		exit(0);
	}

	double *cpy_Theta = new double[t];
	double *cpy_Theta_mi = new double[t];
	double *cpy_Theta_2 = new double[t];
	double *cpy_Theta_mi_2 = new double[t];
	// do not touch t; 

	for (i = 0; i < npara; i++){
		for (j = 0; j < npara; j++){
			if (i == j){
				memcpy(cpy_Theta, Theta, t* sizeof(double));
				cpy_Theta[i] = Theta[i] + tol_1;
				memcpy(cpy_Theta_mi, Theta, t* sizeof(double));
				cpy_Theta_mi[i] = Theta[i] - tol_1;
				out(i, j) = (-ln_lca(dat, cpy_Theta_mi, t, no_ind, Nsym, n_k) + 2 * ln_lca(dat, Theta, t, no_ind, Nsym, n_k) - ln_lca(dat, cpy_Theta, t, no_ind, Nsym, n_k)) / (tol_1*tol_1);
			}
			if (i < j){
				memcpy(cpy_Theta, Theta, t* sizeof(double));
				cpy_Theta[i] = Theta[i] + tol_1;
				cpy_Theta[j] = Theta[j] + tol_2;
				memcpy(cpy_Theta_mi, Theta, t* sizeof(double));
				cpy_Theta_mi[i] = Theta[i] + tol_1;
				cpy_Theta_mi[j] = Theta[j] - tol_2;
				memcpy(cpy_Theta_2, Theta, t* sizeof(double));
				cpy_Theta_2[i] = Theta[i] - tol_1;
				cpy_Theta_2[j] = Theta[j] + tol_2;
				memcpy(cpy_Theta_mi_2, Theta, t* sizeof(double));
				cpy_Theta_mi_2[i] = Theta[i] - tol_1;
				cpy_Theta_mi_2[j] = Theta[j] - tol_2;
				out(i, j) = (ln_lca(dat, cpy_Theta, t, no_ind, Nsym, n_k) - ln_lca(dat, cpy_Theta_mi, t, no_ind, Nsym, n_k) - ln_lca(dat, cpy_Theta_2, t, no_ind, Nsym, n_k) + ln_lca(dat, cpy_Theta_mi_2, t, no_ind, Nsym, n_k)) / (4 * tol_1*tol_2);
			}
			else{
				out(i, j) = out(j, i);
			}

		}
	}

	cout << out << endl;

	//check out 


	delete[] cpy_Theta;
	cpy_Theta = nullptr;

	delete[] cpy_Theta_2;
	cpy_Theta_2 = nullptr;
	delete[] cpy_Theta_mi;
	cpy_Theta_mi = nullptr;

	delete[] cpy_Theta_mi_2;
	cpy_Theta_mi_2 = nullptr;

}
double DIC3_LC(double mean_log_dev, double **p_star, double ***lambda, int NUM_KEEP, int K, int Nsym, int no_ind, struct symptom *dat){
	// based on celeux et al (2006)
	// DIC= \hat{D(\theta)}-D(\hat{\theta})
	// DIC= 2 * (average of log deviance) + 2* log of average deviance
	// we already have average of log deviance from earlier 
	// we need to find the second part
	// step 1: find likelihood
	// step 2: sum over all iterations
	// step 3: divided by all iterations
	// step 4: take log
	// IN:  NUM_KEEP: NUM_ITI-NUM_BURN+1
	int i, j, k, m;
	double DIC;
	double f_y = 0;

	for (i = 0; i < no_ind; i++){
		double temp = 0;
		for (m = 0; m < NUM_KEEP; m++){
			for (k = 0; k < K; k++){
				double tt = log(p_star[m][k]);
				for (j = 0; j < Nsym; j++){
					tt += ln_binomial(lambda[m][k][j], dat[i].sym[j]);
				}
				temp += exp(tt);
			}

		}
		f_y += log(temp) - log(NUM_KEEP);
	}

	DIC = -4 * mean_log_dev + 2 * f_y;
	return DIC;

}

double DIC3_GoM(double mean_log_dev, double ***gik, double ***gamma, int NUM_KEEP, int K, int Nsym, int no_ind, struct symptom *dat){
	// based on celeux et al (2006)
	// DIC= \hat{D(\theta)}-D(\hat{\theta})
	// DIC= 2 * (average of log deviance) + 2* log of average deviance
	// we already have average of log deviance from earlier 
	// we need to find the second part
	// step 1: find likelihood
	// step 2: sum over all iterations
	// step 3: divided by all iterations
	// step 4: take log
	// IN:  NUM_KEEP: NUM_ITI-NUM_BURN+1
	int i, j, k, m;
	double DIC;
	double f_y = 0;
	double u;

	for (i = 0; i < no_ind; i++){
		double temp1 = 0;
		for (m = 0; m < NUM_KEEP; m++){
			double temp = 0;
			for (j = 0; j < Nsym; j++){
				double oo = 0;
				for (k = 0; k < K; k++){
					oo += gik[m][i][k] * gamma[m][k][j];
				}
				if (dat[i].sym[j] == 1){
					u = oo;
				}
				else{
					u = 1 - oo; 
				}
				temp += log(u);
			}
			temp1 += exp(temp);
		}
		f_y += log(temp1) - log(NUM_KEEP);
	}
	DIC = -4 * mean_log_dev + 2 * f_y;
	return DIC;

}