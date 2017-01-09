// Bayesian LCA and GoM
// Part 1 of Method 1 of the BMA paper
// Need to run this code twice,
// 1) K=3
// 2) K=4
// OUTPUT:
// 1) MEMBERSHIP save as Pik.txt
// 2) LAMBDA, save as lambdaik.txt
// 3) weights- DIC and Laplace-Gibbs

#include<stdio.h>
#include <stdlib.h>
#include <string>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cassert>
#include <omp.h>
#include <ctime>
#include<sys/types.h>
#include <random>
#include <thread>
#include "Dist.h"
#include "otherDist.h"

using namespace std;
#define N_sym 12 // define total number of symptoms



void usage()
{
	printf("1) LCA=0 or GoM=1\n2) K\n3) Number of Iterations\n4) Burn-in size\n");
	exit(1);
}

const int get_int(const string& s) {
	stringstream ss(s);
	int ret;
	ss >> ret;
	return ret;
}

double round(double value)
{
	return floor(value * 10e5 + 5e-11) / 10e5;
}

int count_Data(ifstream &filename)// struct weatherdata *enviroData)
{
	int numrec = 0;
	string lines;
	if (filename){
		// ------------ read and store data into enviroData-------------
		while (getline(filename, lines))
		{
			numrec++;
		}

	}
	else{ cout << "Problem in opening file.\n"; exit(0); }

	return numrec;
}


void get_sym_data(ifstream &filename, struct symptom *out){

	int i;
	string lines;

	if (filename){
		i = 0;
		while (getline(filename, lines)){
			stringstream ss(lines);
			vector<std::string> dat;
			string field;
			while (getline(ss, field, '\t')) {
				if (field.empty()){
					field = "-9999";// missing value
				}
				dat.push_back(field);
			}
			out[i].Fam_ID = get_int(dat[0]);
			out[i].ID = get_int(dat[1]);
			out[i].p_id = get_int(dat[2]);
			out[i].m_id = get_int(dat[3]);
			out[i].sex = get_int(dat[4]);
			out[i].t_pheno = get_int(dat[5]);
			out[i].sym[0] = get_int(dat[6]);
			out[i].sym[1] = get_int(dat[7]);
			out[i].sym[2] = get_int(dat[8]);
			out[i].sym[3] = get_int(dat[9]);
			out[i].sym[4] = get_int(dat[10]);
			out[i].sym[5] = get_int(dat[11]);
			out[i].sym[6] = get_int(dat[12]);
			out[i].sym[7] = get_int(dat[13]);
			out[i].sym[8] = get_int(dat[14]);
			out[i].sym[9] = get_int(dat[15]);
			out[i].sym[10] = get_int(dat[16]);
			out[i].sym[11] = get_int(dat[17]);
			i++;
			vector<string>().swap(dat);
		}

	}
	else{ cout << "Problem in opening sympton data.\n"; exit(0); }
}

int main(int argc, char** argv) {

	if (argc != 5)usage();
	int MODEL = get_int(argv[1]);
	int K = get_int(argv[2]); // number of clusters
	int NUM_ITI = get_int(argv[3]);
	int NUM_BURN = get_int(argv[4]);
	ifstream aFile, bFile;
	ofstream myFile;



	// read in symptom data

	int no_ind;
	aFile.open("phenosub.txt", ios::in);
	no_ind = count_Data(aFile);
	cout << "Number of subjects " << no_ind << endl;
	aFile.clear();
	aFile.seekg(0);

	struct symptom *symp_data;
	symp_data = new struct symptom[no_ind]; // check this 
	get_sym_data(aFile, symp_data);
	aFile.close();



	int i, j, k, iti;

	// read in symptom data- for calculating posterior predictive probability

	int no_ind_ppp;
	bFile.open("ppp_pheno.txt", ios::in);
	no_ind_ppp = count_Data(bFile);
	cout << "Number of PPP subjects " << no_ind_ppp << endl;
	bFile.clear();
	bFile.seekg(0);
	struct symptom *ppp_symp_data;
	ppp_symp_data = new struct symptom[no_ind_ppp]; // check this 

	get_sym_data(bFile, ppp_symp_data);
	bFile.close();

	base_generator_type eng_global;
	eng_global.seed(static_cast<unsigned int>(time(0)));
	Context dist;

	int nk = NUM_ITI - NUM_BURN;
	int ctr = 0;


	if (MODEL == 0){
		cout << "LCA MODLE" << endl;
		remove("paras_LCA.txt");
		remove("DIC_LCA.txt");
		remove("pik_LCA.txt");
		myFile.open("paras_LCA.txt", ios::app);
		myFile << "iti" << "\t" << "deviance" << "\t" << "BIC" << "\t" << "ln_PPP" << "\t";
		for (k = 0; k < K; k++){
			myFile << "p_" << k << "\t";
		}
		for (k = 0; k < K; k++){
			for (j = 0; j < N_sym; j++){
				myFile << "lambda_" << k << j << "\t";
			}
		}
		myFile << endl;
		myFile.close();

		myFile.open("pik_LCA.txt", ios::app);
		myFile << "iti" << "\t";
		for (i = 0; i < no_ind; i++){
			for (k = 0; k < K; k++){
				myFile << "p_" << i << "_" << k << "\t";
			}
		}
		myFile << endl;
		myFile.close();



		// LCA
		double *wt = new double[K]; // store the parameters of dirichlet
		double *p = new double[K];// store p_k
		double **pik; // probability belonging to each cluster
		int **zik; // binary 0/1
		double **lambda;
		int no_parameters = K*N_sym + (K - 1);
		double ppp;
		double sum_ln_dev = 0;


		pik = new double *[no_ind];
		zik = new int *[no_ind];
		lambda = new double *[K];
		for (i = 0; i < no_ind; i++){
			pik[i] = new double[K];
			zik[i] = new int[K];
		}
		// set up p_k^m
		// set up lambda_kj^m
		// m is sample size of the posterior samples
		double **p_m = new double *[nk];// store p_k
		double ***lambda_m = new double **[nk];

		for (i = 0; i < nk; i++){
			p_m[i] = new double[K];
			lambda_m[i] = new double *[K];
			for (k = 0; k < K; k++){
				lambda_m[i][k] = new double[N_sym];
			}
		}
		/////// end 

		// setting up initial values
		double *start_pk = new double[K];

		for (k = 0; k < K; k++){
			start_pk[k] = 1;
			lambda[k] = new double[N_sym];
			for (j = 0; j < N_sym; j++){
				lambda[k][j] = dist.uniform_real_dis(0, 1, eng_global);
			}
		}

		dirichlet_sample(K, start_pk, p, eng_global);

		// start iteration

		for (iti = 0; iti < NUM_ITI; iti++){

			for (i = 0; i < no_ind; i++){
				double part_b = 0;
				for (k = 0; k < K; k++){
					pik[i][k] = log(p[k]);
					for (j = 0; j < N_sym; j++){
						pik[i][k] += ln_binomial(lambda[k][j], symp_data[i].sym[j]);
					}
					part_b += exp(pik[i][k]);
				}

				for (k = 0; k < K; k++){
					pik[i][k] -= log(part_b);
					pik[i][k] = round(exp(pik[i][k]));
				}

				multi_boost(1, pik[i], K, zik[i], eng_global);
				if (iti >= NUM_BURN){
					//put p_ik here- of each iteration
					myFile.open("pik_LCA.txt", ios::app);
					myFile << iti << "\t";
					for (i = 0; i < no_ind; i++){
						for (k = 0; k < K; k++){
							myFile << pik[i][k] << "\t";
						}
					}
					myFile << endl;
					myFile.close();
				}
			}
			// update p[k]
			for (k = 0; k < K; k++){
				wt[k] = 1;
				for (i = 0; i < no_ind; i++){
					wt[k] += zik[i][k];
				}
			}

			dirichlet_sample(K, wt, p, eng_global);
			// update lambda[k][j]
			for (k = 0; k < K; k++){
				for (j = 0; j < N_sym; j++){
					int a = 1, b = 1;
					for (i = 0; i < no_ind; i++){
						if (symp_data[i].sym[j] == 1){
							a += zik[i][k];
						}
						else{
							b += zik[i][k];
						}
					}
					lambda[k][j] = dist.beta_dis(a, b, eng_global);
				}
			}

			if (iti >= NUM_BURN){
				// store p_k in p_m & lambda[k][j] in lambda [k][j][m]

				for (k = 0; k < K; k++){
					p_m[ctr][k] = p[k];
					for (j = 0; j < N_sym; j++){
						lambda_m[ctr][k][j] = lambda[k][j];
					}
				}
				ctr++;

				double deviance = 0, de = 0;
				double BIC;
				// write likelihood 
				for (i = 0; i < no_ind; i++){
					double t2 = 0;
					for (k = 0; k < K; k++){
						double tt = 0;
						for (j = 0; j < N_sym; j++){
							tt += ln_binomial(lambda[k][j], symp_data[i].sym[j]);
						}
						tt += log(p[k]);
						t2 += exp(tt);
					}
					deviance += log(t2);
				}

				sum_ln_dev += deviance;
				deviance = -2 * deviance;
				BIC = deviance + (no_parameters)*log(no_ind);


				// posterior predictive probability
				ppp = 0;
				for (i = 0; i < no_ind_ppp; i++){
					double pp2 = 0;
					for (k = 0; k < K; k++){
						double pp1 = 0;
						for (j = 0; j < N_sym; j++){
							pp1 += ln_binomial(lambda[k][j], ppp_symp_data[i].sym[j]);
						}
						pp1 += log(p[k]);
						pp2 += exp(pp1);
					}
					ppp += log(pp2);
				}
				ppp += log(2); // log scale
				cout << " ln_ppp " << ppp << endl;


				myFile.open("paras_LCA.txt", ios::app);
				//myFile << t << "\t" << rho_out[t] << "\t" << choice[t].beta << "\t" << choice[t].detect << "\t" << choice[t].jv_sur << "\t" << choice[t].k << "\t" << choice[t].larv_mor << "\t" << choice[t].x0 << endl;
				myFile << iti << "\t" << deviance << "\t" << BIC << "\t" << ppp << "\t";
				for (k = 0; k < K; k++){
					myFile << p[k] << "\t";
				}
				for (k = 0; k < K; k++){
					for (j = 0; j < N_sym; j++){
						myFile << lambda[k][j] << "\t";
					}
				}
				myFile << endl;
				myFile.close();
			}



		}// end of iteration

		// calculate DIC

		double mean_log_dev = sum_ln_dev / nk;
		double DIC;
		DIC = DIC3_LC(mean_log_dev, p_m, lambda_m, nk, K, N_sym, no_ind, symp_data);
		cout << "DIC LCA " << DIC << endl;
		myFile.open("DIC_LCA.txt", ios::app);
		myFile << DIC << endl;
		myFile.close();

		// calculate Laplace-gibbs
		// do it in R
		// use the posterior mean as theta*
		// esimated the H* use cov.mve 



		// free pointers

		for (i = 0; i < no_ind; i++){
			delete[] pik[i];
			delete[] zik[i];
			pik[i] = nullptr;
			zik[i] = nullptr;
		}
		for (i = 0; i < K; i++){
			delete[] lambda[i];
			lambda[i] = nullptr;
		}
		delete[] p;
		delete[] wt;
		p = nullptr;
		wt = nullptr;

		for (i = 0; i < nk; i++){
			delete[] p_m[i];
			p_m[i] = nullptr;
			for (k = 0; k < K; k++){
				delete[] lambda_m[i][k];
				lambda_m[i][k] = nullptr;
			}
			delete[] lambda_m[i];
			lambda_m[i] = nullptr;
		}
		delete[] p_m;
		p_m = nullptr;
		delete[] lambda_m;
		lambda_m = nullptr;

	}
	else{
		// GoM
		// setting up starting values
		// two model parameters, g_ik and gamma_kj
		cout << "GOM Model " << endl;
		remove("paras_GOM_gamma.txt");
		remove("paras_GOM_gik.txt");
		remove("DIC_GoM.txt");

		myFile.open("paras_GOM_gamma.txt", ios::app);
		myFile << "iti" << "\t" << "deviance" << "\t" << "BIC" << "\t";
		for (k = 0; k < K; k++){
			for (j = 0; j < N_sym; j++){
				myFile << "gamma_" << k << j << "\t";
			}
		}
		for (i = 0; i < no_ind; i++){
			for (k = 0; k < K; k++){
				myFile << "g_" << i << "_" << k << "\t";
			}
		}
		myFile << endl;
		myFile.close();

		myFile.open("paras_GOM_gik.txt", ios::app);
		myFile << "iti" << "\t" << "ID" << "\t";
		for (k = 0; k < K; k++){
			myFile << "g_" << k << "\t";
		}
		myFile << endl;
		myFile.close();

		// store g_ik^m and gammakj^m after burn-in
		double ***gik_m = new double **[nk];
		double ***gammakj_m = new double **[nk];
		int l;
		for (l = 0; l < nk; l++){
			gik_m[l] = new double *[no_ind];
			for (i = 0; i < no_ind; i++){
				gik_m[l][i] = new double[K];
			}
			gammakj_m[l] = new double *[K];
			for (k = 0; k < K; k++){
				gammakj_m[l][k] = new double[N_sym];
			}

		}

		double **gik;
		double **start_gik; // starting weight of dirichlet
		int ***wijk; //omega ik
		double **gammakj;
		double ***kijk;
		double no_parameters_gom = no_ind*(K - 1) + K*N_sym;

		// para
		double *dd = new double[K]; // dirichlet para
		double a1, b1; // beta para
		double BIC_gom;
		double sum_ln_dev = 0;


		gik = new double *[no_ind];
		start_gik = new double *[no_ind];
		wijk = new int **[no_ind];
		kijk = new double **[no_ind];
		gammakj = new double *[K];

		for (i = 0; i < no_ind; i++){
			gik[i] = new double[K];
			start_gik[i] = new double[K];
			kijk[i] = new double*[N_sym];
			wijk[i] = new int*[N_sym];
			for (j = 0; j < N_sym; j++){
				kijk[i][j] = new double[K];
				wijk[i][j] = new int[K];
			}

			for (k = 0; k < K; k++){
				start_gik[i][k] = dist.uniform_int_dis(1, 5, eng_global);
			}
			dirichlet_sample(K, start_gik[i], gik[i], eng_global);
		}

		for (k = 0; k < K; k++){
			gammakj[k] = new double[N_sym];
			for (j = 0; j < N_sym; j++){
				if (k == 0){
					gammakj[k][j] = dist.uniform_real_dis(0, 0.1, eng_global);
				}
				else{
					gammakj[k][j] = dist.uniform_real_dis(0, 1, eng_global);
				}
			}
		}
		// start iteration

		double top, bot;
		int k1;
		for (iti = 0; iti < NUM_ITI; iti++){

			double **temp = new double *[N_sym];

			for (i = 0; i < no_ind; i++){
				/* this is my initial code- this is what works, however this is not what Erosheva did
				for (j = 0; j < N_sym; j++){
					bot = 0;
					for (k = 0; k < K; k++){
						kijk[i][j][k] = log(gik[i][k]) + ln_binomial(gammakj[k][j], symp_data[i].sym[j]);
						bot += exp(kijk[i][j][k]);
					}
					for (k = 0; k < K; k++){
						kijk[i][j][k] = exp(kijk[i][j][k] - log(bot));
					}
					multi_boost(1, kijk[i][j], K, wijk[i][j], eng_global);
				}
				*/

				// this is the correct code- this is what works, change likelihood
				for (j = 0; j < N_sym; j++){
					bot = 0;
					for (k = 0; k < K; k++){
						if (symp_data[i].sym[j] == 1){
							kijk[i][j][k] = gik[i][k] * gammakj[k][j];
						}
						else{
							kijk[i][j][k] = gik[i][k] * (1 - gammakj[k][j]);
						}
						bot += kijk[i][j][k];
					}
					for (k = 0; k < K; k++){
						kijk[i][j][k] = kijk[i][j][k] / bot;
					}
					multi_boost(1, kijk[i][j], K, wijk[i][j], eng_global);
				}

				/* This is what Erosheva did, it does not work- the likelihood is not used!!
				for (j = 0; j < N_sym; j++){
				multi_boost(1, gik[i], K, wijk[i][j], eng_global);
				}
				*/

				for (k = 0; k < K; k++){
					dd[k] = 1;
					for (j = 0; j < N_sym; j++){
						dd[k] += wijk[i][j][k];
					}
				}
				dirichlet_sample(K, dd, gik[i], eng_global);
			}
			// sort of gamma
			for (k = 0; k < K; k++){
				for (j = 0; j < N_sym; j++){
					a1 = 1;
					b1 = 1;
					for (i = 0; i < no_ind; i++){
						if (symp_data[i].sym[j] == 1){
							a1 += wijk[i][j][k];
						}
						else{
							b1 += wijk[i][j][k];
						}
					}
					gammakj[k][j] = dist.beta_dis(a1, b1, eng_global);
				}
			}

			if (iti >= NUM_BURN){
				// store
				for (i = 0; i < no_ind; i++){
					for (k = 0; k < K; k++){
						gik_m[ctr][i][k] = gik[i][k];
					}
				}
				for (k = 0; k < K; k++){
					for (j = 0; j < N_sym; j++){
						gammakj_m[ctr][k][j] = gammakj[k][j];
					}
				}
				ctr++;
				// deviance
				// -2 log(p(y|\hat{theta})
				double deviance = 0;
				double u;


				for (i = 0; i < no_ind; i++){
					for (j = 0; j < N_sym; j++){
						double oo = 0;
						for (k = 0; k < K; k++){
							oo += gik[i][k] * gammakj[k][j];
						}
						if (symp_data[i].sym[j] == 1){
							u = oo;
						}
						else{
							u = 1 - oo;
						}
						deviance += log(u);
					}
				}
				sum_ln_dev += deviance;
				deviance = -2 * deviance;
				BIC_gom = deviance + (no_parameters_gom)*log(no_ind);

			

				myFile.open("paras_GOM_gamma.txt", ios::app);
				//myFile << t << "\t" << rho_out[t] << "\t" << choice[t].beta << "\t" << choice[t].detect << "\t" << choice[t].jv_sur << "\t" << choice[t].k << "\t" << choice[t].larv_mor << "\t" << choice[t].x0 << endl;
				myFile << iti << "\t" << deviance << "\t" << BIC_gom << "\t";

				for (k = 0; k < K; k++){
					for (j = 0; j < N_sym; j++){
						myFile << gammakj[k][j] << "\t";
					}
				}
				for (i = 0; i < no_ind; i++){
					for (k = 0; k < K; k++){
						myFile << gik[i][k] << "\t";
					}
				}
				myFile << endl;
				myFile.close();

				myFile.open("paras_GOM_gik.txt", ios::app);
				for (i = 0; i < no_ind; i++){
					myFile << iti << "\t" << i << "\t";
					for (k = 0; k < K; k++){
						myFile << gik[i][k] << "\t";
					}
					myFile << endl;
				}

				myFile.close();


			}

		}// end iteration

		/*for (l = 0; l < nk; l++){
			for (i = 0; i < no_ind; i++){
			for (k = 0; k < K; k++){
			cout << "gik" << i << "_" << k << " " << gik_m[l][i][k] << "\t";
			}
			}
			for (k = 0; k < K; k++){
			for (j = 0; j < N_sym; j++){
			cout << "gamma" << k << "_" << j << " " << gammakj_m[l][k][j] << "\t" ;
			}
			}
			cout << endl;
			}*/


		// calculate DIC
		double mean_log_dev = sum_ln_dev / nk;
		double DIC;
		DIC = DIC3_GoM(mean_log_dev, gik_m, gammakj_m, nk, K, N_sym, no_ind, symp_data);
		// print results
		cout << " DIC GOM " << DIC << endl;

		myFile.open("DIC_GoM.txt", ios::app);
		myFile << DIC << endl;
		myFile.close();


		for (i = 0; i < no_ind; i++){

			for (j = 0; j < N_sym; j++){
				delete[] wijk[i][j];
				delete[] kijk[i][j];
				wijk[i][j] = nullptr;
				kijk[i][j] = nullptr;
			}
			delete[] gik[i];
			delete[] start_gik[i];
			delete[] wijk[i];
			delete[] kijk[i];
			gik[i] = nullptr;
			start_gik[i] = nullptr;
			wijk[i] = nullptr;
			kijk[i] = nullptr;
		}

		for (k = 0; k < K; k++){
			delete[] gammakj[k];
			gammakj[k] = nullptr;
		}
		delete[] gammakj;
		gammakj = nullptr;
		delete[] start_gik;
		start_gik = nullptr;
		delete[] gik;
		gik = nullptr;
		delete[] wijk;
		wijk = nullptr;
		delete[] kijk;
		kijk = nullptr;

		for (i = 0; i < nk; i++){
			for (k = 0; k < K; k++){
				delete[] gammakj_m[i][k];
				gammakj_m[i][k] = nullptr;
			}
			delete[] gammakj_m[i];
			gammakj_m[i] = nullptr;
			for (l = 0; l < no_ind; l++){
				delete[] gik_m[i][l];
				gik_m[i][l] = nullptr;
			}
			delete[] gik_m[i];
			gik_m[i] = nullptr;
		}
		delete[] gammakj_m;
		gammakj_m = nullptr;
		delete[] gik_m;
		gik_m = nullptr;

	}

	return 0;
}// end of main
