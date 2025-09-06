#include <random>
#include <fstream>
#include <iostream>
#include <valarray>
#include <string>

using namespace std;

int main()
{
	const long int Nrep=100000000;
	const int M_EFA = 8;
	const double P_pass_min = 0.7, P_pass_max = 1;
	const double P_EFA_a_min = 0.8, P_EFA_a_max = 1;
	const int MinEFT = 3, MaxEFT = 9;
	long int T_N_EFP[M_EFA+1] = {};
	long int T_N_EFAA_given_N_EFP[M_EFA+1] = {};
	long int T_EFT_given_N_EFP[M_EFA+1] = {};
	double Mn_EFT_given_N_EFP[M_EFA + 1] = {};
	double T_P_pass_given_N_EFP[M_EFA+1] = {};
	const int NSeg = 100;
	double P_pass_given_N_EFP[M_EFA+1][NSeg + 1] = {};
	
	default_random_engine gen;
	// top level distributions
	uniform_real_distribution<double> P_EFA_trial_dist(P_EFA_a_min, P_EFA_a_max);
	uniform_real_distribution<double> P_pass_dist(P_pass_min, P_pass_max);
	uniform_int_distribution<int> N_EFT_dist(MinEFT, MaxEFT);

	for (int i = 0; i < Nrep; i++)
	{
		double P_EFA_trial = P_EFA_trial_dist(gen);
		binomial_distribution<int> N_EFA_trial_dist(M_EFA, P_EFA_trial);
		int N_EFA_trial = N_EFA_trial_dist(gen);

		double P_pass = P_pass_dist(gen);
		int N_fail[M_EFA];
		int Any_fail[M_EFA];
		int N_EFT[M_EFA];


		int itot_fail = 0;
		int iN_EFP = 0;
		int itot_N_EFT = 0;
		double imean_N_EFT = 0;
		for (int EFA = 0; EFA < N_EFA_trial; EFA++)
		{
			N_EFT[EFA] = N_EFT_dist(gen);
			binomial_distribution<int> N_EFA_fail_dist(N_EFT[EFA], 1 - P_pass);
			//N_EFA_fail_dist.t = N_EFT[EFA];
			N_fail[EFA] = N_EFA_fail_dist(gen);
			itot_fail += N_fail[EFA];
			iN_EFP += (N_fail[EFA] > 0);
			itot_N_EFT += N_EFT[EFA];			
		}
		T_N_EFP[iN_EFP]++;
		T_N_EFAA_given_N_EFP[iN_EFP] += N_EFA_trial;
		T_EFT_given_N_EFP[iN_EFP] += itot_N_EFT;
		if (N_EFA_trial > 0)
			Mn_EFT_given_N_EFP[iN_EFP] += 1.0 * itot_N_EFT / N_EFA_trial;
		else
			Mn_EFT_given_N_EFP[iN_EFP] += 1.0 * (MinEFT + MaxEFT) / 2;
		T_P_pass_given_N_EFP[iN_EFP] += P_pass;

		
		int s =  (P_pass - P_pass_min) / (P_pass_max - P_pass_min)*NSeg;
		double rem =  (P_pass - P_pass_min) / (P_pass_max - P_pass_min)*NSeg-s;
		
		P_pass_given_N_EFP[iN_EFP][s] += 1 - rem;
		P_pass_given_N_EFP[iN_EFP][s + 1] += rem;
	}

	ofstream outp("post.csv");

	outp << "Fails_reported,P_Prob,density\n";
	for (int f = 0; f < (M_EFA+1); f++)
	{
		P_pass_given_N_EFP[f][0] *= 2;
		P_pass_given_N_EFP[f][NSeg] *= 2;
		{
			for (int s = 0; s < (NSeg + 1); s++)
				outp << f <<","<< P_pass_min + (P_pass_max - P_pass_min) * s / NSeg << "," << P_pass_given_N_EFP[f][s] / (1.0*T_N_EFP[f]) << "\n";
		}
	}

	//Output results to files
	outp.close();
	outp.open("rep_errs.csv");
	
	outp << "Reported errors, Est. Prob.\n";
	for (int i = 0; i < (M_EFA+1); i++)
	{
		outp << i << "," <<  1.0*T_N_EFP[i]/Nrep << "\n";
	}
	outp.close();
	outp.open("mean_prob.csv");
	outp << "Reported errors, Exp. Prob. Pass\n";
	for (int i = 0; i < (M_EFA + 1); i++)
	{
		outp << i << "," << T_P_pass_given_N_EFP[i]/T_N_EFP[i]  << "\n";
	}
	outp.close();
	outp.open("exp_testers.csv");
	outp << "Reported errors, Exp N testers\n";
	for (int i = 0; i < (M_EFA + 1); i++)
	{
		outp << i << "," << (1.0)*T_N_EFAA_given_N_EFP[i] / T_N_EFP[i] << "\n";
	}
	outp.close();
	outp.open("exp_tests.csv");
	outp << "Reported errors, Exp N tests\n";
	for (int i = 0; i < (M_EFA + 1); i++)
	{
		outp << i << "," << (1.0)*T_EFT_given_N_EFP[i] / T_N_EFP[i] << "\n";
	}
	outp.close();
	outp.open("mean_tests.csv");
	outp << "Reported errors, Exp mean tests\n";
	for (int i = 0; i < (M_EFA + 1); i++)
	{
		outp << i << "," << Mn_EFT_given_N_EFP[i] / T_N_EFP[i] << "\n";
	}
	outp.close();
}