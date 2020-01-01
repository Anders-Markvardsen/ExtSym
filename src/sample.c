# include <stdlib.h>		// (to include exit(1)) 
# include <stdio.h>
# include <string.h>
# include <math.h>
# include "../lib/nrutil.h"
# include "../lib/nr.h"
# include "../lib/Thompson.h"
# include "data_types.h"   // uses typedef yesno in generate()
# include "read_in_parameters.h"  // max_monte_iteration

# include <float.h>		//DBL_MIN - smallest positive number

/*******************************************************************
*  sample --  function that calculates an expectation         
*             value required in the evaluation                          
*             of the evidence of a space group 
*             symmetry. 
*
*  Inputs --  N       : the dimensionality of the integral
*             Cov     : the likelihood covariance matrix of the
*                       intensities
*             Ave     : the likelihood averages of the 
*                       intensities
*
*             Ave_Int : prior mean value
*
*  Returns -- log_Int : the logarithm of the expectation value
*
*  Credits -- Numerical recipies in C
*
********************************************************************/

void reduce_matrix_and_vector(int from_N, int to_N, double **matrix, double *vector,
							  int *project_vec);  // (from integral_evaluation.c)

void calculate_extra_vector(int N, int no_of_pos_intensity, int *project_vec, 
							double **inv_cov, double *extra_mean_values, 
							double *extra_vector);

double generate(int N, double **sqrt_Inv_M, double *Ave, int *multiplet, 
				double **inv_Cov, double *percent_error, double *monte_ave,
				double *monte_extra_mean_values, char highway_name[]);  

double transfer_to_gaussian(int N, double **Cov, double **inv_Cov, double *Ave);

// temp have decleration

double doublet_integration(double Variance, double Average);


// int seed;               // The seed required for gasdev() (trees.c)
double Ave_Int;         // Average prior intensity        (trees.c)


// *********************************************************
// sample() does some of the set up for the main sampling 
// routine generate().  
//
// *percent_error, double*monte_ave, *extra_mean_values and highway_name[]
// are all for diagnosis purposes
// *********************************************************

double sample(int N, double **Cov_keep, double *Ave, int *multiplet, 
							double **inv_Cov, double *percent_error, double *monte_ave, 
							double *extra_mean_values, char highway_name[]) 
{
	int i, j;                     // Simple indices for matrix calculations
	
	double *d;                   // The diagonal elements of the choldsky
	// decomposition of covariance matrix
	
	double log_Int = 0.0;              // The answer
	
	double **Cov;	// used to avoid Cov_keep gets written over in this function
	
	
	d = dvector(1,N);
	
	Cov = dmatrix(1,N,1,N);
	
	for (i = 1; i <= N; ++i)
		for (j = 1; j <= N; ++j)
			Cov[i][j] = Cov_keep[i][j];
		
		
	// choldc() is a Numerical Recipies function that finds the 
	// choldsky decompostion of an NxN matrix and returns the
	// diagonal elements in d and the lower off diagonal elements
	// in the original matrix. The upper diagonal elements of the
	// original matrix are preserved.
		
	choldc(Cov,N,d);
		
	for(i=1;i<=N;i++)
		Cov[i][i] = d[i];
		
		
	// The answer is calculated by generate()
		
	log_Int += generate(N, Cov, Ave, multiplet, inv_Cov, percent_error, monte_ave,
		extra_mean_values, highway_name); 
		
		
	free_dvector(d,1,N);
	free_dmatrix(Cov, 1, N, 1, N);
		
	return log_Int;
}



// *********************************************************
// generate() samples the Gaussian likelihood probability 
// distribution, defined by the covariance matrix and the 
// amplitudes. It returns the expectation value of the
// prior probability distribution.
//
// The additional inv_Cov in the argument list is used when negative Ave's
// *********************************************************

#define MAX_NUM_HIGHWAY 3

double generate(int N, double **sqrt_M, double *Ave_keep, int *multiplet, 
								double **inv_Cov, double *percent_error, double *monte_ave,
								double *monte_extra_mean_values, char highway_name[]) 
{
	int i,j;                  // Simple indices for matrix calculations. 
	
	long n;                    // The sample number
	
	double *random;           // A vector which stores N normal deviates
	
	double *st_samp;          // A vector which stores the N elements of a sample
	
	double prior;             // The value of the prior
	
	double expct_val[MAX_NUM_HIGHWAY+1];         // The expectation value of the prior 
	
	double prior_sqrd[MAX_NUM_HIGHWAY+1];        // The expectation value of the squared prior
	
	double sum_samp;          // The sum of the elements of a single sample
	
	double prior_err[MAX_NUM_HIGHWAY+1];         // The error on the estimate of the expectation 
	// value of the prior	
	
	yesno has_positive_int_occurred[MAX_NUM_HIGHWAY+1];   // to check when the first full sample is generated
	// with all intensities are positive 	       
	
	double extra_zero_term_neg[MAX_NUM_HIGHWAY+1];		// is non-zero when one of the Ave's are negative
	
	double extra_zero_term_pos[MAX_NUM_HIGHWAY+1];		// is non-zero when one of the Ave's are positive
	
	double store_ratio;				// store a ratio; used when positive intensities
	
	double **extra_mean_values;	// it elements are different from Ave_Int when a Ave is negative 
	
	double *temp_mean_values;   // used when cal new extra_mean_values for positive intensities
	
	double *extra_vector;		// used when positive intensities
	
	int *identity_vector;		// unit vector of length N
	
	int *neg_control_vector;			// has element equal to one  for a negative Ave
	//             equal to zero for a positive Ave
	
	int *pos_control_vector;			// has element equal to one  for a positive Ave
	//             equal to zero for a negative Ave
	
	int no_of_pos_intensity = 0; // count number of positive intensities in pos_control_vector
	
	yesno any_negative_intensities = NO;	// checking if there is any negative intensities
	
	yesno change_scale = NO;		// to check if scale must be changed
	
	double prior_zero_point[MAX_NUM_HIGHWAY+1];  // used to avoid zero prior values only
	
	long lowest_extra_mean;			// used to find the index for which the component of
	// extra_mean_values is smallest
	
	double **Ave;					// because I modify Ave in this function
	
	char name_of_highway[MAX_NUM_HIGHWAY+1][80];
	
	int num_different_highway;		// number of different highways for this integral
	
	yesno finished_monte_carlo;		// boolean when monte carlo loop is done and not done
	
	int highway;					// highway index
	
	int highway_lowest;				// find best highway
	
	int seed  =-1000;	// used in gasdev(&seed) 
	// Interesting the program to not compile if set seed
	// to be a const int
	
	// allocate space 
	
	random  = dvector(1,N);
	st_samp = dvector(1,N);
	identity_vector = ivector(1,N);
	neg_control_vector = ivector(1,N);
	pos_control_vector = ivector(1,N);
	extra_mean_values = dmatrix(1, MAX_NUM_HIGHWAY, 1, N);
	temp_mean_values = dvector(1,N);
	extra_vector = dvector(1,N);
	Ave = dmatrix(1, MAX_NUM_HIGHWAY, 1, N);
	
	
	// initialize before creating highways
	
	for (highway = 1; highway <= MAX_NUM_HIGHWAY; highway++)
	{
		extra_zero_term_neg[highway] = 0.0;
		extra_zero_term_pos[highway] = 0.0;
		prior_zero_point[highway] = 0.0;
		expct_val[highway] = 0.0;
		prior_sqrd[highway] = 0.0;
		prior_err[highway] = 0.0;
		has_positive_int_occurred[highway] = NO;
	}
	
	
	///////////////////////
	// Make plain highway
	///////////////////////
	
	for (i = 1; i <= N; i++)
	{
		Ave[1][i] = Ave_keep[i];
		extra_mean_values[1][i] = 1 / Ave_Int;
		prior_zero_point[1] += Ave[1][i]*extra_mean_values[1][i];
	}
	
	strcpy(name_of_highway[1], "plain");
	num_different_highway = 1;
	
	
	///////////////////////////////////////////
	// Make move-negative-intensities highway
	///////////////////////////////////////////
	
	// check if there is any negative intensities and initializing the identity and control 
	// vectors
	
	for (i = 1; i <= N; i++)
	{
		if (Ave[1][i] < 0.0)
		{
			neg_control_vector[i] = 1;
			pos_control_vector[i] = 0;
			any_negative_intensities = YES;
		}
		else if (Ave[1][i] > 0.0)
		{
			neg_control_vector[i] = 0;
			pos_control_vector[i] = 1;
			no_of_pos_intensity++;     // initialized in declaration statement
		}
		else
		{
			neg_control_vector[i] = 0;
			pos_control_vector[i] = 0;
		}
		identity_vector[i] = 1;
	}
	
	// In case any of the Ave intensities are negative, then this can cause great difficulties in
	// finding a sample within the limits of integration. The mode of the multivariate Gaussian 
	// is Ave. If any of the Ave[i] are negative those Ave's are moved to Ave[i]=0. This way the 
	// probability of generating a sample within the limits of integration will always be 
	// acceptable large. However, this is at the cost of some tiny additions overhead by the 
	// calculation of extra_mean_values[i] = 1/mu - [(Px^{new})^T*C^{-1}]_i, and 
	// extra_zero_term_neg = - (Px^{new})^T*C^{-1}*((I-0.5*P)x^{new}). Where P is in the code the 
	// neg_control_vector, I the identity_vector and x the Ave intensities.
	
	
	if (any_negative_intensities == YES)
	{
		// increase number of different highways by one and updata Ave and mean_extra_values
		
		++num_different_highway;
		
		for (i = 1; i <= N; ++i)
		{
			Ave[2][i] = Ave[1][i];
			extra_mean_values[2][i] = extra_mean_values[1][i];
		}
		
		
		// calculate extra zeror term
		
		for (j = 1; j <= N; ++j)
		{
			random[j] = 0.0;	// used random as dummy vector here
			for (i = 1; i <= N; ++i)
			{
				if (neg_control_vector[i] == 1)
				{
					random[j] += Ave[2][i]*inv_Cov[i][j];
					// random contains (Px^{new})^T*C^{-1} here
				}
			}
			extra_zero_term_neg[2] -= random[j]*
				(identity_vector[j]-0.5*neg_control_vector[j])*Ave[2][j];
		}
		
		
		// calculate extra mean value
		
		for (i = 1; i <= N; ++i)
		{
			extra_mean_values[2][i] = 1/Ave_Int - random[i];
		}
		
		
		// make new Ave
		
		for (i = 1; i <= N; ++i)
		{
			if (neg_control_vector[i] == 1)
			{
				Ave[2][i] = 0.0;
			}
			else
			{
				// leave Ave intensity
			}
		}
		
		// make prior_zero_point
		
		for (i = 1; i <= N; ++i)
		{
			prior_zero_point[2] += Ave[2][i]*extra_mean_values[2][i];
		}
		
		// add name
		
		strcpy(name_of_highway[2], "move neg intensities");
	}
	
	
	///////////////////////////////////////////
	// Make move-positive-intensities highway
	///////////////////////////////////////////
	
	// check if any of the extra_mean_values are negative.
	// Use random[1] as dummy variable and also initialize lowest_extra_mean
	
	random[1] = 0.0;
	lowest_extra_mean = 0; 
	
	for (i = 1; i <= N; ++i)
	{
		if (extra_mean_values[num_different_highway][i] < 0.0)
		{
			if (random[1] > extra_mean_values[num_different_highway][i])
			{
				random[1] = extra_mean_values[num_different_highway][i];
				lowest_extra_mean = i;
			}
		}
	} 
	
	// If any positive intensities
	
	if (no_of_pos_intensity != 0 && lowest_extra_mean == 0) 
	{
		// increase number of different highways by one and updata Ave and mean_extra_values
		
		++num_different_highway;
		
		for (i = 1; i <= N; ++i)
		{
			Ave[num_different_highway][i] = Ave[num_different_highway-1][i];
			extra_mean_values[num_different_highway][i] = extra_mean_values[num_different_highway-1][i];
		}
		extra_zero_term_neg[num_different_highway] = extra_zero_term_neg[num_different_highway-1];
		
		
		
		// calculate I_{1r}^T = \mu_r^T*C_rr returned as the 1*N extra_vector
		
		calculate_extra_vector(N, no_of_pos_intensity, pos_control_vector, inv_Cov, 
			extra_mean_values[num_different_highway], extra_vector); 
		
		
		// check if any of the components of extra_vector are bigger than the equivalent
		// components of the Ave vector. If so reduce the size of the extra_vector 
		// component to be equal to the extra_mean_value component
		
		for (i = 1; i <= N; ++i)
		{
			if (pos_control_vector[i] != 0 && extra_vector[i] > Ave[num_different_highway][i])
			{
				extra_vector[i] = Ave[num_different_highway][i];
			}
		}
		
		
		// finally check if any of the components of mu^T - I_1^T*C^{-1} are negative. If this
		// is the case reduce the I_1[i] such that this is no longer the case
		
		change_scale = NO;
		for (i = 1; i <= N; ++i)
		{
			temp_mean_values[i] = 0.0;
			for (j = 1; j <= N; ++j)
			{
				temp_mean_values[i] += extra_vector[j]*inv_Cov[j][i];
			}
			if (temp_mean_values[i] > extra_mean_values[num_different_highway][i] )
			{
				// first time get to this point initialize store_scale and set change_scale = YES 
				
				if (change_scale == YES)
				{
					if (extra_mean_values[num_different_highway][i] / temp_mean_values[i] < store_ratio)
					{
						store_ratio = extra_mean_values[num_different_highway][i] / temp_mean_values[i];
					}
				}
				else 
				{
					// to initialize store_ratio and set change_scale = YES 
					
					store_ratio = extra_mean_values[num_different_highway][i] / temp_mean_values[i];
					change_scale = YES;
				}	
			}
		}
		
		if (change_scale == YES)
		{
			// readjust the I_1 vector
			
			for (i = 1; i <= N; ++i)
			{
				temp_mean_values[i] *= store_ratio;
				extra_vector[i] *= store_ratio;
			}
			
			
			// double check 
			
			for (i = 1; i <= N; ++i)
			{
				if (temp_mean_values[i] > (extra_mean_values[num_different_highway][i] + 1e-10))
				{
					printf("Error in re-scaling\n\n");
					printf("  Please contact author of program Anders J. Markvardsen,\n");
					printf("  since I would be interested in knowing what caused this error\n");
					printf("  to be displayed. Thank you in advance and I apologize for the\n");
					printf("  inconvenience.\n\n");
					
					exit(1);					
				}	
			}
		}
		
		
		// finally set new extra_mean_value
		
		for (i = 1; i <= N; ++i)
		{
			extra_mean_values[num_different_highway][i] -= temp_mean_values[i];
		}
		
		
		// substrate I_{1r} from the positive intensities
		
		for (i = 1; i <= N; ++i)
		{
			if (pos_control_vector[i] != 0)
			{
				Ave[num_different_highway][i] -= extra_vector[i];
			}
		}
		
		
		// additional extra zero constant
		
		for (j = 1; j <= N; ++j)
		{
			random[j] = 0.0;	// used random as dummy vector here
			for (i = 1; i <= N; ++i)
			{
				random[j] += extra_vector[i]*inv_Cov[i][j];
				// random contains (Px^{new})^T*C^{-1} here
			}
			extra_zero_term_pos[num_different_highway] -= random[j]*
				(0.5*extra_vector[j] + Ave[num_different_highway][j]);
		}
		
		
		// make prior_zero_point
		
		for (i = 1; i <= N; ++i)
		{
			prior_zero_point[num_different_highway] += 
				Ave[num_different_highway][i]*extra_mean_values[num_different_highway][i];
		}
		
		// add name
		
		strcpy(name_of_highway[num_different_highway], "move pos intensities");
	}



	// The main sampling loop


	finished_monte_carlo = NO;


	//num_different_highway= 1;
	seed  =-1000;
	n = 0;
	do
	{
		++n;
		
		// Normal deviates are linearly transformed to the required
		// distribution. Since the Choldsky decompositon is lower diagonal
		// Nx(N+1)/2 multiplications are required instead of the NxN that
		// are required for a decomposition. 
		//
		// The elements of the sample are calculated until one is found
		// to be negative. If a negative element is created the sample is 
		// discarded and the loop began again. The Cholsky decompostion's
		// help here is twofold -- Firstly : if the sample gets half way through
		// and find a negative element, only a quarter of the multiplications
		// have been done. Secondly : the number of normal deviates that 
		// must be calculated is equal to the number of elements of the sample
		// calculated. For general decompostion, N deviates must be calculated
		// before one element of a sample can be calculated
		
		
		for (highway = 1; highway <= num_different_highway; ++highway)
		{
			for(i=1;i<=N;i++)
			{
				random[i] = gasdev(&seed);      // generate normal deviate
				
				st_samp[i] = Ave[highway][i];
				for(j=1;j<=i;j++)
					st_samp[i] += sqrt_M[i][j]*random[j];	// transform to
				// required
				// distribution
				
				if(st_samp[i]<0)  goto loop_end;    // if a negative element is found
				// restart the sample loop
			}
			
			// a-ha an sample where all the intensities in the clump are positive has occurred
			
			has_positive_int_occurred[highway] = YES;
			
			
			for(sum_samp=0,i=1;i<=N;i++) 
				sum_samp+=st_samp[i]*extra_mean_values[highway][i];	// calculate argument of modified 
			// prior 
			
			prior = exp(-sum_samp+prior_zero_point[highway]);     // calculate the value of the
			// modified prior for this sample
			
			
			// modification of prior required if doublets are present
						
			for(i=1;i<=N;i++)
				if(multiplet[i]>=1) 
					prior *= ( pow(st_samp[i]/Ave_Int, multiplet[i]) ) / ( (double) factrl(multiplet[i]) );
					
				
				expct_val[highway]  += prior;                // and add to the running totals
				prior_sqrd[highway] += prior*prior;                      
				
loop_end:             
				
				
				// find the square of the fractional error on the estimate
				
				if (has_positive_int_occurred[highway] == YES && expct_val[highway] >= DBL_MIN)
				{
					prior_err[highway]  = prior_sqrd[highway]/(expct_val[highway]*expct_val[highway]) - (1.0/n);
				}
				
				
				// if the desired accuracy has been reached (after at least n samples)
				// then the sampling has finished
				
				if( (sqrt(prior_err[highway])<=TOL) && (n>=1000) && 
					has_positive_int_occurred[highway]==YES
					&& expct_val[highway] >= DBL_MIN)   
				{
					finished_monte_carlo = YES;
					break;  
				}
				
				// if the accuracy has not been reached
				
				if(n == max_monte_iteration)
				{
					if (has_positive_int_occurred[highway] == YES && expct_val[highway] >= DBL_MIN)
					{
						//printf("\n\n ***********************************************************\n\n");
						//printf(" warning : percentage error = %.3f \t (tolerance = %.3f)\n\n"
						//             , 100*sqrt(prior_err[highway]), 100*TOL);
						//printf(" ***********************************************************\n\n");
					}
					else if (has_positive_int_occurred[highway] == YES && expct_val[highway] < DBL_MIN)
					{
						//printf("\n\n ***********************************************************\n\n");
						//printf(" warning : all prior sample values equal to zero\n");
						//printf(" ***********************************************************\n\n");
						
						prior_err[highway] = 20*20;
					}
					else
					{
						//printf("\n\n ***********************************************************\n\n");
						//printf(" warning : no sample with all intensities positive in clump was");
						//printf(" generated\n");
						//printf(" ***********************************************************\n\n");
						
						prior_err[highway] = 10*10;
					}
					// if in addition highway == num_different_highway then stop monte carlo
					
					if (highway == num_different_highway)
					{
						finished_monte_carlo = YES;
						break;
					}
				}		
		}
} while (finished_monte_carlo == NO);	


// if not successful in Monte carlo then find highway which was best

if (n == max_monte_iteration && highway == num_different_highway)
{ 
	// find highway which has lowest prior_err
	
	highway_lowest = 1;
	random[1] = prior_err[highway_lowest];				// use random[1] as dummy
	for (highway = 2; highway <= num_different_highway; ++highway)
	{
		if (prior_err[highway] < random[1])
		{
			highway_lowest = highway;
			random[1] = prior_err[highway];
		}
	}
	
	// highway which was best is
	
	highway = highway_lowest;
}


// store diagnosis of integral

*percent_error = 100.0 * sqrt(prior_err[highway]);
strcpy(highway_name, name_of_highway[highway]); 

for (i = 1; i<=N; ++i)
{
	monte_ave[i] = Ave[highway][i];
	monte_extra_mean_values[i] = extra_mean_values[highway][i];
}


// deallocate space

free_dvector(random , 1, N);
free_dvector(st_samp, 1, N);
free_ivector(identity_vector, 1, N);
free_ivector(neg_control_vector, 1, N);
free_ivector(pos_control_vector, 1, N);
free_dmatrix(extra_mean_values, 1, MAX_NUM_HIGHWAY, 1, N);
free_dvector(temp_mean_values, 1, N);
free_dvector(extra_vector, 1, N);
free_dmatrix(Ave, 1, MAX_NUM_HIGHWAY, 1, N);

expct_val[highway] /= ( (double) n);        // calculate the expectation value of the modified prior


// calculate the log of the expected value of the prior. 
// log(prior value) = log(modified prior value) + extra_zero_term_neg 

expct_val[highway] = log(expct_val[highway]) + extra_zero_term_neg[highway] + 
extra_zero_term_pos[highway] - prior_zero_point[highway];


return expct_val[highway];      // return the log of the expected value of the prior
}



// ********************************************************************
// A small routine which takes the likelihood variance, average and multiplet of
// an intensity and returns the expectation value of the prior over
// the likelihood, which may be written as <prior>_gaus.
// Notice at the moment is return [mu * <prior>_gaus]
// ********************************************************************

double direct_integration(double Variance, double Average, int multiplet)
{
double answer;                    // the output
double erfc_argument;            // the argument of the erfc function

const double MY_PI = acos(-1);


if (multiplet >= 1)
{
	// This expression is derived from Eq. (18) in https://aip.scitation.org/doi/abs/10.1063/1.1835216
	
	return -0.5*log(2*MY_PI) + multiplet*log(sqrt(Variance)/Ave_Int) - Average*Average/(2*Variance) 
		+ ParCylU(multiplet+0.5, sqrt(Variance)/Ave_Int - Average/sqrt(Variance), 0.000001);
}


// if singlet return

answer =  - log(2) + (Variance/(2*Ave_Int) - Average)/Ave_Int;


// calculate the arguement of the erfc function

erfc_argument = (Variance/Ave_Int - Average)/sqrt(2*Variance);


// erfcc is a Numerical recipies routine that
// calculates the erfc function to an accuracy of better than 1e-7
// over the entire range using Chebyshev fitting. 

// the erfcc c-code function has problems when the argument becomes too large
// e.g. we have erfcc(10.0) = 2.1e-45 and erfcc(20.0) = 5.4e-176.
// Therefore the log(erfc(x)) is calculated directly for large positive argument.
// I have not taken into considereation the case when the argument might be large 
// and negative but hopefully this will not happen.
// Also the below log_erfcc() function below was derived directly from erfcc(). I 
// have no feeling for its accuracy.

if (erfc_argument < 10.0)
{
	answer += log(erfcc(erfc_argument));
}
else
{
	answer += log_erfcc(erfc_argument);
}

return answer;
}



// *********************************************************************
// A small routine which takes the likelihood variance and average of
// an intensity and returns the expectation value of the doublet prior 
// over the likelihood
// *********************************************************************

double doublet_integration(double Variance, double Average)
{
double answer;                    // the output
double erfc_argument;            // the arguement of the effc function
double exp_argument;             // the arguement of the exponent term 
double a;                         // a storage variable
double x;

const double MY_PI = acos(-1);

// put in the prefactor

answer = (Variance/(2*Ave_Int) - Average)/Ave_Int - log(Ave_Int);


// calculate the arguement of the erfc function

a = Average - (Variance/Ave_Int);

erfc_argument = -a/sqrt(2*Variance);

exp_argument  = -erfc_argument*erfc_argument;

// erfcc is a Numerical recipies routine that
// calculates the erfc function to an accuracy of better than 1e-7
// over the entire range using Chebyshev fitting. 

x = sqrt(Variance/(2*MY_PI))*exp(exp_argument);

printf("\n%g\n", 0.5*log(MY_PI/2)  - exp_argument + log(
			 sqrt(2/MY_PI)*exp(exp_argument) - sqrt(2)*erfc_argument*
			 erfcc(erfc_argument)));
printf("\n%g %g\n", 0.5*log(MY_PI/2) + log(sqrt(2/MY_PI)), erfc_argument);


x+= (a/2.0)*erfcc(erfc_argument);

answer += log(x);



return answer;
}


//----------------------------------------------------------------------------
// Modify the Ave's and returns additional constant.
// Return value = +0.5*X^(modify)C^-1X^(modify)-0.5*X^0C^-1X^0 
//----------------------------------------------------------------------------
double transfer_to_gaussian(int N, double **Cov, double **inv_Cov, double *Ave)
{
	long i,j;								// index variables

	double additional_constant = 0.0;		// return constant

	double x;								// dummy variable


	// add -X^0C^-1X^0 to additional_constant

	for (i = 1; i <= N; i++)  
	{
		x = 0.0;
		for (j = 1; j <= N; j++)
		{
			x += inv_Cov[i][j] * Ave[j];
		}
		additional_constant -= Ave[i] * x;
	}
	
	
	// calculate new gaussian averages.  
	// X^(modify) = X^0-I^T*C , where I is here a unit N*1 vector.

	for (i = 1; i <= N; i++)  
	{
		x = 0.0;
		for (j = 1; j <= N; j++)
		{
			x += Cov[j][i];
		}
		
		Ave[i] -= x / Ave_Int;	// Ave is modified here to new values
	}
	

	// add X^(modify)C^-1X^(modify) to additional_constant

	for (i = 1; i <= N; i++)  
	{
		x = 0.0;
		for (j = 1; j <= N; j++)
		{
			x += inv_Cov[i][j] * Ave[j];
		}
		additional_constant += Ave[i] * x;
	}

	// multiply additional_constant by a half

	additional_constant *= 0.5;

	
	// return additional constant

	return additional_constant;
}


//--------------------------------------------------------------------------------------
// Calculate I_{1r}^T = \mu_r^T*C_rr 
//
// Returns I_{1r}^T in extra_vector which N components and no_of_pos_intensity which are
// non-zero.
//--------------------------------------------------------------------------------------

void calculate_extra_vector(int N, int no_of_pos_intensity, int *project_vec, 
							double **inv_cov, double *extra_mean_values, 
							double *extra_vector)
{
	int i, j;

	int pos_counter;		// counts no of positive intensities

	double **reduce_matrix;

	double *reduce_vector;

	
	// allocate space

	reduce_matrix = dmatrix(1, N, 1, N);
	reduce_vector = dvector(1, N);

	
	// initialize reduce_matrix and reduce_vector

	for (i = 1; i <= N; ++i)
	{
		reduce_vector[i] = extra_mean_values[i];
		for (j = 1; j <= N; ++j)
		{
			reduce_matrix[i][j] = inv_cov[i][j];
		}
	}
	
	
	// generate C_rr^{-1} and \mu_r^T

	reduce_matrix_and_vector(N, no_of_pos_intensity, reduce_matrix, reduce_vector, 
		project_vec);


	// generate output vector

	pos_counter = 0;
	for (i = 1; i <= N; ++i)     
	{
		extra_vector[i] = 0.0;
		if (project_vec[i] != 0)
		{
			pos_counter++;
			for (j = 1; j <= no_of_pos_intensity; ++j)
			{
				extra_vector[i] += reduce_vector[j]*reduce_matrix[j][pos_counter];
			}
		}
		else 
		{
			extra_vector[i] = 0.0;
		}
	}

	// check if pos_counter has the right value here

	if (pos_counter != no_of_pos_intensity)
	{
		printf("Error in calculate_extra_vector() in sample.c\n\n");
		exit(1);
	}


	// free space

	free_dmatrix(reduce_matrix, 1, N, 1, N);
	free_dvector(reduce_vector, 1, N);
}
